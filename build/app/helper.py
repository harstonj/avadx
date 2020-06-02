import os
import sys
import shutil
import gzip
import csv
import threading
import time
import re
import subprocess
import errno
import signal
import argparse
import shlex
from functools import wraps
from concurrent.futures import ThreadPoolExecutor
from collections import OrderedDict


FLAGS = []
QUIET = False
ERROR_STATE = False
VALID_FILETYPES = ["fasta", "fa", "fna", "fastq", "fq"]
VALID_ARCHIVETYPES = ["zip", "bzip2", "xz", "tar", "tar.gz", "gz", "bz2"]
VALID_URLS = ("http://", "https://", "ftp://", "ftps://")


flatten = lambda l: [item for sublist in l for item in sublist]


# noinspection
class ClubberTimeoutError(Exception):
    pass


class ReadableDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        prospective_dir = values
        if not os.path.isdir(prospective_dir):
            raise argparse.ArgumentTypeError(f"readable_dir:{prospective_dir} is not a valid path")
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace, self.dest, prospective_dir)
        else:
            raise argparse.ArgumentTypeError(f"readable_dir:{prospective_dir} is not a readable dir")


class Error(Exception):
    """Base class for other exceptions"""
    pass


class InputError(Error):
    """Raised when the input file is not valid"""

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(_, __):
            raise ClubberTimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator


def get_platform():
    if sys.platform == "linux" or sys.platform == "linux2":
        return 'linux'
    elif sys.platform == "darwin":
        return 'osx'
    elif sys.platform == "win32" or sys.platform == "win64":
        return 'win'
    return None


def log(msg, quiet, msgstack=list()):
    if not quiet:
        print(msg)
    else:
        ts = '{:23s}'.format(f'[{time.time()}]')
        msgstack.append(f'{ts} {msg}')


def file_len(fname):
    i = -1
    with open(fname) as f_:
        for i, _ in enumerate(f_):
            pass
    return i + 1


def is_int(x):
    try:
        int(x)
        return True
    except ValueError:
        return False


def fasta_header_count(fname):
    header_count = 0

    is_fasta, _ = check_filetype(fname, [
        "fasta", "fa", "fna",
        "fasta.gz", "fa.gz", "fna.gz",
        "fasta.tar.gz", "fa.tar.gz", "fna.tar.gz"
    ])
    is_fastq, _ = check_filetype(fname, [
        "fastq", "fq",
        "fastq.gz", "fq.gz",
        "fastq.tar.gz", "fq.tar.gz"
    ])

    if is_fasta:
        if fname.endswith('.gz'):
            zgrep = which("zgrep")
            if zgrep:
                cmd = f'{zgrep} -c -e "^>" {fname}'
                zgrep_count = run_command(cmd, shell=True, print_output=False)
                if zgrep_count:
                    zgrep_count_val = zgrep_count[0]
                    if is_int(zgrep_count_val):
                        return int(zgrep_count_val)
        else:
            with open(fname) as f_:
                for l in f_:
                    if l and l[0] == '>':
                        header_count += 1
                        next(f_, None)
    elif is_fastq:
        if fname.endswith('.gz'):
            cat = which("zcat")
        else:
            cat = which("cat")

        wc = which("wc")
        if cat and wc:
            cmd = f'{cat} {fname} | {wc} -l'
            cat_count = run_command(cmd, shell=True, print_output=False)
            if cat_count:
                cat_count_val = cat_count[0]
                if is_int(cat_count_val):
                    return int(cat_count_val)/4

    return header_count


def run_command(command, shell=False, print_output=False, env_exports={}, wait=True, logger=None):
    print_ = logger.info if logger else print
    current_env = os.environ.copy()
    merged_env = {**current_env, **env_exports}
    cmd = shlex.split(command) if type(command) == str else shlex.split(" ".join(command))
    if logger:
        logger.debug(f'Executing: {shlex.join(cmd)}')
    process = subprocess.Popen(
        cmd,
        shell=shell, env=merged_env,
        stdout=subprocess.PIPE if wait else None, stderr=subprocess.STDOUT if wait else None
    )
    if wait:
        stdout = []
        stdout_data, stderr_data = process.communicate()
        for line in stdout_data.splitlines():
            line = line.rstrip().decode('utf8')
            if print_output:
                print_(f'shell> {line}')
            stdout.append(line)
        if process.returncode != 0:
            print_ = logger.error if logger else print
            stderr = []
            stderr_data = "" if not stderr_data else stderr_data
            for line in stderr_data.splitlines():
                line = line.rstrip().decode('utf8')
                stderr.append(line)
            error_msg = stderr if stderr else stderr + stdout
            for msg in error_msg:
                print_(msg)
                if logger and logger.level > 10:
                    break
        return stdout
    else:
        return None


def which(program):
    def is_exe(file_path):
        return os.path.isfile(file_path) and os.access(file_path, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def runningInDocker():
    cgroup = '/proc/self/cgroup'
    if os.path.exists(cgroup):
        with open(cgroup, 'r') as procfile:
            for line in procfile:
                fields = line.strip().split('/')
                if fields[1] == 'docker':
                    return True
    return False


def create_markdown(src_files=['README.md'], other_files=['files/css/stylesheet.css'], out_file='README.pdf', src='markdown', to='pdf', api_url='http://c.docverter.com/convert'):
    if not isinstance(src_files, (list, tuple)):
        src_files = [src_files]
    cmd = [
        'curl', '--form', f'from={src}', '--form',
        f'to={to}', '--form', 'css=stylesheet.css'
    ]
    cmd += [f'--form input_files[]=@{f}' for f in src_files]
    cmd += [f'--form other_files[]=@{f}' for f in other_files]
    if not out_file:
        out_file = f'out.{to}'
    cmd += ['--output', out_file]
    cmd += [api_url]
    run_command(cmd)
    return out_file


def parse_fasta(fasta_in):
    records = {}
    with open(fasta_in, 'r') as fin:
        header, sequence = None,  ''
        for line in fin:
            line = line.strip()
            if line:
                if line[0] == '>':
                    if sequence:
                        records[header] = sequence
                        sequence = ''
                    header = line[1:]
                else:
                    sequence += line

            if sequence:
                next_line = next(fin, None)
                if next_line is None:
                    records[header] = sequence
                else:
                    next_line = next_line.strip()
                    if next_line:
                        if next_line[0] == '>':
                            records[header] = sequence
                            sequence = ''
                            header = next_line[1:]
                        else:
                            sequence += next_line
    return records


def find_in_dir_pattern(root_dir, pattern):
    p = re.compile(pattern)
    found_files = []
    for root, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename[0] == '.':
                continue
            if re.search(p, filename):
                found_files.append(os.path.join(root, filename))
    return found_files


def find_in_dir(root_dir, extensions=None, correct=False):
    found_files = []
    for root, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename[0] == '.':
                continue
            file = os.path.join(root, filename)
            valid, filename_corrected = check_filetype(filename, extensions, correct)
            if valid:
                if filename_corrected:
                    file_old = file
                    file = os.path.join(root, filename_corrected)
                    os.rename(file_old, file)
                found_files.append(file)
    return found_files


def check_filetype(filename, extensions=None, correct=False):
    valid = False
    filename_corrected = None
    if not extensions:
        valid = True
    else:
        for ext in extensions:
            if filename.endswith(f'.{ext}'):
                valid = True
                return valid, filename_corrected
        if correct:
            for ext in extensions:
                if ext in filename:
                    filename_corrected = f'{filename}.corrected.{ext}'
                    valid = True
                    return valid, filename_corrected
    return valid, filename_corrected


def check_valid_input(file_name):
    for file_type in VALID_FILETYPES:
        if file_name.endswith(file_type):
            return True
    for archive_type in VALID_ARCHIVETYPES:
        if file_name.endswith(archive_type):
            return True
    return False


def is_download(inputfiles):
    if inputfiles:
        for inputfile in inputfiles:
            if inputfile.startswith(VALID_URLS) or inputfile.startswith("sra:"):
                return True
    return False


def download_url(url, file_name=None, out_folder=None):
    import ssl
    import http
    from urllib.request import urlopen
    from urllib.error import HTTPError, URLError

    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE

    if not file_name:
        file_name = url.rsplit('/', 1)[-1]
        if not check_valid_input(file_name):
            return None

    if out_folder:
        if not os.path.exists(out_folder):
            os.makedirs(out_folder, exist_ok=True)
        file_name_path = os.path.join(out_folder, file_name)
    else:
        file_name_path = file_name
    try:
        connection_error = False
        print(f'[ urllib ] Retrieving {file_name}...', end='', flush=True)
        with urlopen(url, context=ctx) as u, open(file_name_path, 'wb') as f:
            f.write(u.read())
        print(' done')
    except http.client.IncompleteRead as e:
        print(' canceled')
        print(f"ERROR: Potential incomplete read ({len(e.partial)} bytes read) from url: {url}")
        connection_error = True
    except HTTPError as e:
        print(' canceled')
        if hasattr(e, 'code'):
            print(f'ERROR: The server could not fulfill the request; Error code: {e.code}')
        connection_error = True
    except URLError as e:
        print(' canceled')
        if hasattr(e, 'reason'):
            print(f'ERROR: We failed to reach the server; Reason: {e.reason}')
        connection_error = True
    
    if connection_error:
        return None

    return file_name_path


def download_inputfiles(inputfiles, download_folder=None, env_exports={}):
    inputfiles_downloaded = []
    has_download_folder = True if download_folder else False
    if not has_download_folder:
        download_folder = os.getcwd()
    if inputfiles and len(inputfiles) == 1:
        inputfile = inputfiles[0]
        if inputfile.startswith(VALID_URLS):
            url = inputfile
            download_file_name = None
            url_in = inputfile.split(",")
            if len(url_in) == 2:
                url = url_in[0]
                download_file_name = url_in[1]
            downloaded_file = download_url(url, file_name=download_file_name, out_folder=download_folder)
            if not downloaded_file:
                sys.exit(f"ERROR - Could not download file from url: {inputfile}")
            else:
                inputfiles_downloaded = [os.path.abspath(downloaded_file)]
        elif inputfile.startswith("sra:"):
            sra_acc = inputfile[4:]
            sra_downloaded = sra_downloader(sra_acc, os.path.abspath(download_folder), env_exports)
            if len(sra_downloaded) == 1:
                inputfiles_downloaded = [os.path.abspath(sra_downloaded[0])]
            else:
                print(f"Got SRA with 2-lane format: {sra_acc}")
                for sra_ in sra_downloaded:
                    inputfiles_downloaded+=[os.path.abspath(sra_)]
        for inputfile_ in inputfiles_downloaded:
            if not os.path.exists(inputfile_):
                found = False
                if has_download_folder:
                    inputfile_at_download_folder = os.path.join(download_folder, os.path.basename(inputfile_))
                    if os.path.exists(inputfile_at_download_folder):
                        inputfile_ = inputfile_at_download_folder
                        found = True
                if not found:
                    sys.exit(f"ERROR - Could not read input file: {inputfile_}")
    elif inputfiles:
        if len(inputfiles) != 2:
            sys.exit(f"ERROR - Could not identify 2-lane input files: {inputfiles}")
        lanes_checked = []
        for lane_file in inputfiles:
            if lane_file.startswith(VALID_URLS):
                url = lane_file
                download_file_name = None
                url_in = lane_file.split(",")
                if len(url_in) == 2:
                    url = url_in[0]
                    download_file_name = url_in[1]
                downloaded_file = download_url(url, file_name=download_file_name, out_folder=download_folder)
                if not downloaded_file:
                    sys.exit(f"ERROR - Could not download file from url: {lane_file}")
                else:
                    lanes_checked.append(downloaded_file)
            elif lane_file.startswith("sra:"):
                sra_acc = lane_file[4:]
                sra_downloaded = sra_downloader(sra_acc, os.path.abspath(download_folder))

                if len(sra_downloaded) == 1:
                    lanes_checked.extend(sra_downloaded)
                else:
                    sys.exit(f"ERROR - Got SRA with 2-lane format ({sra_acc}) while processing 2-lane input: {inputfiles}")
            else:
                lanes_checked.append(lane_file)
        r1_ = os.path.abspath(lanes_checked[0])
        r2_ = os.path.abspath(lanes_checked[1])
        if not os.path.exists(r1_) or not os.path.exists(r2_):
            found_r1, found_r2 = (False, False)
            if has_download_folder:
                r1_at_download_folder = os.path.join(download_folder, os.path.basename(r1_))
                if os.path.exists(r1_at_download_folder):
                    r1_ = r1_at_download_folder
                    found_r1 = True
                r2_at_download_folder = os.path.join(download_folder, os.path.basename(r2_))
                if os.path.exists(r2_at_download_folder):
                    r2_ = r2_at_download_folder
                    found_r2 = True
            if not (found_r1 and found_r2):
                sys.exit(f"ERROR - Could not read all input file: {r1_},{r2_}")
        inputfiles_downloaded = [r1_, r2_]

    return inputfiles_downloaded


def pre_execution(method):
    @wraps(method)
    def _impl(self, *method_args, **method_kwargs):
        pass
        return method(self, *method_args, **method_kwargs)
    return _impl


def post_execution(method):
    @wraps(method)
    def _impl(self, *method_args, **method_kwargs):
        return_value = method(self, *method_args, **method_kwargs)
        pass
        return return_value
    return _impl


def retry(exceptions, tries=4, delay=3, backoff=2, logger=None, verbose=False, catchAll=False):
    """Retry calling the decorated function using an exponential backoff.

    adapted from:  http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param exceptions: the exception(s) to check. may be a tuple of
        exceptions to check.
    :type exceptions: Exception type, exception instance, or tuple containing
        any number of both (eg. IOError, IOError(errno.ECOMM), (IOError,), or
        (ValueError, IOError(errno.ECOMM))
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param silent: If set then no logging will be attempted
    :type silent: bool
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    :param verbose: if set, print verbose manages
    :type verbose: bool
    :param catchAll: if set, also catch last retry
    :type catchAll: bool
    """
    try:
        len(exceptions)
    except TypeError:
        exceptions = (exceptions,)
    all_exception_types = tuple(set(x if type(x) == type else x.__class__ for x in exceptions))
    exception_types = tuple(x for x in exceptions if type(x) == type)
    exception_instances = tuple(x for x in exceptions if type(x) != type)

    def deco_retry(f):
        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > (1 if not catchAll else 0):
                try:
                    return f(*args, **kwargs)
                except all_exception_types as e:
                    if (not any(x for x in exception_types if isinstance(e, x))
                        and not any(x for x in exception_instances if type(x) == type(e) and x.args == e.args)):
                        raise
                    msg = f'{str(e) if str(e) != "" else repr(e)}, Retrying in {mdelay} seconds...'
                    if verbose:
                        if logger:
                            logger.warning(msg)
                        else:
                            print(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            if mtries == 1:
                return f(*args, **kwargs)
        return f_retry  # true decorator
    return deco_retry
