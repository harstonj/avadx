#!/usr/bin/env python3

import argparse
import re
import os
import sys
import csv
import math
import uuid
import shutil
import random
import threading
import configparser
import multiprocessing
from pathlib import Path
from datetime import datetime
from timeit import default_timer as timer
from concurrent.futures import ThreadPoolExecutor
from .logger import Logger
from .helper import check_config_ini, run_command, flatten, runningInDocker, runningInSingularity
from . import __version__, __releasedate__, name


QUIET = False
VM_MOUNT = Path('/mnt')
SINGULARITY_LIB = Path(os.getenv('SINGULARITY_LIB', str(Path.cwd())))
LOG_LEVEL = 'INFO'
LOG_FILE = None
USE_PROGRESS = True


try:
    import progressbar
    progressbar.streams.wrap_stderr()
except ImportError:
    USE_PROGRESS = False


class AVADx:

    IMAGES = {
        'avadx': 'bromberglab/avadx',
        'R': 'bromberglab/avadx-rscript',
        'bcftools': 'bromberglab/avadx-bcftools',
        'tabix': 'bromberglab/avadx-tabix',
        'annovar': 'bromberglab/avadx-annovar',
        'varidb': 'bromberglab/varidb'
    }

    def __init__(self, pipeline):
        self.log = self.get_logger()
        self.pipeline = pipeline
        self.threadLock = threading.Lock()

    @staticmethod
    def init_vm(daemon):
        _ = AVADx(None)
        if daemon not in ['docker', 'singularity']:
            _.log.warning(f'Unknown daemon: {daemon}')
        if daemon == 'docker':
            try:
                import docker
            except ImportError:
                _.log.error('Docker SDK for Python not installed. Please install with: "pip install docker"')
                sys.exit(1)
            _.log.info('Checking Docker images ...')
            client = docker.from_env()
            for image in AVADx.IMAGES.values():
                _.log.info(f'Processing {image} ...')
                client.images.pull(image)
            _.log.info('Done.')
        elif daemon == 'singularity':
            for image in AVADx.IMAGES.values():
                _.log.info(f'Processing {image} ...')
                run_command(f'singularity pull docker://{image}')
            _.log.info('Done.')

    def get_logger(self):
        logger = Logger(self.__class__.__name__, level=LOG_LEVEL)
        logger.addConsoleHandler()
        log = logger.getLogger()
        return log

    def run(self, func_name, uid, kwargs, extra):
        func = getattr(self, func_name, None)
        if func:
            return func(uid, **kwargs)
        else:
            self.log.warning(f'No method: {func_name}')

    def run_method(self, container, name, uid, kwargs):
        args = kwargs.get(name).pop(0).split()
        daemon_args = kwargs.get(name + "_darg").pop(0)
        taskflag, taskfile, taskprefix = kwargs.get(name + "_tasks").pop(0)
        if taskfile:
            with (taskfile).open() as fin:
                tasklist = [_.strip() for _ in fin.readlines()]
        else:
            tasklist = [None]
        outdir = kwargs.get(name + "_outd").pop(0)
        mounts = kwargs.get(name + "_mounts").pop(0)
        progress = kwargs.get(name + "_progress").pop(0)
        fns = kwargs.get(name + "_fns").pop(0)
        fns_pre, fns_post = fns.get('pre', (None, [])), fns.get('post', (None, []))
        description = kwargs.get(name + "_desc").pop(0)
        level = kwargs.get(name + "_lvl").pop(0)
        wd_results, out_results = kwargs.get(name + "_results").pop(0)
        reports = kwargs.get(name + "_reports").pop(0)
        log_stdout, log_stderr = kwargs.get(name + "_logs").pop(0)
        resources = kwargs.get(name + "_res").pop(0)
        cpu = resources.get('cpu', 1)
        # mem = resources.get('mem', 0)
        if outdir:
            outdir.mkdir(parents=True, exist_ok=True)
        task_cnt = len(tasklist)
        if task_cnt == 1:
            self.run_task(
                1, tasklist[0], task_cnt, level, args, name, description,
                fns_pre, fns_post, taskprefix, taskflag, log_stdout, log_stderr,
                container, daemon_args, uid, mounts, kwargs, reports, wd_results, out_results, None, progress
            )
        else:
            workers_cnt = int(self.pipeline.get_host_cpu() // cpu)
            self.log.info(f'|{level:.2f}| {name}: Spawning {workers_cnt} threads [reserving {cpu} cores each] ...')
            bar_prefix = f'[     INFO ] --- |{level:.2f}| {name}: '
            bar = progressbar.ProgressBar(max_value=task_cnt, prefix=bar_prefix).start() if USE_PROGRESS and task_cnt > 1 else None
            with ThreadPoolExecutor(max_workers=workers_cnt) as executor:
                args = [
                    (
                        tid, task, task_cnt, level, args, name, description,
                        fns_pre, fns_post, taskprefix, taskflag, log_stdout, log_stderr,
                        container, daemon_args, uid, mounts, kwargs, reports, wd_results, out_results, bar, progress
                    ) for tid, task in enumerate(tasklist, 1)
                ]
                bar.update(1) if bar else None
                executor.map(lambda a: self.run_task(*a), args)
            bar.finish() if bar else None

    def run_task(
        self, tid, task, task_cnt, level, args, name, description,
        fns_pre, fns_post, taskprefix, taskflag, log_stdout, log_stderr,
        container, daemon_args, uid, mounts, kwargs, reports, wd_results, out_results, bar, progress
    ):
        args_ = list(args)
        task_info = f' [{tid}/{task_cnt}] ({task})' if task and task_cnt > 1 else ''
        if not bar:
            self.log.info(f'|{level:.2f}| {name}{task_info}: {description}')
        self.log.debug(f'|{level:.2f}| {name}{task_info}: {description} - started {datetime.now()}')
        timer_start = timer()
        if fns_pre[0]:
            fns_pre[0](*fns_pre[1])
        if task:
            task_idxs = [i for i in range(len(args_)) if args_[i].find('$TASK') != -1]
            for task_idx in task_idxs:
                update = f'{taskprefix if taskprefix else ""}{task}'
                args_[task_idx] = update if args_[task_idx] == '$TASK' else args_[task_idx].replace('$TASK', update)
            if taskflag:
                task_idxs.reverse()
                for task_idx_r in task_idxs:
                    args_.insert(task_idx_r, taskflag)
            if log_stdout:
                log_stdout = f'{task}_{log_stdout}'
            if log_stderr:
                log_stderr = f'{task}_{log_stderr}'
        if '-D' in sys.argv or '--Debug' in sys.argv:
            codebase = self.pipeline.kwargs.get('Debug')
            codebase = (codebase if codebase is not None else Path(__file__).parent.parent).absolute()
            mounts += [
                (codebase / 'avadx', Path('/usr/local/lib/python3.8/site-packages/avadx')),
                (codebase / 'python', Path('/app/python/avadx/')),
                (codebase / 'R', Path('/app/R/avadx/')),
            ]
        self.pipeline.run_container(
            container, args=args_, daemon_args=daemon_args, uid=uid,
            mounts=mounts,
            out_folder=kwargs.get('wd'),
            stdout=log_stdout,
            stderr=log_stderr,
            progress=progress
        )
        if fns_post[0]:
            fns_post[0](*fns_post[1])
        self.reports(kwargs.get('wd'), uid, reports)
        self.log.debug(f'|{level:.2f}| {name}{task_info}: took {(timer() - timer_start):.3f} seconds')
        if bar:
            with self.threadLock:
                bar.update(bar.data()['value'] + 1, force=True)
        self.check_results(kwargs.get('wd'), uid, wd_results, out_results)

    def reports(self, wd, uid, reports):
        for report in reports:
            source = wd / str(uid) / 'wd' / Path(report[0])
            target = wd / str(uid) / 'out' / 'reports' / Path(report[1])
            if source.exists():
                if len(report) == 3 and report[2]:
                    if source.is_file():
                        shutil.copyfile(source, target)
                    else:
                        if target.exists():
                            shutil.rmtree(target)
                        shutil.copytree(source, target)
                else:
                    if target.exists() and target.is_dir():
                        shutil.rmtree(target)
                    shutil.move(source, target)

    def check_results(self, wd, uid, wd_results, out_results):

        def validate(path):
            if not path.exists():
                self.log.error(f'Missing result file: {path} - Aborting.')
                sys.exit(1)

        for f in wd_results:
            validate(wd / str(uid) / 'wd' / f)
        for f in out_results:
            validate(wd / str(uid) / 'out' / f)

    def run_info(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'run_info', uid, kwargs)

    def run_preprocess(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'run_preprocess', uid, kwargs)

    def run_retrieve(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'run_retrieve', uid, kwargs)

    def generate_transcripts_protlength(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'generate_transcripts_protlength', uid, kwargs)

    def bcftools(self, uid, **kwargs):
        self.run_method(self.IMAGES['bcftools'], 'bcftools', uid, kwargs)

    def avadx(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'avadx', uid, kwargs)

    def filterVCF_by_ABAD(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'filterVCF_by_ABAD', uid, kwargs)

    def gnomad_ALLabove0_preprocess(self, uid, **kwargs):
        self.run_method(self.IMAGES['tabix'], 'gnomad_ALLabove0_preprocess', uid, kwargs)

    def generate_gnomad_above0(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'generate_gnomad_above0', uid, kwargs)

    def gnomad_ALLabove0_postprocess(self, uid, **kwargs):
        self.run_method(self.IMAGES['tabix'], 'gnomad_ALLabove0_postprocess', uid, kwargs)

    def filterVCF_by_gnomAD(self, uid, **kwargs):
        self.run_method(self.IMAGES['tabix'], 'filterVCF_by_gnomAD', uid, kwargs)

    def stats_quality_pca(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'stats_quality_pca', uid, kwargs)

    def ethnicity_EthSEQ(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'ethnicity_EthSEQ', uid, kwargs)

    def ethnicity_EthSEQ_summary(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'ethnicity_EthSEQ_summary', uid, kwargs)

    def relatedness(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'relatedness', uid, kwargs)

    def annovar(self, uid, **kwargs):
        self.run_method(self.IMAGES['annovar'], 'annovar', uid, kwargs)

    def generate_SNAP_query(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'generate_SNAP_query', uid, kwargs)

    def varidb_prefilter(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'varidb_prefilter', uid, kwargs)

    def varidb(self, uid, **kwargs):
        self.run_method(self.IMAGES['varidb'], 'varidb', uid, kwargs)

    def cal_genescore_make_genescore(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'cal_genescore_make_genescore', uid, kwargs)

    def merge_genescore(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'merge_genescore', uid, kwargs)

    def ava_model(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'ava_model', uid, kwargs)

    def FS_CVgeneOverRep_kfold(self, uid, **kwargs):
        self.run_method(self.IMAGES['R'], 'FS_CVgeneOverRep_kfold', uid, kwargs)

    def run_postprocess(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'run_postprocess', uid, kwargs)

    def run_visualization(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'run_visualization', uid, kwargs)

    def run_prediction(self, uid, **kwargs):
        self.run_method(self.IMAGES['avadx'], 'run_prediction', uid, kwargs)


class Pipeline:

    REGISTRY = {
        'singularity': {
            'CMD_RUN': 'run',
            'CMD_EXEC': 'exec'
        },
        'docker': {
            'CMD_RUN': 'run',
            'CMD_EXEC': 'exec'
        }
    }

    ASSEMBLY_DEFAULT = 'hg19'

    ASSEMBLY_MAPPING = {
        'hg19': ('GRCh37', 'p13', 'GCF_000001405', '25'),
        'hg38': ('GRCh38', 'p13', 'GCF_000001405', '39')
    }

    def __init__(self, actions=[], kwargs={}, uid=None, config_file=None, daemon='docker'):
        self.log = self.get_logger()
        self.actions = actions
        self.kwargs = kwargs
        self.config_file = config_file
        self.daemon = daemon
        self.uid = self.get_uid(uid)
        self.config = self.load_config()
        self.is_vm = self.check_vm()
        self.entrypoint = kwargs.get('entrypoint') if kwargs.get('entrypoint', None) is not None else 1
        self.exitpoint = kwargs.get('exitpoint', None)
        self.skip = kwargs.get('skip', None)
        self.resources = {
            'cpu': None,
            'mem': None,
            'vm.cpu': None,
            'vm.mem': None
        }
        self.seed = self.config.get('avadx', 'random.seed', fallback=None)
        self.set_seed()
        self.prediction = False
        self.prediction_out = None

    def reset(self, actions=[], kwargs=None, uid=None, config_file=None, prediction=False):
        self.actions = actions
        self.kwargs = kwargs if kwargs is not None else self.kwargs
        self.config_file = config_file if config_file is not None else self.config_file
        self.uid = self.get_uid(uid) if uid is not None else self.uid
        self.config = self.load_config() if config_file is not None else self.config
        self.entrypoint = self.kwargs.get('entrypoint') if self.kwargs.get('entrypoint', None) is not None else 1
        self.exitpoint = self.kwargs.get('exitpoint', None)
        self.skip = self.kwargs.get('skip', None)
        self.prediction = prediction

    def get_logger(self):
        logger = Logger(self.__class__.__name__, level=LOG_LEVEL)
        logger.addConsoleHandler()
        log = logger.getLogger()
        return log

    def set_seed(self):
        random.seed(self.seed)

    def get_wd(self):
        return self.kwargs['wd'] / str(self.uid)

    def get_uid(self, uid):
        if uid is None:
            uid = self.config_file.stem
            if Path(self.kwargs['wd'] / uid).exists():
                uid = f'{uid}_{uuid.uuid1()}'
        return uid

    def is_prediction(self):
        return True if self.prediction or self.config.get('DEFAULT', 'predict', fallback='no') == 'yes' else False

    def get_prediction_out(self):
        return self.prediction_out

    def check_vm(self):
        if runningInDocker():
            return 'docker'
        elif runningInSingularity():
            return 'singularity'
        else:
            return None

    def get_vm_resources(self, dry_run=False):
        if dry_run:
            cpu, mem = 1, 1
        else:
            out = self.run_container(AVADx.IMAGES['avadx'], args=['--info'], mkdirs=False)
            daemon, cpu, mem = [_.split()[1] for _ in out[1].replace('[', '').replace(']', '').split(', ')] if out is not None else [None, 1, 0]
        self.resources['vm.cpu'] = int(cpu)
        self.resources['vm.mem'] = int(mem)

    def set_host_cpu(self):
        self.resources['cpu'] = multiprocessing.cpu_count()

    def get_host_cpu(self):
        return self.resources['cpu']

    def get_vm_cpu(self):
        return self.resources['vm.cpu']

    def get_vm_mem(self):
        return math.ceil(self.resources['vm.mem'] / 1024**3)

    def set_vm_cpu(self, cores):
        self.resources['vm.cpu'] = cores

    def set_vm_mem(self, gb):
        self.resources['vm.mem'] = gb * (1024**3)

    def get_splits(self, threads=1):
        return 100 * (((self.get_vm_mem() / threads) - 2) // 4)

    def load_config(self):
        config = configparser.ConfigParser()
        if self.config_file and self.config_file.exists():
            config.read(str(self.config_file))
        else:
            self.config_file = self.kwargs.get('wd') / str(self.uid) / 'out' / 'pipeline_config.ini'
            if self.config_file.exists():
                self.log.info('Using default pipeline config file')
                config.read(str(self.config_file))
            else:
                self.config_file = VM_MOUNT / 'in' / self.config_file.name
                if self.config_file and self.config_file.exists():
                    config.read(str(self.config_file))
        if self.config_file is None or not self.config_file.exists():
            if self.kwargs.get('info', False):
                pass
            else:
                self.log.error('No config file found. Aborting.')
                sys.exit(1)
        data_base_path = config.get('DEFAULT', 'datadir', fallback=None)
        if data_base_path is not None:
            if not Path(data_base_path).is_absolute():
                config.set('DEFAULT', 'datadir', str(self.config_file.parent / data_base_path))
        configoverride = {_.split('=')[0]: _.split('=')[1] for _ in self.kwargs['configoverride']} if self.kwargs['configoverride'] is not None else {}
        for k, v in configoverride.items():
            if config.get('avadx', k, fallback=None) is not None:
                config.set('avadx', k, v)
            else:
                self.log.warning(f'No entry in [avadx] config for "{k}" - no override possible.')
        return config

    def save_run_config(self):
        if not self.prediction:
            with (self.get_wd() / 'out' / 'pipeline_config.ini').open('w') as configfile:
                self.config.write(configfile)

    def save_run_info(self):
        if not self.prediction:
            cpu = self.resources['cpu' if self.is_vm else 'vm.cpu']
            mem = self.resources['mem' if self.is_vm else 'vm.mem']
            with (self.get_wd() / 'out' / 'pipeline_info.txt').open('w') as infofile:
                infofile.write(f' --- AVA,Dx {__version__} {__releasedate__} ---\n\n [Environment]\n + VM.daemon  : {self.daemon:>13}\n + VM.cpu     : {cpu:>13}\n + VM.memory  : {mem:>13}\n + Random.seed: {self.seed:>13}\n\n [Started]\n + {datetime.now()}\n')

    def info(self, quiet=False, run_args=None):
        if self.is_vm:
            try:
                import psutil
                self.resources['cpu'] = psutil.cpu_count()
                self.resources['mem'] = psutil.virtual_memory().total
            except ImportError:
                self.log.warning('Info requires the "psutil" package: pip install psutil')
        else:
            self.get_vm_resources()
        cpu = self.resources['cpu' if self.is_vm else 'vm.cpu']
        mem = self.resources['mem' if self.is_vm else 'vm.mem']
        if not self.is_vm and run_args:
            pipeline = run_all(*run_args)
            for idx, action in enumerate(pipeline.actions):
                level, description = pipeline.kwargs[f'{action}_lvl'].pop(0), pipeline.kwargs[f'{action}_desc'].pop(0)
                print(f'[ {level:<4} ] - {action:<30} * {description}')
            print()
        elif not quiet:
            print(f'AVA,Dx {__version__} {__releasedate__}\n[VM.daemon: {self.daemon}, VM.cpu: {cpu}, VM.memory: {mem}]')
        return ('AVA,Dx', __version__, __releasedate__, self.daemon, cpu, mem)

    def preprocess(self, rerun=False):
        def indices(labels, search):
            return [i for i, x in enumerate(labels) if x == search]

        def folds(samples, foldcnt):
            return [samples[i * foldcnt:(i + 1) * foldcnt] for i in range((len(samples) + foldcnt - 1) // foldcnt)]

        def chrname2single(X=True, Y=True, M=True):
            return [_ for _ in range(1, 23)] + (['X'] if X else []) + (['Y'] if Y else []) + (['M'] if M else [])

        if self.is_vm:
            wd_folder = self.kwargs.get('wd')
            out_folder = wd_folder.parent / 'out'
            samples_path = VM_MOUNT / 'in' / Path(self.config.get('avadx', 'samples')).name
            samplesids_path = wd_folder / 'tmp' / 'sampleids.txt'
            cvscheme_path = wd_folder / 'tmp' / 'cv-scheme.csv'
            chr2num_path = wd_folder / 'tmp' / 'chr_to_number.txt'
            samples_list = wd_folder / 'tmp' / 'all_samples.txt'
        else:
            wd_folder = self.kwargs.get('wd') / str(self.uid) / 'wd'
            out_folder = self.kwargs.get('wd') / str(self.uid) / 'out'
            wd_folder.mkdir(parents=True, exist_ok=True)
            samples = self.check_config('samples', is_file=True, quiet=False)
            if not samples:
                self.log.warning('Required input missing - sample annotation file')
                return
            samples_path = Path(self.config.get('avadx', 'samples'))
            if not samples_path.is_absolute():
                samples_path = self.config_file.parent / samples_path
            samplesids_path = wd_folder / 'tmp' / 'sampleids.txt'
            cvscheme_path = wd_folder / 'tmp' / 'cv-scheme.csv'
            chr2num_path = wd_folder / 'tmp' / 'chr_to_number.txt'
            samples_list = wd_folder / 'tmp' / 'all_samples.txt'
        if samples_path.is_dir() or not samples_path.exists():
            samples_path_new = out_folder / 'samples.csv'
            if not samples_path_new.exists():
                if not self.is_prediction():
                    self.log.warning(
                        f'No samples file found: {samples_path}. Creating file {samples_path_new.name} in output directory using ALL samples and unknown class labels.\n'
                        '-> IMPORTANT: You have to set class labels manually for model training! <-'
                    )
                with samples_list.open('r') as fin, samples_path_new.open('w') as fout:
                    fout.write(
                        '# auto generated samples class labels file\n'
                        '# replace "?" by the respective class label for any given sample\n'
                    )
                    fout.write('sampleid,class\n')
                    for line in fin:
                        fout.write(f'{line.strip()},?\n')
            samples_path = samples_path_new
        with samples_path.open() as fin, samplesids_path.open('w') as fout_ids, cvscheme_path.open('w') as fout_cv, chr2num_path.open('w') as fout_chr:
            samples, labels, col3 = [], [], []
            reader = csv.reader(filter(lambda row: row[0] != '#', fin))
            header = next(reader)
            header_ncol = len(header)
            has_sampleid = 'sampleid' in header
            has_class = 'class' in header
            has_groups = len(header) == 3 and 'group' in header
            has_folds = len(header) == 3 and 'fold' in header
            for column, check in {'sampleid': has_sampleid, 'class': has_class}.items():
                if not check:
                    self.log.error(f'|1.02| Required column "{column}" not found in sample file header: {header}. Please check valid content/format of: {samples_path.name}')
                    return
            s_idx, c_idx = header.index('sampleid'), header.index('class')
            f_idx = header.index('fold') if has_folds else None
            g_idx = header.index('group') if has_groups else None
            for row in reader:
                if len(row) != header_ncol:
                    self.log.warning(f'|1.02| Format error in sample file, skipping: {",".join(row)}')
                    continue
                samples += [row[s_idx]]
                labels += [row[c_idx]]
                if has_folds:
                    col3 += [row[f_idx]]
                elif has_groups:
                    col3 += [row[g_idx]]
            if not samples:
                self.log.error(f'|1.02| No samples recognized, please check valid content/format of: {samples_path.name}')
                return
            fout_ids.writelines([f'{_}\n' for _ in samples])
            cv_folds_config = int(self.config.get('avadx', 'cv.folds'))
            if has_groups:
                auto_folds = []
                groups = col3
                grouped_samples = {g: indices(groups, g) for g in set(groups)}
                cv_folds_max = len(grouped_samples)
                cv_folds = min(cv_folds_config, cv_folds_max) if cv_folds_config > 1 else cv_folds_max
                if cv_folds_config > cv_folds_max:
                    self.log.warning(
                        f'Too few groups ({cv_folds_max}) for specified {cv_folds_config}-fold cross-validation'
                        f' - proceeding with {cv_folds}-fold split')
                grouped_samples_shuffled = list(grouped_samples)
                random.shuffle(grouped_samples_shuffled)
                group_folds = folds(grouped_samples_shuffled, len(grouped_samples_shuffled) // cv_folds)
                for gfold in group_folds:
                    auto_folds += [flatten([grouped_samples.get(g) for g in gfold])]
                auto_folds = adjust_splits(auto_folds, cv_folds)
                if len(auto_folds) != cv_folds:
                    self.log.warning(f'Auto-generated cross-validation folds ({len(auto_folds)}) differ from specified folds ({cv_folds})')
                grouped_rows = {}
                for idx, fold in enumerate(auto_folds, 1):
                    for sidx in fold:
                        grouped_rows[samples[sidx]] = (idx, labels[sidx])
                fout_cv.writelines([f'{sampleid},{grouped_rows[sampleid][0]},{grouped_rows[sampleid][1]}\n' for sampleid in samples])
                split_type = 'auto-generated group'
                split_description = f'{len(auto_folds)}-fold split' if cv_folds < cv_folds_max else 'leave-one-out splits'
            elif has_folds:
                folds = col3
                for sidx, fold in enumerate(folds):
                    fout_cv.write(f'{samples[sidx]},{fold},{labels[sidx]}\n')
                cv_folds = len(set(folds))
                split_type = 'user-specified fold'
                split_description = f'{cv_folds}-fold split' if cv_folds < len(folds) else 'leave-one-out splits'
            else:
                cv_folds_max = len(samples)
                cv_folds = min(cv_folds_config, cv_folds_max) if cv_folds_config > 1 else cv_folds_max
                if cv_folds_config > cv_folds_max:
                    if not self.is_prediction():
                        self.log.warning(
                            f'Too few samples ({cv_folds_max}) for specified {cv_folds_config}-fold cross-validation'
                            f' - proceeding with {cv_folds}-fold split')
                samples_indices = [_ for _ in range(cv_folds_max)]
                random.shuffle(samples_indices)
                auto_folds = folds(samples_indices, len(samples_indices) // cv_folds)
                auto_folds = adjust_splits(auto_folds, cv_folds)
                if len(auto_folds) != cv_folds:
                    self.log.warning(f'Auto-generated cross-validation folds ({len(auto_folds)}) differ from specified folds ({cv_folds})')
                auto_folds_rows = {}
                for idx, fold in enumerate(auto_folds, 1):
                    for sidx in fold:
                        auto_folds_rows[samples[sidx]] = (idx, labels[sidx])
                fout_cv.writelines([f'{sampleid},{auto_folds_rows[sampleid][0]},{auto_folds_rows[sampleid][1]}\n' for sampleid in samples])
                split_type = 'auto-generated sample'
                split_description = f'{len(auto_folds)}-fold split' if cv_folds < cv_folds_max else 'leave-one-out splits'
            if not rerun and not self.is_prediction():
                self.log.info(f'|1.02| Using {split_type} based cross-validation scheme for {split_description}')
            chr_rename = chrname2single()
            fout_chr.writelines([f'chr{_} {_}\n' for _ in chr_rename])

    def retrieve(self, target, outfolder=None, outfile=None):
        import requests
        if self.is_vm:
            config_datadir_orig = self.config.get('DEFAULT', 'datadir', fallback=None)
            wd_folder = outfolder if outfolder else self.kwargs.get('wd')
            self.config.set('DEFAULT', 'datadir', str(VM_MOUNT / 'data'))
        else:
            wd_folder = outfolder if outfolder else self.kwargs.get('wd') / str(self.uid) / 'wd'
            wd_folder.mkdir(parents=True, exist_ok=True)

        if target == 'cpdb':
            save_as = Path(outfile) if outfile else Path(self.config.get('DEFAULT', 'avadx.data')) / 'CPDB_pathways_genesymbol.tab'
            url = 'http://cpdb.molgen.mpg.de/CPDB/getPathwayGenes?idtype=hgnc-symbol'
            r = requests.get(url, allow_redirects=True)
            open(str(save_as), 'wb').write(r.content)
        elif target == 'varidb':
            avadx_data_path = Path(self.config.get('DEFAULT', 'avadx.data', fallback=self.kwargs.get('wd')))
            save_as = Path(wd_folder) / 'master.zip'
            url = 'https://bitbucket.org/bromberglab/avadx-lfs/get/master.zip'
            r = requests.get(url, allow_redirects=True)
            open(str(save_as), 'wb').write(r.content)
            run_command(['unzip', '-o', '-j', str(save_as), '-d', str(avadx_data_path)])
            run_command(['7z', 'e', str(avadx_data_path / 'varidb.db.7z'), '-aoa', f'-o{avadx_data_path}'])
            md5sum_downloaded = run_command(['md5sum', str(avadx_data_path / 'varidb.db')])
            md5sum_varidb = md5sum_downloaded[0].split()[0] if md5sum_downloaded else ''
            with (avadx_data_path / 'varidb.md5').open() as fin:
                md5_check = fin.readline().strip()
                if md5sum_varidb != md5_check:
                    self.log.warning(f'|0.14| md5 hash not identical - is: {md5sum_varidb} reference: {md5_check}')
            os.remove(str(save_as))
            os.remove(str(avadx_data_path / 'varidb.db.7z'))
        elif target == 'refseq':
            hgref = self.config.get('avadx', 'hgref')
            hgref_mapped = self.ASSEMBLY_MAPPING.get(hgref, None)
            if not hgref_mapped:
                hgref_mapped = self.ASSEMBLY_MAPPING.get(self.ASSEMBLY_DEFAULT)
                self.log.warning(f'|0.17| could not map hgref {hgref} - falling back to {self.ASSEMBLY_DEFAULT}')
            refseq_data_path = Path(self.config.get('DEFAULT', 'refseq.data', fallback=self.kwargs.get('wd')))
            url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/'
            assembly_base = f'{hgref_mapped[2]}.{hgref_mapped[3]}_{hgref_mapped[0]}.{hgref_mapped[1]}'
            url_sequences = f'{url}/{assembly_base}/{assembly_base}_protein.faa.gz'
            sequences_save_as = Path(outfile).with_suffix('.sequences' + Path(outfile).suffix) if outfile else refseq_data_path / f'{hgref_mapped[0]}.{hgref_mapped[1]}_protein.faa.gz'
            url_stats = f'{url}/{assembly_base}/{assembly_base}_feature_table.txt.gz'
            stats_save_as = Path(outfile).with_suffix('.stats' + Path(outfile).suffix) if outfile else refseq_data_path / f'{hgref_mapped[0]}.{hgref_mapped[1]}_feature_table.txt.gz'
            r = requests.get(url_sequences, allow_redirects=True)
            open(str(sequences_save_as), 'wb').write(r.content)
            run_command(['gunzip', '-f', str(sequences_save_as)])
            r = requests.get(url_stats, allow_redirects=True)
            open(str(stats_save_as), 'wb').write(r.content)
            run_command(['gunzip', '-f', str(stats_save_as)])
        else:
            self.log.warning(f'|0.14| Could not retrieve {target} - target unknown')

        if self.is_vm:
            self.config.set('DEFAULT', 'datadir', config_datadir_orig)

    def postprocess(self):

        def cleanup(wd_folder, item):
            base = item[0]
            if item[1] is None:
                path = wd_folder / base
                if path.exists():
                    if path.is_file():
                        os.remove(path)
                    else:
                        shutil.rmtree(path)
            else:
                path = wd_folder / base
                for pattern in item[1]:
                    # TODO implement
                    pass

        if self.is_vm:
            wd_folder = self.kwargs.get('wd')
        else:
            wd_folder = self.kwargs.get('wd') / str(self.uid) / 'wd'

        wd_cleanup = [
            ('tmp', None),
            ('vcf', None),
            ('annovar', ('*.log', '*.avinput'))
        ]

        for item in wd_cleanup:
            cleanup(wd_folder, item)

    def visualization(self, info):
        from .plots import Figure
        vis, dataset = info.split(':')
        if self.is_vm:
            wd_folder = self.kwargs.get('wd')
            out_folder = wd_folder.parent / 'out'
        else:
            wd_folder = self.kwargs.get('wd') / str(self.uid) / 'wd'
            out_folder = self.kwargs.get('wd') / str(self.uid) / 'out'
        genescore_normalize = True if self.config.get('avadx', 'genescore.standardize', fallback='no') == 'yes' else False
        fig = Figure(vis, dataset, genescore_normalize, wd_folder, out_folder)
        fig.create()

    def predict(self, uid, kwargs, extra, config, daemon):
        model_path, features_path = None, None
        # wd_folder = self.kwargs.get('wd') / str(self.uid) / 'wd'
        out_folder = self.kwargs.get('wd') / str(self.uid) / 'out'
        model_pipeline = out_folder / 'results' / 'model.joblib'
        features_pipeline = out_folder / 'results' / 'model_features.txt'
        prediction_folder = out_folder / 'predictions'
        prediction_folder.mkdir(parents=True, exist_ok=True)
        if kwargs['model'] is not None:
            if not kwargs['model'].exists():
                self.log.warn(f'Could not find user model file: {kwargs["model"].absolute()}')
            else:
                self.log.info(f'Loading user model file: {kwargs["model"].absolute()}')
                model_path = kwargs["model"]
                features_path = kwargs["model"].parent / f'{kwargs["model"].stem}_features.txt'
        if model_path is None:
            if model_pipeline.exists():
                self.log.info('Loading pipeline model file...')
                model_path = model_pipeline
                features_path = features_pipeline
        if model_path is None:
            self.log.warning('Please apply the AVA,Dx pipeline first to train a model. Aborting.')
            return

        input_samples = kwargs['sampleprediction']
        if not input_samples.exists():
            self.log.error(f'Could not find sample prediction file: {input_samples.absolute()}. Aborting.')
            return
        elif input_samples.suffix == '.csv':
            self.log.info(f'Loading genescores from file: {input_samples.absolute()}')
            features_samples = prediction_folder / f'{input_samples.stem}_extern.csv'
            shutil.copy(input_samples, features_samples)
        else:
            features_samples = prediction_folder / f'{input_samples.stem}_genescores.csv'
        extract_features = not features_samples.exists()
        if extract_features:
            features_uid = input_samples.stem
            features_wd = self.kwargs.get('wd') / str(self.uid) / 'wd' / 'prediction'
            features_config = prediction_folder / f'{input_samples.stem}_config.ini'
            with features_config.open('w') as configfile:
                configp = configparser.ConfigParser()
                if not config.exists():
                    config = self.kwargs.get('wd') / str(self.uid) / 'out' / 'pipeline_config.ini'
                    if config.exists():
                        self.log.info('Using default pipeline config file.')
                    else:
                        self.log.error('No config file found. Aborting.')
                        sys.exit(1)
                configp.read(str(config))
                configp.set('avadx', 'vcf', str(input_samples.absolute()))
                configp.set('avadx', 'samples', '')
                configp.set('avadx', 'filter.MR.enabled', 'no')
                configp.set('avadx', 'qc.call.MR', '1')
                configp.set('DEFAULT', 'predict', 'yes')
                configp.write(configfile)
            features_kwargs = {k: v for k, v in kwargs.items()}
            features_kwargs['wd'] = features_wd
            exitpoint_main = features_kwargs['exitpoint']
            features_kwargs['skip'] = [5, 7]
            features_kwargs['exitpoint'] = None
            genescore_normalize = True if configp.get('avadx', 'genescore.standardize', fallback='no') == 'yes' else False
            genescores_suffix = 'normalized' if genescore_normalize else 'raw'
            features_generated = features_wd / features_uid / 'out' / 'genescores' / f'GeneScoreTable_{genescores_suffix}.csv'
            extra += [f'pred_id={input_samples.stem}', f'features={features_path if features_path.exists() else "None"}', f'model={model_path}', f'genescores={features_generated}']
            self.reset(kwargs=features_kwargs, uid=features_uid, config_file=features_config)
            self.prediction_out = prediction_folder.absolute()
            run_all_p(self, extra)
            shutil.copy(features_generated, features_samples)
            genescores_prediction = features_wd / features_uid / 'out' / 'predictions' / f'{input_samples.stem}_genescores_prediction.csv'
            shutil.copy(genescores_prediction, prediction_folder)
            predictions_csv = features_wd / features_uid / 'out' / 'predictions' / f'{input_samples.stem}_predictions.csv'
            predictions_png = features_wd / features_uid / 'out' / 'predictions' / f'{input_samples.stem}_predictions.png'
            shutil.copy(predictions_csv, prediction_folder)
            shutil.copy(predictions_png, prediction_folder)
            if exitpoint_main is None and predictions_csv.exists():
                shutil.rmtree(self.get_wd())
        else:
            self.log.info(f'Using existing genescores from: {features_samples}')
            extra += [f'pred_id={input_samples.stem}', f'features={features_path if features_path.exists() else "None"}', f'model={model_path}', f'genescores={features_samples}']
            kwargs_prediction = {k: v for k, v in kwargs.items()}
            kwargs_prediction['entrypoint'], kwargs_prediction['exitpoint'] = 8, None
            self.reset(kwargs=kwargs_prediction, prediction=True)
            self.prediction_out = self.get_wd() / 'out' / 'predictions'
            run_all_p(self, extra)

    def add_action(
            self, action, level=None, description='', args='', daemon_args={},
            tasks=(None, None, None), fns={'pre': (None, []), 'post': (None, [])},
            outdir=None, mounts=[], progress=False, results=([], []), reports=[], logs=(None, None), resources={'cpu': 1, 'mem': 0}):
        if level is None or (level >= self.entrypoint and (self.exitpoint is None or level < self.exitpoint)):
            if self.skip is not None and (level >= self.skip[0] and level <= self.skip[1]):
                return
            self.actions += [action]
            self.kwargs[action] = self.kwargs.get(action, []) + [args]
            self.kwargs[f'{action}_lvl'] = self.kwargs.get(f'{action}_lvl', []) + [level]
            self.kwargs[f'{action}_desc'] = self.kwargs.get(f'{action}_desc', []) + [description]
            self.kwargs[f'{action}_darg'] = self.kwargs.get(f'{action}_darg', []) + [daemon_args]
            self.kwargs[f'{action}_tasks'] = self.kwargs.get(f'{action}_tasks', []) + [tasks]
            self.kwargs[f'{action}_fns'] = self.kwargs.get(f'{action}_fns', []) + [fns]
            self.kwargs[f'{action}_outd'] = self.kwargs.get(f'{action}_outd', []) + [outdir]
            self.kwargs[f'{action}_mounts'] = self.kwargs.get(f'{action}_mounts', []) + [mounts]
            self.kwargs[f'{action}_progress'] = self.kwargs.get(f'{action}_progress', []) + [progress]
            self.kwargs[f'{action}_results'] = self.kwargs.get(f'{action}_results', []) + [results]
            self.kwargs[f'{action}_reports'] = self.kwargs.get(f'{action}_reports', []) + [reports]
            self.kwargs[f'{action}_logs'] = self.kwargs.get(f'{action}_logs', []) + [logs]
            self.kwargs[f'{action}_res'] = self.kwargs.get(f'{action}_res', []) + [resources]

    @staticmethod
    def format_step(step):
        pos1, remaining = f'{step:.2f}'.rsplit(".", 1)[0], list(f'{step:.2f}'.rsplit(".", 1)[-1])
        return f'{pos1}_{remaining[0]}{"-" + remaining[1] if remaining[1] != "0" else ""}'

    def add_stats_report(
        self, step, data, refers_to=None, query='stats --threads <THREADS> -v -s -', mounts=[], save_as=None,
        report_name='stats_report.log', report_description='generate stats report', check_exists=None, keep=False, enabled=True
    ):
        if enabled:
            query = query.replace('<THREADS>', str(self.get_vm_cpu()))
            source = data if mounts else f'$WD/{data}'
            target = save_as if save_as is not None else f'tmp/{report_name}'
            pre_check = f'[[ -f "$WD/{check_exists}" && -s "$WD/{check_exists}" ]] && ' if check_exists else ''
            self.add_action(
                'bcftools', step,
                report_description,
                f'-c \'{pre_check}bcftools {query} {source} > $WD/{target}\'',
                daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
                mounts=mounts,
                reports=[(target, f'{Pipeline.format_step(step) if refers_to is None else Pipeline.format_step(refers_to)}-{report_name}', keep)]
            )

    def check_config(self, name, section='avadx', flag='', is_file=False, default=None, quiet=False):
        option_value = 'NA'
        formatted = ''
        skip = ''
        flag_checked = f'{flag} ' if flag else flag
        if self.config.has_option(section, name):
            config_datadir_orig = self.config.get('DEFAULT', 'datadir', fallback=None)
            self.config.set('DEFAULT', 'datadir', str(VM_MOUNT / 'data'))
            option_value = self.config.get(section, f'{name}')
            self.config.set('DEFAULT', 'datadir', config_datadir_orig)
            if not is_file:
                formatted = f'{flag_checked}{option_value}'
            else:
                if Path(option_value).exists() or (VM_MOUNT / 'in' / option_value).is_file():
                    formatted = f'{flag_checked}{str(VM_MOUNT / "in")}/{option_value}'
                elif not Path(option_value).is_absolute():
                    option_value_abs = self.config_file.parent / option_value
                    if Path(option_value_abs).exists():
                        formatted = f'{flag_checked}{str(VM_MOUNT / "in")}/{option_value}'
                    else:
                        skip = 'file not found'
                else:
                    skip = 'file not found'
        else:
            skip = 'argument not defined in config'
        if not formatted and default:
            if not is_file:
                formatted = f'{flag_checked}{default}'
            else:
                formatted = f'{flag_checked}$WD/{default}'
        if skip and not quiet:
            if not default:
                self.log.warning(f'Ignoring optional argument: {flag_checked}{option_value} -> {skip}')
            else:
                self.log.debug(f'Using optional argument defaults: {flag_checked}{option_value} -> {skip}')
        return formatted

    def check_cv_scheme(self, csv_file):
        self.preprocess(rerun=True)
        with csv_file.open() as fin:
            reader = csv.reader(filter(lambda row: row[0] != '#', fin))
            class_labels = set()
            for row in reader:
                class_labels.add(row[2])
        if not self.is_prediction() and (len(class_labels) != 2 or '?' in class_labels):
            self.log.error(f'|5.10| No binary class labels (0/1) for Cross-Validation and Model Training. Labels found: {class_labels}')
            sys.exit(1)

    def run_container(self, container, args=[], daemon_args={}, uid=uuid.uuid1(), mounts=[], out_folder: Path = Path.cwd(), stdout=None, stderr=None, mkdirs=True, progress=False):
        wd_folder = (out_folder if out_folder else Path.cwd()) / str(uid) / 'wd'
        if mkdirs:
            if out_folder:
                (out_folder / str(uid) / 'out').mkdir(parents=True, exist_ok=True)
                (out_folder / str(uid) / 'out' / 'reports').mkdir(parents=True, exist_ok=True)
            wd_folder.mkdir(parents=True, exist_ok=True)
            (wd_folder / 'tmp').mkdir(parents=True, exist_ok=True)
            (wd_folder / 'vcf').mkdir(parents=True, exist_ok=True)
            (wd_folder / 'annovar').mkdir(parents=True, exist_ok=True)
        data_folder = Path(self.config.get('DEFAULT', 'datadir', fallback=wd_folder)).absolute()
        config_datadir_orig = self.config.get('DEFAULT', 'datadir', fallback=str(wd_folder))
        self.config.set('DEFAULT', 'datadir', str(data_folder))
        if self.daemon not in ['docker', 'singularity']:
            self.log.warning(f'Unknown daemon: {self.daemon}')
        if self.daemon == 'docker':
            wd = VM_MOUNT / 'out' / str(uid) / 'wd'
            out = VM_MOUNT / 'out' / str(uid) / 'out'
            data = VM_MOUNT / 'data'
        elif self.daemon == 'singularity':
            wd = VM_MOUNT / 'out' / str(uid) / 'wd'
            out = VM_MOUNT / 'out' / str(uid) / 'out'
            data = VM_MOUNT / 'data'
        else:
            wd = out_folder / str(uid) / 'wd'
            out = out_folder / str(uid) / 'out'
            data = data_folder

        cmd_base, cmd_suffix, args_parsed, daemon_args_parsed, env_exports_parsed = self.parse_args(args, daemon_args)

        if self.daemon == "docker":
            bind_mounts = []
            for m in mounts:
                if m[0] is None and m[1] is None:
                    self.log.debug(f'Invalid bind mount: {m}')
                elif m[0] is None or m[1] is None:
                    self.log.warning(f'Invalid bind mount: {m}')
                else:
                    bind_mounts += ['-v', f'{m[0]}:{m[1]}']
            cmd_base += ['--tty', '--rm'] \
                + (['-v', f'{out_folder.absolute()}:{VM_MOUNT / "out"}'] if out_folder.absolute().exists() else []) \
                + (['-v', f'{data_folder}:{data}'] if data_folder.exists() else []) \
                + bind_mounts + daemon_args_parsed
            cmd_base += [
                container,
            ] + cmd_suffix
        elif self.daemon == "singularity":
            bind_mounts = []
            for m in mounts:
                if m[0] is None and m[1] is None:
                    self.log.debug(f'Invalid bind mount: {m}')
                elif m[0] is None or m[1] is None:
                    self.log.warning(f'Invalid bind mount: {m}')
                else:
                    bind_mounts += ['-B', f'{m[0]}:{m[1]}']
            cmd_base += [] \
                + (['-B', f'{out_folder.absolute()}:{VM_MOUNT / "out"}'] if out_folder.absolute().exists() else []) \
                + (['-B', f'{data_folder}:{data}'] if data_folder.exists() else []) \
                + bind_mounts + daemon_args_parsed
            cmd_base += [
                '--containall',
                '--pwd', '/app',
                str(SINGULARITY_LIB / f'{container.split("/")[1]}_latest.sif'),
            ] + cmd_suffix
        else:
            cmd_base += [
                container,
            ] + cmd_suffix

        args_parsed = [
            a
            .replace('$WD', str(wd))
            .replace('$OUT', str(out))
            for a in args_parsed
        ]
        self.config.set('DEFAULT', 'datadir', str(data))
        config_pattern = r'config\[([^\]]*)\]'
        args_parsed = [
            re.sub(
                config_pattern,
                lambda x: self.config.get(x.group(1).split('.', 1)[0], x.group(1).split('.', 1)[1], fallback=''),
                a
            )
            for a in args_parsed
        ]
        self.config.set('DEFAULT', 'datadir', config_datadir_orig)
        cmd = cmd_base + [_.replace('//', '/') if _.startswith('/') else _ for _ in args_parsed]
        out = run_command(cmd, env_exports=env_exports_parsed, logger=self.log, poll=(LOG_LEVEL <= 10) or progress, progress=progress)
        if stdout and out:
            if stdout == 'print':
                print('\n'.join(out))
            elif stdout == 'log':
                self.log.info(' *** '.join(out))
            else:
                with (wd_folder / stdout).open('w') as fout:
                    fout.writelines([f'{_}\n' for _ in out])
        return out

    def parse_args(self, args, daemon_args):
        if self.daemon == 'docker':
            cmd_base = [
                self.daemon,
                self.REGISTRY[self.daemon]['CMD_RUN']
            ]
            return cmd_base, [], args, [arg for arg in daemon_args.get('docker', [])], {}
        elif self.daemon == "singularity":
            parsed_daemon_args, cmd_suffix, env_exports = [], [], {}
            cmd_base = [
                self.daemon,
                self.REGISTRY[self.daemon]['CMD_RUN']
            ]
            escape_args = True
            for arg in daemon_args.get('singularity', []):
                if arg.startswith('exec:'):
                    cmd_base[-1] = self.REGISTRY[self.daemon]['CMD_EXEC']
                    cmd_suffix += [arg[5:]]
                    escape_args = False
                elif arg.startswith('env:'):
                    key, val = arg[4:].split('=')
                    env_exports[f'SINGULARITYENV_{key}'] = val
                else:
                    pass
            parsed_args = [_.replace('"', r'\"') for _ in args] if escape_args else args
            return cmd_base, cmd_suffix, parsed_args, parsed_daemon_args, env_exports


def get_extra(x, y, fallback=None):
    res = [_.split('=')[1] for _ in x if _.startswith(y)]
    if len(res) == 1:
        return res[0]
    elif res:
        return res
    return fallback


def adjust_splits(auto_folds, cv_folds):
    while len(auto_folds) > 1 and len(auto_folds) > cv_folds:
        auto_folds = auto_folds[0:-2] + [auto_folds[-2] + auto_folds[-1]]
    return auto_folds


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog=f"{name}",
        description=(
            "AVA,Dx pipeline\n \n"
            "AVA,Dx (Analysis of Variation for Association with Disease)\n"
            "is a computational method for defining the functional role of\n"
            "DNA variation in complex diseases. AVA,Dx uses exonic variants\n"
            "from whole exome or genome sequencing data to extract a\n"
            "*disease* signal and predict *disease* status."
            "\n\n"
            "Version: %s [%s]\n" % (__version__, __releasedate__)
        ),
        epilog=(
            f"If you use *{name}* in published research, please cite:\n \n"
            "Wang Y., Miller M., ... Bromberg, Y. (2019).\n"
            "Identifying Crohn's disease signal from variome analysis.\n"
            "Genome medicine, 11(1), 59. [doi:10.1186/s13073-019-0670-6]\n"
            "(https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkx1209/4670955)\n \n"
            f"{name} is maintained by Maximilian Miller and Yana Bromberg.\n"
            "Feel free to contact us for support at services@bromberglab.org.\n \n"
            "This project is licensed under [NPOSL-3.0](http://opensource.org/licenses/NPOSL-3.0)\n \n"
            "Public web service: https://services.bromberglab.org/avadx\n \n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('config', nargs='?', type=Path, default=Path('avadx.ini'), action=check_config_ini(Path(__file__).parent))
    parser.add_argument('-c', '--configoverride', action='append')
    parser.add_argument('-a', '--action', action='append')
    parser.add_argument('-u', '--uid', type=str)
    parser.add_argument('-w', '--wd', type=Path, default=Path.cwd(),
                        help='working directory')
    parser.add_argument('-d', '--daemon', type=str, default='docker', choices=['docker', 'singularity'],
                        help='container engine')
    parser.add_argument('-i', '--info', action='store_true',
                        help='print pipeline info')
    parser.add_argument('-U', '--update', type=str, choices=['all', 'data', 'vm'],
                        help='update pipeline databases/datasources')
    parser.add_argument('-I', '--init', action='store_true',
                        help='init pipeline - retrieve all required databases/datasources')
    parser.add_argument('-p', '--preprocess', action='store_true',
                        help='run input pre-processing')
    parser.add_argument('-P', '--postprocess', action='store_true',
                        help='run pipeline post-processing')
    parser.add_argument('-V', '--visualize', type=str,
                        help='create visualizations / result plots')
    parser.add_argument('-s', '--sampleprediction', type=Path,
                        help='path to samples file for which predictions should be generated')
    parser.add_argument('-M', '--model', type=Path,
                        help='path to joblib model file other then the default when using -s/--sampleprediction')
    parser.add_argument('-r', '--retrieve', type=str,
                        help='retrieve data source')
    parser.add_argument('-S', '--steps', type=str,
                        help='steps to run, e.g. 1 (start from 1) or 1-2 (start from 1, stop at 2)')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='set verbosity level (v: ERROR, vv: WARNING, vvv: info, vvvv: debug); default: log level INFO')
    parser.add_argument('-D', '--Debug', type=Path, nargs='?',
                        help='path of local avadx codebase to mount to containers')
    parser.add_argument('-L', '--logfile', type=Path, const=Path('pipeline.log'),
                        nargs='?', help='redirect logs to file')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='if flag is set console output is logged to file')

    parser.lvl = 1
    namespace, extra = parser.parse_known_args()
    parser_actions = []
    for action in namespace.action if namespace.action else []:
        if action not in parser_actions:
            parser.add_argument(f'--{action}', action='append', default=[])
            parser_actions += [action]
    parser.lvl = 2
    namespace, extra = parser.parse_known_args()
    if namespace.verbose == 0:
        namespace.verbose = 20
    else:
        namespace.verbose -= 1
        namespace.verbose = max(10, 40 - min(10 * namespace.verbose, 40)) if namespace.verbose >= 0 else 0
    global LOG_LEVEL, LOG_FILE, QUIET
    LOG_LEVEL = namespace.verbose
    LOG_FILE = namespace.logfile
    QUIET = namespace.quiet
    del namespace.verbose
    del namespace.logfile
    if namespace.steps:
        namespace.entrypoint = float(namespace.steps.split('-')[0]) if namespace.steps.split('-')[0] != '' else None
        namespace.exitpoint = float(namespace.steps.split('-')[1]) if len(namespace.steps.split('-')) == 2 else None
    else:
        namespace.entrypoint, namespace.exitpoint = None, None
    del namespace.steps

    return namespace, extra


def main(pipeline, extra):
    app = AVADx(pipeline=pipeline)
    if pipeline.kwargs.get('wd'):
        pipeline.kwargs.get('wd').mkdir(parents=True, exist_ok=True)
    timer_start = timer()
    for action in pipeline.actions:
        app.run(action, pipeline.uid, pipeline.kwargs, extra)
    app.log.info(f'Runtime: {(timer() - timer_start):.3f} seconds')
    if pipeline.is_prediction():
        app.log.info(f'Predictions: {pipeline.get_prediction_out()}')
    else:
        app.log.info(f'Results: {pipeline.get_wd().absolute()}')


def get_mounts(pipeline, *config, exit_on_error=False, mount_as=None, show_warning=True):
    mounts = []
    for section, option in config:
        cfg = pipeline.config.get(section, option, fallback=None)
        if cfg:
            cfg_path = Path(cfg)
            if not cfg_path.is_absolute():
                cfg_path = pipeline.config_file.parent / cfg_path
            if cfg_path.exists():
                if mount_as is None:
                    mount_path = cfg_path.name
                else:
                    mount_path = mount_as if Path(mount_as).is_absolute() else mount_as
                mounts += [(cfg_path.absolute(), VM_MOUNT / 'in' / mount_path)]
            else:
                if show_warning:
                    pipeline.log.warning(f'Could not mount {cfg_path}. [Path not found]')
                if exit_on_error:
                    sys.exit(1)
                else:
                    mounts += [(None, None)]
        else:
            mounts += [(None, None)]
    return mounts


def run_init(uid, kwargs, extra, config, daemon):
    kwargs['entrypoint'] = 0
    kwargs['exitpoint'] = 1
    pipeline = run_all(uid, kwargs, extra, config, daemon)
    shutil.rmtree(pipeline.get_wd())


def run_all(uid, kwargs, extra, config, daemon, dry_run=False):
    pipeline = Pipeline(kwargs=kwargs, uid=uid, config_file=config, daemon=daemon)
    return run_all_p(pipeline, extra, dry_run)


def run_all_p(pipeline, extra, dry_run=False):
    pipeline.set_host_cpu()
    pipeline.get_vm_resources(dry_run)
    VM_CPU = pipeline.get_vm_cpu()
    VM_MEM = pipeline.get_vm_mem()
    cpu_limit = int(pipeline.config.get('avadx', 'resources.cpu', fallback=0))
    mem_limit = int(pipeline.config.get('avadx', 'resources.mem', fallback=0))
    if cpu_limit > 0:
        VM_CPU_new = min(cpu_limit, VM_CPU)
        if VM_CPU != VM_CPU_new:
            VM_CPU = VM_CPU_new
            pipeline.set_vm_cpu(VM_CPU)
    if mem_limit > 0:
        VM_MEM_new = min(mem_limit, VM_MEM)
        if VM_MEM != VM_MEM_new:
            VM_MEM = VM_MEM_new
            pipeline.set_vm_mem(VM_MEM)
    if not dry_run:
        (pipeline.get_wd() / 'out').mkdir(parents=True, exist_ok=True)
        pipeline.save_run_config()
        pipeline.save_run_info()
        CFG = VM_MOUNT / 'out' / str(pipeline.uid) / 'out' / 'pipeline_config.ini'
    else:
        CFG = VM_MOUNT / 'in' / 'avadx.ini'
    WD = pipeline.kwargs['wd'] / str(pipeline.uid) / 'wd'
    OUT = pipeline.kwargs['wd'] / str(pipeline.uid) / 'out'
    hgref = pipeline.config.get('avadx', 'hgref', fallback='hg19')
    hgref_mapped = pipeline.ASSEMBLY_MAPPING[hgref]
    gnomAD_filter = True if pipeline.config.get('avadx', 'filter.gnomad.enabled', fallback='yes') == 'yes' else False
    vqsr_PASS_filter = True if pipeline.config.get('avadx', 'filter.vqsrPASS.enabled', fallback='yes') == 'yes' else False
    site_quality_filter = True if pipeline.config.get('avadx', 'filter.sitequality.enabled', fallback='yes') == 'yes' else False
    ABAD_filter = True if pipeline.config.get('avadx', 'filter.ABAD.enabled', fallback='yes') == 'yes' else False
    MR_filter = True if pipeline.config.get('avadx', 'filter.MR.enabled', fallback='yes') == 'yes' else False
    sex_filter = True if pipeline.config.get('avadx', 'filter.sex.enabled', fallback='yes') == 'yes' else False
    mito_filter = True if pipeline.config.get('avadx', 'filter.mito.enabled', fallback='yes') == 'yes' else False
    analyze_indels = True if pipeline.config.get('avadx', 'analyze.indels', fallback='no') == 'yes' else False
    ethseq_splits_cfg = int(pipeline.config.get('avadx', 'ethseq.split', fallback=0))
    ethseq_threads = int(pipeline.get_host_cpu() // VM_CPU)
    ethseq_splits = ethseq_splits_cfg if ethseq_splits_cfg > 0 else pipeline.get_splits(threads=ethseq_threads)
    R_MAX_VSIZE_ethseq = f'{int(VM_MEM/ethseq_threads)}Gb'
    outliers_available = True if pipeline.check_config('outliers', is_file=True, quiet=True) else False
    outliers_break = True if pipeline.config.get('avadx', 'outliers.break', fallback='no') == 'yes' else False
    genescorefn_available = True if pipeline.check_config('genescore.fn', is_file=True, quiet=True) else False
    variantscorefn_available = True if pipeline.check_config('variantscore.fn', is_file=True, quiet=True) else False
    genescore_normalize = True if pipeline.config.get('avadx', 'genescore.standardize', fallback='no') == 'yes' else False
    fselectionclass_available = True if pipeline.check_config('fselection.class', is_file=True, quiet=True) else False
    modelclass_available = True if pipeline.check_config('model.class', is_file=True, quiet=True) else False
    featurelist_available = True if pipeline.check_config('featurelist', is_file=True, quiet=True) else False
    create_filter_reports = True if pipeline.config.get('avadx', 'filter.reports', fallback='yes') == 'yes' else False
    is_init = pipeline.kwargs['init']
    pipeline.entrypoint = 2.40 if outliers_available else pipeline.entrypoint
    if pipeline.kwargs['entrypoint'] is not None:
        pipeline.entrypoint = pipeline.kwargs['entrypoint']
    if outliers_break and not outliers_available:
        pipeline.exitpoint = 2.40
    if pipeline.kwargs['exitpoint'] is not None:
        pipeline.exitpoint = pipeline.kwargs['exitpoint']

    header = f'    AVA,Dx - {__version__} {__releasedate__}    '
    divider = f' {"".join(["_" for _ in range(0,len(header))])}'
    print(
        f'\n {header} \n{divider}\n\n'
        f'  * VM      : {pipeline.daemon}\n'
        f'  * Cores   : {pipeline.get_vm_cpu()}\n'
        f'  * Memory  : {pipeline.get_vm_mem()}\n'
        f'{divider}\n'
    )

    # 0     Init (downloads & data source preprocessing)  -------------------------------------------------------------------- #
    gnomad_base_path = Path(pipeline.config.get("DEFAULT", "annovar.humandb", fallback=WD))
    avadx_base_path = Path(pipeline.config.get("DEFAULT", "avadx.data", fallback=WD))
    refseq_base_path = Path(pipeline.config.get("DEFAULT", "refseq.data", fallback=WD))
    models_base_path = Path(pipeline.config.get("DEFAULT", "ethseq.models", fallback=WD))

    # 0.11  Retrieve gnomad_exome database
    pipeline.add_action(
        'annovar', 0.11,
        f'verify/download database: {hgref}_gnomad_exome',
        f'-c \'annotate_variation.pl -thread {VM_CPU} -buildver {hgref} -downdb -webfrom annovar gnomad_exome config[DEFAULT.annovar.humandb]/; '
        'mv config[DEFAULT.annovar.humandb]/annovar_downdb.log config[DEFAULT.annovar.humandb]/annovar_downdb_gnomad_exome.log\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(gnomad_base_path)
    )

    # 0.12  Retrieve gnomad_genome database
    pipeline.add_action(
        'annovar', 0.12,
        f'verify/download database: {hgref}_gnomad_genome',
        f'-c \'annotate_variation.pl -thread {VM_CPU} -buildver {hgref} -downdb -webfrom annovar gnomad_genome config[DEFAULT.annovar.humandb]/; '
        'mv config[DEFAULT.annovar.humandb]/annovar_downdb.log config[DEFAULT.annovar.humandb]/annovar_downdb_gnomad_genome.log\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(gnomad_base_path)
    )

    # 0.121 preprocess: Generate gnomad_exome_allAFabove0 / Generate gnomad_genome_allAFabove0
    pipeline.add_action(
        'gnomad_ALLabove0_preprocess', 0.121,
        'preprocess: split gnomad database files',
        f'-c \'/app/bash/avadx/gnomad_ALLabove0_preprocess.sh config[DEFAULT.annovar.humandb]/{hgref}_gnomad_exome.txt '
        f'config[DEFAULT.annovar.humandb]/{hgref}_gnomad_genome.txt config[DEFAULT.avadx.data]\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(avadx_base_path)
    )

    # 0.122 main: Generate gnomad_exome_allAFabove0 / Generate gnomad_genome_allAFabove0
    pipeline.add_action(
        'generate_gnomad_above0', 0.122,
        'filter gnomad database above0 (exome & genome)',
        f'/app/R/avadx/generate_gnomad_above0.R config[DEFAULT.avadx.data] {hgref}',
        outdir=(avadx_base_path)
    )

    # 0.123 postprocess: Generate gnomad_exome_allAFabove0 / Generate gnomad_genome_allAFabove0
    pipeline.add_action(
        'gnomad_ALLabove0_postprocess', 0.123,
        'postprocess: bgzip and tabix for above0 gnomad database files',
        f'-c \'/app/bash/avadx/gnomad_ALLabove0_postprocess.sh config[DEFAULT.avadx.data]/{hgref}_gnomad_exome_allAFabove0.txt '
        f'config[DEFAULT.avadx.data]/{hgref}_gnomad_genome_allAFabove0.txt config[DEFAULT.avadx.data]\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(avadx_base_path)
    )

    # 0.13  Retrieve refGene database
    pipeline.add_action(
        'annovar', 0.13,
        f'verify/download database: {hgref}_refGene',
        f'-c \'annotate_variation.pl -thread {VM_CPU} -buildver {hgref} -downdb -webfrom annovar refGene config[DEFAULT.annovar.humandb]/; '
        'mv config[DEFAULT.annovar.humandb]/annovar_downdb.log config[DEFAULT.annovar.humandb]/annovar_downdb_refGene.log\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(refseq_base_path)
    )

    # 0.14  Retrieve varidb database
    pipeline.add_action(
        'run_retrieve', 0.14,
        'verify/download database: varidb',
        f'{CFG} --retrieve varidb --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=[(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')],
        outdir=(avadx_base_path),
        logs=('print', None)
    )

    # 0.15  Retrieve CPDB pathway mapping
    pipeline.add_action(
        'run_retrieve', 0.15,
        'verify/download database: CPDB pathway mapping',
        f'{CFG} --retrieve cpdb --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=[(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')],
        outdir=(avadx_base_path),
        logs=('print', None)
    )

    # 0.16  Retrieve EthSEQ Model
    ethseq_model_name = 'Exonic.All.Model.gds'
    ethseq_model_baseurl = 'https://github.com/cibiobcg/EthSEQ_Data/raw/master/EthSEQ_Models/'
    pipeline.add_action(
        'avadx', 0.16,
        f'verify/download EthSEQ Model: {ethseq_model_name}',
        f'-c \'[[ -f config[DEFAULT.ethseq.models]/{ethseq_model_name} ]] || wget -O config[DEFAULT.ethseq.models]/{ethseq_model_name} {ethseq_model_baseurl}{ethseq_model_name}\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(models_base_path)
    )

    # 0.17  Retrieve prot_seqs.fa
    pipeline.add_action(
        'run_retrieve', 0.17,
        'verify/download reference sequences/stats: refseq',
        f'{CFG} --retrieve refseq --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=[(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')],
        outdir=(avadx_base_path),
        logs=('print', None)
    )

    # 0.18  Generate Transcript-ProtLength.csv
    pipeline.add_action(
        'generate_transcripts_protlength', 0.18,
        'generate reference proteins (refseq) stats',
        f'/app/R/avadx/generate_refseq_stats.R config[DEFAULT.refseq.data]/{hgref_mapped[0]}.{hgref_mapped[1]}_feature_table.txt '
        f'config[DEFAULT.refseq.data]/{hgref_mapped[0]}.{hgref_mapped[1]}_protein.faa config[DEFAULT.avadx.data]',
        outdir=(avadx_base_path)
    )

    # 1   Preprocess -------------------------------------------------------------------------- #
    # 1.0.0 POTENTIAL PREPROCESSING SLOT

    # 1.0.1 Extract list of samples
    mounts_step1_0_1 = get_mounts(pipeline, ('avadx', 'vcf'), exit_on_error=False if is_init or pipeline.is_prediction() else True, show_warning=False if pipeline.is_prediction() else True)
    step1_0_1_out = 'tmp/all_samples.txt'
    pipeline.add_stats_report(
        1.01, mounts_step1_0_1[0][1], query='query -l', mounts=mounts_step1_0_1, save_as=step1_0_1_out,
        report_name='all_samples.log', report_description='extract sample list', keep=True
    )

    # 1.0.2 Run pre-processing
    mounts_step1_0_2 = get_mounts(pipeline, ('avadx', 'samples'), show_warning=not pipeline.prediction) + [(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')]
    pipeline.add_action(
        'run_preprocess', 1.02,
        'AVA,Dx pipeline preprocess',
        f'{CFG} --preprocess --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=mounts_step1_0_2,
        outdir=(avadx_base_path),
        logs=('print', None)
    )

    # 1.0.3 Generate stats report
    mounts_step1_0_1 = get_mounts(pipeline, ('avadx', 'vcf'), exit_on_error=False if is_init or pipeline.is_prediction() else True, show_warning=False if pipeline.is_prediction() else True)
    pipeline.add_stats_report(1.03, mounts_step1_0_1[0][1], refers_to=1.00, mounts=mounts_step1_0_1, enabled=create_filter_reports)

    # 1.1-7 Variant QC -------------------------------------------------------------------------- #

    # 1.1   Extract individuals of interest (diseased and healthy individuals of interest).
    step1_1_out = 'vcf/1_1.vcf.gz'
    mounts_step1_1 = get_mounts(pipeline, ('avadx', 'vcf'), exit_on_error=False if is_init or pipeline.is_prediction() else True, show_warning=False if pipeline.is_prediction() else True)
    pipeline.add_action(
        'bcftools', 1.10,
        'filter for individuals of interest ',
        f'view --threads {VM_CPU} -S $WD/tmp/sampleids.txt {mounts_step1_1[0][1]} -Oz -o $WD/{step1_1_out}',
        mounts=mounts_step1_1,
        results=([step1_1_out], [])
    )

    # 1.1.1 Generate stats report
    pipeline.add_stats_report(1.11, step1_1_out, refers_to=1.10, check_exists='tmp/sampleids.txt', enabled=create_filter_reports)

    # 1.2   POTENTIAL PREPROCESSING SLOT

    # 1.3   OPTIONAL - vqsr PASS filter: Remove variant sites which did not pass the VQSR standard.
    if vqsr_PASS_filter:
        step1_3_out = 'vcf/1_3.vcf.gz'
        pipeline.add_action(
            'bcftools', 1.30,
            'filter variant sites < VQSR standard',
            f'filter --threads {VM_CPU} -i \'FILTER="PASS"\' $WD/{step1_1_out} -Oz -o $WD/{step1_3_out}'
        )
    step1_3_out = step1_3_out if vqsr_PASS_filter else step1_1_out

    # 1.3.1 Generate stats report
    if vqsr_PASS_filter:
        pipeline.add_stats_report(1.31, step1_3_out, refers_to=1.30, enabled=create_filter_reports)

    # 1.4   OPTIONAL - Remove variant sites by site-wise quality.
    #       Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150.
    if site_quality_filter:
        step1_4_out = 'vcf/1_4.vcf.gz'
        pipeline.add_action(
            'bcftools', 1.40,
            'filter variant sites by site-wise quality',
            f'view --threads {VM_CPU} -i \'QUAL>config[avadx.qc.site.quality]&'
            'AVG(FMT/DP)<=config[avadx.qc.site.mean_dp_upper]&'
            'AVG(FMT/DP)>=config[avadx.qc.site.mean_dp_lower]\' '
            f'$WD/{step1_3_out} -Oz -o $WD/{step1_4_out}'
        )
    step1_4_out = step1_4_out if site_quality_filter else step1_3_out

    # 1.4.1 Generate stats report
    if site_quality_filter:
        pipeline.add_stats_report(1.41, step1_4_out, refers_to=1.40, enabled=create_filter_reports)

    # 1.5   OPTIONAL - Missing Rate (MR) Filter:
    #       only variants with MR > threshold are retained
    if MR_filter:
        step1_5_out = 'vcf/1_5.vcf.gz'
        pipeline.add_action(
            'bcftools', 1.5,
            'filter by missing rate (MR)',
            f'view --threads {VM_CPU} -i \'F_MISSING<=config[avadx.qc.call.MR]\' $WD/{step1_4_out} -Oz -o $WD/{step1_5_out}'
        )
    step1_5_out = step1_5_out if MR_filter else step1_4_out

    # 1.5.1 Generate stats report
    if MR_filter:
        pipeline.add_stats_report(1.51, step1_5_out, refers_to=1.50, enabled=create_filter_reports)

    # 1.6   OPTIONAL - Check individual call quality. In filterVCF_by_ABAD.py:
    #       good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4;
    #       bad individual GTs are converted into missing "./.";
    #       low call rate is determined as a call rate < 80%,
    #       i.e. missing rate >= 20%. Variant sites with a low call rate are removed.
    if ABAD_filter:
        step1_6_out = 'vcf/1_6.vcf.gz'
        pipeline.add_action(
            'filterVCF_by_ABAD', 1.60,
            'check individual call quality by allele balance and missing rate',
            f'/app/python/avadx/filterVCF_by_ABAD.py $WD/{step1_5_out} $WD/{step1_6_out} '
            'config[avadx.qc.call.AB_low] config[avadx.qc.call.AB_high] '
            'config[avadx.qc.call.DP] config[avadx.qc.call.GQ] config[avadx.qc.call.MR]',
            daemon_args={'docker': ['--entrypoint=python'], 'singularity': ['exec:python']},
            reports=[(f'{step1_6_out}.log', '1_6-call_quality_filter.log')]
        )
    step1_6_out = step1_6_out if ABAD_filter else step1_5_out

    # 1.7   Convert the chromosome annotation if the chromosomes are recorded as "chr1" instead of "1":
    step1_7_out = 'vcf/1_7.vcf.gz'
    pipeline.add_action(
        'bcftools', 1.7,
        'convert the chromosome annotations',
        f'annotate --threads {VM_CPU} --rename-chrs $WD/tmp/chr_to_number.txt $WD/{step1_6_out} -Oz -o $WD/{step1_7_out}'
    )
    step1_out = step1_7_out

    # 2     Individual QC ----------------------------------------------------------------------- #

    # 2.2   Quality check - Check quality outliers by examining nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing:

    # 2.1.0 Generate stats report
    step2_1_0_out = 'tmp/qc_filters.stats.txt'
    pipeline.add_stats_report(2.10, step1_out, refers_to=1.60, save_as=step2_1_0_out, keep=True)

    # 2.1.1 Draw individual quality figure:
    step2_1_1_out = 'tmp/quality_control_samples_PCA.pdf'
    pipeline.add_action(
        'stats_quality_pca', 2.11,
        'draw individual quality figures',
        f'/app/R/avadx/stats_quality_pca.R -f $WD/{step2_1_0_out} -o $WD/{step2_1_1_out}',
        reports=[(step2_1_1_out, '2_1-1-quality_control_samples_PCA.pdf')]
    )

    # 2.2   Ethnicity check - Annotate ethnicity with EthSEQ R package:

    # 2.2.1 OPTIONAL: If the number of individuals exceeds certain number, "memory exhausted" error may occur.
    #  Manually divide input VCF into chunks of individuals and run EthSEQ separately for each chunk:
    step2_2_1_out = 'tmp/EthSEQ_split_ids.txt'
    step2_2_1_outfolder = 'tmp/EthSEQ_splits'
    step2_2_1_splits = 'tmp/EthSEQ_splits.txt'
    pipeline.add_action(
        'bcftools', 2.21,
        'create samples batches',
        f'-c \'bcftools query -l $WD/{step1_out} > $WD/{step2_2_1_out}; '
        f'SPLIT=$(SAMPLES=$(wc -l < $WD/{step2_2_1_out}); echo $((SAMPLES <= {ethseq_splits} ? 0 : {ethseq_splits}))); '
        f'[[ $SPLIT -gt 0 ]] && split -d -l $SPLIT $WD/{step2_2_1_out} $WD/{step2_2_1_outfolder}/split_ || cp -f $WD/{step2_2_1_out} $WD/{step2_2_1_outfolder}/split_00; '
        f'(cd $WD/{step2_2_1_outfolder} && ls -f -1 split_*) > $WD/{step2_2_1_splits}\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(WD / step2_2_1_outfolder)
    )

    # 2.2.2 OPTIONAL: Clean VCF format for EthSEQ input for all splits:
    step2_2_2_splits = 'tmp/EthSEQ_splits'
    pipeline.add_action(
        'bcftools', 2.22,
        'EthSEQ preprocessing (VCF cleanup)',
        f'-c \'bcftools view --threads {VM_CPU} -S $TASK $WD/{step1_out} | '
        f'bcftools annotate --threads {VM_CPU} --remove "ID,INFO,FORMAT" | '
        f'sed /^##/d | gzip --stdout > $WD/{step2_2_2_splits}/source_$(basename $TASK)_EthSEQinput.vcf.gz\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/'),
        outdir=(WD / step2_2_2_splits),
        resources={'cpu': VM_CPU}
    )

    # 2.2.3 Run EthSEQ:
    #       Increase memory for larger datasets
    step2_2_3_outfolder = 'EthSEQ_reports'
    pipeline.add_action(
        'ethnicity_EthSEQ', 2.23,
        'run EthSEQ',
        f'-c \'mkdir -p $WD/{step2_2_3_outfolder}/$(basename $TASK) && '
        'Rscript /app/R/avadx/ethnicity_EthSEQ.R '
        f'$WD/{step2_2_2_splits}/source_$(basename $TASK)_EthSEQinput.vcf.gz '
        f'config[DEFAULT.ethseq.models] $WD/{step2_2_3_outfolder}/$(basename $TASK) {VM_CPU}\'',
        daemon_args={'docker': ['--entrypoint=bash', f'--env=R_MAX_VSIZE={R_MAX_VSIZE_ethseq}'], 'singularity': ['exec:/bin/bash', f'env:R_MAX_VSIZE={R_MAX_VSIZE_ethseq}']},
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/'),
        resources={'cpu': VM_CPU}
    )

    # 2.2.4 Collect ethnicity prediction reports
    pipeline.add_action(
        'avadx', 2.24,
        'collect EthSEQ reports',
        f'$WD/{step2_2_3_outfolder} -mindepth 2 -regex \'.*\(txt\|pdf\)$\' -exec bash -c \'cp $1 $OUT/reports/2_2-4-EthSEQ_$(basename $(dirname $1))_$(basename $1)\' _ {{}} \;',  # noqa: W605
        daemon_args={'docker': ['--entrypoint=find'], 'singularity': ['exec:find']}
    )

    # 2.2.5 Summarize and visualize ethnicity prediction reports
    mounts_visualization = [(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')]
    pipeline.add_action(
        'run_visualization', 2.25,
        'summarize/visualize EthSEQ reports',
        f'{CFG} --visualize piechart:ethseq --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=mounts_visualization,
        logs=('print', None)
    )

    # 2.3   Relatedness check:
    #       Check relatedness within datasets withe the SNPRelate R package.
    #       A default kinship > 0.3 is considered to be related.
    step2_3_outfolder = 'SNPRelate'
    step2_3_out = 'tmp/relatedness.gds'
    pipeline.add_action(
        'relatedness', 2.30,
        'check relatedness using SNPRelate',
        f'/app/R/avadx/relatedness.R -i $WD/{step1_out} -g $WD/{step2_3_out} -c config[avadx.kinship] -o $WD/{step2_3_outfolder} -t {VM_CPU}',
        outdir=(WD / step2_3_outfolder),
        reports=[(Path(step2_3_outfolder) / 'IBD_related.txt', '2_3-IBD_related.txt'), (Path(step2_3_outfolder) / 'IBD.txt', '2_3-IBD.txt')]
    )

    # 2.4   Remove individual outliers
    #       Outlier individual IDs should be combined from the above PCA,
    #       ethnicity annotation, and relatedness calculation
    #       to a file outliers.txt (one ID per row).
    step2_4_out = 'vcf/2_4.vcf.gz'
    if outliers_available:
        mounts_step2_4 = get_mounts(pipeline, ('avadx', 'outliers'), exit_on_error=False if is_init or pipeline.is_prediction() else True)
        pipeline.add_action(
            'bcftools', 2.40,
            'summarize outliers',
            f'view --threads {VM_CPU} -S ^{mounts_step2_4[0][1]} $WD/{step1_out} -Oz -o $WD/{step2_4_out}',
            mounts=mounts_step2_4
        )
    step2_4_out = step2_4_out if outliers_available else step1_out

    # 2.4.1 Generate stats report
    if outliers_available:
        pipeline.add_stats_report(2.41, step2_4_out, refers_to=2.40, enabled=create_filter_reports)

    # 2.5   Additional filters
    #       - Remove sex chromosomes flag X,Y (false)
    #       - Remove mito chromosome flag M (false)
    run_2_5_filter = any([sex_filter, mito_filter])
    if run_2_5_filter:
        exclude = (['X', 'Y'] if sex_filter else []) + (['M', 'MT'] if mito_filter else [])
        regions_flt = f'-t ^{",".join(exclude)} ' if exclude else ' '
        step2_5_out = 'vcf/2_5.vcf.gz'
        pipeline.add_action(
            'bcftools', 2.5,
            'run remaining filters',
            f'view --threads {VM_CPU} '
            f'{regions_flt}'
            f'$WD/{step2_4_out} -Oz -o $WD/{step2_5_out}'
        )
    step2_5_out = step2_5_out if run_2_5_filter else step2_4_out

    # 2.6   Split SNV and InDel calls to separated files because they use different QC thresholds.
    #       Current AVA,Dx works mainly with SNPs. InDels need another set of standards for QC.

    # 2.6.1 snps
    step2_6_1_out = 'vcf/2_6-1-snps.vcf.gz'
    pipeline.add_action(
        'bcftools', 2.61,
        'filter snps',
        f'view --threads {VM_CPU} --types snps $WD/{step2_5_out} -Oz -o $WD/{step2_6_1_out}'
    )

    # 2.6.2 indels
    if analyze_indels:
        step2_6_2_out = 'vcf/2_6-2-indels.vcf.gz'
        pipeline.add_action(
            'bcftools', 2.62,
            'filter indels',
            f'view --threads {VM_CPU} --types indels $WD/{step2_5_out} -Oz -o $WD/{step2_6_2_out}'
        )
    step2_indels_out = step2_6_2_out if analyze_indels else None

    # 2.7   OPTIONAL - gnomAD filter: filtering out variants that were not recorded in the gnomAD database.
    #       The gnomAD reference used here is the ANNOVAR gnomAD filexx_gnomad_exome.txt and hgxx_gnomad_genome.txt.
    #       Note that tabix is required for indexing to run this script.
    step2_7_out = 'vcf/2_7.vcf.gz'
    if gnomAD_filter:
        pipeline.add_action(
            'filterVCF_by_gnomAD', 2.7,
            'filter variants missing in gnomAD database',
            f'avadx.filterVCF_by_gnomAD $WD $WD/{step2_6_1_out} $WD/{step2_7_out} '
            f'config[DEFAULT.avadx.data]/{hgref}_gnomad_exome_allAFabove0.txt.gz '
            f'config[DEFAULT.avadx.data]/{hgref}_gnomad_genome_allAFabove0.txt.gz'
        )
    step2_snps_out = step2_7_out if gnomAD_filter else step2_6_1_out

    # 2.7.1 Generate stats report
    pipeline.add_stats_report(2.71, step2_snps_out, refers_to=2.70, enabled=create_filter_reports)

    # 3     Query/Calculate SNAP scores for all variants ------------------------------------------ #

    # 3.1   Get all variant annotations with ANNOVAR for cleaned VCF:
    step3_1_out = 'annovar/3_1.avinput'
    pipeline.add_action(
        'annovar', 3.10,
        'convert VCF file to ANNOVAR input format',
        f'convert2annovar.pl -format vcf4old $WD/{step2_snps_out} -outfile $WD/{step3_1_out}',
        reports=[('annovar/convert2annovar.log', '3_1-convert2annovar.log')],
        logs=('annovar/convert2annovar.log', None)
    )

    # 3.2   Annotate using <hgref> RefSeq:
    step3_2_out = 'annovar/3_1.avinput.exonic_variant_function'
    pipeline.add_action(
        'annovar', 3.20,
        f'annotate using {hgref} RefSeq',
        f'annotate_variation.pl -thread {VM_CPU} -buildver {hgref} $WD/{step3_1_out} config[DEFAULT.annovar.humandb]/',
        logs=('annovar/annotate_variation.log', None)
    )

    # 3.3   Create varidb query
    step3_3_outfolder = 'varidb'
    pipeline.add_action(
        'generate_SNAP_query', 3.30,
        'generate SNAP scores query',
        f'/app/R/avadx/check_missing_SNAP.R $WD/{step3_2_out} config[DEFAULT.avadx.data] $WD/{step3_3_outfolder}',
        outdir=(WD / step3_3_outfolder),
        reports=[
            (Path('varidb') / 'annovar_stats.csv', '3_3-annovar_stats.csv'),
            (Path('varidb') / 'SNAP_extract_vars.log', '3_3-extract_variants.log')
        ]
    )

    # 3.4   Run varidb pre-filter
    step3_4_outfolder = 'varidb'
    pipeline.add_action(
        'varidb_prefilter', 3.40,
        'check variants for SNAP length limit and correct sequence mapping',
        f'/app/python/avadx/filter_vars.py $WD/{step3_4_outfolder}',
        daemon_args={'docker': ['--entrypoint=python'], 'singularity': ['exec:python']},
        reports=[(Path(step3_4_outfolder) / 'varidb_prefilter_report.log', '3_4-varidb_prefilter.log')]
    )

    # 3.5   Query varidb for SNAP mutations
    step3_5_out = 'varidb/varidb_query_result.csv'
    pipeline.add_action(
        'varidb', 3.50,
        'query SNAP variants from varidb',
        '-D config[DEFAULT.avadx.data]/varidb.db '
        f'-Q $WD/{step3_3_outfolder}/varidb_query_nonsyn_filtered.ids '
        f'-f $WD/{step3_3_outfolder}/varidb_query_nonsyn_filtered.fa '
        f'-o $WD/{step3_5_out} '
        '-R $WD/varidb/varidb_query_report.txt '
        '-C query variant score -S tab -s',
        reports=[
            (Path('varidb') / 'varidb_query_report.txt', '3_5-varidb_report.log'),
            (Path('varidb') / 'varidb_query_result_failed.csv', '3_5-varidb_failed.tsv'),
            (Path('varidb') / 'varidb_query_report_info.txt', '3_5-scoring_functions_info.log')
        ]
    )

    # 4     Gene score calculation ---------------------------------------------------------------- #

    # 4.1   Convert cleaned VCF to individual ANNOVAR annotation files by:
    step4_1_outfolder = 'annovar/annotations'
    step4_1_out = 'tmp/annovar_annotations.txt'
    pipeline.add_action(
        'annovar', 4.10,
        'convert cleaned VCF to single annovar annotation files',
        f'-c \'convert2annovar.pl -format vcf4 $WD/{step2_snps_out} -outfile $WD/{step4_1_outfolder}/sample -allsample; '
        f'(cd $WD/{step4_1_outfolder} && ls -f -1 sample.*.avinput) > $WD/{step4_1_out}\'',
        daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
        outdir=(WD / step4_1_outfolder)
    )

    if analyze_indels:
        step4_1_1_outfolder = 'annovar/annotations_indels'
        step4_1_1_out = 'tmp/annovar_annotations_indels.txt'
        pipeline.add_action(
            'annovar', 4.11,
            'convert cleaned VCF to single annovar annotation files (indels)',
            f'-c \'convert2annovar.pl -format vcf4 $WD/{step2_indels_out} -outfile $WD/{step4_1_1_outfolder}/sample -allsample; '
            f'(cd $WD/{step4_1_1_outfolder} && ls -f -1 sample.*.avinput) > $WD/{step4_1_1_out}\'',
            daemon_args={'docker': ['--entrypoint=bash'], 'singularity': ['exec:/bin/bash']},
            outdir=(WD / step4_1_1_outfolder)
        )

    # 4.2   Annotate all sample*.avinput files
    pipeline.add_action(
        'annovar', 4.20,
        'generate annovar annotation',
        f'annotate_variation.pl -thread {VM_CPU} -build {hgref} $TASK config[DEFAULT.annovar.humandb]/',
        tasks=(None, WD / step4_1_out, f'$WD/{step4_1_outfolder}/'),
        resources={'cpu': VM_CPU}
    )

    if analyze_indels:
        pipeline.add_action(
            'annovar', 4.21,
            'generate annovar annotation (indels)',
            f'annotate_variation.pl -thread {VM_CPU} -build {hgref} $TASK config[DEFAULT.annovar.humandb]/',
            tasks=(None, WD / step4_1_1_out, f'$WD/{step4_1_1_outfolder}/'),
            resources={'cpu': VM_CPU}
        )

    # 4.3   Compute gene scores
    step4_3_outfolder = 'genescores'
    genescore_file = Path(pipeline.check_config('genescore.fn', quiet=dry_run)).suffix == '.py'
    genescorefn_mnt = get_mounts(pipeline, ('avadx', 'genescore.fn'), exit_on_error=False if is_init or pipeline.is_prediction() else genescorefn_available) if genescore_file else []
    variantscore_file = Path(pipeline.check_config('variantscore.fn', quiet=dry_run)).suffix == '.py'
    variantscorefn_mnt = get_mounts(pipeline, ('avadx', 'variantscore.fn'), exit_on_error=False if is_init or pipeline.is_prediction() else variantscorefn_available) if variantscore_file else []
    mounts_step4_3 = genescorefn_mnt + variantscorefn_mnt
    pipeline.add_action(
        'cal_genescore_make_genescore', 4.30,
        'calculate gene score',
        '/app/python/avadx/genescore.py -a $TASK.exonic_variant_function '
        f'-s $WD/{step3_5_out} -m config[DEFAULT.avadx.data]/Transcript-ProtLength_cleaned.csv '
        f'-g {genescorefn_mnt[0][1] if genescore_file else "config[avadx.genescore.fn]"} '
        f'-v {variantscorefn_mnt[0][1] if variantscore_file else "config[avadx.variantscore.fn]"} '
        f'-t config[avadx.varidb.predictors] -n config[avadx.genescore.standardize] {"--indels " if analyze_indels else ""}'
        f'-o $WD/{step4_3_outfolder}',
        tasks=(None, WD / step4_1_out, f'$WD/{step4_1_outfolder}/'),
        daemon_args={'docker': ['--entrypoint=python'], 'singularity': ['exec:python']},
        mounts=mounts_step4_3,
        outdir=(WD / step4_3_outfolder),
        resources={'cpu': 1}
    )

    # 4.4   Merge gene scores:
    #  Assuming there are 500 individuals in the dataset, 500 resulting files will be generated (e.g. sample.S001.gs).
    #  Merge them into a data frame where a row is an individual and a column is a gene (protein):
    step4_4_outfolder = 'genescores'
    pipeline.add_action(
        'merge_genescore', 4.40,
        'merge single gene score result files',
        f'/app/python/avadx/genescore.py -M -m config[DEFAULT.avadx.data]/Transcript-ProtLength_cleaned.csv '
        f'-c $WD/tmp/cv-scheme.csv -g {genescorefn_mnt[0][1] if genescore_file else "config[avadx.genescore.fn]"} '
        f'-v {variantscorefn_mnt[0][1] if variantscore_file else "config[avadx.variantscore.fn]"} '
        f'-n config[avadx.genescore.standardize] -r $WD/{step4_3_outfolder} -o $OUT/{step4_4_outfolder} ',
        fns={'pre': (pipeline.check_cv_scheme, [WD / 'tmp' / 'cv-scheme.csv'])},
        daemon_args={'docker': ['--entrypoint=python'], 'singularity': ['exec:python']},
        mounts=mounts_step4_3,
        outdir=(OUT / step4_4_outfolder)
    )

    # 5   Feature selection (FS) and model building --------------------------------------------- #
    #  AVA,Dx by default uses K-S (Kolmogorov–Smirnov) test for FS, random forest for model building,
    #  and 10-fold cross validation to test the predictability of the top-ranking genes.
    #  Other FS methods and machine learning models are also included.
    #  User needs to provide a cross-validation scheme file. For example, we split Tourette dataset (e.g. yale-1)
    #  into 10 folds and the individuals from the same family enter the same fold,
    #  so that to compare sick v.s. healthy instead of differentiating families.

    # 5.1   Cross-validation:
    step5_1_outfolder = 'results'
    fselection_file = Path(pipeline.check_config('fselection.class', quiet=dry_run)).suffix == '.py'
    fselectionclass_mnt = get_mounts(pipeline, ('avadx', 'fselection.class'), exit_on_error=False if is_init or pipeline.is_prediction() else fselectionclass_available, mount_as='/app/python/avadx/feature_selections/fselection_avadx.py') if fselection_file else []
    model_file = Path(pipeline.check_config('model.class', quiet=dry_run)).suffix == '.py'
    modelclass_mnt = get_mounts(pipeline, ('avadx', 'model.class'), exit_on_error=False if is_init or pipeline.is_prediction() else modelclass_available, mount_as='/app/python/avadx/models/model_avadx.py') if model_file else []
    featurelist_mnt = get_mounts(pipeline, ('avadx', 'featurelist'), exit_on_error=False if is_init or pipeline.is_prediction() else featurelist_available) if featurelist_available else []
    use_featurelist = f' -F {featurelist_mnt[0][1]}' if featurelist_mnt else ''
    mounts_step5_1 = fselectionclass_mnt + modelclass_mnt + featurelist_mnt
    pipeline.add_action(
        'ava_model', 5.10,
        'perform model cross-validation',
        f'/app/python/avadx/model.py -g $OUT/{step4_4_outfolder}/GeneScoreTable_{"normalized" if genescore_normalize else "raw"}.csv '
        '-f config[avadx.cv.featureselection] -m config[avadx.cv.model] -c $WD/tmp/cv-scheme.csv '
        '-p config[DEFAULT.avadx.data]/Transcript-ProtLength_cleaned.csv -v config[avadx.cv.varcutoff] -V config[avadx.cv.sklearnvariance] '
        '-P config[avadx.cv.ks.pvalcutoff] -G config[avadx.cv.topgenes] -S config[avadx.cv.steps] -T config[avadx.model.topgenes] '
        f'-o $OUT/{step5_1_outfolder} -w $WD/{step5_1_outfolder} -C {VM_CPU}{use_featurelist} {"-q" if pipeline.daemon == "singularity" else ""}',
        fns={'pre': (pipeline.check_cv_scheme, [WD / 'tmp' / 'cv-scheme.csv'])},
        daemon_args={'docker': ['--entrypoint=python'], 'singularity': ['exec:python']},
        mounts=mounts_step5_1,
        progress=True,
        outdir=(OUT / step5_1_outfolder)
    )

    # 5.2   Check pathway over-representation:
    step5_2_outfolder = 'pathways'
    pipeline.add_action(
        'FS_CVgeneOverRep_kfold', 5.20,
        'check pathway over-representation',
        f'/app/R/avadx/FS-CVgeneOverRep-kfold.R -f $WD/results/fselection_performance.csv '
        f'-t $WD/{step5_1_outfolder}/crossval_genes.csv -b $OUT/{step4_4_outfolder}/GeneScoreTable_{"normalized" if genescore_normalize else "raw"}_variation_filtered.csv '
        '-n config[avadx.pathways.topgenes] -d config[DEFAULT.avadx.data]/CPDB_pathways_genesymbol.tab '
        f'-o $OUT/{step5_2_outfolder}',
        outdir=(OUT / step5_2_outfolder)
    )

    # 6   Visualizations ------------------------------------------------------------------------ #
    mounts_visualization = [(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')]
    pipeline.add_action(
        'run_visualization', 6.00,
        'AVA,Dx visualizations',
        f'{CFG} --visualize heatmap_clust_dend:genescores --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=mounts_visualization,
        logs=('print', None)
    )

    # 7   Cleanup ------------------------------------------------------------------------------- #
    mounts_postprocess = [(pipeline.config_file.absolute(), VM_MOUNT / 'in' / 'avadx.ini')]
    pipeline.add_action(
        'run_postprocess', 7.00,
        'AVA,Dx pipeline postprocess',
        f'{CFG} --postprocess --wd $WD {("-" if LOG_LEVEL > 0 else "") + "v"*((50-LOG_LEVEL)//10)}',
        mounts=mounts_postprocess,
        logs=('print', None)
    )

    # 8   Prediction ---------------------------------------------------------------------------- #
    if pipeline.is_prediction() or dry_run:
        predictions_outfolder = 'predictions'
        predictions_pred_id, predictions_model, predictions_genescores, predictions_features = get_extra(extra, 'pred_id', 'null'), get_extra(extra, 'model', False), get_extra(extra, 'genescores', False), get_extra(extra, 'features', False)
        mounts_predictions = []
        mounts_predictions += [(predictions_model, VM_MOUNT / 'in' / 'model.joblib')] if predictions_model else []
        mounts_predictions += [(predictions_genescores, VM_MOUNT / 'in' / 'genescores.csv')] if predictions_genescores else []
        pred_features = Path(predictions_features).exists() if predictions_features else False
        mounts_predictions += [(predictions_features, VM_MOUNT / 'in' / 'features.txt')] if pred_features else []
        pipeline.add_action(
            'run_prediction', 8.00,
            'AVA,Dx predictions',
            f'/app/python/avadx/model.py -M {predictions_pred_id} {mounts_predictions[0][1] if predictions_model else "None"} {mounts_predictions[1][1] if predictions_genescores else "None"} {mounts_predictions[2][1] if pred_features else "None"} '
            f'{variantscorefn_mnt[0][1] if variantscore_file else "config[avadx.variantscore.fn]"} {genescorefn_mnt[0][1] if genescore_file else "config[avadx.genescore.fn]"} '
            f'-o $OUT/{predictions_outfolder} -C {VM_CPU}',
            daemon_args={'docker': ['--entrypoint=python'], 'singularity': ['exec:python']},
            mounts=mounts_predictions + mounts_step4_3,
            logs=('print', None),
            outdir=(OUT / predictions_outfolder)
        )

    # RUN --------------------------------------------------------------------------------------- #
    if not dry_run:
        main(pipeline, extra)

    return pipeline


def init():
    namespace, extra = parse_arguments()
    actions = namespace.action
    config = namespace.config
    daemon = namespace.daemon
    uid = namespace.uid
    del namespace.action
    del namespace.config
    del namespace.daemon
    del namespace.uid
    if namespace.info:
        Pipeline(actions, kwargs=vars(namespace), config_file=config, daemon=daemon) \
            .info(quiet=False, run_args=[uid, vars(namespace), extra, config, daemon, True])
    elif namespace.init:
        AVADx.init_vm(daemon)
        run_init(uid, vars(namespace), extra, config, daemon)
    elif namespace.update:
        namespace.init = True
        if namespace.update == 'vm':
            AVADx.init_vm(daemon)
        elif namespace.update == 'data':
            run_init(uid, vars(namespace), extra, config, daemon)
        else:
            AVADx.init_vm(daemon)
            run_init(uid, vars(namespace), extra, config, daemon)
    elif namespace.preprocess:
        pipeline = Pipeline(actions, kwargs=vars(namespace), uid=uid, config_file=config, daemon=daemon)
        pipeline.preprocess()
    elif namespace.retrieve:
        pipeline = Pipeline(actions, kwargs=vars(namespace), uid=uid, config_file=config, daemon=daemon)
        pipeline.retrieve(namespace.retrieve)
    elif namespace.postprocess:
        pipeline = Pipeline(actions, kwargs=vars(namespace), uid=uid, config_file=config, daemon=daemon)
        pipeline.postprocess()
    elif namespace.visualize:
        namespace.init = True
        pipeline = Pipeline(actions, kwargs=vars(namespace), uid=uid, config_file=config, daemon=daemon)
        pipeline.visualization(namespace.visualize)
    elif namespace.sampleprediction:
        namespace.init = False
        pipeline = Pipeline(actions, kwargs=vars(namespace), uid=uid, config_file=config, daemon=daemon)
        pipeline.predict(uid, vars(namespace), extra, config, daemon)
    elif actions is None:
        run_all(uid, vars(namespace), extra, config, daemon)
    else:
        pipeline = Pipeline(actions, kwargs=vars(namespace), uid=uid, config_file=config, daemon=daemon)
        main(pipeline, extra)


if __name__ == "__main__":
    init()
