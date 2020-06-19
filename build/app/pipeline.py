#!/usr/bin/env python3

import argparse
import re
import csv
import uuid
import shlex
import random
import configparser
import multiprocessing
from pathlib import Path
from datetime import datetime
from timeit import default_timer as timer
from .logger import Logger
from .helper import run_command, flatten, runningInDocker
from . import __version__, __releasedate__


QUIET = False
DOCKER_MNT = Path('/mnt')
LOG_LEVEL = 'INFO'
LOG_FILE = None


class AVADxMeta:
    def __init__(self, pipeline):
        self.log = self.get_logger()
        self.pipeline = pipeline

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
        args = shlex.split(kwargs.get(name).pop(0))
        daemon_args = shlex.split(kwargs.get(name + "_darg").pop(0))
        taskflag, taskfile, taskprefix = kwargs.get(name + "_tasks").pop(0)
        if taskfile:
            with (taskfile).open() as fin:
                tasklist = [_.strip() for _ in fin.readlines()]
        else:
            tasklist = [None]
        outdir = kwargs.get(name + "_outd").pop(0)
        mounts = kwargs.get(name + "_mounts").pop(0)
        fns_pre, fns_post = kwargs.get(name + "_fns").pop(0)
        description = kwargs.get(name + "_desc").pop(0)
        level = kwargs.get(name + "_lvl").pop(0)
        log_stdout, log_stderr = kwargs.get(name + "_logs").pop(0)
        if outdir:
            outdir.mkdir(parents=True, exist_ok=True)
        for tid, task in enumerate(tasklist, 1):
            args_ = list(args)
            task_info = f' [{tid}/{len(tasklist)}] ({task})' if task else ''
            self.log.info(f'|{level:.2f}| {name}{task_info}: {description}')
            self.log.debug(f'|{level:.2f}| {name}{task_info}: started {datetime.now()}')
            timer_start = timer()
            if fns_pre:
                fns_pre()
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
            self.pipeline.run_container(
                container, args=args_, daemon_args=daemon_args, uid=uid,
                mounts=mounts,
                out_folder=kwargs.get('wd'),
                stdout=log_stdout,
                stderr=log_stderr
            )
            if fns_post:
                fns_post()
            self.log.info(f'|{level:.2f}| {name}{task_info}: took {(timer() - timer_start):.3f} seconds')

    def run_preprocess(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-meta', 'run_preprocess', uid, kwargs)

    def bcftools(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-bcftools', 'bcftools', uid, kwargs)

    def avadx(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-meta', 'avadx', uid, kwargs)

    def filterVCF_by_ABAD(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-meta', 'filterVCF_by_ABAD', uid, kwargs)

    def filterVCF_by_gnomAD(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-tabix', 'filterVCF_by_gnomAD', uid, kwargs)

    def stats_quality_pca(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'stats_quality_pca', uid, kwargs)

    def ethnicity_EthSEQ(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'ethnicity_EthSEQ', uid, kwargs)

    def ethnicity_EthSEQ_summary(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'ethnicity_EthSEQ_summary', uid, kwargs)

    def relatedness(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'relatedness', uid, kwargs)

    def annovar(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-annovar', 'annovar', uid, kwargs)

    def check_missing_SNAP(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'check_missing_SNAP', uid, kwargs)

    def varidb(self, uid, **kwargs):
        self.run_method('bromberglab/varidb', 'varidb', uid, kwargs)

    def cal_genescore_make_genescore(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'cal_genescore_make_genescore', uid, kwargs)

    def merge_genescore(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'merge_genescore', uid, kwargs)

    def FS_CVperf_kfold(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'FS_CVperf_kfold', uid, kwargs)


class Pipeline:

    REGISTRY = {
        'singularity': {
            'CMD_RUN': 'run'
        },
        'docker': {
            'CMD_RUN': 'run'
        }
    }

    def __init__(self, actions=[], kwargs={}, uid=uuid.uuid1(), config_file=None, daemon='docker'):
        self.log = self.get_logger()
        self.actions = actions
        self.kwargs = kwargs
        self.uid = uid
        self.config_file = config_file
        self.config = self.load_config()
        self.daemon = daemon
        self.is_docker_vm = runningInDocker()
        self.init_daemon()
        self.entrypoint = 0
        self.exitpoint = None

    def get_logger(self):
        logger = Logger(self.__class__.__name__, level=LOG_LEVEL)
        logger.addConsoleHandler()
        log = logger.getLogger()
        return log

    def init_daemon(self):
        if self.daemon == 'docker':
            # import docker
            pass

    def load_config(self):
        config = configparser.ConfigParser()
        if self.config_file and self.config_file.exists():
            config.read(str(self.config_file))
        else:
            config_file_mnt = DOCKER_MNT / 'in' / self.config_file.name
            if config_file_mnt and config_file_mnt.exists():
                config.read(str(config_file_mnt))
        return config

    def info(self):
        print(f'AVA,Dx {__version__} {__releasedate__}')

    def preprocess(self):
        if self.is_docker_vm:
            wd_folder = self.kwargs.get('wd')
            samples_path = DOCKER_MNT / 'in' / self.config.get('avadx', 'samples')
            samplesids_path = wd_folder / 'sampleids.txt'
            cvscheme_path = wd_folder / 'cv-scheme.csv'
        else:
            wd_folder = self.kwargs.get('wd') / str(self.uid) / 'wd'
            wd_folder.mkdir(parents=True, exist_ok=True)
            samples = self.check_config('samples', is_file=True, quiet=False)
            if not samples:
                self.log.warning('Required input missing - sample annotation file')
                return
            samples_path = Path(self.config.get('avadx', 'samples'))
            samplesids_path = wd_folder / 'sampleids.txt'
            cvscheme_path = wd_folder / 'cv-scheme.csv'
            if not samples_path.exists():
                self.log.warning(f'Could not read required input: {samples_path}')
                return
        indices = lambda labels, search: [i for i, x in enumerate(labels) if x == search]
        folds = lambda samples, foldcnt: [samples[i * foldcnt:(i + 1) * foldcnt] for i in range((len(samples) + foldcnt - 1) // foldcnt)]
        with samples_path.open() as fin, samplesids_path.open('w') as fout_ids, cvscheme_path.open('w') as fout_cv:
            samples, labels, col3 = [], [], []
            reader = csv.reader(fin)
            header = next(reader)
            has_groups = len(header) == 3 and header[2] == 'group'
            has_folds = len(header) == 3 and header[2] == 'fold'
            for row in reader:
                samples += [row[0]]
                labels += [row[1]]
                if has_groups or has_folds:
                    col3 += [row[2]]
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
                        + f' - proceeding with {cv_folds}-fold split')
                grouped_samples_shuffled = list(grouped_samples)
                random.shuffle(grouped_samples_shuffled)
                group_folds = folds(grouped_samples_shuffled, len(grouped_samples_shuffled)//cv_folds)
                for gfold in group_folds:
                    auto_folds += [flatten([grouped_samples.get(g) for g in gfold])]
                auto_folds = adjust_splits(auto_folds, cv_folds)
                if len(auto_folds) != cv_folds:
                    self.log.warn(f'Auto-generated cross-validation folds ({len(auto_folds)}) differ from specified folds ({cv_folds})')
                for idx, fold in enumerate(auto_folds, 1):
                    for sidx in fold:
                        fout_cv.write(f'{samples[sidx]},{idx},{labels[sidx]}\n')
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
                    self.log.warning(
                        f'Too few samples ({cv_folds_max}) for specified {cv_folds_config}-fold cross-validation'
                        + f' - proceeding with {cv_folds}-fold split')
                samples_indices = [_ for _ in range(6)]
                random.shuffle(samples_indices)
                auto_folds = folds(samples_indices, len(samples_indices)//cv_folds)
                auto_folds = adjust_splits(auto_folds, cv_folds)
                if len(auto_folds) != cv_folds:
                    self.log.warn(f'Auto-generated cross-validation folds ({len(auto_folds)}) differ from specified folds ({cv_folds})')
                for idx, fold in enumerate(auto_folds, 1):
                    for sidx in fold:
                        fout_cv.write(f'{samples[sidx]},{idx},{labels[sidx]}\n')
                split_type = 'auto-generated sample'
                split_description = f'{len(auto_folds)}-fold split' if cv_folds < cv_folds_max else 'leave-one-out splits'
            self.log.info(f'|0.10| Using {split_type} based cross-validation scheme for {split_description}')

    def add_action(
            self, action, level=None, description='', args='', daemon_args='',
            tasks=(None, None, None), fns=(None, None), outdir=None, mounts=[], logs=(None, None)
            ):
        if level is None or (level >= self.entrypoint and (self.exitpoint is None or level < self.exitpoint)):
            self.actions += [action]
            self.kwargs[action] = self.kwargs.get(action, []) + [args]
            self.kwargs[f'{action}_lvl'] = self.kwargs.get(f'{action}_lvl', []) + [level]
            self.kwargs[f'{action}_desc'] = self.kwargs.get(f'{action}_desc', []) + [description]
            self.kwargs[f'{action}_darg'] = self.kwargs.get(f'{action}_darg', []) + [daemon_args]
            self.kwargs[f'{action}_tasks'] = self.kwargs.get(f'{action}_tasks', []) + [tasks]
            self.kwargs[f'{action}_fns'] = self.kwargs.get(f'{action}_fns', []) + [fns]
            self.kwargs[f'{action}_outd'] = self.kwargs.get(f'{action}_outd', []) + [outdir]
            self.kwargs[f'{action}_mounts'] = self.kwargs.get(f'{action}_mounts', []) + [mounts]
            self.kwargs[f'{action}_logs'] = self.kwargs.get(f'{action}_logs', []) + [logs]

    def check_config(self, name, section='avadx', flag='', is_file=False, default=None, quiet=False):
        option_value = 'NA'
        formatted = ''
        skip = ''
        flag_checked = f'{flag} ' if flag else flag
        if self.config.has_option(section, name):
            config_datadir_orig = self.config.get('DEFAULT', 'datadir', fallback=None)
            self.config.set('DEFAULT', 'datadir', str(DOCKER_MNT / 'data'))
            option_value = self.config.get(section, f'{name}')
            self.config.set('DEFAULT', 'datadir', config_datadir_orig)
            if not is_file:
                formatted = f'{flag_checked}{option_value}'
            elif Path(option_value).exists() or (DOCKER_MNT / 'in' / option_value).is_file():
                formatted = f'{flag_checked}{str(DOCKER_MNT / "in")}/{option_value}'
            elif is_file and not (Path(option_value).exists() or (DOCKER_MNT / 'in' / option_value).is_file()):
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

    def run_container(self, container, args=[], daemon_args=[], uid=uuid.uuid1(), mounts=[], out_folder:Path=Path.cwd(), stdout=None, stderr=None):
        if out_folder:
            (out_folder / str(uid) / 'out').mkdir(parents=True, exist_ok=True)
        wd_folder = out_folder / str(uid) / 'wd'
        wd_folder.mkdir(parents=True, exist_ok=True)
        data_folder = Path(self.config.get("DEFAULT", "datadir", fallback=wd_folder)).absolute()
        config_datadir_orig = self.config.get('DEFAULT', 'datadir')
        self.config.set('DEFAULT', 'datadir', str(data_folder))
        if self.daemon == 'docker':
            wd = DOCKER_MNT / 'out' / str(uid) / 'wd'
            out = DOCKER_MNT / 'out' / str(uid) / 'out'
            data = DOCKER_MNT / 'data'
        else:
            wd = out_folder / str(uid) / 'wd'
            out = out_folder / str(uid) / 'out'

        cmd_base = [
            self.daemon,
            self.REGISTRY[self.daemon]['CMD_RUN']
        ]
        if self.daemon == "docker":
            bind_mounts = []
            for m in mounts:
                bind_mounts += ['-v', f'{m[0]}:{m[1]}']
            cmd_base += [
                '--rm',
                '-v', f'{out_folder.absolute()}:{DOCKER_MNT / "out"}',
                '-v', f'{data_folder}:{data}',
            ] + bind_mounts + daemon_args
        cmd_base += [
            container,
        ]
        args_parsed = [
            a
             .replace('$WD', str(wd))
             .replace('$OUT', str(out))
            for a in args
        ]
        self.config.set('DEFAULT', 'datadir', str(data))
        config_pattern = r'config\[([^\]]*)\]' #r'config\[(.*)\]'
        args_parsed = [
            re.sub(
                config_pattern,
                lambda x: self.config.get(x.group(1).split('.', 1)[0], x.group(1).split('.', 1)[1], fallback=''),
                a
            )
            for a in args_parsed
        ]
        self.config.set('DEFAULT', 'datadir', config_datadir_orig)
        cmd = cmd_base + [_.replace('//', '/') for _ in args_parsed]
        out = run_command(cmd, logger=self.log)
        if stdout:
            if stdout == 'print':
                print('\n'.join(out))
            else:
                with (wd_folder / stdout).open('w') as fout:
                    fout.writelines([f'{_}\n' for _ in out])


get_extra = lambda x, y: [_.split('=')[1] for _ in x if _.startswith(y)]


def adjust_splits(auto_folds, cv_folds):
    last_elem = -1
    while len(auto_folds) > cv_folds:
        tmp = []
        for idx, elem in enumerate(auto_folds):
            if idx == last_elem:
                tmp += [elem]
            elif idx == last_elem + 1 and idx + 1 < len(auto_folds):
                tmp += [elem + auto_folds[idx + 1]]
            elif idx > last_elem + 2:
                tmp += [elem]
        last_elem += 1
        auto_folds = tmp
    return auto_folds


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog=f"{__name__}",
        description="AVA,Dx pipeline\n \n" \
            + "description\n" \
            + "XXX\n" \
            + "XXX\n \n"
            + "Public web service: https://services.bromberglab.org/avadx-meta\n \n" \
            + "Version: %s [%s]\n" % (__version__, __releasedate__),
        epilog=f"If you use *{__name__}* in published research, please cite:\n \n" \
            + "Wang Y., Miller M., ... Bromberg, Y. (2019).\n" \
            + "Identifying Crohn's disease signal from variome analysis.\n" \
            + "Genome medicine, 11(1), 59. [doi:10.1186/s13073-019-0670-6]\n" \
            + "(https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkx1209/4670955)\n \n" \
            + f"{__name__} is developed by Yanran Wang and Maximilian Miller.\n" \
            + "Feel free to contact us for support at services@bromberglab.org.\n \n" \
        "This project is licensed under [NPOSL-3.0](http://opensource.org/licenses/NPOSL-3.0)\n \n" \
            + "Test: XXX\n \n",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('config', nargs='?', type=Path, default=Path('pipeline.ini'))
    parser.add_argument('-a', '--action', action='append')
    parser.add_argument('-w', '--wd', type=Path, default=Path.cwd(),
                        help='working directory')
    parser.add_argument('-c', '--cpu', type=int, default=multiprocessing.cpu_count(),
                        help='max cpus per thread; default: all available cores')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='number of threads; default: 1')
    parser.add_argument('-d', '--daemon', type=str, default='docker', choices=['docker', 'singularity'],
                        help='container engine')
    parser.add_argument('-i', '--info', action='store_true',
                        help='print pipeline info')
    parser.add_argument('-p', '--preprocess', action='store_true',
                        help='run input preprocessing')
    parser.add_argument('-e', '--entrypoint', type=float, default=0.0,
                        help='(re)start pipeline at specified entrypoint')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='set verbosity level; default: log level INFO')
    parser.add_argument('-L', '--logfile', type=Path, const=Path('pipeline.log'),
                        nargs='?', help='redirect logs to file')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help='if flag is set console output is logged to file')

    namespace, extra = parser.parse_known_args()
    parser_actions = []
    for action in namespace.action if namespace.action else []:
        if action not in parser_actions:
            parser.add_argument(f'--{action}',  action='append', default=[])
            parser_actions += [action]
    namespace, extra = parser.parse_known_args()
    namespace.verbose = 40 - min(10*namespace.verbose, 40) if namespace.verbose >= 0 else 0
    global LOG_LEVEL, LOG_FILE
    LOG_LEVEL = namespace.verbose
    LOG_FILE = namespace.logfile
    del namespace.verbose
    del namespace.logfile

    return namespace, extra


def main(pipeline, extra):
    app = AVADxMeta(pipeline=pipeline)
    uid = get_extra(extra, '--uid')
    if uid:
        pipeline.uid = uid[0]
    if pipeline.kwargs.get('wd'):
        pipeline.kwargs.get('wd').mkdir(parents=True, exist_ok=True)
    timer_start = timer()
    for action in pipeline.actions:
        app.run(action, pipeline.uid, pipeline.kwargs, extra)
    app.log.info(f'Total runtime: {(timer() - timer_start):.3f} seconds')


def get_mounts(pipeline, *config):
    mounts = []
    for section, option in config:
        cfg = pipeline.config.get(section, option, fallback=None)
        if cfg:
            cfg_path = Path(cfg)
            if cfg_path.exists():
                mounts += [(cfg_path.absolute(), DOCKER_MNT / 'in' / cfg_path)]
            else:
                pipeline.log.debug(f'Could not mount {cfg_path}. Path not found.')
    return mounts


def run_all(kwargs, extra, config, daemon):
    pipeline = Pipeline(kwargs=kwargs, config_file=config, daemon=daemon)
    uid = get_extra(extra, '--uid')
    if uid:
        pipeline.uid = uid[0]
    WD = kwargs['wd'] / str(pipeline.uid) / 'wd'
    OUT = kwargs['wd'] / str(pipeline.uid) / 'out'
    hgref = pipeline.config.get('avadx', 'hgref')
    gnomADfilter = True if pipeline.config.get('avadx', 'gnomadfilter.enabled', fallback='no') == 'yes' else False
    outliers_available = True if pipeline.check_config('outliers', is_file=True, quiet=True) else False
    outliers_break = True if pipeline.config.get('avadx', 'outliers.break', fallback='no') == 'yes' else False
    pipeline.entrypoint = 2.40 if outliers_available else 0
    if kwargs['entrypoint']:
        pipeline.entrypoint = kwargs['entrypoint']
    if outliers_break and not outliers_available:
        pipeline.exitpoint = 2.40

    # 0     Downloads & Info -------------------------------------------------------------------- #

    # 0.1   Print Info
    mounts_preprocess = get_mounts(pipeline, ('avadx', 'samples')) + [(pipeline.config_file.absolute(), DOCKER_MNT / 'in' / 'pipeline.ini')]
    pipeline.add_action(
        'run_preprocess', 0.10,
        f'AVA,Dx pipeline preprocess',
        f'app.pipeline --preprocess --wd $WD {"-v "*(LOG_LEVEL//10)}',
        mounts=mounts_preprocess,
        logs=('print', None)
    )

    # 0.2   Retrieve gnomad_exome database
    pipeline.add_action(
        'annovar', 0.20,
        f'verify/download database: {hgref}_gnomad_exome',
        f'-c \\"[[ -f config[DEFAULT.annovar.humandb]/{hgref}_gnomad_exome.txt ]] || (annotate_variation.pl -buildver {hgref} -downdb -webfrom annovar gnomad_exome config[DEFAULT.annovar.humandb]/; '
        + f'mv config[DEFAULT.annovar.humandb]/annovar_downdb.log config[DEFAULT.annovar.humandb]/annovar_downdb_gnomad_exome.log)\\"',
        '--entrypoint=bash'
    )

    # 0.3   Retrieve refGene database
    pipeline.add_action(
        'annovar', 0.30,
        f'verify/download database: {hgref}_refGene',
        f'-c \\"[[ -f config[DEFAULT.annovar.humandb]/{hgref}_refGene.txt ]] || (annotate_variation.pl -buildver {hgref} -downdb -webfrom annovar refGene config[DEFAULT.annovar.humandb]/; '
        + f'mv config[DEFAULT.annovar.humandb]/annovar_downdb.log config[DEFAULT.annovar.humandb]/annovar_downdb_refGene.log)\\"',
        '--entrypoint=bash'
    )

    # 1     Variant QC -------------------------------------------------------------------------- #
    
    # 1.1   Extract individuals of interest (diseased and healthy individuals of interest).
    step1_1_in = 'config[avadx.vcf]'
    step1_1_out = 'source_samp.vcf.gz'
    mounts_step1_1 = get_mounts(pipeline, ('avadx', 'vcf'))
    pipeline.add_action(
        'bcftools', 1.10,
        'filter for individuals of interest ',
        f'view -S $WD/sampleids.txt {mounts_step1_1[0][1]} -Oz -o $WD/{step1_1_out}',
        mounts=mounts_step1_1
    )

    # 1.2   Remove variant sites which did not pass the VQSR standard.
    step1_2_out = 'source_samp_pass.vcf.gz'
    pipeline.add_action(
        'bcftools', 1.20,
        'filter variant sites < VQSR standard',
        f'filter -i \'FILTER=\\"PASS\\"\' $WD/{step1_1_out} -Oz -o $WD/{step1_2_out}'
    )

    # 1.3   Split SNV and InDel calls to separated files because they use different QC thresholds.
    #       Current AVA,Dx works mainly with SNPs. InDels need another set of standards for QC.

    # 1.3.1 snps
    step1_3_1_out = 'source_samp_pass_snps.vcf.gz'
    pipeline.add_action(
        'bcftools', 1.31,
        'split SNV and InDel calls',
        f'view --types snps $WD/{step1_2_out} -Oz -o $WD/{step1_3_1_out}'
    )
    # 1.3.2 indels
    step1_3_2_out = 'source_samp_pass_indels.vcf.gz'
    pipeline.add_action(
        'bcftools', 1.32,
        'split SNV and InDel calls',
        f'view --types indels $WD/{step1_2_out} -Oz -o $WD/{step1_3_2_out}'
    )

    # 1.4   Remove variant sites by site-wise quality.
    #       Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150.
    step1_4_out = 'source_samp_pass_snps_site-v.vcf.gz'
    pipeline.add_action(
        'bcftools', 1.40,
        'filter variant sites by site-wise quality',
        f'view -i \\"QUAL>config[avadx.qc.site.quality] & '
        + f'AVG(FMT/DP)<=config[avadx.qc.site.mean_dp_upper] & '
        + f'AVG(FMT/DP)>=config[avadx.qc.site.mean_dp_lower]\\" '
        + f'$WD/{step1_3_1_out} -Oz -o $WD/{step1_4_out}'
    )

    # 1.5   Check individual call quality. In filterVCF_by_ABAD.py:
    #       good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4;
    #       bad individual GTs are converted into missing "./.";
    #       low call rate is determined as a call rate < 80%,
    #       i.e. missing rate >= 20%. Variant sites with a low call rate are removed.
    step1_5_out = 'source_samp_pass_snps_site-v_gt-v.vcf.gz'
    pipeline.add_action(
        'filterVCF_by_ABAD', 1.50,
        'check individual call quality',
        f'avadx.filterVCF_by_ABAD $WD/{step1_4_out} $WD/{step1_5_out} '
        + f'config[avadx.qc.call.AB_low] config[avadx.qc.call.AB_high] '
        + f'config[avadx.qc.call.DP] config[avadx.qc.call.GQ] config[avadx.qc.call.MR]'
    )

    # 1.6   OPTIONAL - gnomAD filter: filtering out variants that were not recorded in the gnomAD database.
    #       The gnomAD reference used here is the ANNOVAR gnomAD filexx_gnomad_exome.txt and hgxx_gnomad_genome.txt.
    #       Check the input path of the two reference files before running the script.
    #       Note that tabix is required for indexing to run this script.
    step1_6_out = 'source_samp_pass_snps_site-v_gt-v_rmchr_gnomad.vcf.gz'
    if gnomADfilter:
        # 1.6.1 Convert the chromosome annotation if the chromosomes are recorded as "chr1" instead of "1":
        step1_6_1_out = 'source_samp_pass_snps_site-v_gt-v_rmchr.vcf.gz'
        pipeline.add_action(
            'bcftools', 1.61,
            'convert the chromosome annotations',
            f'annotate --rename-chrs config[avadx.chr2number] $WD/{step1_5_out} -Oz -o $WD/{step1_6_1_out}'
        )

        # 1.6.2 Then remove variants that are not in gnomAD database:
        pipeline.add_action(
            'filterVCF_by_gnomAD', 1.62,
            'filter variants missing in gnomAD database',
            f'avadx.filterVCF_by_gnomAD $WD/{step1_6_1_out} $WD/{step1_6_out} '
            + f'config[DEFAULT.annovar.gnomad]/{hgref}_config[avadx.gnomadfilter.exome] '
            + f'config[DEFAULT.annovar.gnomad]/{hgref}_config[avadx.gnomadfilter.genome]'
        )
    step1_out = step1_6_out if gnomADfilter else step1_5_out

    # 2     Individual QC ----------------------------------------------------------------------- #
    
    # 2.2   Quality check - Check quality outliers by examine nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing:

    # 2.1.0 Output quality metrics after variant QC:
    step2_1_0_out = 'source_samp_pass_snps_site-v_gt-v_rmchr_gnomad.stats.txt'
    pipeline.add_action(
        'bcftools', 2.10,
        'generate quality metrics after variant QC',
        f'-c \\"bcftools stats -v -s - $WD/{step1_out} > $WD/{step2_1_0_out}\\"',
        '--entrypoint=bash'
    )

    # 2.1.1 Draw individual quality figure:
    step2_1_1_out = 'source_samp_pass_snps_site-v_gt-v_stats.PSC.PCA.pdf'
    pipeline.add_action(
        'stats_quality_pca', 2.11,
        'draw individual quality figures',
        f'/app/R/avadx/stats_quality_pca.R -f $WD/{step2_1_0_out} -o $OUT/{step2_1_1_out}'
    )

    # 2.2   Ethnicity check - Annotate ethnicity with EthSEQ R package:

    # 2.2.1 OPTIONAL: If the number of individuals exceeds certain number, "memory exhausted" error may occur.
    #  Manually divide input VCF into chunks of individuals and run EthSEQ separately for each chunk:
    splits_cfg = int(pipeline.config.get('avadx', 'ethseq.split', fallback=0))
    splits = splits_cfg if splits_cfg > 0 else int(10e12)
    step2_2_1_out = 'sample_list.txt'
    step2_2_1_outfolder = 'splits'
    step2_2_1_splits = 'splits.txt'
    pipeline.add_action(
        'bcftools', 2.21,
        'extract sample list',
        f'-c \\"bcftools query -l $WD/{step1_out} > $WD/{step2_2_1_out}; '
        + f'SPLIT=$(SAMPLES=$(wc -l < $WD/{step2_2_1_out}); echo $((SAMPLES <= {splits} ? 0 : {splits}))); '
        + f'[[ $SPLIT -gt 0 ]] && split -d -l $SPLIT $WD/{step2_2_1_out} $WD/{step2_2_1_outfolder}/xx || cp -f $WD/{step2_2_1_out} $WD/{step2_2_1_outfolder}/xx00; '
        + f'(cd $WD/{step2_2_1_outfolder} && ls -f -1 xx*) > $WD/{step2_2_1_splits}\\"',
        '--entrypoint=bash',
        outdir=(WD / step2_2_1_outfolder)
    )

    # 2.2.2 OPTIONAL: Clean VCF format for EthSEQ input for all splits:
    step2_2_2_splits = 'splits_EthSEQ'
    pipeline.add_action(
        'bcftools', 2.22,
        'EthSEQ preprocessing (VCF cleanup)',
        f'-c \\"bcftools view -S $TASK $WD/{step1_out} | '
        + f'bcftools annotate --remove \\"ID,INFO,FORMAT\\" | '
        + f'bcftools view --no-header -Oz -o $WD/{step2_2_2_splits}/source_$(basename $TASK)_EthSEQinput.vcf.gz\\"',
        '--entrypoint=bash',
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/'),
        outdir=(WD / step2_2_2_splits)
    )

    # 2.2.3 Run EthSEQ:
    #       Note: "export R_MAX_VSIZE=32000000000" can be used
    #       to increase memory before running below for larger datasets
    step2_2_3_outfolder = 'EthSEQ_reports'
    pipeline.add_action(
        'ethnicity_EthSEQ', 2.23,
        'run EthSEQ',
        f'-c \\"mkdir -p $WD/{step2_2_3_outfolder}/$(basename $TASK) && '
        + f'Rscript /app/R/avadx/ethnicity_EthSEQ.R '
        + f'$WD/{step2_2_2_splits}/source_$(basename $TASK)_EthSEQinput.vcf.gz '
        + f'config[DEFAULT.ethseq.models] $WD/{step2_2_3_outfolder}/$(basename $TASK)\\"',
        '--entrypoint=bash --env=R_MAX_VSIZE=32000000000',
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/')
    )

    # 2.2.4 Ethnicity prediction summary
    step2_2_4_outfolder = 'EthSEQ_summary'
    pipeline.add_action(
        'ethnicity_EthSEQ_summary', 2.24,
        'generate EthSEQ summaries',
        f'-c \\"mkdir -p $WD/{step2_2_4_outfolder}/$(basename $TASK) && '
        + f'Rscript /app/R/avadx/ethnicity_EthSEQ_summary.R $WD/{step2_2_3_outfolder}/$(basename $TASK)/Report.txt $WD/{step2_2_1_outfolder}/$(basename $TASK) $WD/{step2_2_4_outfolder}/$(basename $TASK)\\"',
        '--entrypoint=bash',
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/')
    )

    # 2.2.5 Merge ethnicity prediction summaries
    pipeline.add_action(
        'avadx', 2.25,
        'merge EthSEQ summaries',
        f'$WD/{step2_2_4_outfolder} -mindepth 2 -name "sampleID_*" -exec bash -c \\"cat $1 >> $WD/{step2_2_4_outfolder}/$(basename $1)\\" _ {{}} \;',
        '--entrypoint=find'
    )

    # 2.3   Relatedness check:
    #       Check relatedness within datasets withe the SNPRelate R package.
    #       A default kinship > 0.3 is considered to be related.
    step2_3_outfolder = 'SNPRelate'
    step2_3_out = 'source_samp_pass_snps_site-v_gt-v_rmchr_gnomad.gds'
    pipeline.add_action(
        'relatedness', 2.30,
        'check relatedness using SNPRelate',
        f'/app/R/avadx/relatedness.R -i $WD/{step1_out} -g $WD/{step2_3_out} -c config[avadx.kinship] -o $WD/{step2_3_outfolder}',
        outdir=(WD / step2_3_outfolder)
    )

    # 2.4   Remove individual outliers
    #       Outlier individual IDs should be combined from the above PCA,
    #       ethnicity annotation, and relatedness calculation
    #       to a file outliers.txt (one ID per row).
    step2_4_out = 'source_samp_pass_snps_site-v_gt-v_rmchr_gnomad_ind-cleaned.vcf.gz'
    if outliers_available:
        outliers_file = pipeline.config.get('avadx', 'outliers')
        mounts_step2_4 = get_mounts(pipeline, ('avadx', 'outliers'))
        pipeline.add_action(
            'bcftools', 2.40,
            'summarize outliers',
            f'view -S ^{mounts_step2_4[0][1]} $WD/{step1_out} -Oz -o $WD/{step2_4_out}',
            mounts=mounts_step2_4
        )
    step2_out = step2_4_out if outliers_available else step1_out

    # 3     Query/Calculate SNAP scores for all variants ------------------------------------------ #
    
    # 3.1   Get all variant annotations with ANNOVAR for cleaned VCF:
    step3_1_out = 'source_samp_pass_snps_site-v_gt-v_rmchr_gnomad_ind-cleaned.avinput'
    pipeline.add_action(
        'annovar', 3.10,
        'convert VCF file to ANNOVAR input format',
        f'convert2annovar.pl -format vcf4old $WD/{step2_out} -outfile $WD/{step3_1_out}',
        logs=('convert2annovar.log', None)
    )

    # 3.2   Annotate using <hgref> RefSeq:
    step3_2_out = 'source_samp_pass_snps_site-v_gt-v_rmchr_gnomad_ind-cleaned.avinput.exonic_variant_function'
    pipeline.add_action(
        'annovar', 3.20,
        f'annotate using {hgref} RefSeq',
        f'annotate_variation.pl -buildver {hgref} $WD/{step3_1_out} config[DEFAULT.annovar.humandb]/',
        logs=('annotate_variation.log', None)
    )

    # 3.3   Create varidb query
    step3_3_outfolder = 'check_missing_SNAP'
    pipeline.add_action(
        'check_missing_SNAP', 3.30,
        'check missing SNAP predictions',
        f'/app/R/avadx/check_missing_SNAP.R $WD/{step3_2_out} config[DEFAULT.avadx.data] $WD/{step3_3_outfolder}',
        outdir=(WD / step3_3_outfolder)
    )

    # 3.4   Query varidb for SNAP mutations
    step3_4_out = 'varidb_query_result.csv'
    pipeline.add_action(
        'varidb', 3.40,
        'query SNAP variants from varidb',
        f'-D config[DEFAULT.avadx.data]/varidb.db '
        + f'-l $WD/{step3_3_outfolder}/varidb_query.ids '
        + f'-f $WD/{step3_3_outfolder}/varidb_query.fa '
        + f'-o $WD/{step3_4_out} '
        + f'-R $WD/varidb_query_report.txt '
        + f'-C query variant score -S"\t" -H -s'
    )

    # 3.5   Reformat varidb results to Mutations.mutOut file
    # TODO

    # 4     Gene score calculation ---------------------------------------------------------------- #
    
    # 4.1   Convert cleaned VCF to individual ANNOVAR annotation files by:
    step4_1_outfolder = 'annovar_annotations'
    step4_1_out = 'annovar_annotations.txt'
    pipeline.add_action(
        'annovar', 4.10,
        'convert cleaned VCF to single annovar annotation files',
        f'-c \\"convert2annovar.pl -format vcf4 $WD/{step2_out} -outfile $WD/{step4_1_outfolder}/sample -allsample; '
        + f'(cd $WD/{step4_1_outfolder} && ls -f -1 sample.*.avinput) > $WD/{step4_1_out}\\"',
        '--entrypoint=bash',
        outdir=(WD / step4_1_outfolder)
    )

    # 4.2 - Next, annotate all sample*.avinput files:
    pipeline.add_action(
        'annovar', 4.20,
        'generate annovar annotation',
        f'annotate_variation.pl -build {hgref} $TASK config[DEFAULT.annovar.humandb]/',
        tasks=(None, WD / step4_1_out, f'$WD/{step4_1_outfolder}/')
    )

    # 4.3 - Then, calculate gene score:
    step4_3_outfolder = 'gene_scores'
    pipeline.add_action(
        'cal_genescore_make_genescore', 4.30,
        'calculate gene score',
        f'/app/R/avadx/cal_genescore_make_genescore.R -f $TASK.exonic_variant_function '
        + f'-s $WD/{step3_4_out} -l config[avadx.trsprotlengthcleaned] '
        + f'-m config[avadx.gscoremethod] -n config[avadx.normalizeby] -o $WD/{step4_3_outfolder}',
        tasks=(None, WD / step4_1_out, f'$WD/{step4_1_outfolder}/'),
        outdir=(WD / step4_3_outfolder)
    )

    # 4.4 - Merge gene scores:
    #  Assuming there are 500 individuals in the dataset, 500 resulting files will be generated (e.g. sample.S001.gs).
    #  Merge them into a data frame where a row is an individual and a column is a gene (protein):
    step4_4_outfolder = 'genescores'
    pipeline.add_action(
        'merge_genescore', 4.40,
        'merge single gene score result files',
        f'/app/R/avadx/merge_genescore.R -f $WD/{step4_3_outfolder} -o $OUT/{step4_4_outfolder}',
        outdir=(OUT / step4_4_outfolder)
    )

    # 5   Feature selection (FS) and model building --------------------------------------------- #
    #  AVA,Dx by default uses K-S (Kolmogorov–Smirnov) test for FS, random forest for model building,
    #  and 10-fold cross validation to test the predictability of the top-ranking genes.
    #  Other FS methods and machine learning models are also included.
    #  User needs to provide a cross-validation scheme file. For example, we split Tourette dataset (e.g. yale-1)
    #  into 10 folds and the individuals from the same family enter the same fold,
    #  so that to compare sick v.s. healthy instead of differentiating families.

    # 5.1 - Cross-validation:
    step5_1_outfolder = 'crossvalidation'
    pipeline.add_action(
        'FS_CVperf_kfold', 5.10,
        'perform model cross-validation',
        f'/app/R/avadx/FS-CVperf-kfold.R -f $OUT/{step4_4_outfolder}/GeneScoreTable_normed.txt '
        + f'-m config[avadx.cv.featureselection] -M config[avadx.cv.model] -s $WD/cv-scheme.csv '
        + f'-k config[avadx.cv.folds] -l config[avadx.trsprotlengthcleaned] -t config[avadx.cv.steps] '
        + f'-n config[avadx.cv.topgenes] -v config[avadx.cv.varcutoff] -o $WD/{step5_1_outfolder}',
        outdir=(WD / step5_1_outfolder)
    )

    # RUN --------------------------------------------------------------------------------------- #
    main(pipeline, extra)


def init():
    namespace, extra = parse_arguments()
    actions = namespace.action
    config = namespace.config
    daemon = namespace.daemon
    del namespace.action
    del namespace.config
    del namespace.daemon
    if namespace.info:
        Pipeline(actions, kwargs=vars(namespace), config_file=config, daemon=daemon).info()
    elif namespace.preprocess:
        pipeline = Pipeline(actions, kwargs=vars(namespace), config_file=config, daemon=daemon)
        uid = get_extra(extra, '--uid')
        if uid:
            pipeline.uid = uid[0]
        pipeline.preprocess()
    elif actions == None:
        run_all(vars(namespace), extra, config, daemon)
    else:
        pipeline = Pipeline(actions, kwargs=vars(namespace), config_file=config, daemon=daemon)
        main(pipeline, extra)


if __name__ == "__main__":
    init()
