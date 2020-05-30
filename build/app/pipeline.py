#!/usr/bin/env python3

import argparse
import re
import copy
import uuid
import shlex
import subprocess
import configparser
import multiprocessing
from pathlib import Path
from datetime import datetime
from timeit import default_timer as timer
from .logger import Logger
from .helper import run_command
from . import __version__, __releasedate__


QUIET = False
DOCKER_MNT = Path('/mnt')
LOG_LEVEL = 'INFO'
LOG_FILE = None


class AVADxMeta:
    def __init__(self, config=None, daemon:str='docker'):
        self.log = self.get_logger()
        self.daemon = daemon
        self.config = config
        self.REGISTRY = {
            'singularity': {
                'CMD_RUN': 'run'
            },
            'docker': {
                'CMD_RUN': 'run'
            }
        }

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
            self.log.warn(f'No method: {func_name}')

    def run_container(self, container, args=[], uid=uuid.uuid1(), in_folder:Path=Path.cwd(), out_folder:Path=Path.cwd()):
        config = copy.copy(self.config)
        if in_folder:
            in_folder.mkdir(parents=True, exist_ok=True)
        if out_folder:
            (out_folder / str(uid) / 'out').mkdir(parents=True, exist_ok=True)
        wd_folder = out_folder / str(uid) / 'wd'
        wd_folder.mkdir(parents=True, exist_ok=True)
        data_folder = Path(config.get("DEFAULT", "datadir", fallback=wd_folder)).absolute()
        if self.daemon == 'docker':
            in_ = DOCKER_MNT / 'in'
            wd = DOCKER_MNT / 'out' / str(uid) / 'wd'
            out = DOCKER_MNT / 'out' / str(uid) / 'out'
            data = DOCKER_MNT / 'data'
            config.set('DEFAULT', 'datadir', str(data.absolute()))
        else:
            wd = out_folder / str(uid) / 'wd'
            out = out_folder / str(uid) / 'out'

        cmd_base = [
            self.daemon,
            self.REGISTRY[self.daemon]['CMD_RUN']
        ]
        if self.daemon == "docker":
            cmd_base += [
                '--rm',
                '-v', f'{in_folder.absolute()}:{in_}',
                '-v', f'{out_folder.absolute()}:{out}',
                '-v', f'{wd_folder.absolute()}:{wd}',
                '-v', f'{data_folder}:{data}',
            ]
        cmd_base += [
            container,
        ]
        args_parsed = [
            a
             .replace('$WD', str(wd))
             .replace('$IN', str(in_))
             .replace('$OUT', str(out))
            for a in args
        ]
        config_pattern = r'config\[(.*)\]'
        args_parsed = [
            re.sub(
                config_pattern,
                lambda x: config.get(x.group(1).split('.', 1)[0], x.group(1).split('.', 1)[1], fallback=''),
                a
            )
            for a in args_parsed
        ]
        cmd = cmd_base + args_parsed
        self.log.debug(f'Executing: {" ".join(cmd)}')

        run_command(cmd)

    def run_method(self, container, name, uid, kwargs):
        self.log.info(f'Start: {name} [{datetime.now()}]')
        timer_start = timer()
        args = shlex.split(kwargs.get(name).pop(0))
        self.run_container(container, args=args, uid=uid, out_folder=kwargs.get('wd'))
        self.log.info(f'  End: bcftools [{(timer() - timer_start):.3f} seconds]')

    def bfctools(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-bcftools', 'bfctools', uid, kwargs)

    def filterVCF_by_ABAD(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-meta', 'filterVCF_by_ABAD', uid, kwargs)

    def filterVCF_by_gnomAD(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-tabix', 'filterVCF_by_gnomAD', uid, kwargs)

class Pipeline():
    def __init__(self, actions=[], kwargs={}, uid=uuid.uuid1(), config_file=None):
        self.actions = actions
        self.kwargs = kwargs
        self.uid = uid
        self.config = self.load_config(config_file)

    def load_config(self, config_file):
        config = configparser.ConfigParser()
        if config_file and config_file.exists():
            config.read(str(config_file))
        return config

    def add_action(self, action, args=''):
        self.actions += [action]
        self.kwargs[action] = self.kwargs.get(action, []) + [args]


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
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help="number of threads; default: 1")
    parser.add_argument('-w', '--wd', type=Path, default=Path.cwd(),
                        help="working directory")
    parser.add_argument('-c', '--cpu', type=int, default=multiprocessing.cpu_count(),
                        help="max cpus per thread; default: all available cores")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help="set verbosity level; default: log level INFO")
    parser.add_argument('-L', '--logfile', type=Path, const=Path('pipeline.log'),
                        nargs='?', help='redirect logs to file')
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="if flag is set console output is logged to file")

    namespace, extra = parser.parse_known_args()
    parser_actions = []
    for action in namespace.action if namespace.action else []:
        if action not in parser_actions:
            parser.add_argument(f'--{action}',  action='append', default=[])
            parser_actions += [action]
    namespace, extra = parser.parse_known_args()
    namespace.verbose = 40 - (10*namespace.verbose) if namespace.verbose > 0 else 0
    global LOG_LEVEL, LOG_FILE
    LOG_LEVEL = namespace.verbose
    LOG_FILE = namespace.logfile
    del namespace.verbose
    del namespace.logfile

    return namespace, extra


def main(pipeline, extra):
    app = AVADxMeta(pipeline.config)
    if pipeline.kwargs.get('wd'):
        pipeline.kwargs.get('wd').mkdir(parents=True, exist_ok=True)
    timer_start = timer()
    for action in pipeline.actions:
        app.run(action, pipeline.uid, pipeline.kwargs, extra)
    app.log.info(f'Total runtime: {(timer() - timer_start):.3f} seconds')


def run_all(kwargs, extra):
    pipeline = Pipeline(kwargs=kwargs, config_file=kwargs['config'])
    get_extra = lambda x, y: [_.split('=')[1] for _ in x if _.startswith(y)]
    uid = get_extra(extra, '--uid')
    if uid:
        pipeline.uid = uid[0]

    # # 1   Variant QC ---------------------------------------------------------------------------- #
    # # 1.1 Extract individuals of interest (diseased and healthy individuals of interest).
    # step1_1_in = '$IN/test/input_vcf.vcf'
    # step1_1_out = '$WD/source_s-selected.vcf.gz'
    # pipeline.add_action(
    #     'bfctools',
    #     f'view {step1_1_in} -Oz -o {step1_1_out}'
    # )

    # # 1.2 Remove variant sites which did not pass the VQSR standard.
    # step1_2_out = '$WD/source_s-selected_v-PASS.vcf.gz'
    # pipeline.add_action(
    #     'bfctools',
    #     f'filter -i \'FILTER="PASS"\' {step1_1_out} -Oz -o {step1_2_out}'
    # )

    # # 1.3 Split SNV and InDel calls to separated files because they use different QC thresholds.
    # #    Current AVA,Dx works mainly with SNPs. InDels need another set of standards for QC.
    # step1_3_out = '$WD/source_s-selected_v-PASS_snps.vcf.gz'
    # pipeline.add_action(
    #     'bfctools',
    #     f'view --types snps {step1_2_out} -Oz -o {step1_3_out}' # -S sampleID.txt
    # )

    # # 1.4 Remove variant sites by site-wise quality.
    # #  Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150.
    # step1_4_out = '$WD/source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz'
    # pipeline.add_action(
    #     'bfctools',
    #     f'view -i "QUAL>30 & AVG(FMT/DP)<=150 & AVG(FMT/DP)>=6" {step1_3_out} -Oz -o {step1_4_out}'
    # )

    # # 1.5 Check individual call quality. In filterVCF_by_ABAD.py:
    # #  good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4;
    # #  bad individual GTs are converted into missing "./.";
    # #  low call rate is determined as a call rate < 80%,
    # #  i.e. missing rate >= 20%. Variant sites with a low call rate are removed.
    # step1_5_out = '$WD/source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz'
    # pipeline.add_action(
    #     'filterVCF_by_ABAD',
    #     f'avadx.filterVCF_by_ABAD {step1_4_out} {step1_5_out}'
    # )

    # # 1.6 Lastly, gnomAD filter: filtering out variants that were not recorded in the gnomAD database.
    # #  The gnomAD reference used here is the ANNOVAR gnomAD file hg19_gnomad_exome.txt and hg19_gnomad_genome.txt.
    # #  Check the input path of the two reference files before running the script.
    # #  Note that tabix is required for indexing to run this script.
    
    # # 1.6.1 Conver the chromosome annotation if the chromosomes are recorded as "chr1" instead of "1":
    # step1_6_1_out = '$WD/input_rmchr.vcf.gz'
    # pipeline.add_action(
    #     'bfctools',
    #     f'annotate {step1_5_out} -Oz -o {step1_6_1_out}' # --rename-chrs chr_to_number.txt
    # )
    step1_6_1_out = '$WD/input_rmchr.vcf.gz'

    # 1.6.2 Then remove variants that are not in gnomAD database:
    step1_6_2_out = '$WD/output.vcf.gz'
    pipeline.add_action(
        'filterVCF_by_gnomAD',
        f'avadx.filterVCF_by_gnomAD {step1_6_1_out} {step1_6_2_out} config[annovar.humandb]'
    )

    # 2   Individual QC ------------------------------------------------------------------------- #

    # run the pipeline
    main(pipeline, extra)


def init():
    namespace, extra = parse_arguments()
    actions = namespace.action
    del namespace.action
    if actions == None:
        run_all(vars(namespace), extra)
    else:
        pipeline = Pipeline(actions, vars(namespace))
        main(pipeline, extra)


if __name__ == "__main__":
    init()
