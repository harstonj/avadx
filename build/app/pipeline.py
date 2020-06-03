#!/usr/bin/env python3

import argparse
import re
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
            self.log.warn(f'No method: {func_name}')

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
        fns_pre, fns_post =  kwargs.get(name + "_fns").pop(0)
        if outdir:
            outdir.mkdir(parents=True, exist_ok=True)
        for tid, task in enumerate(tasklist, 1):
            task_info = f' [{tid}/{len(tasklist)}]' if task else ''
            self.log.info(f'{name}{task_info}: {kwargs.get(name + "_desc").pop(0)}')
            self.log.debug(f'{name}{task_info}: started {datetime.now()}')
            timer_start = timer()
            if fns_pre:
                fns_pre()
            if task:
                task_idxs = [i for i in range(len(args)) if args[i].find('$TASK') != -1 ]
                for task_idx in task_idxs:
                    update = f'{taskprefix if taskprefix else ""}{task}'
                    args[task_idx] = update if args[task_idx] == '$TASK' else args[task_idx].replace('$TASK', update)
                if taskflag:
                    task_idxs.reverse()
                    for task_idx_r in task_idxs:
                        args.insert(task_idx_r, taskflag)
            self.pipeline.run_container(
                container, args=args, daemon_args=daemon_args, uid=uid,
                in_folder=Path(self.pipeline.config.get('DEFAULT', 'inputdir', fallback=Path.cwd())),
                out_folder=kwargs.get('wd')
            )
            if fns_post:
                fns_post()
            self.log.info(f'{name}{task_info}: took {(timer() - timer_start):.3f} seconds')

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

    def exonic_variant_function2snap_query(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'exonic_variant_function2snap_query', uid, kwargs)

    def check_missing_SNAP(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'check_missing_SNAP', uid, kwargs)

    def cal_genescore_make_genescore(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'cal_genescore_make_genescore', uid, kwargs)

    def FS_CVperf_kfold(self, uid, **kwargs):
        self.run_method('bromberglab/avadx-rscript', 'FS_CVperf_kfold', uid, kwargs)

class Pipeline():

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
        self.config = self.load_config(config_file)
        self.daemon = daemon
        self.init_daemon()

    def get_logger(self):
        logger = Logger(self.__class__.__name__, level=LOG_LEVEL)
        logger.addConsoleHandler()
        log = logger.getLogger()
        return log

    def init_daemon(self):
        if self.daemon == 'docker':
            # import docker
            pass

    def load_config(self, config_file):
        config = configparser.ConfigParser()
        if config_file and config_file.exists():
            config.read(str(config_file))
        return config

    def add_action(self, action, description='', args='', daemon_args='', tasks=(None, None, None), fns=(None, None), outdir=None):
        self.actions += [action]
        self.kwargs[action] = self.kwargs.get(action, []) + [args]
        self.kwargs[f'{action}_desc'] = self.kwargs.get(f'{action}_desc', []) + [description]
        self.kwargs[f'{action}_darg'] = self.kwargs.get(f'{action}_darg', []) + [daemon_args]
        self.kwargs[f'{action}_tasks'] = self.kwargs.get(f'{action}_tasks', []) + [tasks]
        self.kwargs[f'{action}_fns'] = self.kwargs.get(f'{action}_fns', []) + [fns]
        self.kwargs[f'{action}_outd'] = self.kwargs.get(f'{action}_outd', []) + [outdir]

    def check_optional(self, name, flag='', is_file=False):
        formatted = ''
        if self.config.has_option('inputs', f'optional.{name}'):
            config_datadir_orig = self.config.get('DEFAULT', 'datadir', fallback=None)
            config_inputdir_orig = self.config.get('DEFAULT', 'inputdir', fallback=None)
            self.config.set('DEFAULT', 'datadir', str(DOCKER_MNT / 'data'))
            self.config.set('DEFAULT', 'inputdir', str(DOCKER_MNT / 'in'))
            option_value = self.config.get('inputs', f'optional.{name}')
            self.config.set('DEFAULT', 'datadir', config_datadir_orig)
            self.config.set('DEFAULT', 'inputdir', config_inputdir_orig)
            if not is_file or Path(option_value).exists():
                formatted = f'{flag} {option_value}'
            elif is_file and not Path(option_value).exists():
                self.log.warn(f'Ignoring optional argument: {flag} {option_value}. Reason: file not found')
        return formatted

    def run_container(self, container, args=[], daemon_args=[], uid=uuid.uuid1(), in_folder:Path=Path.cwd(), out_folder:Path=Path.cwd()):
        if in_folder:
            in_folder.mkdir(parents=True, exist_ok=True)
        if out_folder:
            (out_folder / str(uid) / 'out').mkdir(parents=True, exist_ok=True)
        wd_folder = out_folder / str(uid) / 'wd'
        wd_folder.mkdir(parents=True, exist_ok=True)
        data_folder = Path(self.config.get("DEFAULT", "datadir", fallback=wd_folder)).absolute()
        config_datadir_orig = self.config.get('DEFAULT', 'datadir')
        self.config.set('DEFAULT', 'datadir', str(data_folder))
        input_folder = Path(self.config.get("DEFAULT", "inputdir", fallback=wd_folder)).absolute()
        config_inputdir_orig = self.config.get('DEFAULT', 'inputdir')
        self.config.set('DEFAULT', 'inputdir', str(input_folder))
        if self.daemon == 'docker':
            in_ = DOCKER_MNT / 'in'
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
            cmd_base += [
                '--rm',
                '-v', f'{in_folder.absolute()}:{in_}',
                '-v', f'{out_folder.absolute()}:{DOCKER_MNT / "out"}',
                '-v', f'{data_folder}:{data}',
            ] + daemon_args
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
        self.config.set('DEFAULT', 'datadir', str(data))
        self.config.set('DEFAULT', 'inputdir', str(in_))
        config_pattern = r'config\[(.*)\]'
        args_parsed = [
            re.sub(
                config_pattern,
                lambda x: self.config.get(x.group(1).split('.', 1)[0], x.group(1).split('.', 1)[1], fallback=''),
                a
            )
            for a in args_parsed
        ]
        self.config.set('DEFAULT', 'datadir', config_datadir_orig)
        self.config.set('DEFAULT', 'inputdir', config_inputdir_orig)
        cmd = cmd_base + args_parsed
        run_command(cmd, logger=self.log)


get_extra = lambda x, y: [_.split('=')[1] for _ in x if _.startswith(y)]


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
                        help="working directory")
    parser.add_argument('-c', '--cpu', type=int, default=multiprocessing.cpu_count(),
                        help="max cpus per thread; default: all available cores")
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help="number of threads; default: 1")
    parser.add_argument('-d', '--daemon', type=str, default='docker', choices=['docker', 'singularity'],
                        help="container engine")
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


def run_all(kwargs, extra, config, daemon):
    pipeline = Pipeline(kwargs=kwargs, config_file=config, daemon=daemon)
    uid = get_extra(extra, '--uid')
    if uid:
        pipeline.uid = uid[0]
    WD = kwargs['wd'] / str(pipeline.uid) / 'wd'
    OUT = kwargs['wd'] / str(pipeline.uid) / 'out'

    # 1   Variant QC ---------------------------------------------------------------------------- #
    # 1.1 Extract individuals of interest (diseased and healthy individuals of interest).
    step1_1_in = 'config[inputs.vcf]'
    step1_1_out = 'source_s-selected.vcf.gz'
    pipeline.add_action(
        'bcftools',
        'filter for individuals of interest ',
        f'view {pipeline.check_optional("sampleids", flag="-S", is_file=True)} {step1_1_in} -Oz -o $WD/{step1_1_out}'
    )

    # 1.2 Remove variant sites which did not pass the VQSR standard.
    step1_2_out = 'source_s-selected_v-PASS.vcf.gz'
    pipeline.add_action(
        'bcftools',
        'filter variant sites < VQSR standard',
        f'filter -i \'FILTER=\\"PASS\\"\' $WD/{step1_1_out} -Oz -o $WD/{step1_2_out}'
    )

    # 1.3 Split SNV and InDel calls to separated files because they use different QC thresholds.
    #    Current AVA,Dx works mainly with SNPs. InDels need another set of standards for QC.
    step1_3_out = 'source_s-selected_v-PASS_snps.vcf.gz'
    pipeline.add_action(
        'bcftools',
        'split SNV and InDel calls',
        f'view --types snps $WD/{step1_2_out} -Oz -o $WD/{step1_3_out}'
    )

    # 1.4 Remove variant sites by site-wise quality.
    #  Good site-wise qualities are: QUAL > 30, mean DP > 6, mean DP < 150.
    step1_4_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150.vcf.gz'
    pipeline.add_action(
        'bcftools',
        'filter variant sites by site-wise quality',
        f'view -i \\"QUAL>30 & AVG(FMT/DP)<=150 & AVG(FMT/DP)>=6\\" $WD/{step1_3_out} -Oz -o $WD/{step1_4_out}'
    )

    # 1.5 Check individual call quality. In filterVCF_by_ABAD.py:
    #  good individual call qualities are: AB > 0.3 and AB < 0.7, GQ > 15, DP > 4;
    #  bad individual GTs are converted into missing "./.";
    #  low call rate is determined as a call rate < 80%,
    #  i.e. missing rate >= 20%. Variant sites with a low call rate are removed.
    step1_5_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.vcf.gz'
    pipeline.add_action(
        'filterVCF_by_ABAD',
        'check individual call quality',
        f'avadx.filterVCF_by_ABAD $WD/{step1_4_out} $WD/{step1_5_out}'
    )

    # 1.6 Lastly, gnomAD filter: filtering out variants that were not recorded in the gnomAD database.
    #  The gnomAD reference used here is the ANNOVAR gnomAD file hg19_gnomad_exome.txt and hg19_gnomad_genome.txt.
    #  Check the input path of the two reference files before running the script.
    #  Note that tabix is required for indexing to run this script.
    
    # 1.6.1 Convert the chromosome annotation if the chromosomes are recorded as "chr1" instead of "1":
    step1_6_1_out = 'input_rmchr.vcf.gz'
    pipeline.add_action(
        'bcftools',
        'convert the chromosome annotations',
        f'annotate --rename-chrs config[avadx.chr2number] $WD/{step1_5_out} -Oz -o $WD/{step1_6_1_out}'
    )

    # 1.6.2 Then remove variants that are not in gnomAD database:
    step1_6_2_out = 'output.vcf.gz'
    pipeline.add_action(
        'filterVCF_by_gnomAD',
        'filter variants missing in gnomAD database',
        f'avadx.filterVCF_by_gnomAD $WD/{step1_6_1_out} $WD/{step1_6_2_out} config[DEFAULT.annovar.humandb]'
    )

    # 2   Individual QC ------------------------------------------------------------------------- #
    # 2.2 Quality check - Check quality outliers by examine nRefHom, nNonRefHom, nHets, nTransitions, nTransversions, average depth, nSingletons, and nMissing:
    # 2.1.0 Output quality metrics after variant QC:
    step2_1_0_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.stats.txt'
    pipeline.add_action(
        'bcftools',
        'generate quality metrics after variant QC',
        f'-c \\"bcftools stats -v -s - $WD/{step1_5_out} > $WD/{step2_1_0_out}\\"',
        '--entrypoint=bash'
    )

    # TODO: broken
    # 2.1.1 Draw individual quality figure:
    pipeline.add_action(
        'stats_quality_pca',
        'draw individual quality figures',
        f'/app/R/avadx/stats_quality_pca.R -f $WD/{step2_1_0_out}'
    )

    # 2.2 Ethnicity check - Annotate ethnicity with EthSEQ R package:
    # 2.2.1 - OPTIONAL: If the number of individuals exceeds certain number, "memory exhausted" error may occur.
    #  Manually divide input VCF into chunks of individuals and run EthSEQ separately for each chunk:
    splits_cfg = int(pipeline.config.get('avadx', 'optional.split', fallback=0))
    splits = splits_cfg if splits_cfg > 0 else int(10e12)
    step2_2_1_out = 'sample_list.txt'
    step2_2_1_outfolder = 'splits'
    step2_2_1_splits = 'splits.txt'
    pipeline.add_action(
        'bcftools',
        'extract sample list',
        f'-c \\"bcftools query -l $WD/{step1_5_out} > $WD/{step2_2_1_out}; '
        + f'SPLIT=$(SAMPLES=$(wc -l < $WD/{step2_2_1_out}); echo $((SAMPLES <= {splits} ? 0 : {splits}))); '
        + f'[[ $SPLIT -gt 0 ]] && split -d -l $SPLIT $WD/{step2_2_1_out} $WD/{step2_2_1_outfolder}/xx || cp -f $WD/{step2_2_1_out} $WD/{step2_2_1_outfolder}/xx00; '
        + f'(cd $WD/{step2_2_1_outfolder} && ls -f -1 xx*) > $WD/{step2_2_1_splits}\\"',
        '--entrypoint=bash',
        outdir=(WD / step2_2_1_outfolder)
    )

    # 2.2.2 - OPTIONAL: Clean VCF format for EthSEQ input for all splits:
    step2_2_2_splits = 'splits_EthSEQ'
    pipeline.add_action(
        'bcftools',
        'EthSEQ preprocessing (VCF cleanup)',
        f'-c \\"bcftools view -S $TASK $WD/{step1_5_out} | '
        + f'bcftools annotate --remove \\"ID,INFO,FORMAT\\" | '
        + f'bcftools view --no-header -Oz -o $WD/{step2_2_2_splits}/source_$(basename $TASK)_EthSEQinput.vcf.gz\\"',
        '--entrypoint=bash',
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/'),
        outdir=(WD / step2_2_2_splits)
    )

    # 2.2.3 - Run EthSEQ:
    # "export R_MAX_VSIZE=32000000000" can be used to increase memory before running below for larger datasets
    step2_2_3_outfolder = 'EthSEQ_reports'
    pipeline.add_action(
        'ethnicity_EthSEQ',
        'run EthSEQ',
        f'-c \\"mkdir -p $WD/{step2_2_3_outfolder}/$(basename $TASK) && '
        + f'Rscript /app/R/avadx/ethnicity_EthSEQ.R $WD/{step2_2_2_splits}/source_$(basename $TASK)_EthSEQinput.vcf.gz $WD/{step2_2_3_outfolder}/$(basename $TASK)\\"',
        '--entrypoint=bash --env=R_MAX_VSIZE=32000000000',
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/')
    )

    # 2.2.4 - Ethnicity prediction summary
    step2_2_4_outfolder = 'EthSEQ_summary'
    pipeline.add_action(
        'ethnicity_EthSEQ_summary',
        'generate EthSEQ summaries',
        f'-c \\"mkdir -p $WD/{step2_2_4_outfolder}/$(basename $TASK) && '
        + f'Rscript /app/R/avadx/ethnicity_EthSEQ_summary.R $WD/{step2_2_3_outfolder}/$(basename $TASK)/Report.txt $WD/{step2_2_1_out} $WD/{step2_2_4_outfolder}/$(basename $TASK)\\"',
        '--entrypoint=bash',
        tasks=(None, WD / step2_2_1_splits, f'$WD/{step2_2_1_outfolder}/')
    )

    # 2.2.4 - Merge ethnicity prediction summaries
    pipeline.add_action(
        'avadx',
        'merge EthSEQ summaries',
        f'$WD/{step2_2_4_outfolder} -mindepth 2 -name "sampleID_*" -exec bash -c \\"cat $1 >> $WD/{step2_2_4_outfolder}/$(basename $1)\\" _ {{}} \;',
        '--entrypoint=find'
    )

    # 2.3 - Relatedness check:
    #  Check relatedness within datasets withe the SNPRelate R package.
    #  A default kinship > 0.3 is considered to be related.
    step2_3_outfolder = 'SNPRelate'
    step2_3_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc.gds'
    pipeline.add_action(
        'relatedness',
        'check relatedness using SNPRelate',
        f'/app/R/avadx/relatedness.R -i $WD/{step1_5_out} -g $WD/{step2_3_out} -c 0.3 -o $WD/{step2_3_outfolder}',
        outdir=(WD / step2_3_outfolder)
    )

    # TODO: manual intervention?
    # 2.4 - Remove individual outliers
    # 2.4.1 Outlier individual IDs should be combined from the above PCA, ethnicity annotation, and relatedness calculation
    #  to a file outliers.txt (one ID per row).
    step2_4_1_out = 'outliers.txt'
    step2_4_1_fnpre = lambda: (WD / step2_4_1_out).touch()

    # 2.4.2 Then remove individual outliers by:
    step2_4_2_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.vcf.gz'
    pipeline.add_action(
        'bcftools',
        'summarize outliers',
        f'view -S ^$WD/{step2_4_1_out} $WD/{step1_5_out} -Oz -o $WD/{step2_4_2_out}',
        fns=(step2_4_1_fnpre, None)
    )

    # 3   Query/Calculate SNAP scores for all variants ------------------------------------------ #
    # 3.1 Get all variant annotations with ANNOVAR for cleaned VCF:
    step3_1_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.avinput'
    pipeline.add_action(
        'annovar',
        'convert VCF file to ANNOVAR input format',
        f'convert2annovar.pl -format vcf4old $WD/{step2_4_2_out} -outfile $WD/{step3_1_out}'
    )

    # 3.2.1 Retrieve gnomad_exome database
    pipeline.add_action(
        'annovar',
        'download database: hg19_gnomad_exome',
        f'-c \\"[[ -f config[DEFAULT.annovar.humandb]/hg19_gnomad_exome.txt ]] || annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad_exome config[DEFAULT.annovar.humandb]/\\"',
        '--entrypoint=bash'
    )

    # 3.2.2 Retrieve refGene database
    pipeline.add_action(
        'annovar',
        'download database: hg19_refGene',
        f'-c \\"[[ -f config[DEFAULT.annovar.humandb]/hg19_refGene.txt ]] || annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene config[DEFAULT.annovar.humandb]/\\"',
        '--entrypoint=bash'
    )

    # 3.2.3 Annotate using hg19 RefSeq:
    step3_2_3_out = 'source_s-selected_v-PASS_snps_site-v-Q30-minavgDP6-maxavgDP150_gt-v-DP4-AB37-GQ15-MR20perc_ind-cleaned.avinput.exonic_variant_function'
    pipeline.add_action(
        'annovar',
        'annotate using hg19 RefSeq',
        f'annotate_variation.pl -buildver hg19 $WD/{step3_1_out} config[DEFAULT.annovar.humandb]/'
    )

    # 3.3 Make snapfun.db query file for missing SNAPs
    # Then, extract all variants from *.exonic_variant_function to query snap scores from snapfun.db.
    # The db folder already contains a file (Mutations.mutOut) of pre-calculated SNAP scores for variants from previous studies.
    # Below steps will generate query file for variants which are not included in Mutations.mu
    step3_3_out = 'query_result.txt'
    pipeline.add_action(
        'exonic_variant_function2snap_query',
        'extract variants to query snapfun.db',
        f'/app/R/avadx/exonic_variant_function2snap_query.R -f $WD/{step3_2_3_out} -m config[avadx.refseq2uniprot] -M config[avadx.mutations] -o $WD/{step3_3_out}'
    )

    # 3.4 Make SNAP input files for missing SNAPs
    step3_4_outfolder = 'check_missing_SNAP'
    pipeline.add_action(
        'check_missing_SNAP',
        'check missing SNAP predictions',
        f'/app/R/avadx/check_missing_SNAP.R $WD/{step3_2_3_out} config[DEFAULT.avadx.data] $WD/{step3_4_outfolder}',
        outdir=(WD / step3_4_outfolder)
    )

    # 4   Gene score calculation ---------------------------------------------------------------- #
    # 4.1 - Convert cleaned VCF to individual ANNOVAR annotation files by:
    step4_1_outfolder = 'annovar_annotations'
    step4_1_out = 'annovar_annotations.txt'
    pipeline.add_action(
        'annovar',
        'convert cleaned VCF to single annovar annotation files',
        f'-c \\"convert2annovar.pl -format vcf4 $WD/{step2_4_2_out} -outfile $WD/{step4_1_outfolder}/sample -allsample; '
        + f'(cd $WD/{step4_1_outfolder} && ls -f -1 sample.*.avinput) > $WD/{step4_1_out}\\"',
        '--entrypoint=bash',
        outdir=(WD / step4_1_outfolder)
    )

    # 4.2 - Next, annotate all sample*.avinput files:
    pipeline.add_action(
        'annovar',
        'generate annovar annotation',
        f'annotate_variation.pl -build hg19 $TASK config[DEFAULT.annovar.humandb]/',
        tasks=(None, WD / step4_1_out, f'$WD/{step4_1_outfolder}/')
    )

    # 4.3 - Then, calculate gene score:
    step4_3_outfolder = 'gene_scores'
    pipeline.add_action(
        'cal_genescore_make_genescore',
        'calculate gene score',
        f'/app/R/avadx/cal_genescore_make_genescore.R -f $TASK.exonic_variant_function '
        + f'-s config[avadx.mutations] -l config[avadx.trsprotlengthcleaned] '
        + f'-m config[avadx.gscoremethod] -n config[avadx.normalizeby] -o $WD/{step4_3_outfolder}',
        tasks=(None, WD / step4_1_out, f'$WD/{step4_1_outfolder}/'),
        outdir=(WD / step4_3_outfolder)
    )

    # 4.4 - Merge gene scores:
    #  Assuming there are 500 individuals in the dataset, 500 resulting files will be generated (e.g. sample.S001.gs).
    #  Merge them into a data frame where a row is an individual and a column is a gene (protein):
    step4_4_outfolder = 'genescores'
    pipeline.add_action(
        'exonic_variant_function2snap_query',
        'extract variants to query snapfun.db',
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
        'FS_CVperf_kfold',
        'perform model cross-validation',
        f'/app/R/avadx/FS-CVperf-kfold.R -f $OUT/{step4_4_outfolder}/GeneScoreTable_normed.txt '
        + f'-m config[avadx.cv.featureselection] -M config[avadx.cv.model] -s config[avadx.cv.scheme] '
        + f'-k config[avadx.cv.folds] -l config[avadx.trsprotlengthcleaned] -t config[avadx.cv.steps] '
        + f'-n config[avadx.cv.topgenes] -v config[avadx.cv.varcutoff] -o $WD/{step5_1_outfolder}',
        outdir=(WD / step5_1_outfolder)
    )

    # ------------------------------------------------------------------------------------------- #
    # run the pipeline
    main(pipeline, extra)


def init():
    namespace, extra = parse_arguments()
    actions = namespace.action
    config = namespace.config
    daemon = namespace.daemon
    del namespace.action
    del namespace.config
    del namespace.daemon
    if actions == None:
        run_all(vars(namespace), extra, config, daemon)
    else:
        pipeline = Pipeline(actions, kwargs=vars(namespace), config=config, daemon=daemon)
        main(pipeline, extra)


if __name__ == "__main__":
    init()
