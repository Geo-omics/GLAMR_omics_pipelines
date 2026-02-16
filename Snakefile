from collections import Counter
from contextlib import ExitStack
import json
import os
import re
import subprocess
import snakemake.io
from glob import glob
import pandas as pd
from tempfile import NamedTemporaryFile
import time

import pypelib.amplicon.dispatch
import pypelib.amplicon.guess_target
import pypelib.amplicon.hmm_check_asvs
import pypelib.amplicon.hmm_summarize
import pypelib.amplicon.remove_primers
import pypelib.amplicon.tabulate_targets
from pypelib.post import post_production
import pypelib.raw_reads
from pypelib.raw_reads import parse_runinfo
import pypelib.sra
from pypelib.utils import load_stats, logme, save_error_file


configfile: "config.yaml"
config.setdefault('binning_bam_dir', '/nosuchdirectory')

report: "code/report/workflow.rst"

shell.prefix('printf "Job executed on: ${{HOSTNAME}}\n" && printf "SLURM job id: ${{SLURM_JOB_ID}}\n\n"; ')

current_dir = os.getcwd()
#humann_ref_dir = "/home/kiledal/scratch_gdick1/GVHD/data/reference/humann" # for running on Great Lakes
humann_ref_dir = "/geomicro/data2/kiledal/projects/GVHD/data/reference/humann"

# Set which rules can be run on cluster head node
localrules: make_rulegraph, link_reads_w_sample_names

# Post-production: no-op unless checkout_file and/or version_file are configured

onerror: post_production(log[0], config, rules)
onsuccess: post_production(log[0], config, rules)

rule make_rulegraph:
    output:
        "rulegraph.pdf",
        "rulegraph.png"
    shell:
        """
        snakemake metaG_annotation --rulegraph --dry-run | dot -Tpdf > rulegraph.pdf
        snakemake metaG_annotation --rulegraph --dry-run | dot -Tpng > rulegraph.png
        """

rule make_rulegraph_bins:
    output:
        pdf = "rulegraph_bins.pdf",
        png = "rulegraph_bins.png"
    shell:
        """
        snakemake run_drep_sep metaG_annotation run_humann_fastp run_sourmash annotate_bins data/sample_data/bracken_counts.tsv --rulegraph --dry-run | dot -Tpdf > {output.pdf}
        snakemake run_drep_sep metaG_annotation run_humann_fastp run_sourmash annotate_bins data/sample_data/bracken_counts.tsv --rulegraph --dry-run | dot -Tpng > {output.png}
        """


def parse_accession(file):
    """ parse_input function for rule get_reads """
    file.seek(0)  # we're called multiple times
    if value := file.readline().strip():
        return value
    else:
        raise RuntimeError(f'bad accession file? {value=}')


rule get_sra_metadata:
    input: accession = "data/omics/{sample_type}/{sample}/reads/accession"
    output: "data/omics/{sample_type}/{sample}/reads/sra_metadata.xml"
    params:
        accn = parse_input(input.accession, parse_accession),
        ncbi_api_key = os.environ.get('NCBI_API_KEY', '')
    log: "logs/get_sra_metadata/{sample_type}.{sample}.err"
    resources: time_min = 5, heavy_network = 1, cpus = 1
    run:
        with save_error_file(log[0]):
            with pypelib.sra.get_entry_raw(params.accn, multi=True, slow=True) as data:
                with open(output[0], 'wb') as ofile:
                    ofile.write(data.read())


rule get_reads_prep:
    input:
        sra_metadata = lambda wc: rules.get_sra_metadata.output if exists(rules.get_sra_metadata.input.accession.format(**wc)) else []
    output:
        runinfo0 = "data/omics/{sample_type}/{sample}/reads/runinfo0.json",
        runinfo = "data/omics/{sample_type}/{sample}/reads/runinfo.json"
    params:
        ncbi_api_key = os.environ.get('NCBI_API_KEY', '')
    conda: "config/conda_yaml/kingfisher.yaml"
    log: "logs/get_reads/{sample_type}-{sample}-prep.log"
    resources: time_min = 5, heavy_network = 1, cpus = 1
    run:
        with logme(log):
            if input.sra_metadata:
                # SRA dataset
                data = pypelib.sra.get_entry(file=str(input.sra_metadata), multi=True, slow=True)
                with open(rules.get_sra_metadata.input.accession.format(**wildcards)) as ifile:
                    accn = parse_accession(ifile)
                with save_error_file(log[0]):
                    expack = pypelib.sra.get_experiment(accn, data, sample_type=wildcards.sample_type)

                info = pypelib.sra.compile_srr_info(expack)
                with open(output.runinfo0, 'w') as ofile:
                    json.dump(info, ofile, indent=4)
                    ofile.write('\n')
                srr_accn = expack['RUN_SET']['RUN']['accession']
                accn_str = accn + ' => ' + srr_accn

                print(f'Accession for {wildcards.sample}: {accn_str}')
                shell("""
                    . code/shell_prelude {log}

                    [[ -n {params.ncbi_api_key:q} ]] && export NCBI_API_KEY={params.ncbi_api_key:q}
                    [[ -v NCBI_API_KEY ]] || echo "[NOTICE] environment variable NCBI_API_KEY is not set"

                    [[ -n "${{KINGFISHER_SLEEP:-}}" ]] && sleep $((RANDOM % 30))
                    kingfisher annotate -r {srr_accn:q} -a -f json -o {output.runinfo}
                """)
                with open(output.runinfo) as ifile:
                    runinfo = parse_runinfo(ifile)
                if runinfo['number_of_runs_for_sample'] > 1:
                    raise RuntimeError(
                        f'Multiple runs per sample! {wildcards=} {accn=}\n'
                        f'{runinfo=}'
                    )
            else:
                # non-SRA dataset

                # Compile runinfo from existing raw reads files.  We will
                # meaningfully populate "layout" and "spots".  There will be
                # placeholders (None) for keys "srr_accession" and "run" because
                # some rules' params will asked for them.  NOTE: Care must be taken
                # that these None values are not used unwittingly in those rules.
                fwd_reads = rules.get_reads_paired.output.fwd_reads.format(**wildcards)
                rev_reads = rules.get_reads_paired.output.rev_reads.format(**wildcards)
                sgl_reads = rules.get_reads_single.output.single_reads.format(**wildcards)
                if exists(fwd_reads) and exists(rev_reads):
                    layout = 'PAIRED'
                    check_file = fwd_reads
                elif exists(sgl_reads):
                    layout = 'SINGLE'
                    check_file = sgl_reads
                else:
                    reads_dir = Path(fwd_reads).parent
                    raise RuntimeError(
                        f'Non-SRA data / manual mode: raw reads file(s) not found '
                        f'in {reads_dir}'
                    )

                with NamedTemporaryFile('rt') as statsfile:
                    pypelib.raw_reads.make_stats(check_file, statsfile.name)
                    stats = load_stats(statsfile)

                info = dict(
                    spots=stats[Path(check_file).name]['num_seqs'],
                    run=None,
                )

                with open(output.runinfo0, 'w') as ofile:
                    info['layout'] = layout
                    json.dump(info, ofile, indent=4)
                    ofile.write('\n')
                with open(output.runinfo, 'w') as ofile:
                    info['library_layout'] = info.pop('layout')
                    # dump list with one dict, simulating kingfisher annotate
                    json.dump([info], ofile, indent=4)
                    ofile.write('\n')

rule summarize_sra_errors:
    output: 'sra_error_summary.txt'
    run:
        unique_pref = re.compile(r'^SR[^:]+: ')  # rm accn to allow grouping
        stats = Counter()
        for err in Path('data/omics/').glob('*/samp_*/reads/sra_error.json'):
            # 1. check if error file is stale
            runinfos = [
                err.parent / Path(info_file).name
                for info_file in rules.get_reads_prep.output
            ]
            if all(i.is_file() for i in runinfos):
                err_mtime = err.stat().st_mtime
                if all(i.stat().st_mtime > err_mtime for i in runinfos):
                    print(f'[Outdated] {err}')

            with open(err) as ifile:
                data = json.load(ifile)
            name = data['Exception']
            try:
                mesg = data['args'][0]
            except IndexError:
                megs = ''  # no args
            mesg = unique_pref.sub('', mesg)
            stats[(name, mesg)] += 1
        rows = [
            f'{count:>5} {name:<25} {msg.strip()}'
            for (name, msg), count in stats.most_common()
        ]
        print(*rows[:20], sep='\n')
        with open(output[0], 'w') as ofile:
            for row in rows:
                ofile.write(f'{row}\n')

rule about_multipack_sra_errors:
    output: 'sra_multipack_errors.json'
    run:
        out = {}
        for errfile in Path('data/omics/').glob('*/samp_*/reads/sra_error.json'):
            with open(errfile) as ifile:
                sample = errfile.parent.parent.name
                errinfo = json.load(ifile)
                for i in errinfo['args']:
                    if 'multiple matching experiment packages' in i:
                        break
                else:
                    # is something else
                    continue
                epacks = errinfo['args'][-1]['EXPERIMENT_PACKAGE_list']
                studies = [i['STUDY']['accession'] for i in epacks]
                if len(set(studies)) != 1:
                    out[sample] = f'[EE] multiple studies: {studies}'
                    continue
                samples = [i['SAMPLE']['accession'] for i in epacks]
                if len(set(samples)) != 1:
                    out[sample] = f'[EE] multiple samples: {samples}'
                    continue
                data = []
                for i in epacks:
                    data.append({
                        'experiment': [
                            i['EXPERIMENT']['accession'],
                            i['EXPERIMENT']['TITLE'],
                        ],
                        'design_description': i['EXPERIMENT']['DESIGN']['DESIGN_DESCRIPTION'],
                        'library': [
                            i['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY'],
                            i['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_SOURCE'],
                            i['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_SELECTION'],
                            i['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_LAYOUT'],
                        ],
                        'run': i['RUN_SET']['RUN']['accession'],
                    })
            out[sample] = data
        with open(output[0], 'w') as ofile:
            json.dump(out, ofile, indent=4)


rule get_reads:
    input:
        runinfo = ancient(rules.get_reads_prep.output.runinfo)
    output:
        download_dir = temp(directory("data/omics/{sample_type}/{sample}/reads/tmp_downloads"))
    params:
        read_dir = subpath(input.runinfo, parent=True),
        srr_accn = parse_input(input.runinfo, parse_runinfo, key='run'),
        ncbi_api_key = os.environ.get('NCBI_API_KEY', '')
    conda: "config/conda_yaml/kingfisher.yaml"
    log: "logs/get_reads/{sample_type}-{sample}-download.log"
    resources: time_min = 5000, heavy_network = 1, cpus = 8
    shell:
        """
        . code/shell_prelude {log}

        [[ -n {params.ncbi_api_key:q} ]] && export NCBI_API_KEY={params.ncbi_api_key:q}
        [[ -v NCBI_API_KEY ]] || echo "[NOTICE] environment variable NCBI_API_KEY is not set"

        [[ -n "${{KINGFISHER_SLEEP:-}}" ]] && sleep $((RANDOM % 30))

        # Adding methods in same order as in kingfisher documentation:
        command -v prefetch && have_sra_toolkit=true || have_sra_toolkit=false
        meths=()
        command -v ascp && meths=(ena-ascp) || true
        $have_sra_toolkit && meths=("${{meths[@]}}" ena-ftp prefetch) || true
        command -v aria2c && $have_sra_toolkit && meths=("${{meths[@]}}" aws-http) || true
        $have_sra_toolkit && meths=("${{meths[@]}}" aws-cp) || true
        echo "Available methods: ${{meths[*]}}"

        echo "Using accession {params.srr_accn:q}"
        kingfisher get \
            --download-threads {resources.cpus} \
            --extraction-threads {resources.cpus} \
            --hide-download-progress \
            -r {params.srr_accn:q} \
            -m ${{meths[@]}}  \
            --output-format-possibilities fastq.gz \
            --output-directory {output.download_dir} \
            --check-md5sums
        """

checkpoint get_runinfo:
    input: rules.get_reads_prep.output.runinfo
    output: "data/omics/{sample_type}/{sample}/reads/runinfo1.tsv"
    resources: time_min = 1, cpus = 1
    shell: "ln -- {input} {output}"

def get_raw_reads_files(wc, format=False):
    """
    Input function to get the raw reads files

    Wildcards need to contain sample_type and sample
    """
    runinfo = checkpoints.get_runinfo.get(**wc).output[0]
    match parse_runinfo(runinfo, 'library_layout'):
        case 'PAIRED':
            ret = rules.get_reads_paired.output
        case 'SINGLE':
            ret = rules.get_reads_single.output
        case other:
            raise ValueError(f'unsupported library layout: {other}')
    if format:
        return [i.format(**wc) for i in ret]
    else:
        return ret

rule raw_reads_stats:
    """
    A stand-in rule, is mostly here to define the name of the stats file.

    The normal way to produce this file via the rules get_reads_paired or
    get_reads_single rules.
    """
    input: get_raw_reads_files
    output:
        stats = update("data/omics/{sample_type}/{sample}/reads/stats.tsv"),
    conda: "config/conda_yaml/seqkit.yaml"
    resources: time_min = 1, cpus = 1
    run: pypelib.raw_reads.make_stats(input, output.stats, keep_existing=True)

def get_download_dir(wc):
    """
    Input function that may return empty or the temp downloads directory, which
    is the output for get_reads, depending on whether raw reads need to be
    downloaded.  If the raw reads exist they will get checked against the stats
    file.  If no stats file exists it will be created from the existing raw reads
    files (and the test should obviously succeed then.)  A failed test will result
    in a new download.
    """
    raw_reads = partial(get_raw_reads_files, format=True)(wc)
    for i in raw_reads:
        if not exists(i):
            print(f'Raw reads file missing: {i}')
            break
    else:
        # all raw reads exist
        stats_file = rules.raw_reads_stats.output.stats.format(**wc)
        pypelib.raw_reads.make_stats(raw_reads, stats_file, keep_existing=True)
        runinfo = checkpoints.get_runinfo.get(**wc).output[0]
        num_spots = parse_runinfo(runinfo, key='spots')
        try:
            pypelib.raw_reads.check(stats=stats_file, num_spots=num_spots)
        except Exception as e:
            print(f'Raw reads check failed: {e}')
        else:
            # do not download again
            return []
    # download (again)
    return rules.get_reads.output.download_dir.format(**wc)

rule get_reads_paired:
    input:
        runinfo = rules.get_reads_prep.output.runinfo,
        download_dir = get_download_dir,
    output:
        fwd_reads = update("data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz"),
        rev_reads = update("data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz"),
    params:
        reads_dir = subpath(input.runinfo, parent=True),
        srr_accession = parse_input(input.runinfo, parse_runinfo, key='run'),
        num_spots = parse_input(input.runinfo, parse_runinfo, key='spots'),
        stats_file = rules.raw_reads_stats.output.stats,
    conda: "config/conda_yaml/seqkit.yaml"
    log: "logs/get_reads/{sample_type}-{sample}-paired.log"
    resources: time_min = 5, heavy_network = 0, cpus = 1
    run: pypelib.raw_reads.post_download_paired(input, output, params, log=log)


rule get_reads_single:
    input:
        runinfo = rules.get_reads_prep.output.runinfo,
        download_dir = get_download_dir,
    output:
        single_reads = update("data/omics/{sample_type}/{sample}/reads/raw_single_reads.fastq.gz"),
    params:
        reads_dir = subpath(input.runinfo, parent=True),
        srr_accession = parse_input(input.runinfo, parse_runinfo, key='run'),
        num_spots = parse_input(input.runinfo, parse_runinfo, key='spots'),
        stats_file = rules.raw_reads_stats.output.stats,
    conda: "config/conda_yaml/seqkit.yaml"
    log: "logs/get_reads/{sample_type}-{sample}-single.log"
    resources: time_min = 5, heavy_network = 0, cpus = 1
    run: pypelib.raw_reads.post_download_single(input, output, params, log=log)


rule clumpify:
    input: 
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz"
    output: 
        done = touch("data/omics/{sample_type}/{sample}/reads/done.touch")
    params:
        clumped_fwd_reads = "data/omics/{sample_type}/{sample}/reads/clumped_raw_fwd_reads.fastq.gz",
        clumped_rev_reads = "data/omics/{sample_type}/{sample}/reads/clumped_raw_rev_reads.fastq.gz",
    conda: "config/conda_yaml/main.yaml"
    benchmark:
        "benchmarks/clumpify/{sample_type}-{sample}.txt"
    log: "logs/clumpify/{sample_type}-{sample}_initial.log"
    resources: cpus=16, time_min=2880, mem_mb = lambda wildcards, attempt: attempt * 175000 # standard
    #resources: partition = "largemem", cpus=16, time_min=2880, mem_mb = 500000 # coassembly or large samples
    shell:
        """
        #BBtools use more memory than given, reduce amount given by 20% to stay within job specs.
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)

        clumpify.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.fwd_reads} \
            in2={input.rev_reads} \
            out1={params.clumped_fwd_reads} \
            out2={params.clumped_rev_reads} \
            groups=24 \
            zl=9 pigz \
            t={resources.cpus} \
            2>&1 | tee {log} &&

        rm {input.fwd_reads} {input.rev_reads} &&
        mv {params.clumped_fwd_reads} {input.fwd_reads} &&
        mv {params.clumped_rev_reads} {input.rev_reads}
        """


rule deduplicate:
    input: 
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        clumpify = "data/omics/{sample_type}/{sample}/reads/done.touch"
    output: 
        dedup_interleaved = temp("data/omics/{sample_type}/{sample}/reads/dedup_interleaved.fastq.gz"),
        dedup_reads_fwd = "data/omics/{sample_type}/{sample}/reads/dedup_reads_fwd.fastq.gz",
        dedup_reads_rev = "data/omics/{sample_type}/{sample}/reads/dedup_reads_rev.fastq.gz"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/dedup/{sample_type}-{sample}_dedup.log"
    # resources: cpus=36, time_min=2880,
    #     mem_mb = lambda wildcards, attempt: attempt * 170000,
    #     #partition = "largemem"
    resources: 
        partition = "largemem",
        cpus = 24, 
        time_min = 7200,
        mem_mb = 1000000
    shell:
        """
        #BBtools use more memory than given, amount given by 20% to stay within job specs.
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.6)/1" | bc)

        dedupe.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.fwd_reads} \
            in2={input.rev_reads} \
            out={output.dedup_interleaved} \
            t={resources.cpus} \
            2>&1 | tee {log}

        # Dedup only outputs interleaved files, this just converts back to paired
        reformat.sh in={output.dedup_interleaved} \
            out1={output.dedup_reads_fwd} \
            out2={output.dedup_reads_rev} 2>&1 | tee -a {log}
        """


# rule de_interleave:
#     input: 
#         reads = "import/staging/jgi_2022/all_sample_filtered_reads/{sample}_interleaved.fastq.gz"
#     output: 
#         reads_fwd = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
#         reads_rev = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz"
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/de_interleave/{sample_type}-{sample}.log"
#     # resources: cpus=36, time_min=2880,
#     #     mem_mb = lambda wildcards, attempt: attempt * 170000,
#     #     #partition = "largemem"
#     resources: 
#         #partition = "largemem",
#         cpus = 24, 
#         time_min = 7200,
#         mem_mb = 120000
#     shell:
#         """
#         #BBtools use more memory than given, amount given by 20% to stay within job specs.
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.6)/1" | bc)

#         # Dedup only outputs interleaved files, this just converts back to paired
#         reformat.sh in={input.reads} \
#             out1={output.reads_fwd} \
#             out2={output.reads_rev} 2>&1 | tee -a {log}
#         """

#ruleorder: remove_contaminants_fastp > de_interleave


rule fastp:
    input: 
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        clumped = "data/omics/{sample_type}/{sample}/reads/done.touch"
    output: 
        tmp_fwd = temp("data/omics/{sample_type}/{sample}/reads/temp_fastp_fwd_reads.fastq.gz"),
        tmp_rev = temp("data/omics/{sample_type}/{sample}/reads/temp_fastp_rev_reads.fastq.gz"),
        fwd_reads = temp("data/omics/{sample_type}/{sample}/reads/fastp_fwd_reads.fastq.gz"),
        rev_reads = temp("data/omics/{sample_type}/{sample}/reads/fastp_rev_reads.fastq.gz"),
        html_dedup = "data/omics/{sample_type}/{sample}/reads/qc/fastp_dedup.html",
        json_dedup = "data/omics/{sample_type}/{sample}/reads/qc/fastp_dedup.json",
        html = "data/omics/{sample_type}/{sample}/reads/qc/fastp.html",
        json = "data/omics/{sample_type}/{sample}/reads/qc/fastp.json"
    conda: "config/conda_yaml/fastp.yaml"
    benchmark:
        "benchmarks/fastp/{sample_type}-{sample}.txt"
    log: "logs/fastp/{sample_type}-{sample}_fastp.log"
    resources: cpus=16, mem_mb = lambda wildcards, attempt: attempt * 60000,
        time_min=2880
    shell:
        """
        # First deduplicate
        fastp \
            -i {input.fwd_reads} -I {input.rev_reads} \
            -o {output.tmp_fwd} -O {output.tmp_rev} \
            -h {output.html_dedup} -j {output.json_dedup} \
            --thread {resources.cpus} \
            -z 3 \
            --dedup \
            --dup_calc_accuracy 6 2>&1 | tee {log}
        
        # Then trim and remove adapters
        fastp \
            -i {output.tmp_fwd} -I {output.tmp_rev} \
            -o {output.fwd_reads} -O {output.rev_reads} \
            -h {output.html} -j {output.json} \
            --thread {resources.cpus} \
            -z 9 \
            --length_required 50 \
            --n_base_limit 5 \
            --low_complexity_filter --complexity_threshold 7 \
            --detect_adapter_for_pe \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size=4 \
            --cut_mean_quality 20 \
            --overrepresentation_analysis 2>&1 | tee -a {log}
        """

rule fastp_no_dedup:
    input: 
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        clumped = "data/omics/{sample_type}/{sample}/reads/done.touch"
    output: 
        fwd_reads = temp("data/omics/{sample_type}/{sample}/reads/fastp_no_dedup_fwd_reads.fastq.gz"),
        rev_reads = temp("data/omics/{sample_type}/{sample}/reads/fastp_no_dedup_rev_reads.fastq.gz"),
        html = "data/omics/{sample_type}/{sample}/reads/qc/fastp_no_dedup.html",
        json = "data/omics/{sample_type}/{sample}/reads/qc/fastp_no_dedup.json"
    conda: "config/conda_yaml/fastp.yaml"
    benchmark:
        "benchmarks/fastp_no_dedup/{sample_type}-{sample}.txt"
    log: "logs/fastp_no_dedup/{sample_type}-{sample}_fastp.log"
    resources: cpus=16, mem_mb = lambda wildcards, attempt: attempt * 60000,
        time_min=2880
    shell:
        """
        # Filter, trim, and remove adapters
        fastp \
            -i {input.fwd_reads} -I {input.rev_reads} \
            -o {output.fwd_reads} -O {output.rev_reads} \
            -h {output.html} -j {output.json} \
            --thread {resources.cpus} \
            -z 9 \
            --length_required 50 \
            --n_base_limit 5 \
            --low_complexity_filter --complexity_threshold 7 \
            --detect_adapter_for_pe \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size=4 \
            --cut_mean_quality 20 \
            --overrepresentation_analysis 2>&1 | tee -a {log}
        """

rule fastqc_fastp:
    input:
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/fastp_fwd_reads.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/fastp_rev_reads.fastq.gz",
    output:
        #fwd_report = "data/omics/{sample_type}/{sample}/reads/fastqc_raw/fwd_fastqc.html",
        #rev_report = "data/omics/{sample_type}/{sample}/reads/fastqc_raw/rev_fastqc.html"
        touch("data/omics/{sample_type}/{sample}/reads/fastqc_fastp/.done")
    conda:
          "config/conda_yaml/fastqc.yaml"
    resources: time_min = 7200, cpus = 24, mem_mb = 60000
    shell:
        """
        mkdir -p data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_fastp
        fastqc -o data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_fastp -t {resources.cpus} {input.fwd_reads}
        fastqc -o data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_fastp -t {resources.cpus} {input.rev_reads}
        """

rule fastqc_raw:
    input:
        clumpify_done = "data/omics/{sample_type}/{sample}/reads/done.touch",
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
    output:
        #fwd_report = "data/omics/{sample_type}/{sample}/reads/fastqc_raw/fwd_fastqc.html",
        #rev_report = "data/omics/{sample_type}/{sample}/reads/fastqc_raw/rev_fastqc.html"
        touch("data/omics/{sample_type}/{sample}/reads/fastqc_raw/.done")
    conda:
          "config/conda_yaml/fastqc.yaml"
    resources: time_min = 7200, cpus = 24, mem_mb = 60000
    shell:
        """
        mkdir -p data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_raw
        fastqc -o data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_raw -t {resources.cpus} {input.fwd_reads}
        fastqc -o data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_raw -t {resources.cpus} {input.rev_reads}
        """

rule multiqc:
    input: 
        "data/omics/{sample_type}/{sample}/reads/fastqc_fastp/.done",
        "data/omics/{sample_type}/{sample}/reads/fastqc_decontam/.done",
        "data/omics/{sample_type}/{sample}/reads/fastqc_raw/.done",
        #"data/omics/{sample_type}/{sample}/reads/fastqc_teal_decon/.done"
    output: 
        multiqc_dir = directory("data/omics/{sample_type}/{sample}/reads/qc/multiqc")
    conda: "config/conda_yaml/multiqc.yaml"
    benchmark:
        "benchmarks/multiqc/{sample_type}-{sample}.txt"
    log: "logs/multiqc/{sample_type}-{sample}.log"
    resources: cpus=1, mem_mb = 20000, time_min=2880
    shell:
        """
        multiqc --interactive -d data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_* -o {output.multiqc_dir}
        """

rule get_contaminants:
    output: 
        human_genome = "data/reference/contaminants/human.fa.gz",
        spike_ins = "data/reference/contaminants/spike-ins.fa"
    resources: cpus = 1
    shell:
        """
        wget -O {output.human_genome} http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.p13.genome.fa.gz
        wget -O {output.spike_ins} https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis/raw/master/spike-ins.fa
        """

rule bb_index:
    input:
        human_genome = ancient(rules.get_contaminants.output.human_genome)
    output:
        "data/reference/contaminants/ref/genome/1/summary.txt",
        index = directory("data/reference/contaminants/ref/")
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/bbmap_index.log"
    benchmark:
        "benchmarks/bb_index.txt"
    resources: cpus = 8, mem_mb = 50000
    shell:
        """        
        bbmap.sh \
            ref={input.human_genome} \
            path={params.bbmap_index_path} \
            t={resources.cpus} \
            2>&1 | tee {log}
        """

rule remove_contaminants:
    input:
        dedup_reads_fwd = rules.deduplicate.output.dedup_reads_fwd,
        dedup_reads_rev = rules.deduplicate.output.dedup_reads_rev,
        human_genome = rules.get_contaminants.output.human_genome,
        spike_ins = rules.get_contaminants.output.spike_ins,
        adapters = "data/reference/contaminants/adapters.fa",
        bbmap_index = "data/reference/contaminants/ref"
    output:
        trimmed_fwd = "data/omics/{sample_type}/{sample}/reads/trimmed_fwd_reads.fastq.gz",
        trimmed_rev = "data/omics/{sample_type}/{sample}/reads/trimmed_rev_reads.fastq.gz",
        phix_rm_fwd = "data/omics/{sample_type}/{sample}/reads/phix_fwd_reads.fastq.gz",
        phix_rm_rev = "data/omics/{sample_type}/{sample}/reads/phix_rev_reads.fastq.gz",
        decon_fwd = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads.fastq.gz",
        decon_rev = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads.fastq.gz",
        cleaned_fwd = "data/omics/{sample_type}/{sample}/reads/cleaned_fwd_reads.fastq.gz",
        cleaned_rev = "data/omics/{sample_type}/{sample}/reads/cleaned_rev_reads.fastq.gz"
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/read_qc/{sample_type}-{sample}.log"
    benchmark:
        "benchmarks/remove_contaminants/{sample_type}-{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 120000, time_min = 2880
    shell:
        """       
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
        echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

        bbduk.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.dedup_reads_fwd} \
            in2={input.dedup_reads_rev} \
            out1={output.trimmed_fwd} \
            out2={output.trimmed_rev} \
            t={resources.cpus} \
            minlen=50 \
            qtrim=rl \
            trimq=15 \
            ref={input.adapters} \
            path={params.bbmap_index_path} \
            ktrim=r k=23 mink=11 hdist=1 \
            2>&1 | tee {log}
        
        echo "\n\n***doing spike-in removal***\n\n" >> {log}
        
        bbduk.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={output.trimmed_fwd} \
            in2={output.trimmed_rev} \
            outu1={output.phix_rm_fwd} \
            outu2={output.phix_rm_rev} \
            t={resources.cpus} k=31 hdist=1 \
            ref={input.spike_ins} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}
        
        echo "\n\n***doing remove contaminants***\n\n" >> {log}
        bbmap.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={output.phix_rm_fwd} \
            in2={output.phix_rm_rev} \
            outu1={output.decon_fwd} \
            outu2={output.decon_rev} \
            t={resources.cpus} fast=t \
            ref={input.human_genome} \
            path={params.bbmap_index_path} \
            -Xmx{resources.mem_mb}m  \
            2>&1 | tee -a {log}
        
        echo "\n\n***Running RemovePolyPairs.pl***\n\n" | tee -a {log}
        perl code/RemovePolyPairs.pl {output.decon_fwd} {output.decon_rev} 50 {output.cleaned_fwd} {output.cleaned_rev} 2>&1 | tee -a {log}
        """

rule count_reads:
    input:
        raw_reads_fwd = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        raw_reads_rev = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        deduped_reads_fwd = "data/omics/{sample_type}/{sample}/reads/dedup_reads_fwd.fastq.gz",
        deduped_reads_rev = "data/omics/{sample_type}/{sample}/reads/dedup_reads_rev.fastq.gz",
        qual_filt_and_trimmed_fwd = "data/omics/{sample_type}/{sample}/reads/trimmed_fwd_reads.fastq.gz",
        qual_filt_and_trimmed_rev = "data/omics/{sample_type}/{sample}/reads/trimmed_rev_reads.fastq.gz",
        decon_reads_fwd = "data/omics/{sample_type}/{sample}/reads/cleaned_fwd_reads.fastq.gz",
        decon_reads_rev = "data/omics/{sample_type}/{sample}/reads/cleaned_rev_reads.fastq.gz"
    output:
        "data/omics/{sample_type}/{sample}/reads/{sample}_read_count.tsv"
    shell:
        """
        printf "read_state\tfwd_read_count\trev_read_count\n" > {output}
        printf "raw_reads\t$(($(zcat {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output}
        printf "deduped_reads\t$(($(zcat {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output}
        printf "filt_and_trimmed_reads\t$(($(zcat {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(zcat {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output}
        printf "decon_reads\t$(($(zcat {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
        """

rule remove_contaminants_fastp:
    input:
        dedup_reads_fwd = rules.fastp.output.fwd_reads,
        dedup_reads_rev = rules.fastp.output.rev_reads,
        human_genome = ancient(rules.get_contaminants.output.human_genome),
        spike_ins = ancient(rules.get_contaminants.output.spike_ins),
        bbmap_index = ancient("data/reference/contaminants/ref")
    output:
        phix_rm_fwd = temp("data/omics/{sample_type}/{sample}/reads/phix_fwd_reads_fastp.fastq.gz"),
        phix_rm_rev = temp("data/omics/{sample_type}/{sample}/reads/phix_rev_reads_fastp.fastq.gz"),
        decon_fwd = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        decon_rev = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/remove_contaminants_fastp/{sample_type}-{sample}.log"
    benchmark:
        "benchmarks/remove_contaminants_fastp/{sample_type}-{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 120000, time_min = 2880
    shell:
        """
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
        echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

        bbduk.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.dedup_reads_fwd} \
            in2={input.dedup_reads_rev} \
            outu1={output.phix_rm_fwd} \
            outu2={output.phix_rm_rev} \
            t={resources.cpus} k=31 hdist=1 \
            ref={input.spike_ins} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}
        
        echo "\n\n***doing remove contaminants***\n\n" >> {log}
        bbmap.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={output.phix_rm_fwd} \
            in2={output.phix_rm_rev} \
            outu1={output.decon_fwd} \
            outu2={output.decon_rev} \
            t={resources.cpus} fast=t \
            ref={input.human_genome} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}
        """

rule remove_contaminants_fastp_no_dedup:
    input:
        dedup_reads_fwd = rules.fastp_no_dedup.output.fwd_reads,
        dedup_reads_rev = rules.fastp_no_dedup.output.rev_reads,
        human_genome = ancient(rules.get_contaminants.output.human_genome),
        spike_ins = ancient(rules.get_contaminants.output.spike_ins),
        bbmap_index = ancient("data/reference/contaminants/ref")
    output:
        phix_rm_fwd = temp("data/omics/{sample_type}/{sample}/reads/phix_fwd_reads_fastp_no_dedup.fastq.gz"),
        phix_rm_rev = temp("data/omics/{sample_type}/{sample}/reads/phix_rev_reads_fastp_no_dedup.fastq.gz"),
        decon_fwd = "data/omics/{sample_type}/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz",
        decon_rev = "data/omics/{sample_type}/{sample}/reads/noDedup_decon_rev_reads_fastp.fastq.gz"
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/remove_contaminants_fastp_no_dedup/{sample_type}-{sample}.log"
    benchmark:
        "benchmarks/remove_contaminants_fastp_no_dedup/{sample_type}-{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 120000, time_min = 2880
    shell:
        """
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
        echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

        bbduk.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.dedup_reads_fwd} \
            in2={input.dedup_reads_rev} \
            outu1={output.phix_rm_fwd} \
            outu2={output.phix_rm_rev} \
            t={resources.cpus} k=31 hdist=1 \
            ref={input.spike_ins} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}
        
        echo "\n\n***doing remove contaminants***\n\n" >> {log}
        bbmap.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={output.phix_rm_fwd} \
            in2={output.phix_rm_rev} \
            outu1={output.decon_fwd} \
            outu2={output.decon_rev} \
            t={resources.cpus} fast=t \
            ref={input.human_genome} \
            path={params.bbmap_index_path} \
            2>&1 | tee -a {log}
        """

rule ribodetector:
    input:
        decon_fwd = "data/omics/{sample_type}/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz",
        decon_rev = "data/omics/{sample_type}/{sample}/reads/noDedup_decon_rev_reads_fastp.fastq.gz"
    output:
        ribo_removed_fwd = "data/omics/{sample_type}/{sample}/reads/ribodeplete_noDedup_decon_fwd_reads_fastp.fastq.gz",
        ribo_removed_rev = "data/omics/{sample_type}/{sample}/reads/ribodeplete_noDedup_decon_rev_reads_fastp.fastq.gz",
        ribo_fwd = "data/omics/{sample_type}/{sample}/reads/ribo_only_fwd_reads_fastp.fastq.gz",
        ribo_rev = "data/omics/{sample_type}/{sample}/reads/ribo_only_rev_reads_fastp.fastq.gz"
    params:
    conda: "config/conda_yaml/ribodetector.yaml"
    log: "logs/ribodetector/{sample_type}-{sample}.log"
    benchmark:
        "benchmarks/ribodetector/{sample_type}-{sample}.txt"
    resources: cpus = 24, mem_mb = 120000, time_min = 2880, gpu_mem_gb = 12, partition = "gpu", gpu = 1
    shell:
        """
        ribodetector \
            -t {resouces.cpus} \
            -i {input.decon_fwd} {input.decon_rev} \
            -m {resources.gpu_mem_gb} \
            -e rrna \
            --chunk_size 256 \
            -r {output.ribo_fwd} {output.ribo_rev} \
            -o {ouptput.ribo_removed_fwd} {ouptput.ribo_removed_rev} 2>&1 | tee -a {log}
        """

rule make_read_blastdb:
    input: 
        decon_fwd = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        decon_rev = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        concat_reads = temp("data/omics/metagenomes/{sample}/reads/concat_decon_reads_fastp.fasta"),
        db_index_made = touch("data/omics/metagenomes/{sample}/reads/.read_blast_index_created")
    conda: "config/conda_yaml/seqtk.yaml"
    log:
        "logs/makeblastdb_reads/{sample}.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        seqtk seq -a {input.decon_fwd} > {output.concat_reads}
        seqtk seq -a {input.decon_rev} >> {output.concat_reads}

        makeblastdb -in {output.concat_reads} -dbtype nucl -logfile {log}
        """


rule blast_nuc:
    input:
        blast_db = rules.make_read_blastdb.output.concat_reads,
        blast_db_index = rules.make_read_blastdb.output.db_index_made,
        gene = "data/reference/blast_queries/{query}.fasta"
    output:
        blast_res = "data/omics/metagenomes/{sample}/BLAST/{query}__{sample}.blastn"
    log:
        "logs/BLAST/{query}__{sample}.log"
    resources: cpus=32, mem_mb=5000, time_min=120
    shell:
        """
        blastn -query {input.gene} \
            -db {input.blast_db} \
            -out {output.blast_res} \
            -outfmt '6 std qcovs stitle qseq sseq' \
            -num_threads {resources.cpus}

        # removed additional output columns: -outfmt '6 std qcovs stitle' \
        # also removed database size standardization: -dbsize 1000000 \
        """

rule fastqc_decontam:
    input:
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        touch("data/omics/{sample_type}/{sample}/reads/fastqc_decontam/.done")
    conda:
          "config/conda_yaml/fastqc.yaml"
    resources: cpus = 4, mem_mb = 16000, time_min = 2880
    shell:
        """
        mkdir -p data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_decontam
        fastqc -o data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_decontam -t {resources.cpus} {input.fwd_reads}
        fastqc -o data/omics/{wildcards.sample_type}/{wildcards.sample}/reads/fastqc_decontam -t {resources.cpus} {input.rev_reads}
        """

rule bbnorm:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        fwd_norm = temp("data/omics/{sample_type}/{sample}/reads/bbnorm_fwd_reads.fastq.gz"),
        rev_norm = temp("data/omics/{sample_type}/{sample}/reads/bbnorm_rev_reads.fastq.gz")
    params: "target=100 mindepth=2 bits=16 prefilter ecc=t"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/bbnorm/{sample_type}/{sample}.log"
    benchmark:
        "benchmarks/bbnorm/{sample_type}/{sample}.txt"
    resources: cpus = 36, mem_mb = lambda wildcards, attempt: attempt * 80000, time_min = 2880
    shell:
        """        
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
        echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"
        
        bbnorm.sh \
            -Xmx${{bbmap_mem}}m -eoom \
            in1={input.fwd_reads} \
            in2={input.rev_reads} \
            out1={output.fwd_norm} \
            out2={output.rev_norm} \
            t={resources.cpus} \
            {params} \
            1>>{log} 2>&1
        """

rule count_reads_fastp:
    input:
        raw_reads_fwd = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        raw_reads_rev = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        deduped_reads_fwd = "data/omics/{sample_type}/{sample}/reads/fastp_fwd_reads.fastq.gz",
        deduped_reads_rev = "data/omics/{sample_type}/{sample}/reads/fastp_rev_reads.fastq.gz",
        qual_filt_and_trimmed_fwd = rules.fastp.output.fwd_reads,
        qual_filt_and_trimmed_rev = rules.fastp.output.rev_reads,
        decon_reads_fwd = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        decon_reads_rev = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        #bbnorm_reads_fwd = rules.bbnorm.output.fwd_norm,
        #bbnorm_reads_rev = rules.bbnorm.output.rev_norm,
    output:
        "data/omics/{sample_type}/{sample}/reads/{sample}_read_count_fastp.tsv"
    resources: cpus=4
    shell:
        """
        printf "read_state\tfwd_read_count\trev_read_count\n" > {output} &&
        printf "raw_reads\t$(($(pigz -dc -p {resources.cpus} {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output} &&
        printf "deduped_reads\t$(($(pigz -dc -p {resources.cpus} {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output} &&
        printf "filt_and_trimmed_reads\t$(($(pigz -dc -p {resources.cpus} {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output} &&
        printf "decon_reads\t$(($(pigz -dc -p {resources.cpus} {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
        """
        #printf "bbnorm_reads\t$(($(pigz -dc -p {resources.cpus} {input.bbnorm_reads_fwd} | wc -l) / 4 ))\t$(($(pigz -dc -p {resources.cpus} {input.bbnorm_reads_rev} | wc -l) / 4 ))\n" >> {output}


rule assemble_metaspades:
    input:
        fwd_reads = rules.bbnorm.output.fwd_norm,
        rev_reads = rules.bbnorm.output.rev_norm
    output:
        assembly_dir = directory("data/omics/{sample_type}/{sample}/assembly/metaspades"),
        contigs = "data/omics/{sample_type}/{sample}/assembly/metaspades/contigs.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/metaspades/{sample_type}_{sample}.log"
    benchmark: "benchmarks/metaspades/{sample_type}_{sample}.txt"
    #resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 120000
    #resources: cpus = 24, time_min=20000, mem_mb = 500000, partition = "largemem"
    resources: cpus = 36, time_min=20000, mem_mb = 170000
    #resources: cpus = 64, time_min=20000, mem_mb = 500000
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus}

        metaspades.py \
            -t {resources.cpus} \
            --memory $(({resources.mem_mb}/1024)) \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -o {output.assembly_dir} 2>&1 | tee {log}
        """

rule assemble_metaspades_noNORM:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        assembly_dir = directory("data/omics/{sample_type}/{sample}/assembly/metaspades_noNORM"),
        contigs = "data/omics/{sample_type}/{sample}/assembly/metaspades_noNORM/contigs.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/metaspades_noNORM/{sample_type}_{sample}.log"
    benchmark: "benchmarks/metaspades_noNORM/{sample_type}_{sample}.txt"
    #resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 120000
    #resources: cpus = 24, time_min=20000, mem_mb = 500000, partition = "largemem"
    #resources: cpus = 36, time_min=20000, mem_mb = 170000
    resources: cpus = 64, time_min=20000, mem_mb = 500000
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus}

        metaspades.py \
            -t {resources.cpus} \
            --memory $(({resources.mem_mb}/1024)) \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -o {output.assembly_dir} 2>&1 | tee {log}
        """

rule assemble_biosyntheticSPAdes:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        assembly_dir = directory("data/omics/{sample_type}/{sample}/assembly/biosyntheticSPAdes"),
        contigs = "data/omics/{sample_type}/{sample}/assembly/biosyntheticSPAdes/scaffolds.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/biosyntheticSPAdes/{sample_type}_{sample}.log"
    benchmark: "benchmarks/biosyntheticSPAdes/{sample_type}_{sample}.txt"
    #resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 120000
    resources: cpus = 24, time_min=20000, mem_mb = 500000, partition = "largemem"
    #resources: cpus = 36, time_min=20000, mem_mb = 170000
    #resources: cpus = 64, time_min=20000, mem_mb = 500000
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus}

        metaspades.py \
            -t {resources.cpus} \
            --bio \
            --memory $(({resources.mem_mb}/1024)) \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -o {output.assembly_dir} 2>&1 | tee {log}
        """

rule assemble_RNAspades:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        contigs = "data/omics/{sample_type}/{sample}/assembly/RNAspades/transcripts.fasta"
    params:
        assembly_dir = "data/omics/{sample_type}/{sample}/assembly/RNAspades"
    conda: "config/conda_yaml/spades4.yaml"
    log: "logs/assembly/RNAspades/{sample_type}_{sample}.log"
    benchmark: "benchmarks/RNAspades/{sample_type}_{sample}.txt"
    #resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 120000
    #resources: cpus = 24, time_min=20000, mem_mb = 500000, partition = "largemem"
    #resources: cpus = 36, time_min=20000, mem_mb = 170000
    resources: cpus = 64, time_min=20000, mem_mb = 500000
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus}

        rnaspades.py \
            -t {resources.cpus} \
            --memory $(({resources.mem_mb}/1024)) \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -o {params.assembly_dir} 2>&1 | tee {log}
        """

rule assemble_biosyntheticSPAdes_100x:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        assembly_dir = directory("data/omics/{sample_type}/{sample}/assembly/biosyntheticSPAdes_100x"),
        contigs = "data/omics/{sample_type}/{sample}/assembly/biosyntheticSPAdes_100x/scaffolds.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/metaspades_100x/{sample_type}_{sample}.log"
    benchmark: "benchmarks/metaspades_100x/{sample_type}_{sample}.txt"
    #resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 120000
    resources: cpus = 24, time_min=20000, mem_mb = 500000, partition = "largemem"
    #resources: cpus = 36, time_min=20000, mem_mb = 170000
    #resources: cpus = 64, time_min=20000, mem_mb = 500000
    shell:
        """
        export OMP_NUM_THREADS={resources.cpus}

        metaspades.py \
            -t {resources.cpus} \
            --bio \
            --memory $(({resources.mem_mb}/1024)) \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -o {output.assembly_dir} 2>&1 | tee {log}
        """

rule rename_metaspades_contigs:
    input: 
        script = "code/rename_contigs.R",
        contigs = "data/projects/{project}/{sample_type}/{sample}/assembly/metaspades_noNORM/contigs.fasta",
        #assembly_done = "data/omics/{sample_type}/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/projects/{project}/{sample_type}/{sample}/assembly/metaspades_noNORM/contigs.renamed.fasta",
        contig_info = "data/projects/{project}/{sample_type}/{sample}/assembly/metaspades_noNORM/contigs_info.tsv"
        #done = touch("data/omics/{sample_type}/{sample}/assembly/megahit/.contigs_renamed")
    container: "docker://eandersk/r_microbiome"
    resources: cpus = 1, time_min=200, mem_mb = 50000
    shell:
        """
        pwd
        
        ./{input.script} \
            -i {input.contigs} \
            -o {output.contigs} \
            -s {output.contig_info} \
            -p {wildcards.sample}
        """

rule rename_RNAspades_contigs:
    input: 
        script = "code/rename_contigs.R",
        contigs = "data/omics/{sample_type}/{sample}/assembly/RNAspades/transcripts.fasta",
        #assembly_done = "data/omics/{sample_type}/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/omics/{sample_type}/{sample}/assembly/RNAspades/contigs.renamed.fasta",
        contig_info = "data/omics/{sample_type}/{sample}/assembly/RNAspades/contigs_info.tsv"
        #done = touch("data/omics/{sample_type}/{sample}/assembly/megahit/.contigs_renamed")
    container: "docker://eandersk/r_microbiome"
    resources: cpus = 1, time_min=200, mem_mb = 50000
    shell:
        """
        pwd
        
        ./{input.script} \
            -i {input.contigs} \
            -o {output.contigs} \
            -s {output.contig_info} \
            -p {wildcards.sample}
        """


rule assemble_megahit:
    input:
        fwd_reads = rules.bbnorm.output.fwd_norm,
        rev_reads = rules.bbnorm.output.rev_norm
    output:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit/final.contigs.fa",
        #touch("data/omics/{sample_type}/{sample}/assembly/megahit/.done")
    params:
        assembly_dir = "data/omics/{sample_type}/{sample}/assembly/megahit"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit/{sample_type}_{sample}.log"
    benchmark: "benchmarks/megahit/{sample_type}_{sample}.txt"
    resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        rm -r {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
        megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {params.assembly_dir} > {log}
        """

rule assemble_megahit_noNORM:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.fa",
        #touch("data/omics/{sample_type}/{sample}/assembly/megahit/.done")
    params:
        assembly_dir = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit_noNORM/{sample_type}_{sample}.log"
    benchmark: "benchmarks/megahit_noNORM/{sample_type}_{sample}.txt"
    resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 150000
    shell:
        """
        rm -r {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
        megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {params.assembly_dir} > {log}
        """


rule rename_contigs:
    input: 
        script = ancient("code/rename_contigs.R"),
        #assembly_dir = "data/omics/{sample_type}/{sample}/assembly/megahit",
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.fa",
        #assembly_done = "data/omics/{sample_type}/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        contig_info = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/contigs_info.tsv"
        #done = touch("data/omics/{sample_type}/{sample}/assembly/megahit/.contigs_renamed")
    container: "docker://eandersk/r_microbiome"
    resources: cpus = 1, time_min=200, mem_mb = 50000
    shell:
        """
        pwd
        
        ./{input.script} \
            -i {input.contigs} \
            -o {output.contigs} \
            -s {output.contig_info} \
            -p {wildcards.sample}
        """

rule quast_megahit:
    input:
        megahit_contigs = rules.rename_contigs.output.contigs
    output:
        report = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/quast/report.tsv"
    params:
        out_dir = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/quast"
    log: "logs/assembly/quast_megahit/{sample_type}_{sample}.log"
    benchmark:
        "benchmarks/assembly/quast_megahit/{sample_type}_{sample}.txt"
    conda:
        "config/conda_yaml/quast.yaml"
    resources:
        cpus = 1, mem_mb = 20000
    shell:
        """
        quast.py {input.megahit_contigs} -o {params.out_dir} 2>&1 | tee {log}
        """

# rule concat_reads_for_COassembly:
#     input:
#         fwd_reads = expand(rules.remove_contaminants_fastp.output.decon_fwd, sample = jgi_samples, sample_type = "metagenomes"),
#         rev_reads = expand(rules.remove_contaminants_fastp.output.decon_rev, sample = jgi_samples, sample_type = "metagenomes")
#     output:
#         concat_fwd = temp("tmp/fwd_concat.fastq"),
#         concat_rev = temp("tmp/rev_concat.fastq")
#     resources: cpus = 1, time_min=20000, mem_mb = 4000
#     shell:
#         """
#         zcat {input.fwd_reads} > {output.concat_fwd}
#         zcat {input.rev_reads} > {output.concat_rev}
#         """

rule COassemble_megahit:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM_coassembly/final.contigs.fa",
        #touch("data/omics/{sample_type}/{sample}/assembly/megahit/.done")
    params:
        assembly_dir = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM_coassembly"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit/{sample_type}_{sample}_coassembly.log"
    benchmark: "benchmarks/megahit/{sample_type}_{sample}_coassembly.txt"
    #resources: cpus = 92, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 150000
    resources: cpus = 32, time_min=20000, mem_mb = lambda wildcards, attempt: attempt * 1200000, partition = "largemem"
    shell:
        """
        rm -rf {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
        megahit \
            -t {resources.cpus} \
            --min-count 2 \
            --k-list 31,37,47,57,67,77,87,99 \
            -m $(({resources.mem_mb} * 1000000)) \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            -o {params.assembly_dir} > {log}
        """

rule merge_assemblies:
    input:
        metaspades_contigs = rules.assemble_metaspades.output.contigs,
        megahit_contigs = rules.rename_contigs.output.contigs,
        merge_contigs_script = "code/Strain-Level_Metagenome_Analysis/Merge_Contigs.pl"
    output:
        concat_contigs = "data/omics/metagenomes/{sample}/assembly/contigs_concat.fa",
        dedup1 = "data/omics/metagenomes/{sample}/assembly/{sample}_dedup1.fa",
        dedup2 = "data/omics/metagenomes/{sample}/assembly/{sample}_dedup2.fa",
        dedup3 = "data/omics/metagenomes/{sample}/assembly/{sample}_dedup3.fa",
        dedup4 = "data/omics/metagenomes/{sample}/assembly/{sample}_dedup4.fa",
        dedup4_dot = "data/omics/metagenomes/{sample}/assembly/{sample}_graph4.dot",
        dedup5 = "data/omics/metagenomes/{sample}/assembly/{sample}_MERGED_CONTIGS.fa",
        dedup6 = "data/omics/metagenomes/{sample}/assembly/{sample}_MCDD.fa"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/merge/{sample}_merge_assemblies.log"
    benchmark: "benchmarks/assembly/merge/{sample}.txt"
    resources: cpus = 8, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        cat {input.metaspades_contigs} {input.megahit_contigs} > {output.concat_contigs}
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.concat_contigs} out={output.dedup1} tuc mid=99 minscaf=200 rnc=f ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f sort=length absorbcontainment=t mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 2>&1 | tee {log}
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup1} out={output.dedup2} tuc mid=99 minscaf=200 rnc=t ngn=t fo c pc=t fmj=t rc=t cc=t fcc=t mst=f sort=length absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 2>&1 | tee -a {log}
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup2} out={output.dedup3} tuc mid=99 minscaf=200 rnc=t ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f ordered=t absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 2>&1 | tee -a {log}
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup3} out={output.dedup4} tuc mid=99 minscaf=200 rnc=f ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f ordered=t absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t dot={output.dedup4_dot} t={resources.cpus} 2>&1 | tee -a {log}
        perl {input.merge_contigs_script} data/omics/metagenomes/{wildcards.sample}/assembly/{wildcards.sample} 99 2>&1 | tee -a {log}
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup5} out={output.dedup6} t={resources.cpus} tuc mid=99 minscaf=200 overwrite=f 2>&1 | tee -a {log}
        """

rule get_MEC:
    output: 
        mec_dir = directory("code/MEC"),
        mec_script = "code/MEC/src/mec.py"
    shell:
        """
        rm -rf {output.mec_dir}
        cd code
        git clone https://github.com/bioinfomaticsCSU/MEC.git
        chmod +x ../{output.mec_script}
        """

rule correct_contigs:
    input:
        assembly = rules.merge_assemblies.output.dedup6,
        fwd_reads = rules.remove_contaminants.output.cleaned_fwd,
        rev_reads = rules.remove_contaminants.output.cleaned_rev,
        mec_script = "code/MEC/src/mec.py"
    output:
        read_mapping = temp("data/omics/metagenomes/{sample}/assembly/{sample}_cleaned_reads_unsorted.sam"),
        read_mapping_sorted = temp("data/omics/metagenomes/{sample}/assembly/{sample}_cleaned_reads.bam")
    params:
        assembly_dir = "data/omics/metagenomes/{sample}/assembly"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/correct_contigs/{sample}_correct_contigs.log"
    benchmark: "benchmarks/assembly/correct_contigs/{sample}.txt"
    resources: cpus = 4, mem_mb = lambda wildcards, attempt: attempt * 100000, time_min = 2880
    shell:
        """ 
        use_mem=$(({resources.mem_mb} - 4000)) # bbmap uses more memory than it's told
        mem_per_thread=$(( ({resources.mem_mb} - 2000) / {resources.cpus}))

        #ALIGN READS TO CONTIGS
        bbmap.sh -da -Xmx${{use_mem}}m -eoom ref={input.assembly} path={params.assembly_dir} t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 idfilter=0.999 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} out={output.read_mapping} 2>&1 | tee {log}
        samtools view -bShu {output.read_mapping} | samtools sort -m ${{mem_per_thread}}M -@ {resources.cpus} -o {output.read_mapping_sorted} 2>&1 | tee -a {log}
        samtools index {output.read_mapping_sorted} 2>&1 | tee -a {log}
        """

rule MEC:
    input:
        assembly = rules.merge_assemblies.output.dedup6,
        read_mapping_sorted = rules.correct_contigs.output.read_mapping_sorted,
        mec_script = rules.get_MEC.output.mec_script
    output:        
        corrected_assembly = "data/omics/metagenomes/{sample}/assembly/{sample}_merged_and_corrected_assembly.fa"
    conda: "config/conda_yaml/mec.yaml"
    log: "logs/assembly/MEC/{sample}.log"
    benchmark: "benchmarks/assembly/MEC/{sample}.txt"
    resources: cpus = 1, mem_mb = lambda wildcards, attempt: attempt * 40000, time_min = 2880
    shell:
        """
        python {input.mec_script} \
            -bam {input.read_mapping_sorted} \
            -i {input.assembly} \
            -o {output.corrected_assembly} 2>&1 | tee {log}
        """

rule quast:
    input:
        combined_contigs = rules.MEC.output.corrected_assembly,
        metaspades_contigs = rules.assemble_metaspades.output.contigs,
        megahit_contigs = rules.rename_contigs.output.contigs
    output: directory("data/omics/metagenomes/{sample}/assembly/quast")
    log: "logs/assembly/quast/{sample}.log"
    benchmark:
        "benchmarks/assembly/quast/{sample}.txt"
    conda:
        "config/conda_yaml/quast.yaml"
    resources:
        cpus = 1, mem_mb = 20000
    shell:
        """
        quast.py {input.combined_contigs} {input.metaspades_contigs} {input.megahit_contigs} -o {output} 2>&1 | tee {log}
        """

rule prodigal:
    input:
        assembly = rules.rename_contigs.output.contigs
    output:
        proteins = "data/omics/{sample_type}/{sample}/proteins/{sample}_PROTEINS.faa",
        genes = "data/omics/{sample_type}/{sample}/genes/{sample}_GENES.fna",
        gbk = "data/omics/{sample_type}/{sample}/genes/{sample}_GENES.gbk"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/prodigal/{sample_type}-{sample}.log"
    benchmark: "benchmarks/prodigal/{sample_type}-{sample}.txt"
    resources: cpus = 1, mem_mb = lambda wildcards, attempt: attempt * 16000, time_min = 2880
    shell:
        """
        prodigal \
            -p meta \
            -i {input.assembly} \
            -a {output.proteins} \
            -o {output.gbk} \
            -d {output.genes} 2>&1 | tee {log}
        """

rule prodigal_metaT:
    input:
        assembly = "data/omics/metatranscriptomes/{sample}/assembly/RNAspades/contigs.renamed.fasta"
    output:
        proteins = "data/omics/metatranscriptomes/{sample}/proteins/{sample}_PROTEINS.faa",
        genes = "data/omics/metatranscriptomes/{sample}/genes/{sample}_GENES.fna",
        gbk = "data/omics/metatranscriptomes/{sample}/genes/{sample}_GENES.gbk"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/prodigal_metaT/{sample}.log"
    benchmark: "benchmarks/prodigal_metaT/{sample}.txt"
    resources: cpus = 1, mem_mb = lambda wildcards, attempt: attempt * 16000, time_min = 2880
    shell:
        """
        prodigal \
            -p meta \
            -i {input.assembly} \
            -a {output.proteins} \
            -o {output.gbk} \
            -d {output.genes} 2>&1 | tee {log}
        """

ruleorder: prodigal_metaT > prodigal

rule calc_gene_abundance:
    input:
        genes = rules.prodigal.output.genes,
        proteins = rules.prodigal.output.proteins,
        assembly = rules.rename_contigs.output.contigs,
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        reads_vs_genes_rpkm = "data/omics/{sample_type}/{sample}/genes/{sample}_READSvsGENES.rpkm",
        reads_vs_contigs_rpkm = "data/omics/{sample_type}/{sample}/assembly/{sample}_READSvsCONTIGS.rpkm",
        reads_vs_assembly_sam_gz = "data/omics/{sample_type}/{sample}/assembly/{sample}_READSvsCONTIGS.sam.gz"
    params:
        reads_vs_assembly_sam = "data/omics/{sample_type}/{sample}/assembly/{sample}_READSvsCONTIGS.sam"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/calc_gene_abundance/{sample_type}-{sample}.log"
    benchmark: "benchmarks/calc_gene_abundance/{sample_type}-{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 64000, time_min = 2880
    shell:
        """
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.genes}) ref={input.genes} rpkm={output.reads_vs_genes_rpkm} 2>&1 | tee -a {log}
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.assembly}) ref={input.assembly} rpkm={output.reads_vs_contigs_rpkm} 32bit=t outm={params.reads_vs_assembly_sam} 2>&1 | tee -a {log}
        pigz -9 -p {resources.cpus} {params.reads_vs_assembly_sam}
        """

rule calc_gene_abundance_metaT:
    input:
        genes = "data/omics/metatranscriptomes/{sample}/genes/{sample}_GENES.fna",
        proteins = "data/omics/metatranscriptomes/{sample}/proteins/{sample}_PROTEINS.faa",
        assembly = "data/omics/metatranscriptomes/{sample}/assembly/RNAspades/contigs.renamed.fasta",
        fwd_reads = "data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_rev_reads_fastp.fastq.gz"
    output:
        reads_vs_genes_rpkm = "data/omics/metatranscriptomes/{sample}/genes/{sample}_READSvsGENES.rpkm",
        reads_vs_contigs_rpkm = "data/omics/metatranscriptomes/{sample}/assembly/{sample}_READSvsCONTIGS.rpkm",
        reads_vs_assembly_sam_gz = "data/omics/metatranscriptomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam.gz"
    params:
        reads_vs_assembly_sam = "data/omics/metatranscriptomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/calc_gene_abundance_metaT/{sample}.log"
    benchmark: "benchmarks/calc_gene_abundance_metaT/{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 64000, time_min = 2880
    shell:
        """
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.genes}) ref={input.genes} rpkm={output.reads_vs_genes_rpkm} 2>&1 | tee -a {log}
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.assembly}) ref={input.assembly} rpkm={output.reads_vs_contigs_rpkm} 32bit=t outm={params.reads_vs_assembly_sam} 2>&1 | tee -a {log}
        pigz -9 -p {resources.cpus} {params.reads_vs_assembly_sam}
        """
    
ruleorder: calc_gene_abundance_metaT > calc_gene_abundance

# rule download_uniref:
#     output: 
#         uniref100="data/reference/uniref/uniref100.fasta.gz"
#     #conda: "config/conda_yaml/main.yaml"
#     log: "logs/make_diamond_uniref_db/download_uniref.log"
#     resources: cpus = 1, mem_mb=1000
#     shell:
#         """
#         #cd data/reference/uniref
#         #wget -N https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
        
#         #mkdir data/reference/uniref
#         cp data/reference/uniref100.fasta.gz data/reference/uniref/
#         """

rule make_diamond_uniref_db:
    input: 
        uniref100="data/reference/uniref/uniref100.fasta.gz"
    output: 
        uniref100_diamond="data/reference/uniref/uniref100.dmnd"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/make_diamond_uniref_db/make_diamond_uniref_db.log"
    benchmark: "benchmarks/make_diamond_uniref_db/make_diamond_uniref_db.txt"
    resources: cpus = 32, mem_mb=150000, time_min = 2880
    shell:
        """
        diamond makedb --threads {resources.cpus} --in {input.uniref100} -d data/reference/uniref/uniref100 2>&1 | tee {log}
        """

rule align_to_uniref:
    input:
        diamond_db = rules.make_diamond_uniref_db.output.uniref100_diamond,
        genes = rules.prodigal.output.genes
    output:
        gene_uniref_alignment = "data/omics/{sample_type}/{sample}/{sample}_GENES.m8"
    params:
        "--top 0.5 --query-cover 50 --strand both -f 6 qseqid qlen sseqid slen qstart qend sstart send evalue pident mismatch qcovhsp scovhsp"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/align_to_uniref/{sample_type}-{sample}_align_to_uniref.log"
    benchmark: "benchmarks/align_to_uniref/{sample_type}-{sample}_align_to_uniref.txt"
    resources: cpus = 32, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        diamond blastx \
            -d {input.diamond_db} \
            -q {input.genes} \
            --threads {resources.cpus} \
            -o {output.gene_uniref_alignment} \
            {params} 2>&1 | tee {log}
        """

rule align_to_uniref_mmseqs2:
    input:
        genes = rules.prodigal.output.genes
    output:
        gene_uniref_alignment = "data/omics/{sample_type}/{sample}/{sample}_GENES_mmseqs.m8"
    params:
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/gene_map_to_uniref/{sample}",
        tmp_dir = "/tmp/kiledal/mmseqs2_genes/{sample}",
        sensitivity="--sensitivity 6",  # Sensitivity setting for MMseqs2
        filter="--min-seq-id 0.5 -c 0.5 --cov-mode 0", # Adjust sequence identity and coverage mode
        format="--format-output 'query,qlen,target,tlen,qstart,qend,tstart,tend,evalue,pident,mismatch,alnlen,gapopen,qcov,tcov'"
    conda:  "config/conda_yaml/mmseqs.yaml"
    log: "logs/align_to_uniref_mmseqs2/{sample_type}-{sample}_align_to_uniref.log"
    benchmark: "benchmarks/align_to_uniref_mmseqs2/{sample_type}-{sample}_align_to_uniref.txt"
    resources: cpus = 32, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        mmseqs easy-search \
            {input.genes} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            {params.sensitivity} \
            {params.filter} \
            {params.format} \
            --threads {resources.cpus} \
            --split-memory-limit 100G \
            2>&1 | tee {log}
        """
        
rule contig_abund:
    input:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        f_seq = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output: 
        coverage_full = "data/omics/{sample_type}/{sample}/{sample}_contig_abund.tsv"
    params:
        tmpdir = "tmp/coverm_contig_abund/{sample}"
    conda: "config/conda_yaml/coverm.yaml"
    log: "logs/contig_abund/{sample}-{sample_type}.log"
    benchmark: "benchmarks/contig_abund/{sample}-{sample_type}.txt"
    resources: cpus=24, mem_mb=120000, time_min=2880 # standard assemblies
    #resources: cpus=24, mem_mb=1000000, time_min=2880, partition = "largemem" # coassembly
    shell:
        """
        export TMPDIR={params.tmpdir}
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/coverm_contig_abund/{wildcards.sample}
        mkdir -p $TMPDIR

        # Link reads w/ naming convention prefered by coverM
        fwd_reads=$(dirname {input.f_seq})/{wildcards.sample}_R1.fq.gz
        rev_reads=$(dirname {input.r_seq})/{wildcards.sample}_R2.fq.gz
        ln {input.f_seq} $fwd_reads
        ln {input.r_seq} $rev_reads

        coverm contig \
            -c $fwd_reads $rev_reads \
            -r {input.contigs} \
            -t {resources.cpus} \
            -m mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
            --output-format sparse \
            --output-file {output.coverage_full} 2>&1 | tee {log}

        rm -r $TMPDIR $fwd_reads $rev_reads 
        """

rule contig_abund_metaT:
    input:
        contigs = "data/omics/metatranscriptomes/{sample}/assembly/RNAspades/contigs.renamed.fasta",
        f_seq = "data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz"
    output: 
        coverage_full = "data/omics/metatranscriptomes/{sample}/assembly/RNAspades/{sample}_contig_abund.tsv"
    params:
        tmpdir = "tmp/coverm_contig_abund/{sample}"
    conda: "config/conda_yaml/coverm.yaml"
    log: "logs/contig_abund_metaT/{sample}.log"
    benchmark: "benchmarks/contig_abund_metaT/{sample}.txt"
    resources: cpus=24, mem_mb=120000, time_min=2880 # standard assemblies
    #resources: cpus=24, mem_mb=1000000, time_min=2880, partition = "largemem" # coassembly
    shell:
        """
        export TMPDIR={params.tmpdir}
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/coverm_contig_abund/{wildcards.sample}
        mkdir -p $TMPDIR

        # Link reads w/ naming convention prefered by coverM
        fwd_reads=$(dirname {input.f_seq})/{wildcards.sample}_R1.fq.gz
        rev_reads=$(dirname {input.r_seq})/{wildcards.sample}_R2.fq.gz
        ln {input.f_seq} $fwd_reads
        ln {input.r_seq} $rev_reads

        coverm contig \
            -c $fwd_reads $rev_reads \
            -r {input.contigs} \
            -t {resources.cpus} \
            -m mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
            --output-format sparse \
            --output-file {output.coverage_full} 2>&1 | tee {log}

        rm -r $TMPDIR $fwd_reads $rev_reads 
        """


# Inspect Kraken2 databases, nescessary to produce nicely formatted taxonomy files for downstream analysis
rule kraken_inspect:
    input:
        db = ancient("data/reference/kraken_databases/{database}")
    output: 
        inspect_file = "data/reference/kraken_databases/{database}/inspect.txt"
    conda: "config/conda_yaml/kraken.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        kraken2-inspect --db {input.db} > {output.inspect_file}
        """

rule add_lineage_to_inspect_gtdb:
    input:
        db = ancient("data/reference/kraken_databases/gtdb_r202"),
        inspect_file = "data/reference/kraken_databases/gtdb_r202/inspect.txt"
    output: 
        inspect_w_lineage = "data/reference/kraken_databases/gtdb_r202/inspect_w_lineage.txt"
    conda: "config/conda_yaml/taxonkit.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        taxonkit lineage \
            {input.inspect_file} \
            --data-dir {input.db}/taxonomy \
            -i 5 \
            -o {output.inspect_w_lineage} 2>&1 | tee {log}
        """

rule add_lineage_to_inspect_refseq:
    input:
        db = ancient("data/reference/kraken_databases/refseq"),
        inspect_file = "data/reference/kraken_databases/refseq/inspect.txt"
    output: 
        inspect_w_lineage_unformatted = temp("data/reference/kraken_databases/refseq/unformatted_inspect_w_lineage.txt"),
        inspect_w_lineage = "data/reference/kraken_databases/refseq/inspect_w_lineage.txt"
    conda: "config/conda_yaml/taxonkit.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        taxonkit lineage \
            {input.inspect_file} \
            --data-dir {input.db}/taxonomy \
            -i 5 \
            -o {output.inspect_w_lineage_unformatted} 2>&1 | tee {log}

        taxonkit reformat \
            {output.inspect_w_lineage_unformatted} \
            -i 7 \
            -P \
            -o {output.inspect_w_lineage} 2>&1 | tee -a {log}
        """


rule kraken_database_tax_merge:
    input:
        script = "code/merge_kraken_tax.R",
        gtdb_tax_info = rules.add_lineage_to_inspect_gtdb.output.inspect_w_lineage,
        refseq_tax_info = rules.add_lineage_to_inspect_refseq.output.inspect_w_lineage
    output:
        combined_tax_info = "data/reference/kraken_tax_info_merged.tsv"
    resources: cpus=1, mem_mb=4000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        ./{input.script} -g {input.gtdb_tax_info} -r {input.refseq_tax_info} -o {output.combined_tax_info}
        """



# Run kraken2 with KrakenUniq like functionality to screen for higher confidence results
# First runs with a GTDB database to annotate bacterial and archaeal reads
rule kraken2_gtdb_w_uniq:
    input:
        f_seq = rules.remove_contaminants.output.cleaned_fwd,
        r_seq = rules.remove_contaminants.output.cleaned_rev,
        # f_seq = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        # r_seq = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        db = "data/reference/kraken_databases/gtdb_r202",
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_report.txt",
        out = temp("data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_out.txt"),
        bracken = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_bracken.txt",
        bracken_report = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_brackenMpa.txt",
        unclass_f = temp("data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_unclassified_1.fasta"),
        unclass_r = temp("data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_unclassified_2.fasta"),
        bracken_input = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150,
        unclass_out = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_unclassified#.fasta"
    conda: "config/conda_yaml/kraken.yaml"
    log: "logs/kraken2_gtdb_w_uniq/{sample_type}-{sample}.log"
    benchmark: "benchmarks/kraken2_gtdb_w_uniq/{sample_type}-{sample}.txt"
    resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            --db {input.db} \
            --minimum-hit-groups 3 \
            --unclassified-out {params.unclass_out} \
            --paired {input.f_seq} {input.r_seq} 2>&1 | tee {log}

        echo "Kraken complete, filtering..." | tee -a {log}

        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

        echo "Filtering complete, running bracken..." | tee -a {log}

        bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report} 2>&1 | tee -a {log}

        ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages 2>&1 | tee -a {log}

        echo "Bracken complete. Quitting." | tee -a {log}
        """

# Any reads not annotated with the GTDB database are then annotated with a RefSeq database
rule kraken2_refseq_w_uniq: ##Run kraken2
    input:
        f_seq = rules.kraken2_gtdb_w_uniq.output.unclass_f,
        r_seq = rules.kraken2_gtdb_w_uniq.output.unclass_r,
        db = "data/reference/kraken_databases/refseq",
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/{sample_type}/{sample}/kraken/refseq_{sample}_report.txt",
        out = "data/omics/{sample_type}/{sample}/kraken/refseq_{sample}_out.txt",
        bracken = "data/omics/{sample_type}/{sample}/kraken/refseq_{sample}_bracken.txt",
        bracken_report = "data/omics/{sample_type}/{sample}/kraken/refseq_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/{sample_type}/{sample}/kraken/refseq_{sample}_brackenMpa.txt",
        bracken_input = "data/omics/{sample_type}/{sample}/kraken/refseq_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150
    conda: "config/conda_yaml/kraken.yaml"
    log: "logs/kraken2_refseq_w_uniq/{sample_type}-{sample}.log"
    benchmark: "benchmarks/kraken2_refseq_w_uniq/{sample_type}-{sample}.txt"
    resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            --db {input.db} \
            --minimum-hit-groups 3 \
            --paired {input.f_seq} {input.r_seq} 2>&1 | tee {log}

        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

        bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

        ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages
        """


# Combine the kraken annotations and produce count table
# rule kraken_summarize:
#     input:
#         script = "code/merge_bracken.R",
#         kraken_results = expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_bracken.txt", database = ["refseq","gtdb"], sample = metaG_samples),
#         combined_tax_info = rules.kraken_database_tax_merge.output.combined_tax_info
#     output:
#         counts = "data/sample_data/bracken_counts.tsv",
#         rel_abund = "data/sample_data/bracken_rel_abund.tsv"
#     resources: cpus=1, mem_mb=5000, time_min=60
#     container: "docker://eandersk/r_microbiome"
#     shell:
#         """
#         ./{input.script} --taxonomy={input.combined_tax_info} --counts-out={output.counts} --rel-out={output.rel_abund}
#         """


#####################
### FASTP Kraken ####
#####################

rule kraken2_load_gtdb_DB:
    input:
        db = "data/reference/kraken_databases/gtdb_r202"
    output:
        #db = service("/dev/shm/gtdb_r202"),
        db = temp(directory("/dev/shm/gtdb_r202")),
        #loaded = service("/tmp/gtdb_copied")
        #temp(directory("/dev/shm/gtdb_r202"))
    resources: cpus=1, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        mkdir {output.db}
        cp -r {input.db} /dev/shm/
        touch /tmp/gtdb_copied
        """

rule kraken2_gtdb_w_uniq_fastp:
    input:
        f_seq = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        # f_seq = "data/omics/{sample_type}/{sample}/reads/raw_fwd_reads.fastq.gz",
        # r_seq = "data/omics/{sample_type}/{sample}/reads/raw_rev_reads.fastq.gz",
        db = rules.kraken2_load_gtdb_DB.output.db,
        #db_loaded = rules.kraken2_load_gtdb_DB.output.loaded,
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_report.txt",
        out = temp("data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_out.txt"),
        bracken = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_bracken.txt",
        bracken_report = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_brackenMpa.txt",
        unclass_f = temp("data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_unclassified_1.fasta"),
        unclass_r = temp("data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_unclassified_2.fasta"),
        bracken_input = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150,
        unclass_out = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_unclassified#.fasta"
    conda: "config/conda_yaml/kraken.yaml"
    benchmark: "benchmarks/kraken2_gtdb_w_uniq_fastp/{sample_type}-{sample}.txt"
    resources: cpus=16, mem_mb=25000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --memory-mapping \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            --db {input.db} \
            --minimum-hit-groups 3 \
            --unclassified-out {params.unclass_out} \
            --paired {input.f_seq} {input.r_seq}

        echo "Kraken complete, filtering..."

        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

        echo "Filtering complete, running bracken..."

        bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

        ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages

        echo "Bracken complete. Quitting."
        """

rule kraken2_load_refseq_DB:
    input:
        db = "data/reference/kraken_databases/refseq"
    output:
        #db = service("/dev/shm/refseq"),
        db = temp(directory("/dev/shm/refseq")),
        #loaded = service("/tmp/refseq_copied")
        #temp(directory("/dev/shm/refseq"))
    resources: cpus=1, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        mkdir {output.db}
        cp -r {input.db} /dev/shm/
        touch /tmp/refseq_copied
        """

# Any reads not annotated with the GTDB database are then annotated with a RefSeq database
rule kraken2_refseq_w_uniq_fastp: ##Run kraken2
    input:
        f_seq = rules.kraken2_gtdb_w_uniq_fastp.output.unclass_f,
        r_seq = rules.kraken2_gtdb_w_uniq_fastp.output.unclass_r,
        db = rules.kraken2_load_refseq_DB.output.db,
        #db_loaded = rules.kraken2_load_refseq_DB.output.loaded,
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_report.txt",
        out = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_out.txt",
        bracken = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_bracken.txt",
        bracken_report = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_brackenMpa.txt",
        bracken_input = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150
    benchmark: "benchmarks/kraken2_refseq_w_uniq_fastp/{sample_type}-{sample}.txt"
    conda: "config/conda_yaml/kraken.yaml"
    resources: cpus=16, mem_mb=25000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --memory-mapping \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            --db {input.db} \
            --minimum-hit-groups 3 \
            --paired {input.f_seq} {input.r_seq}

        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > {output.bracken_input}

        bracken -d {input.db} -i {output.bracken_input} -o {output.bracken} -w {output.bracken_report}

        ./{input.kreport2mpa} -r {output.bracken_report} -o {output.bracken_mpa} --percentages
        """


rule bracken_metacodeR:
        input:
            script = "code/plot_metacoder_single_sample.R",
            bracken_refseq = "data/omics/{sample_type}/{sample}/kraken_fastp/refseq_{sample}_bracken.txt",
            bracken_gtdb = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_bracken.txt",
            tax_ref = "data/reference/kraken_tax_info_merged.tsv"
        output: "data/omics/{sample_type}/{sample}/kraken_fastp/{sample}_braken_metacodeR.pdf"
        resources: cpus=1, mem_mb=8000, time_min=60
        container: "docker://eandersk/r_microbiome"
        priority: 3
        shell:
            """
            {input.script} \
                --abund-refseq={input.bracken_refseq} \
                --abund-gtdb={input.bracken_gtdb} \
                --sample={wildcards.sample} \
                --tax-ref={input.tax_ref} \
                --output={output}
            """

rule contig_abund_metacodeR:
        input:
            contig_abund = "data/omics/{sample_type}/{sample}/{sample}_lca_abund_summarized.tsv"
        output: "data/omics/{sample_type}/{sample}/{sample}_lca_abund_metacoder.pdf"
        resources: cpus=1, mem_mb=8000, time_min=60
        container: "docker://eandersk/r_microbiome"
        priority: 3
        shell:
            """
            code/plot_contig_abund_uniref_LCA_single_sample.R \
                --abund={input.contig_abund} \
                --sample={wildcards.sample} \
                --output={output}
            """

rule reads_unirefLCA_mmseqs:
    input:
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        report = "data/omics/{sample_type}/{sample}/{sample}_report"
    conda:  "config/conda_yaml/mmseqs.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/{sample}",
        tmp_dir = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/tmp",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        #tmp_dir = "/tmp/kiledal/mmseqs2/{sample}",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/reads_unirefLCA_mmseqs/{sample_type}-{sample}.txt"
    log: "logs/reads_unirefLCA_mmseqs/{sample_type}-{sample}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        #mem_mb = 160000, cpus=48, time_min=20000
        cpus=32, mem_mb=170000, time_min=19440
    shell:
        """
        export TMPDIR={params.tmp_dir}
        
        # Different tmp dir if running on lab servers
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/mmseqs/{wildcards.sample}/tmp

        mkdir -p $TMPDIR

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        # Previous copied reads to temp space, but no longer think this is nescessary. Should be deleted.
        #cp {input.fwd_reads} {params.tmp_fwd_reads}
        #cp {input.rev_reads} {params.tmp_rev_reads}

        #mmseqs touchdb {params.unirefDB} # loads  the database into memory, only use in large memory systems / nodes
        
        mmseqs \
            easy-taxonomy \
            {input.fwd_reads} {input.rev_reads} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            $TMPDIR \
            --lca-mode 3 \
            --orf-filter 1 \
            --orf-filter-s 3.5 \
            -s 4 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit $((({resources.mem_mb} * 8/10) / 1024))G \
            2>&1 | tee -a {log}

            #{params.tmp_fwd_reads} {params.tmp_rev_reads} \
            # --db-load-mode 2 # this loads the database into memory, was an attempt to speed up on Great Lakes but limited memory & scratch space were limiting factors

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r {params.tmp_fwd_reads} {params.tmp_rev_reads} $TMPDIR 2>&1 | tee -a {log}
        rm -r $TMPDIR 2>&1 | tee -a {log}
        printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """

rule add_lineage_to_unirefLCAtax:
    input:
        db = ancient("data/reference/ncbi_tax"),
        report = "data/omics/{sample_type}/{sample}/{sample}_report"
    output: 
        inspect_w_lineage = "data/omics/{sample_type}/{sample}/{sample}_report_w_full_lineage",
        inspect_w_7_lev_lineage = "data/omics/{sample_type}/{sample}/{sample}_report_w_standardized_lineage"
    conda: "config/conda_yaml/taxonkit.yaml"
    resources: cpus=1, mem_mb=10000, time_min=5440, mem_gb = 10
    log: "logs/unirefLCA_mmseqs_add_lineage/{sample_type}-{sample}.log"
    shell:
        """
        taxonkit lineage \
            {input.report} \
            --data-dir {input.db} \
            -i 5 \
            -o {output.inspect_w_lineage} 2>&1 | tee {log}

        taxonkit reformat \
            {output.inspect_w_lineage} \
            --data-dir {input.db} \
            -i 7 \
            -P \
            -o {output.inspect_w_7_lev_lineage} 2>&1 | tee -a {log}
        """
    
rule contig_unirefLCA_mmseqs:
    input:
        contigs = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa"
    output:
        report = "data/omics/{sample_type}/{sample}/{sample}_contig_report",
        lca = "data/omics/{sample_type}/{sample}/{sample}_contig_lca.tsv"
    conda:  "config/conda_yaml/mmseqs.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/{sample}_contig",
        tmp_dir = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/tmp",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        #tmp_dir = "/tmp/kiledal/mmseqs2/{sample}_contigs",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/contig_unirefLCA_mmseqs/{sample_type}-{sample}.txt"
    log: "logs/contig_unirefLCA_mmseqs/{sample_type}-{sample}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        mem_mb = 160000, cpus=32, time_min=7200
        #mem_mb = 160000, cpus=32, time_min=500
    shell:
        """
        export TMPDIR={params.tmp_dir}
        
        # Different tmp dir if running on lab servers
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/mmseqs2/{wildcards.sample}_contigs
        [[ "${{HOSTNAME}}" == "alpena" ]] && export TMPDIR=/old-geomicro/kiledal-extra/mmseqs_contig/$USER/mmseqs2/{wildcards.sample}_contigs
        
        mkdir -p $TMPDIR

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #mmseqs touchdb {params.unirefDB} # loads the database into memory, only use in large memory systems / nodes

        mmseqs \
            easy-taxonomy \
            {input.contigs} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            $TMPDIR \
            --lca-mode 3 \
            --orf-filter 1 \
            --orf-filter-s 3.5 \
            -s 4 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit $(echo "scale=0;({resources.mem_mb}*0.8)/1024" | bc)G \
            --db-load-mode 1 \
            2>&1 | tee -a {log}

            # --db-load-mode 2 # this loads the database into memory, was an attempt to speed up on Great Lakes but limited memory & scratch space were limiting factors

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r {params.tmp_fwd_reads} {params.tmp_rev_reads} $TMPDIR 2>&1 | tee -a {log}
        rm -r $TMPDIR 2>&1 | tee -a {log}
        printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """

rule tax_abund_summary_from_contigs:
    input: 
        mmseqs_report = "data/omics/{sample_type}/{sample}/{sample}_contig_report",
        script = "code/tax_abund_from_contigs.R",
        contig_abund = "data/omics/{sample_type}/{sample}/{sample}_contig_abund.tsv",
        lca = "data/omics/{sample_type}/{sample}/{sample}_contig_lca.tsv"
        #assembly_done = "data/omics/{sample_type}/{sample}/assembly/megahit/.done"
    output:
        abund_summary = "data/omics/{sample_type}/{sample}/{sample}_lca_abund_summarized.tsv"
    params:
        taxonkit_path = "code/dependencies/taxonkit",
        taxdump = "data/reference/ncbi_tax"
    benchmark: "benchmarks/tax_abund_summary_from_contigs/{sample_type}-{sample}.txt"
    container: "docker://eandersk/r_microbiome"
    resources: cpus = 24, time_min=1000, mem_mb = 150000
    priority: 2
    shell:
        """
        pwd #check that the proper working directory is being used
        
        ./{input.script} \
            -l {input.lca} \
            -r {input.contig_abund} \
            -o {output.abund_summary} \
            -c {resources.cpus} \
            -t {params.taxonkit_path} \
            -d {params.taxdump}
        """


find_mmseq_reports = """
for filepath in data/omics/metagenomes/*/; do
  sample=$(basename $filepath)
  if [ -f "$filepath/${sample}_report" ]; then
        echo -n $sample" "
  fi
done
"""

rule kofam_scan:
    input:
        genes = "data/omics/{sample_type}/{sample}/proteins/{sample}_PROTEINS.faa",
        profile = "data/reference/kegg/kofamscan/profiles",
        ko_list = "data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/omics/{sample_type}/{sample}/kofam_scan/{sample}_kofam_results.txt"
    conda: "config/conda_yaml/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofamscan/{sample_type}-{sample}.txt"
    log: "logs/kofamscan/{sample_type}-{sample}.log"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --format=detail-tsv \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.sample}_kofamscan \
            --ko-list {input.ko_list} {input.genes} | tee {log}
        """
        
rule annotate_contigs:
    input:
        genes = rules.prodigal.output.genes,
        gene_alignment = rules.align_to_uniref.output.gene_uniref_alignment,
        gene_rpkm = rules.calc_gene_abundance.output.reads_vs_genes_rpkm,
        contig_rpkm = rules.calc_gene_abundance.output.reads_vs_contigs_rpkm,
        #contigs = rules.MEC.output.corrected_assembly,
        contigs = rules.rename_contigs.output.contigs,
        script = "code/Strain-Level_Metagenome_Analysis/AnnotateContigs.pl",
        UMRAD = "data/reference/UMRAD"
    output:
        done = touch("data/omics/{sample_type}/{sample}/.annotation_done"),
    params:
        annotation_dir = directory("data/omics/{sample_type}/{sample}/annotation")
    conda: "config/conda_yaml/main.yaml"
    resources: cpus=1, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 170000
    log: "logs/annotate_contigs/{sample_type}-{sample}.log"
    benchmark: "benchmarks/annotate_contigs/{sample_type}-{sample}.log"
    shell:
        """
        rm -r {params.annotation_dir}
        mkdir -p {params.annotation_dir}
        ln {input.genes} {input.gene_alignment} {input.gene_rpkm} {input.contig_rpkm} {input.contigs} {params.annotation_dir}/
        mv  {params.annotation_dir}/$(basename {input.contigs}) {params.annotation_dir}/{wildcards.sample}_MCDD.fa
        proj_dir=$(pwd)

        cp {input.script} {params.annotation_dir}/
        ln {input.UMRAD}/* {params.annotation_dir}/

        cd {params.annotation_dir}
        perl AnnotateContigs.pl -s {wildcards.sample} 2>&1 | tee $proj_dir/{log}
        
        #-d ${{proj_dir}}/{input.UMRAD}
        #> ${{proj_dir}}/{log}
        """


rule download_gtdbtk_refs:
    output:
        dir = directory("/geomicro/data2/kiledal/references/gtdbtk"),
        tar = "/geomicro/data2/kiledal/references/gtdbtk/gtdbtk_data.tar.gz"
    resources: cpus=1, mem_mb=8000, time_min=2880, mem_gb = 8
    shell:
        """
        mkdir -p {output.dir}
        cd {output.dir}
        #wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
        cp ../backup_gtdbtk_data.tar.gz $PWD/gtdbtk_data.tar.gz
        tar xzf gtdbtk_data.tar.gz
        """


rule drep_mag_coverage:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev,
        drep_done = "data/projects/{project}/{sample_type}/{sample}/bins/.drep_done"
    output: "data/projects/{project}/{sample_type}/{sample}/bins/coverage_drep_bins.tsv"
    params:
        drep_bins = "data/projects/{project}/{sample_type}/{sample}/bins/drep/dereplicated_genomes"
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=24, mem_mb=150000, time_min=2880
    shell:
        """
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/
        
        coverm genome \
            -t {resources.cpus} \
            -m relative_abundance mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
            --output-format sparse \
            --min-covered-fraction 0 \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            --genome-fasta-files {params.drep_bins}/*.fa \
            -o {output}
        """

rule humann:
    input:
        # Sample data #
        f_seq = rules.remove_contaminants.output.cleaned_fwd,
        r_seq = rules.remove_contaminants.output.cleaned_rev,
        bracken_mpa = "data/omics/{sample_type}/{sample}/kraken/gtdb_{sample}_brackenMpa.txt",
        # Reference database #
        NUC_DB = "data/reference/humann/genome_reps_filt_annot.fna.gz",
        PROT_DB = "data/reference/humann/protein_database/uniref90_201901b.dmnd",
        NUC_fol = "data/reference/humann/",
        PROT_fol = "data/reference/humann/protein_database/"
    output:
        humann_output = directory("data/omics/{sample_type}/{sample}/humann"),
        concat_unzipped_reads = temp("data/omics/{sample_type}/{sample}/reads/for_humann.fastq")
    params:
        mem_use = "maximum"
    conda: "config/conda_yaml/humann.yaml"
    log: "logs/humann/{sample_type}-{sample}.log"
    benchmark: "benchmarks/humann/{sample_type}-{sample}.txt"
    #resources: cpus=36, mem_mb=175000, time_min=20160
    resources: cpus=36, time_min=20160, partition = "largemem", mem_mb = lambda wildcards, attempt: attempt * 250000
    shell:
        """
        #Humann needs non-compressed fastqs, and forward and reverse files should be concatenated
        gunzip -c {input.f_seq} > {output.concat_unzipped_reads}
        gunzip -c {input.r_seq} >> {output.concat_unzipped_reads}

        echo "Combined and uncompressed fastq files.\n"

        #proj_dir=$PWD

        humann3 --bypass-nucleotide-index \
            --threads {resources.cpus} \
            --memory-use {params.mem_use} \
            --nucleotide-database {input.NUC_fol} \
            --protein-database {input.PROT_fol} \
            --taxonomic-profile {input.bracken_mpa} \
            --input {output.concat_unzipped_reads} \
            --output-basename {wildcards.sample}_humann \
            --output {output.humann_output}/ > {log}
        """

rule humann_fastp:
    input:
        # Sample data #
        f_seq = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        bracken_mpa = "data/omics/{sample_type}/{sample}/kraken_fastp/gtdb_{sample}_brackenMpa.txt",
        # Reference database #
        NUC_DB = humann_ref_dir + "/genome_reps_filt_annot.fna.gz",
        PROT_DB = humann_ref_dir + "/protein_database/uniref90_201901b.dmnd",
        NUC_fol = humann_ref_dir + "/",
        PROT_fol = humann_ref_dir +"/protein_database/"
    output:
        humann_output = directory("data/omics/{sample_type}/{sample}/humann_fastp"),
        concat_unzipped_reads = temp("data/omics/{sample_type}/{sample}/reads/for_humann_fastp.fastq")
    params:
        mem_use = "maximum"
    conda: "config/conda_yaml/humann.yaml"
    log: "logs/humann/{sample_type}-{sample}_fastp.log"
    benchmark: "benchmarks/humann/{sample_type}-{sample}_fastp.txt"
    #resources: cpus=36, mem_mb=175000, time_min=20160
    #resources: cpus=36, time_min=20160, partition = "largemem", mem_mb = lambda wildcards, attempt: attempt * 250000
    resources: cpus=36, time_min=4320, mem_mb = 175000
    shell:
        """
        #Humann needs non-compressed fastqs, and forward and reverse files should be concatenated
        gunzip -c {input.f_seq} > {output.concat_unzipped_reads}
        gunzip -c {input.r_seq} >> {output.concat_unzipped_reads}

        echo "Combined and uncompressed fastq files.\n"

        #proj_dir=$PWD

        humann3 --bypass-nucleotide-index \
            --threads {resources.cpus} \
            --memory-use {params.mem_use} \
            --nucleotide-database {input.NUC_fol} \
            --protein-database {input.PROT_fol} \
            --taxonomic-profile {input.bracken_mpa} \
            --input {output.concat_unzipped_reads} \
            --output-basename {wildcards.sample}_humann \
            --output {output.humann_output}/ > {log}

        rm -r {output.humann_output}/{wildcards.sample}_humann_humann_temp
        """


rule sourmash_sketch:
    input:
       f_seq = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
       r_seq = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
       sig = "data/omics/{sample_type}/{sample}/sourmash/{sample}.sig"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash_sketch/{sample_type}-{sample}.log"
    benchmark: "benchmarks/sourmash_sketch/{sample_type}-{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 20000
    shell:
        """
        sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge {wildcards.sample} -o {output.sig} {input.f_seq} {input.r_seq} 2>&1 | tee {log}
        """

rule sourmash_gather:
    input:
       sig = "data/omics/{sample_type}/{sample}/sourmash/{sample}.sig",
       #gtdb_refDB = "data/reference/sourmash/gtdb-rs214-k31.zip",
       gtdb_refDB = "data/reference/sourmash/gtdb-rs207.dna.k31.zip",
       #genbank_fungi = "data/reference/sourmash/genbank-2022.03-fungi-k31.zip",
       #genbank_protozoa = "data/reference/sourmash/genbank-2022.03-protozoa-k31.zip",
       #genbank_viral = "data/reference/sourmash/genbank-2022.03-viral-k31.zip",
       microcystis_refDB = "data/reference/sourmash/Microcystis_sigs.sig",
       taxDB = "data/reference/sourmash/gtdb_and_Microcystis_tax.db"
    output:
       reps = "data/omics/{sample_type}/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.csv",
       tax = "data/omics/{sample_type}/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.with-lineages.csv"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash/{sample_type}-{sample}.log"
    benchmark: "benchmarks/sourmash/{sample_type}-{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 80000
    shell:
        """
        sourmash gather {input.sig} {input.gtdb_refDB} {input.microcystis_refDB} -o {output.reps} 2>&1 | tee -a {log}

        sourmash tax annotate -g {output.reps} -t {input.taxDB} 2>&1 | tee -a {log}

        mv $(basename {output.tax}) data/omics/{wildcards.sample_type}/{wildcards.sample}/sourmash/
        """

rule sourmash_gather_Microcystis:
    input:
       sig = "data/omics/{sample_type}/{sample}/sourmash/{sample}.sig",
       refDB = "data/reference/sourmash/Microcystis_sigs.sig"
    output:
       reps = "data/omics/{sample_type}/{sample}/sourmash/{sample}_gather_Microcystis.csv",
       #tax = "data/omics/{sample_type}/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.with-lineages.csv"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash/{sample_type}-{sample}.log"
    benchmark: "benchmarks/sourmash/{sample_type}-{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 20000
    shell:
        """
        sourmash gather {input.sig} {input.refDB} -o {output.reps} 2>&1 | tee -a {log}

        #mv $(basename {output.reps}) data/omics/{wildcards.sample_type}/{wildcards.sample}/sourmash/
        """

## New binning


rule link_reads_w_sample_names:
    input: 
        reads = "data/omics/{sample_type}/{sample}/reads/decon_{dir}_reads_fastp.fastq.gz"
        #"data/projects/{project}/{sample_type}/{sample}/reads/decon_{dir}_reads_fastp.fastq.gz"
    output: 
        touch = touch("data/omics/{sample_type}/{sample}/reads/.linked_{dir}_w_sample_name"),
        linked_reads = "data/omics/{sample_type}/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz"
    resources: cpus=1, mem_mb = 500
    shell:
        """
        mkdir -p $(dirname {output.linked_reads})
        cd $(dirname {output.linked_reads})
        ln -sf ../$(basename {input.reads}) $(basename {output.linked_reads})
        """


rule map_to_contigs:
    input:
        #expand("data/projects/{project}/{sample_type}/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz", sample = metaG_samples, dir = ["fwd", "rev"], project = "2022_geomicro_JGI_CSP"),
        #expand("data/projects/{project}/{sample_type}/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz", sample = metaG_samples, dir = ["fwd", "rev"], project = current_project, sample_type = "metagenomes"),
        #reads_linked = "data/omics/{sample_type}/{sample}/reads/.linked_w_sample_name",
        contigs = rules.rename_contigs.output.contigs
    output: 
        #bam_dir = temp(directory("data/projects/{project}/{sample_type}/{sample}/bins/bam"))
        bam_dir = temp(directory(config["binning_bam_dir"]))
    params:
        mapper = "minimap2-sr",
        #mapper = "strobealign"
    conda: "config/conda_yaml/coverm.yaml"
    # benchmark: "benchmarks/map_to_contigs/{sample_type}-{project}__{sample}.txt"
    resources: cpus=32, mem_mb=150000, time_min=7200, disk_mb=500000, scratch_disk_mb = 1000000
    shell:
        """
        #ulimit -u 1000000 # Increase the number of open files to avoid errors

        # Was running out of space on GL when using the local /tmp drives on each node, so use the /scratch space instead
        export TMPDIR={output.bam_dir}/coverm_tmp
        
        # Different tmp dir if running on lab servers
        #[[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/
        
        # Ensure directory exists
        mkdir -p $TMPDIR

        coverm make -c data/projects/{wildcards.project}/metagenomes/*/reads/fastp_decon/*.fastq.gz \
            -r {input.contigs} \
            --discard-unmapped \
            -t {resources.cpus} \
            --mapper {params.mapper} \
            -o {output.bam_dir} 

        #mkdir -p {output.bam_dir}
        #rsync -av -P kiledal@greatlakes-xfer.arc-ts.umich.edu:/scratch/gdick_root/gdick2/kiledal/GLAMR/binning/bams/{wildcards.project}/{wildcards.sample_type}/{wildcards.sample} /scratch/gdick_root/gdick2/kiledal/GLAMR/binning/bams/{wildcards.project}/{wildcards.sample_type}
        """


rule contig_coverage:
    input:
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam"
        bam_dir = config["binning_bam_dir"]
    output: 
        coverage = "data/projects/{project}/{sample_type}/{sample}/bins/contig_coverage.tsv",
        coverage_full = "data/projects/{project}/{sample_type}/{sample}/bins/contig_coverage_all_stats.tsv",
        coverage_metabat = "data/projects/{project}/{sample_type}/{sample}/bins/metabat_style_contig_coverage.tsv"
    params:
        tmpdir = "tmp/coverm_contig_coverage/{sample}"
    benchmark: "benchmarks/contig_coverage/{sample_type}-{project}__{sample}.txt"
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=24, mem_mb=120000, time_min=10000 # standard assemblies
    #resources: cpus=24, mem_mb=1000000, time_min=14400, partition = "largemem" # coassembly
    #resources: cpus=24, mem_mb=1500000, time_min=14400, partition = "largemem" # XL coassembly
    priority: 2
    shell:
        """
        export TMPDIR={params.tmpdir}
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER
        mkdir -p $TMPDIR

        coverm contig \
            -b {input.bam_dir}/*.bam \
            -t {resources.cpus} \
            --output-file {output.coverage}

        coverm contig \
            -b {input.bam_dir}/*.bam \
            -t {resources.cpus} \
            -m mean trimmed_mean covered_bases variance length count reads_per_base rpkm tpm \
            --output-file {output.coverage_full}

        coverm contig \
            -b {input.bam_dir}/*.bam \
            -t {resources.cpus} \
            --methods metabat \
            --output-file {output.coverage_metabat}

        rm -r {params.tmpdir}
        """


rule index_contig_coverage:
    input:
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam"
        bam_dir = config["binning_bam_dir"]
    output: 
        index_done = touch("data/projects/{project}/{sample_type}/{sample}/bins/.bam_indexed")
    params:
        #bam = "data/projects/{project}/{sample_type}/{sample}/bins/bam/final.contigs.renamed.fa.decon_fwd_reads_fastp.fastq.gz.bam"
        bam = config["binning_bam_dir"] + "/final.contigs.renamed.fa.decon_fwd_reads_fastp.fastq.gz.bam"
    benchmark: "benchmarks/index_contig_coverage/{sample_type}-{project}__{sample}.txt"
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=4, mem_mb=60000, time_min=2880
    priority: 2
    shell:
        """
        parallel -j {resources.cpus} samtools index -@ 1 ::: {input.bam_dir}/*.bam
        """

rule concoct:
    input:
        contigs = rules.rename_contigs.output.contigs,
        bam_index = rules.index_contig_coverage.output.index_done,
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam"
        bam_dir = config["binning_bam_dir"] # This is the default
    output:
        cut_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/CONCOCT/cut_contigs_10K.fa",
        cut_contigs_bed = "data/projects/{project}/{sample_type}/{sample}/bins/CONCOCT/contigs_10K.bed",
        cut_coverage = "data/projects/{project}/{sample_type}/{sample}/bins/CONCOCT/cut_coverage_table.tsv"
    params:
        outdir = "data/projects/{project}/{sample_type}/{sample}/bins/CONCOCT/output",
        #bam = "data/projects/{project}/{sample_type}/{sample}/bins/bam/*.bam"
        bam = config["binning_bam_dir"] + "/*.bam" # This is the default
        #bam = "/ssd/GLAMR/binning/bams/bams/{project}/{sample_type}/{sample}/*.bam" #changed just for Paul's coassembly
    benchmark: "benchmarks/concoct/{sample_type}-{project}__{sample}.txt"
    conda: "config/conda_yaml/concoct.yaml"
    resources: cpus=24, mem_mb=150000, time_min=10080, mem_gb = 50 # standard samples
    #resources: cpus=24, mem_mb=170000, time_min=10080, mem_gb = 50 # coassembly
    #resources: cpus=32, mem_mb=1200000, time_min=18720, partition = "largemem" # XLcoassembly
    priority: 3
    shell:
        """
        export OMP_THREAD_LIMIT={resources.cpus}

        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output.cut_contigs_bed} > {output.cut_contigs}
        
        concoct_coverage_table.py {output.cut_contigs_bed} {params.bam} > {output.cut_coverage}

        concoct --threads {resources.cpus} --composition_file {output.cut_contigs} --coverage_file {output.cut_coverage} -b {params.outdir}/
        
        merge_cutup_clustering.py {params.outdir}/clustering_gt1000.csv > {params.outdir}/clustering_merged.csv

        mkdir -p {params.outdir}/fasta_bins
        extract_fasta_bins.py {input.contigs} {params.outdir}/clustering_merged.csv --output_path {params.outdir}/fasta_bins
        """


rule metabat2:
    input:
        contigs = rules.rename_contigs.output.contigs,
        bam_index = rules.index_contig_coverage.output.index_done,
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam",
        bam_dir = config["binning_bam_dir"],
        coverm_depth = "data/projects/{project}/{sample_type}/{sample}/bins/metabat_style_contig_coverage.tsv"
    output:
        #depth = "data/omics/{sample_type}/{sample}/bins/jgi_depth_summary.txt",
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/METABAT2/.done")
    params:
        bin_name = directory("data/projects/{project}/{sample_type}/{sample}/bins/METABAT2/metabat2")
    benchmark: "benchmarks/metabat2/{sample_type}-{project}__{sample}.txt"
    singularity: "docker://metabat/metabat"
    resources: cpus=16, mem_mb=20000, time_min=2880 # standard samples
    #resources: cpus=36, mem_mb=150000, time_min=5880 # coassembly
    priority: 3
    shell:
        """
        pwd 
        cd {current_dir}
        pwd

        metabat2 -i {input.contigs} -a {input.coverm_depth} -o {params.bin_name} -m 2000 -t {resources.cpus} --unbinned
        """

        #jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam_dir}/*.bam

rule maxbin2_coverage:
    input:
        script = "code/create_maxbin_coverage.R",
        coverm_depth = "data/projects/{project}/{sample_type}/{sample}/bins/metabat_style_contig_coverage.tsv"
    output:
        depths_file = "data/projects/{project}/{sample_type}/{sample}/bins/maxbin/depths.txt"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb=50000, time_min=1000 # standard samples
    #resources: cpus=1, mem_mb=120000, time_min=2000 # coassembly
    priority: 3
    shell:
        """
        cd {current_dir}
        pwd

        ./{input.script} {input.coverm_depth}
        """

rule maxbin2:
    input:
        contigs = rules.rename_contigs.output.contigs,
        depth = rules.maxbin2_coverage.output.depths_file
    output:
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/maxbin/.done")
    params:
        bin_dir = "data/projects/{project}/{sample_type}/{sample}/bins/maxbin/maxbin"
    benchmark: "benchmarks/maxbin/{sample_type}-{project}__{sample}.txt"
    conda: "config/conda_yaml/maxbin.yaml"
    resources: cpus=16, mem_mb=20000, time_min=10080 # standard samples
    #resources: cpus=16, mem_mb=80000, time_min=20130 # coassembly
    priority: 3
    shell:
        """
        pwd 
        cd {current_dir}
        pwd
        
        run_MaxBin.pl -contig {input.contigs} \
            -markerset 107 \
            -thread {resources.cpus} \
            -min_contig_length 2000\
	        -out {params.bin_dir} \
	        -abund_list {input.depth}
        """

rule semibin_ref_download:
    output: directory("data/reference/semibin/gtdb")
    params:
    conda: "config/conda_yaml/semibin.yaml"
    resources: cpus =1, mem_mb = 2000, time_min=2880
    shell:
        """
        SemiBin download_GTDB --reference-db {output}
        """


rule semibin:
    input:
        #"data/reference/semibin/.gtdb_downloaded",
        ref_db = "data/reference/semibin/gtdb",
        contigs = rules.rename_contigs.output.contigs,
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam"
        bam_dir = config["binning_bam_dir"]
    output:
        out_dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/semibin"),
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/semibin/.done")
    params:
    conda: "config/conda_yaml/semibin.yaml"
    benchmark: "benchmarks/semibin/{sample_type}-{project}__{sample}.txt"
    log: "logs/semibin/{sample_type}-{project}__{sample}.log"
    resources: cpus=16, mem_mb=170000, time_min=2880, mem_gb = 50 # standard samples
    #resources: cpus=32, mem_mb=1250000, time_min=2880, partition = "largemem" # coassembly
    #resources: cpus=32, mem_mb=700000, time_min=2880, partition = "largemem" # coassembly on our servers
    priority: 3
    shell:
        r"""
        WORK_DIR=$PWD

        SemiBin \
            single_easy_bin \
            --threads {resources.cpus} \
            --reference-db {input.ref_db}/gtdb \
            -i {input.contigs} \
            -b {input.bam_dir}/*.bam \
            -o {output.out_dir} | tee {log}

            #  --environment mouse_gut ## can only be used for single sample binning

        # Fix permissions on outputs (by default only readable by owner)
        find {output.out_dir} -type f -exec chmod g+r,o+r {{}} \; -o -type d -exec chmod g+rx,o+rx {{}} \;
        """

rule semibin2:
    input:
        #"data/reference/semibin/.gtdb_downloaded",
        ref_db = "data/reference/semibin/gtdb",
        contigs = rules.rename_contigs.output.contigs,
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam"
        bam_dir = config["binning_bam_dir"]
    output:
        out_dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/semibin2"),
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/semibin2/.done")
    params:
    conda: "config/conda_yaml/semibin2.yaml"
    benchmark: "benchmarks/semibin2/{sample_type}-{project}__{sample}.txt"
    log: "logs/semibin2/{sample_type}-{project}__{sample}.log"
    resources: cpus=36, mem_mb=170000, time_min=2880 # standard samples
    #resources: cpus=32, mem_mb=1250000, time_min=2880, partition = "largemem" # coassembly
    #resources: cpus=32, mem_mb=700000, time_min=2880, partition = "largemem" # coassembly on our servers
    priority: 3
    shell:
        r"""
        WORK_DIR=$PWD

        SemiBin2 single_easy_bin \
            --threads {resources.cpus} \
            --compression none \
            -i {input.contigs} \
            -b {input.bam_dir}/*.bam \
            -o {output.out_dir} | tee {log}

        # Fix permissions on outputs (by default only readable by owner)
        find {output.out_dir} -type f -exec chmod g+r,o+r {{}} \; -o -type d -exec chmod g+rx,o+rx {{}} \;
        """

rule VAMB:
    input:
        contigs = rules.rename_contigs.output.contigs,
        #coverm_depth = "data/projects/{project}/{sample_type}/{sample}/bins/metabat_style_contig_coverage.tsv"
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam"
        bam_dir = config["binning_bam_dir"]
    output:
        outdir = directory("data/projects/{project}/{sample_type}/{sample}/bins/VAMB")
    #conda: "config/conda_yaml/VAMB.yaml"
    params:
        #bam_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bam" # Should be switched to input after fixing samples that need to be fixed
        bam_dir = config["binning_bam_dir"]
    conda: "vamb"
    #conda: "/home/kiledal/miniconda3/envs/vamb"
    benchmark: "benchmarks/VAMB/{sample_type}-{project}__{sample}.txt"
    log: "logs/VAMB/{sample_type}-{project}__{sample}.log"
    resources: cpus=1, mem_mb=40000, time_min=5000, partition = "gpu", gpu = 1 # standard samples
    #resources: cpus=1, mem_mb=120000, time_min=14400, partition = "gpu", gpu = 1 # coassembly
    priority: 3
    shell:
        """
        vamb -o _ \
            --outdir {output.outdir} \
            --fasta {input.contigs} \
            --bamfiles {params.bam_dir}/*.bam \
            --minfasta 200000 \
            --model vae-aae \
            --cuda
        """

rule format_coverage_for_metadecoder:
    input:
        script = "code/make_metadecoder_coverage.R",
        coverage = "data/projects/{project}/{sample_type}/{sample}/bins/contig_coverage.tsv",
        contigs = rules.rename_contigs.output.contigs
    output: "data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/coverage.tsv"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb = 50000, time_min=360 # standard samples
    #resources: cpus=1, mem_mb = 140000, time_min=2000 # coassembly
    priority: 3
    shell:
        """
        pwd && cd {current_dir} && pwd

        ./{input.script} --coverage={input.coverage} --contigs={input.contigs} --out={output}
        """
        
rule metadecoder:
    input:
        contigs = rules.rename_contigs.output.contigs,
        coverage = "data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/coverage.tsv"
    output:
        touch("data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/.done"),
        seed = "data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/seed.txt",
        bins_dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/bins")
    params:
        out_prefix = "data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/bins/{sample}"
    #conda: "config/conda_yaml/VAMB.yaml"
    conda: "metadecoder"
    shadow: "minimal"
    benchmark: "benchmarks/metadecoder/{sample_type}-{project}__{sample}.txt"
    log: "logs/metadecoder/{sample_type}-{project}__{sample}.log"
    resources: cpus=1, mem_mb=150000, time_min=10080, partition = "gpu", gpu = 1 # standard samples
    #resources: cpus=1, mem_mb=160000, time_min=15080, partition = "gpu", gpu = 1 # coassembly
    priority: 3
    shell:
        """
        mkdir -p {output.bins_dir}
        
        metadecoder seed --threads {resources.cpus} -f {input.contigs} -o {output.seed}

        metadecoder cluster -f {input.contigs} -c {input.coverage} -s {output.seed}  -o {params.out_prefix} | tee {log}
        """


rule standardize_bins:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/CONCOCT/cut_contigs_10K.fa",
        #"data/projects/{project}/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        "data/projects/{project}/{sample_type}/{sample}/bins/semibin",
        "data/projects/{project}/{sample_type}/{sample}/bins/METABAT2/.done/",
        "data/projects/{project}/{sample_type}/{sample}/bins/maxbin/.done", 
        "data/projects/{project}/{sample_type}/{sample}/bins/VAMB",
        "data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/.done",
        script = ancient("code/standardize_bins.R"),
        #contig_info = "data/projects/{project}/{sample_type}/{sample}/assembly/megahit_noNORM/contigs_info.tsv"
        assembly = rules.rename_contigs.output.contigs,
        contig_info = rules.rename_contigs.output.contig_info
    output: 
        contig_bin_mapping = "data/projects/{project}/{sample_type}/{sample}/bins/contig_bins.rds",
        bin_info = "data/projects/{project}/{sample_type}/{sample}/bins/bins_prelim_info.tsv",
        bins_linked = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/.bins_linked"
    params: 
        sample = "{sample}",
        sample_dir = "data/projects/{project}/{sample_type}/{sample}"
    benchmark: "benchmarks/standardize_bins/{sample_type}-{project}__{sample}.txt"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb = 50000, time_min=360
    priority: 4
    shell:
        """
        pwd && cd {current_dir} && pwd

        ./{input.script} --sample_dir={params.sample_dir} --contig_info={input.contig_info}
        """

rule checkm_new_per_sample:
    input: "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/.bins_linked"
    output:
        results = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/checkm.txt"
    params:
        in_dir = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins",
        out_dir = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/checkm"
    benchmark: "benchmarks/checkm_new_per_sample/{sample_type}-{project}__{sample}.txt"
    conda: "config/conda_yaml/checkm.yaml"
    resources: cpus=16, mem_mb=80000, time_min=7200
    priority: 4
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
        """

rule checkm2:
    input: "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/.bins_linked"
    output:
        results = directory("data/projects/{project}/{sample_type}/{sample}/bins/checkm2")
    params:
        in_dir = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins",
        database = "data/reference/checkm2/CheckM2_database/uniref100.KO.1.dmnd",
        #out_dir = "data/projects/{project}/{sample_type}/{sample}/bins/checkm2"
    benchmark: "benchmarks/checkm2/{sample_type}-{project}__{sample}.txt"
    conda: "config/conda_yaml/checkm2.yaml"
    resources: cpus=16, mem_mb=80000, time_min=7200
    priority: 4
    shell:
        """
        # Database can be dowloaded with checkm2 database --download --path {{path}}

        checkm2 predict \
            --threads {resources.cpus} \
            --input {params.in_dir} \
            --database_path {params.database} \
            --output-directory {output.results}
        """


rule make_das_and_drep_inputs:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/.bins_linked",
        "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/checkm.txt",
        contig_bin_mapping = "data/projects/{project}/{sample_type}/{sample}/bins/contig_bins.rds",
        script = "code/make_das_and_drep_inputs.R"
    output: 
        drep_bin_info = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/genome_info.csv",
        drep_bins_linked = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/.bins_linked"
    params: 
        sample_dir = "data/projects/{project}/{sample_type}/{sample}",
        metabat2_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/metabat2_contigs.tsv",
        maxbin_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/maxbin_contigs.tsv",
        concoct_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/concoct_contigs.tsv",
        metadecoder_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/metadecoder_contigs.tsv",
        semibin_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/semibin_contigs.tsv",
        VAMB_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/VAMB_contigs.tsv"
    singularity: "docker://eandersk/r_microbiome"
    benchmark: "benchmarks/make_das_and_drep_inputs/{sample_type}-{project}__{sample}.txt"
    resources: cpus=1, mem_mb = 50000, time_min=1440
    priority: 4
    shell:
        """
        pwd && cd {current_dir} && pwd

        ./{input.script} --sample_dir={params.sample_dir}
        """

rule checkm_new:
    #input: "data/omics/{sample_type}/metagenome_bins/raw_combined_bins"
    output:
        #dir = temp(directory("data/omics/{sample_type}/metagenome_bins/raw_combined_bins")),
        results = "data/projects/{project}/{sample_type}/metagenome_bins/raw_combined_bins/checkm.txt"
    params:
        in_dir = "data/projects/{project}/{sample_type}/metagenome_bins/raw_combined_bins",
        out_dir = "data/projects/{project}/{sample_type}/metagenome_bins/raw_combined_bins/checkm"
    benchmark: "benchmarks/checkm_new/{sample_type}-{project}.txt"
    conda: "config/conda_yaml/checkm.yaml"
    resources: cpus=24, mem_mb=120000, time_min=2880
    priority: 4
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
        """


checkpoint dastool_new:
    input:
        contigs = rules.rename_contigs.output.contigs,
        checkm_res = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/checkm.txt"
    params:
        bin_folder = "data/omics/{sample_type}/{sample}/bins",
        das_prefix = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/output/{sample}",
        metabat2_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/metabat2_contigs.tsv",
        maxbin_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/maxbin_contigs.tsv",
        concoct_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/concoct_contigs.tsv",
        metadecoder_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/metadecoder_contigs.tsv",
        semibin_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/semibin_contigs.tsv",
        VAMB_contigs = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/VAMB_contigs.tsv"
    output: 
        #summary = "data/projects/{project}/{sample_type}/{sample}/bins/DASTool/_DASTool_summary.txt",
        das_done = touch("data/projects/{project}/{sample_type}/{sample}/bins/das_tool/.done"),
        #das_bins_folder = directory("data/omics/{sample_type}/{sample}/bins/das_tool/output/_DASTool_bins")
    conda: "config/conda_yaml/das_tool.yaml"
    benchmark: "benchmarks/dastool/{sample_type}-{project}__{sample}.txt"
    log: "logs/dastool/{sample_type}-{project}__{sample}.log"
    resources: cpus=8, mem_mb=50000, time_min=2880, mem_gb = 50
    priority: 4
    shell:
        """
        mkdir -p $(dirname {params.das_prefix})
        
        DAS_Tool \
            -i {params.metabat2_contigs},{params.maxbin_contigs},{params.concoct_contigs},{params.metadecoder_contigs},{params.semibin_contigs},{params.VAMB_contigs} \
            -c {input.contigs} \
            -o {params.das_prefix} \
            -l metabat2,maxbin,concoct,metadecoder,semibin,VAMB \
            --threads {resources.cpus} \
            --write_bins \
             | tee {log}
        """


rule drep_new:
    input: 
        "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/.bins_linked"
    output:
        touch("data/projects/{project}/{sample_type}/{sample}/bins/.drep_done")
    params:
        bins_linked = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/.bins_linked",
        input_bins = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/*.fa",
        genome_info = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/genome_info.csv", # Will need to be moved to inputs 
        main_dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/drep/"),
        MAGs= directory("data/projects/{project}/{sample_type}/metagenome_bins/derep/dereplicated_genomes")
    conda: "config/conda_yaml/drep.yaml"
    benchmark: "benchmarks/drep/{sample_type}-{project}__{sample}.txt"
    resources: cpus=8, mem_mb=150000, time_min=2880
    priority: 4
    shell:
        """
        rm -rf {params.main_dir} # Clear any old drep output
        
        dRep dereplicate \
            {params.main_dir} \
            -p {resources.cpus} \
            --contamination 50 \
            --completeness 30 \
            -pa 0.9 \
            -sa 0.99 \
            --length 10000 \
            --genomeInfo {params.genome_info} \
            -g {params.input_bins}
        """        

rule GTDB:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/.bins_linked",
        #"/home/kiledal/geomicro_home/references/.done_gtdb_refs_downloaded"
    params:
        input_bin_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep",
        refs = "/nfs/turbo/lsa-Erie/GVHD/data/reference/GTDBtk/release207_v2",
        out_dir = "data/projects/{project}/{sample_type}/{sample}/bins/GTDB",
        pplacer_cpus = 1
    output:
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/.done_GTDB")
    conda: "config/conda_yaml/gtdbtk.yaml"
    benchmark: "benchmarks/GTDB/{sample_type}-{project}__{sample}.txt"
    log: "logs/GTDB/{sample_type}-{project}__{sample}.log"
    resources: cpus=24, mem_mb=120000, time_min=2880
    priority: 4
    shell:
        """
        export GTDBTK_DATA_PATH={params.refs}

        # When upgrading version of GTDB, set this paramater for mashDB location
        #--mash_db $GTDBTK_DATA_PATH/mash_db

        gtdbtk classify_wf \
            --extension fa \
            --genome_dir {params.input_bin_dir} \
            --out_dir {params.out_dir} \
            --cpus {resources.cpus} \
            --pplacer_cpus {params.pplacer_cpus}
        """

rule GTDB_to_NCBI:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/.done_GTDB"
        #"/home/kiledal/geomicro_home/references/.done_gtdb_refs_downloaded"
    params:
        refs = "/nfs/turbo/lsa-Erie/GVHD/data/reference/GTDBtk/release207_v2",
        out_dir = "data/projects/{project}/{sample_type}/{sample}/bins/GTDB",
    output: "data/projects/{project}/{sample_type}/{sample}/bins/GTDB/gtdb_to_ncbi_taxonmy.tsv"
    conda: "config/conda_yaml/gtdbtk.yaml"
    benchmark: "benchmarks/GTDB_to_NCBI/{sample_type}-{project}__{sample}.txt"
    log: "logs/GTDB_to_NCBI/{sample_type}-{project}__{sample}.log"
    resources: cpus=1, mem_mb=4000, time_min=240
    priority: 4
    shell:
        """
        export GTDBTK_DATA_PATH={params.refs}

        python code/GTDBtk_scripts/gtdb_to_ncbi_majority_vote.py \
            --gtdbtk_output_dir {params.out_dir} \
            --output_file {params.out_dir}/gtdb_to_ncbi_taxonmy.tsv \
            --bac120_metadata_file {params.refs}/bac120_metadata_r207.tar.gz \
            --ar53_metadata_file {params.refs}/ar53_metadata_r207.tar.gz
        """

rule GTDB_versioned:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/.bins_linked",
        refs = "data/reference/GTDBtk/{database_version}"
    params:
        input_bin_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep",
        out_dir = "data/projects/{project}/{sample_type}/{sample}/bins/GTDB_{database_version}",
        pplacer_cpus = 1
    output:
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/.done_GTDB_{database_version}")
    conda: "config/conda_yaml/gtdbtk_2.4.0.yaml"
    benchmark: "benchmarks/GTDB/{sample_type}-{project}__{sample}_database-{database_version}.txt"
    log: "logs/GTDB/{sample_type}-{project}__{sample}_database-{database_version}.log"
    resources: cpus=16, mem_mb=100000, time_min=2880
    shell:
        """
        export GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf \
            --extension fa \
            --genome_dir {params.input_bin_dir} \
            --out_dir {params.out_dir} \
            --cpus {resources.cpus} \
            --pplacer_cpus {params.pplacer_cpus} \
            --mash_db $GTDBTK_DATA_PATH/mash_db
        """


rule gunc_GTDB_db_download:
    output: directory("data/reference/gunc_gtdb")
    resources: cpus=1, time_min=2880
    conda: "config/conda_yaml/gunc.yaml"
    shell:
        """
        mkdir -p {output}
        gunc download_db -db gtdb {output}
        """

rule gunc:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep/.bins_linked",
        ref_file = "data/reference/gunc_gtdb/gunc_db_gtdb95.dmnd",
        ref_dir = "data/reference/gunc_gtdb"
    output: 
        done = touch("data/projects/{project}/{sample_type}/{sample}/bins/.done_gunc")
    params: 
        bin_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bins_for_drep",
        out_dir = "data/projects/{project}/{sample_type}/{sample}/bins/gunc"
    resources: cpus = 24, mem_mb = 120000, time_min = 2880
    conda: "config/conda_yaml/gunc.yaml"
    benchmark: "benchmarks/gunc/{sample_type}-{project}__{sample}.txt"
    log: "logs/gunc/{sample_type}-{project}__{sample}.log"
    priority: 4
    shell:
        """
        mkdir -p {params.out_dir}
        gunc run \
            --input_dir {params.bin_dir} \
            -r {input.ref_file} \
            --threads {resources.cpus} \
            --temp_dir /tmp \
            --out_dir {params.out_dir}
        """


rule kofam_scan_bins:
    input:
        genes = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{bin}.faa",
        profile = "data/reference/kegg/kofamscan/profiles",
        ko_list = "data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/projects/{project}/{sample_type}/{sample}/bins/kofamscan/{bin}_kofam_results.txt"
    conda: "config/conda_yaml/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofamscan/{project}_{sample}_{sample_type}-{bin}.txt"
    log: "logs/kofamscan/{project}_{sample}_{sample_type}-{bin}.log"
    resources: cpus=12, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 50000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --format=detail-tsv \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.bin}_kofamscan \
            --ko-list {input.ko_list} {input.genes}
        """

rule traitar:
    input:
        #bin_genes = expand("data/omics/{sample_type}/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/{sample_type}/metagenome_bins/{BIN,[^/]+}.fa").BIN),
        pfam_db = "data/reference/traitar_pfamDB",
        #bin_dir = "data/omics/{sample_type}/metagenome_bins",
        sample_file = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar_sample_list.tsv"
    params:
        bin_dir = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes",
        gene_dir = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal",
        out_dir = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar"
    output: 
        #directory("data/omics/{sample_type}/metagenome_bins/traitar")
        touch("data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar/.done")
    benchmark: "logs/traitar_coassembly/{sample_type}-benchmark.txt"
    container: "library://a_gihawi/traitar3/traitar3"
    resources: cpus = 16, mem_mb = 170000, time_min=20000 #, partition = "largemem"
    shell:
        """
        #mkdir -p {params.out_dir}

        pwd
        cd ~/scratch_gdick1/GVHD/
        pwd

        traitar phenotype --overwrite -c {resources.cpus} /db {params.gene_dir} {input.sample_file} from_genes {params.out_dir}
        """

rule prodigal_mags_DREP:
    input:
        bin = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{bin}.fa"
    output: 
        proteins = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa",
        genes = "data/omics/{sample_type}/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.gff"
    conda: "config/conda_yaml/main.yaml"
    resources: cpus = 1, mem_mb = 10000
    shell:
        """
        prodigal -p meta -i {input.bin} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
        """

rule prodigal_MAGs_generic:
    input:
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/{genome}.fa"
    output: 
        proteins = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.faa",
        genes = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.gff"
    conda: "config/conda_yaml/main.yaml"
    resources: cpus = 1, mem_mb = 10000
    shell:
        """
        prodigal -p meta -i {input.genome} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
        """

rule prodigal_MAGs:
    input:
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/output/{sample}_DASTool_bins/{genome}.fa"
    output: 
        proteins = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.faa",
        genes = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.gff"
    conda: "config/conda_yaml/main.yaml"
    resources: cpus = 1, mem_mb = 10000
    shell:
        """
        prodigal -p meta -i {input.genome} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
        """

ruleorder: prodigal_MAGs_generic > prodigal_MAGs

rule bakta_das_bins:
    input:
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/das_tool/output/{sample}_DASTool_bins/{genome}.fa"
    output:
        dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}")
    params:
        db = "data/reference/bakta/db6"
    conda: "config/conda_yaml/bakta.yaml"
    log: "logs/bakta/{sample_type}-{project}__{sample}__{genome}.tsv"
    benchmark: "benchmarks/bakta/{sample_type}-{project}__{sample}__{genome}.tsv"
    resources: cpus=8, mem_mb=32000, time_min=5000, 
    shell:
        """
        bakta --db {params.db} \
            --keep-contig-headers \
            --output {output.dir} \
            --threads {resources.cpus} \
            {input.genome} | tee {log}
        """

ruleorder: bakta_generic > bakta_das_bins

rule bakta_generic:
    input:
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/all_raw_bins/{genome}.fa"
    output:
        dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}")
    params:
        db = "data/reference/bakta/db6"
    conda: "config/conda_yaml/bakta.yaml"
    log: "logs/bakta/{sample_type}-{project}__{sample}__{genome}.tsv"
    benchmark: "benchmarks/bakta/{sample_type}-{project}__{sample}__{genome}.tsv"
    resources: cpus=8, mem_mb=32000, time_min=5000, 
    shell:
        """
        bakta --db {params.db} \
            --keep-contig-headers \
            --output {output.dir} \
            --threads {resources.cpus} \
            {input.genome} | tee {log}
        """

rule antismash:
    input:
        #dir = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}",
        gff = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.gff",
        #genome = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}/{genome}.gbff",
        bakta_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}"
    output:
        dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/antismash/{genome}")
    params:
        db = "data/reference/antismash",
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}/{genome}.gbff"
    conda: "config/conda_yaml/antismash.yaml"
    log: "logs/antismash/{sample_type}-{project}__{sample}__{genome}.tsv"
    benchmark: "benchmarks/antismash/{sample_type}-{project}__{sample}__{genome}.tsv"
    resources: cpus=8, mem_mb=10000, time_min=5000, 
    shell:
        """
        antismash \
            --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees \
            --databases {params.db} \
            --genefinding-tool none \
            --output-dir {output.dir} \
            --cpus {resources.cpus} \
            {params.genome} | tee {log}
        """

rule antismash7:
    input:
        #dir = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}",
        gff = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.gff",
        #genome = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}/{genome}.gbff",
        bakta_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}"
    output:
        dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/antismash7/{genome}")
    params:
        db = "/local/databases",
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}/{genome}.gbff"
    singularity: "docker://antismash/standalone"
    log: "logs/antismash7/{sample_type}-{project}__{sample}__{genome}.tsv"
    benchmark: "benchmarks/antismash7/{sample_type}-{project}__{sample}__{genome}.tsv"
    resources: cpus=16, mem_mb=10000, time_min=5000, 
    shell:
        """
        pwd 
        cd {current_dir}
        pwd

        antismash \
            --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --cc-mibig --tfbs \
            -t bacteria \
            --databases {params.db} \
            --genefinding-tool none \
            --output-dir {output.dir} \
            --cpus {resources.cpus} \
            {params.genome} | tee {log}
        """

rule antismash8_db_download:
    output:
        database_dir = directory("data/reference/antismash8")
    conda: "config/conda_yaml/antismash8.yaml"
    log: "logs/antismash8_db_download.txt"
    benchmark: "benchmarks/antismash8_db_download.txt"
    resources: cpus=4, mem_mb=10000, time_min=10000
    shell:
        """
        pwd 
        cd {current_dir}
        pwd

        download-antismash-databases --database-dir {output.database_dir} | tee {log}

        #antismash --prepare-data --databases {output.database_dir}
        """

rule antismash8_db_prep:
    input:
        database_dir = "data/reference/antismash8"
    output:
        database_prep = touch("data/reference/.antismash8_dbPrep")
    conda: "config/conda_yaml/antismash8.yaml"
    log: "logs/antismash8_db_prep.txt"
    benchmark: "benchmarks/antismash8_db_prep.txt"
    resources: cpus=4, mem_mb=10000, time_min=10000
    shell:
        """
        antismash --prepare-data --databases {input.database_dir} | tee {log}
        """

rule antismash8:
    input:
        #dir = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}",
        gff = "data/projects/{project}/{sample_type}/{sample}/bins/prodigal/{genome}.gff",
        #genome = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}/{genome}.gbff",
        bakta_dir = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}"
    output:
        dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/antismash8/{genome}")
    params:
        db = "data/reference/antismash8",
        genome = "data/projects/{project}/{sample_type}/{sample}/bins/bakta/{genome}/{genome}.gbff"
    conda: "config/conda_yaml/antismash8.yaml"
    log: "logs/antismash8/{sample_type}-{project}__{sample}__{genome}.tsv"
    benchmark: "benchmarks/antismash8/{sample_type}-{project}__{sample}__{genome}.tsv"
    resources: cpus=16, mem_mb=10000, time_min=5000, 
    shell:
        """
        antismash \
            --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --cc-mibig --tfbs \
            -t bacteria \
            --databases {params.db} \
            --genefinding-tool none \
            --output-dir {output.dir} \
            --cpus {resources.cpus} \
            {params.genome} | tee {log}
        """

rule antismash_summary:
    input:
        "data/projects/{project}/{sample_type}/{sample}/bins/antismash{version}/{genome}"
    output:
        counts = "data/projects/{project}/{sample_type}/{sample}/bins/antismash{version}/summaries/{genome}-counts.tsv",
        region_summary = "data/projects/{project}/{sample_type}/{sample}/bins/antismash{version}/summaries/{genome}-region_summary.tsv"
    params:
        count_script = "code/multismash/workflow/scripts/count_regions_single.py",
        summarize_script = "code/multismash/workflow/scripts/tabulate_regions_single.py"
    log: "logs/antismash_summary/AS{version}-{sample_type}-{project}__{sample}__{genome}.tsv"
    benchmark: "benchmarks/antismash_summary/AS{version}-{sample_type}-{project}__{sample}__{genome}.tsv"
    resources: cpus=1, mem_mb=4000, time_min=120 
    shell:
        """
        python3 {params.count_script} {input} {output.counts} | tee {log}
        python3 {params.summarize_script} {input} {output.region_summary} | tee -a {log}
        """

rule get_bigscape_db:
    output: directory("data/reference/bigscape")
    conda: "config/conda_yaml/bigscape.yaml"
    resources:
    shell:
        """
        mkdir -p {output}
        cd {output}
        wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && gunzip Pfam-A.hmm.gz
        hmmpress Pfam-A.hmm
        """

rule setup_bigscape:
    output: directory("data/reference/BiG-SCAPE")
    resources: cpus=1, mem_mb=4000, time_min=100
    conda: "config/conda_yaml/bigscape.yaml"
    log: "logs/setup_bigscape.txt"
    shell:
        """
        cd data/reference
        git clone -b release_v1 https://github.com/medema-group/BiG-SCAPE.git | tee {current_dir}/{log}
        
        chmod +x BiG-SCAPE/*py | tee -a {current_dir}/{log}
        chmod a+w BiG-SCAPE/domains_color_file.tsv | tee -a {current_dir}/{log}
        chmod a+w BiG-SCAPE/Annotated_MIBiG_reference/ | tee -a {current_dir}/{log}
        """

rule bigscape:
    input: 
        antismash_gbk_dir = "data/projects/{project}/{sample_type}/coassembly/bins/bigscape/input_gbks",
        db = "data/reference/bigscape",
        bigscape_repo = "data/reference/BiG-SCAPE"
    output:
        #dir = directory("data/projects/{project}/{sample_type}/coassembly/bins/bigscape/output")
        done = touch("data/projects/{project}/{sample_type}/coassembly/bins/bigscape/.done")
    params:
        dir = "data/projects/{project}/{sample_type}/coassembly/bins/bigscape/output"
    conda: "config/conda_yaml/bigscape.yaml"
    #singularity: "docker://eandersk/big-scape"
    log: "logs/bigscape/{sample_type}-{project}.log"
    benchmark: "benchmarks/bigscape/{sample_type}-{project}.tsv"
    resources: cpus=24, mem_mb=100000, time_min=20000
    shell:
        """
        PATH={input.bigscape_repo}:$PATH

        bigscape.py \
            --inputdir {input.antismash_gbk_dir} \
            --outputdir {params.dir} \
            --pfam_dir {input.db} \
            --cores {resources.cpus} | tee {log}
        """

rule bakta_assembly:
    input:
        genome = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        prodigal = rules.prodigal.output.gbk
    output:
        dir = directory("data/omics/{sample_type}/{sample}/bakta_assembly")
    params:
        db = "data/reference/bakta/db6"
    conda: "config/conda_yaml/bakta.yaml"
    log: "logs/bakta_assembly/{sample_type}-{sample}.tsv"
    benchmark: "benchmarks/bakta_assembly/{sample_type}-{sample}.tsv"
    resources: cpus=32, mem_mb=120000, time_min=10000
    priority: 4 
    shell:
        """
        bakta --db {params.db} \
            --keep-contig-headers \
            --min-contig-length 1000 \
            --meta \
            --skip-trna \
            --skip-tmrna \
            --skip-crispr \
            --skip-pseudo \
            --skip-plot \
            --regions {input.prodigal} \
            --output {output.dir} \
            --threads {resources.cpus} \
            {input.genome} | tee {log}
        """

# rule antismash_assembly:
#     input:
#         #"data/omics/{sample_type}/{sample}/bakta_assembly",
#         genome = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
#         #gbk = "data/omics/{sample_type}/{sample}/genes/{sample}_GENES.gbk"
#     output:
#         dir = directory("data/omics/{sample_type}/{sample}/antismash_assembly")
#     params:
#         db = "data/reference/antismash",
#         #genome = "data/omics/{sample_type}/{sample}/bakta_assembly/{sample}.gbff"
#     conda: "config/conda_yaml/antismash.yaml"
#     log: "logs/antismash_assembly/{sample_type}-{sample}.tsv"
#     benchmark: "benchmarks/antismash_assembly/{sample_type}-{sample}.tsv"
#     resources: cpus=8, mem_mb=90000, time_min=20000
#     shell:
#         """
#         antismash \
#             --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees \
#             --databases {params.db} \
#             --genefinding-tool prodigal-m \
#             --output-dir {output.dir} \
#             --cpus {resources.cpus} \
#             {input.genome} | tee {log}

#             #--genefinding-tool none \
#         """

rule antismash8_assembly:
    input:
        "data/omics/{sample_type}/{sample}/bakta_assembly",
        genome = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        #gbk = "data/omics/{sample_type}/{sample}/genes/{sample}_GENES.gbk"
    output:
        dir = directory("data/omics/{sample_type}/{sample}/antismash8_assembly")
    params:
        db = "data/reference/antismash8_mibig4",
        genome = "data/omics/{sample_type}/{sample}/bakta_assembly/final.contigs.renamed.gbff"
    conda: "config/conda_yaml/antismash8.yaml"
    log: "logs/antismash8_assembly/{sample_type}__{sample}.tsv"
    benchmark: "benchmarks/antismash8_assembly/{sample_type}-{sample}.tsv"
    resources: cpus=1, mem_mb=100000, time_min=14400
    priority: 5
    shell:
        """
        antismash \
            --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --cc-mibig --tfbs \
            -t bacteria \
            --databases {params.db} \
            --output-dir {output.dir} \
            --cpus {resources.cpus} \
            --genefinding-tool none \
            --no-abort-on-invalid-records \
            {params.genome} | tee {log}

            #--genefinding-tool error \
            #--genefinding-gff3  \
            #--genefinding-tool prodigal-m \
            
        """

rule antismash_assembly_summary:
    input:
        "data/omics/{sample_type}/{sample}/antismash{version}_assembly"
    output:
        done = touch("data/omics/{sample_type}/{sample}/.done-antismash{version}_assembly_summary")
    params:
        count_script = "code/multismash/workflow/scripts/count_regions_single.py",
        summarize_script = "code/multismash/workflow/scripts/tabulate_regions_single.py",
        counts = "data/omics/{sample_type}/{sample}/antismash{version}_assembly/summaries/counts.tsv",
        region_summary = "data/omics/{sample_type}/{sample}/antismash{version}_assembly/summaries/region_summary.tsv"
    log: "logs/antismash_assembly_summary/AS{version}-{sample_type}-{sample}.tsv"
    benchmark: "benchmarks/antismash_assembly_summary/AS{version}-{sample_type}-{sample}.tsv"
    resources: cpus=1, mem_mb=16000, time_min=120 
    priority: 6
    shell:
        """
        mkdir -p $(dirname {params.counts})

        python3 {params.count_script} {input} {params.counts} | tee {log}
        python3 {params.summarize_script} {input} {params.region_summary} | tee -a {log}
        """

rule ref_read_mapping:
    input:
        f_reads = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_reads = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        ref = "data/reference/blast_queries/{ref_seqs}.fasta"
    output:
        temp_bam = temp("data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped_temp.bam"),
        sam = temp("data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped.sam"),
        unsorted_bam = temp("data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped_unsorted.bam"),
        bam = "data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped.bam"
    conda: "config/conda_yaml/minimap2.yaml"
    log: "logs/ref_read_mapping/{sample_type}-{sample}.{ref_seqs}.log"
    benchmark: "benchmarks/ref_read_mapping/{sample_type}-{sample}.{ref_seqs}.tsv"
    resources: cpus=8
    shell:
        """
        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --secondary=no \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS {output.sam} > {output.temp_bam}
        
        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule ref_read_mapping_pileup:
    input:
        bam = "data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_mapped.bam",
        ref = "data/reference/blast_queries/{ref_seqs}.fasta"
    output:
        pileup = "data/omics/{sample_type}/{sample}/ref_read_mapping/{ref_seqs}_pileup.txt"
    conda: "config/conda_yaml/minimap2.yaml"
    log: "logs/ref_read_mapping_pileup/{sample_type}-{sample}.{ref_seqs}.log"
    benchmark: "benchmarks/ref_read_mapping_pileup/{sample_type}-{sample}.{ref_seqs}.tsv"
    resources: cpus=1
    shell:
        """
        samtools mpileup -f {input.ref} -o {output.pileup} {input.bam}
        """

rule contig_search:
    input:
        query="data/reference/blast_queries/{ref_seqs}.fasta",
        contigs="data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa"
    output:
        "data/omics/{sample_type}/{sample}/contig_search/{ref_seqs}.m8"
    params:
        tmp_results = "/tmp/{sample_type}/{sample}/contig_search/tmp_results_{ref_seqs}"
    conda: "config/conda_yaml/mmseqs.yaml"
    benchmark: "benchmarks/contig_search/{sample_type}_{sample}--{ref_seqs}.txt"
    log: "logs/contig_search/{sample_type}_{sample}--{ref_seqs}.log"
    resources: cpus=24, mem_mb=120000, time_min=20000
    shell:
        """
        mkdir -p {params.tmp_results}

        mmseqs easy-search {input.query} {input.contigs} {output} {params.tmp_results} \
        --threads {resources.cpus} \
        --exhaustive-search \
        --search-type 3
        """

rule map_reads_to_microcystis_markers:
    input:
        f_reads = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_reads = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        ref = "data/reference/microcystis_markers/{marker}.fasta"
    output:
        temp_bam = temp("data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_mapped_temp.bam"),
        sam = temp("data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_mapped.sam"),
        bam = "data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_mapped.bam",
        unsorted_bam = temp("data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_mapped_unsorted.bam")
    params:
        min_cover = 25,
        min_id = 95
    conda: "config/conda_yaml/bwa.yaml"
    benchmark: "benchmarks/map_reads_to_microcystis_markers/{sample_type}_{sample}--{marker}.txt"
    log: "logs/map_reads_to_microcystis_markers/{sample_type}_{sample}--{marker}.log"
    resources: cpus=16
    shell:
        """
        mkdir -p $(dirname {output.bam})

        minimap2 \
            -ax sr \
            -t {resources.cpus} \
            --secondary=no \
            --sam-hit-only \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS {output.sam} > {output.temp_bam}
        
        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover {params.min_cover} \
            --minId {params.min_id}
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule summarize_marker_mapping:
    input:
        bam = "data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_mapped.bam",
        read_counts = "data/omics/{sample_type}/{sample}/reads/{sample}_read_count_fastp.tsv",
        marker_info = "data/reference/microcystis_markers/info/20240709_groupings.tsv"
    output:
        marker_summary = "data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_summary.tsv",
        clade_summary = "data/omics/{sample_type}/{sample}/microcystis_markers/{sample}--{marker}_clade-summary.tsv"
    params:
    resources: cpus=1, mem_mb=5000, time_min=60
    benchmark: "benchmarks/summarize_marker_mapping/{sample_type}_{sample}--{marker}.txt"
    log: "logs/summarize_marker_mapping/{sample_type}_{sample}--{marker}.log"
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        code/summarize_marker_gene_read_mapping.R \
            -i {input.bam} \
            --clade-summary {output.clade_summary} \
            --marker-summary {output.marker_summary} \
            --prefix {wildcards.sample} \
            --info {input.marker_info} \
            --read-counts {input.read_counts}
        """


rule amplicon_hmm:
    input: "data/omics/amplicons/{sample}/reads/raw_{direc}_reads.fastq.gz"
    output:
        hmm_tbl = "data/omics/amplicons/{sample}/detect_region/tbl_{direc}.txt",
        full_out = "data/omics/amplicons/{sample}/detect_region/full_{direc}.txt",
    params:
        amplicon_hmm_db ="data/reference/hmm_amplicons/combined.hmm",
        outdir = subpath(input[0], parent=True),
    resources: cpus=4, mem_mb=20000, time_min=500
    benchmark: "benchmarks/amplicon_hmm/{sample}_{direc}.txt"
    log: "logs/amplicon_hmm/{sample}_{direc}.log"
    conda: "config/conda_yaml/hmmer.yaml"
    shell:
        """
        . code/shell_prelude {log}
        mkdir -p -- {params.outdir}

        seqkit head -n 1000 {input} |
            seqkit fq2fa |
            nhmmscan --cpu {resources.cpus} --tblout {output.hmm_tbl} {params.amplicon_hmm_db} - > {output.full_out}
        """

rule amplicon_hmm_summarize_r:
    """
    The original R-implemented hmm summary

    This is not used downstream but its output may be of interest for
    debugging.  So needs to be triggered manually.
    """
    input: rules.amplicon_hmm.output.hmm_tbl
    output: "data/omics/amplicons/{sample}/detect_region/{direc}_summary.tsv"
    resources: cpus=1, mem_mb=10000, time_min=500
    benchmark: "benchmarks/amplicon_hmm_summarize_a/{sample}_{direc}.txt"
    log: "logs/amplicon_hmm_summarize_a/{sample}_{direc}.log"
    container: "docker://eandersk/r_microbiome"
    shell: "./code/summarize_hmmscan.R --input {input} --output {output} |& tee {log}"

rule amplicon_hmm_summarize:
    input: rules.amplicon_hmm.output.hmm_tbl
    output: "data/omics/amplicons/{sample}/detect_region/{direc}_summary.json"
    resources: cpus=1, mem_mb=100, time_min=1
    benchmark: "benchmarks/amplicon_hmm_summarize/{sample}_{direc}.txt"
    log: "logs/amplicon_hmm_summarize/{sample}_{direc}.log"
    run: pypelib.amplicon.hmm_summarize.main(input[0], output[0], log=log)


def get_hmm_summaries(wc):
    """
    Input function to get the hmm summaries

    Wildcards need to contain sample_type and sample
    """
    runinfo = checkpoints.get_runinfo.get(**wc).output[0]
    match parse_runinfo(runinfo, 'library_layout'):
        case 'PAIRED':
            # This MUST have fwd first and rev second, this is what guess_target expects!
            return [
                str(rules.amplicon_hmm_summarize.output).format(direc='fwd', sample=wc['sample']),
                str(rules.amplicon_hmm_summarize.output).format(direc='rev', sample=wc['sample']),
            ]
        case 'SINGLE':
            return str(rules.amplicon_hmm_summarize.output).format(direc='single', sample=wc['sample'])
        case _:
            raise ValueError(f'unsupported library layout: {_}')

rule amplicon_guess_target:
    input:
        summaries=get_hmm_summaries,
        stats = rules.raw_reads_stats.output.stats
    output:
        target_info = "data/omics/{sample_type}/{sample}/detect_region/target_info.json"
    resources: cpus=1, mem_mb=100, time_min=1
    log: "logs/amplicon_guess_target/{sample_type}_{sample}.log"
    run: pypelib.amplicon.guess_target.main(input.summaries, input.stats, output.target_info, log=log)


def get_target_info_files(wc):
    """
    input function for amplicon_dispatch

    Collects target info files for each amplicon sample in given dataset.
    """
    assignments=Path(str(rules.amplicon_dispatch.output.assignments).format(dataset=wc.dataset))
    if assignments.exists():
        # check for manual overrides
        skips = [
            sample for sample, group
            in pypelib.amplicon.dispatch.get_assignments(path=assignments).items()
            if group == pypelib.amplicon.dispatch.SKIP
        ]
    else:
        skips = []

    project_samples = Path("data/projects/") / wc.dataset / 'amplicons'
    target_info = rules.amplicon_guess_target.output.target_info
    ret = []
    for i in project_samples.glob("samp_*"):
        if i in skips:
            continue
        ret.append(target_info.format(sample_type='amplicons', sample=i.name))
    return ret

checkpoint amplicon_collect_target_guesses:
    input:
        get_target_info_files
    output:
        target_tab="data/projects/{dataset}/target_info.tsv"
    params:
        project_dir = subpath(output.target_tab, parent=True)
    log: "logs/amplicon_collect_target_guesses/{dataset}.log"
    run: pypelib.amplicon.tabulate_targets.main(input, output=output.target_tab, log=log, **params, **wildcards)

rule remove_primers_pe:
    input:
        target_info=rules.amplicon_guess_target.output.target_info,
        fwd=rules.get_reads_paired.output.fwd_reads,
        rev=rules.get_reads_paired.output.rev_reads,
    output:
        fwd="data/omics/{sample_type}/{sample}/reads/clean.fwd_reads.fastq.gz",
        rev="data/omics/{sample_type}/{sample}/reads/clean.rev_reads.fastq.gz",
    params:
        reads_dir = subpath(output.fwd, parent=True)
    #conda: uses cutadapt installed in the the runtime environment
    log: "logs/remove_primers/{sample_type}-{sample}.log"
    run:
        pypelib.amplicon.remove_primers.main_paired(
            input.target_info,
            input.fwd,
            input.rev,
            output.fwd,
            output.rev,
            log=log,
        )

rule remove_primers_se:
    input:
        target_info=rules.amplicon_guess_target.output.target_info,
        single=rules.get_reads_single.output.single_reads,
    output:
        single="data/omics/{sample_type}/{sample}/reads/clean.single_reads.fastq.gz"
    params:
        reads_dir = subpath(output.single, parent=True)
    #conda: uses cutadapt installed in the the runtime environment
    log: "logs/remove_primers/{sample_type}-{sample}.log"
    run:
        pypelib.amplicon.remove_primers.main_single(
            input.target_info,
            input.single,
            output.single,
            log=log,
        )

def get_dataset_fastq_files(wc):
    target_tab = checkpoints.amplicon_collect_target_guesses.get(**wc).output.target_tab
    return pypelib.amplicon.dispatch.get_fastq_files(target_tab, **wc)

checkpoint amplicon_dispatch:
    input:
        fastqs = get_dataset_fastq_files,
        target_tab = rules.amplicon_collect_target_guesses.output.target_tab
    output:
        assignments=update("data/projects/{dataset}/amplicon_target_assignments.tsv"),
        samples="data/projects/{dataset}/amplicon_sample_info.tsv",
    params:
        project_dir = subpath(output.assignments, parent=True)
    run:
        pypelib.amplicon.dispatch.make(
            input.fastqs,
            input.target_tab,
            params.project_dir,
            out_assignments=output.assignments,
            out_samples=output.samples,
        )


rule hmm_summary_plots:
    input: rules.amplicon_dispatch.output.assignments
    output: pdf = "data/projects/{dataset}/hmm_summary_plots.pdf"
    params:
        project_dir = subpath(output.pdf, parent=True)
    run: pypelib.amplicon.hmm_summarize.multiplot(input[0], output=output.pdf)


rule amplicon_dada2_target:
    input:
        assignments = rules.amplicon_dispatch.output.assignments,
        samples = rules.amplicon_dispatch.output.samples,
        target_tab = rules.amplicon_collect_target_guesses.output.target_tab
    output:
        seqs = "data/projects/{dataset}/dada2.{target_spec}/rep_seqs.fasta",
        abund = "data/projects/{dataset}/dada2.{target_spec}/asv_table.tsv",
    params:
        outdir = subpath(output.abund, parent=True),
        quality_threshold=25,
    log: "logs/dada2/{dataset}_{target_spec}.log"
    container: "docker://eandersk/r_microbiome"
    resources: cpus=16
    shell:
        """
        . code/shell_prelude {log}
        ./code/ampliconTrunc.R \
            --quality {params.quality_threshold} \
            --outdir {params.outdir} \
            --assignments {input.assignments} \
            --samples {input.samples} \
            --targets {input.target_tab} \
            --cpus {resources.cpus} \
            {wildcards.target_spec}
        status=$?
        if [[ $status -ne 0 ]]; then
            mv -v --backup=numbered -T -- {params.outdir} {params.outdir}_ERROR
            exit $?
        fi
        """

rule amplicon_asv_test:
    input: rules.amplicon_dada2_target.output.seqs
    output: "data/projects/{dataset}/dada2.{target_spec}/rep_seqs_hmm.txt"
    params: hmm_db = "data/reference/hmm_amplicons/combined.hmm"
    run:
        pypelib.amplicon.hmm_check_asvs.main(
            params.hmm_db,
            input[0],
            wildcards.target_spec,
            output[0],
        )


def get_dada2_output(wc):
    """
    input function: get dada2 output files for given dataset

    Derives output directory names from the distinct target (from first column)
    in sample listing.

    This will also check for missing target assignments and may raise a
    RuntimeError prompting to manually curate the target assignment file.
    """
    assignment_file = Path(checkpoints.amplicon_dispatch.get(**wc).output.assignments)
    assignments = pypelib.amplicon.dispatch.get_assignments(path=assignment_file)

    wc_str = ' '.join(f'{k}={v}' for k, v in wc.items())
    targets = Counter(assignments.values())
    if bad_count := targets.pop(pypelib.amplicon.dispatch.UNKNOWN, 0):
        raise RuntimeError(
            f'[{wc_str}] There are {bad_count} samples that do not have an amplicon '
            f'target assigned.  Manual curation required for {assignment_file}'
        )
    if skip_count := targets.pop(pypelib.amplicon.dispatch.SKIP, 0):
        print(f'[INFO]: [{wc_str}] {skip_count} samples got excluded from analysis')

    specs = sorted(set(
        pypelib.amplicon.dispatch.target2spec(i)
        for i in targets
    ))
    if specs:
        print(f'{len(specs)} dada2 job(s) for wildcards: {wc_str}:')
        for i in specs:
            print(' ', i)
    else:
        print('[WARNING] no dada2 jobs for wildcard {wc_str} ???')
    files = []
    for i in specs:
        for j in rules.amplicon_dada2_target.output:
            files.append(j.format(dataset=wc.dataset, target_spec=i))
        for j in rules.amplicon_asv_test.output:
            files.append(j.format(dataset=wc.dataset, target_spec=i))
    return files

rule amplicon_pipeline_dataset:
    """ The top rule for the amplicon pipeline """
    input: get_dada2_output
    output: "data/projects/{dataset}/amplicon_pipeline_done"
    shell:
        """
        for i in {input:q}; do echo "$i" >> {output}; done
        """

rule deeparg_ls: #this rule is using LS mode and annotated genes  
    input:
        genes="data/omics/{sample_type}/{sample}/genes/{sample}_GENES.fna"  # For sets 35, 41, 42
    output:
        predictions = "data/omics/{sample_type}/{sample}/deeparg/predictions.out.mapping.ARG"
    params:
        reference="data/reference/deeparg",
        output_prefix="data/omics/{sample_type}/{sample}/deeparg/predictions.out"
    benchmark: "benchmarks/deeparg_ls/{sample_type}_{sample}.txt"
    log: "logs/deeparg_ls/{sample_type}_{sample}.log"
    conda:
        "config/conda_yaml/deeparg.yml"
    resources:
        cpus=1, mem_mb=20000, time_min=10000
    shell:
        """
        deeparg predict \
            --model LS \
            -i {input.genes} \
            -o {params.output_prefix} \
            -d {params.reference} \
            --type nucl \
            --min-prob 0.8 \
            --arg-alignment-identity 30 \
            --arg-alignment-evalue 1e-10 \
            --arg-num-alignments-per-entry 1000 | tee {log}
        """


rule metaEuk:
    input:
        assembly = "data/omics/{sample_type}/{sample}/assembly/RNAspades/contigs.renamed.fasta"
    output: 
        done = touch("data/omics/{sample_type}/{sample}/metaeuk/.done"),
        aa_seqs = "data/omics/{sample_type}/{sample}/metaeuk/{sample}.fas",
        header_map = "data/omics/{sample_type}/{sample}/metaeuk/{sample}.headersMap.tsv",
        tax_res = "data/omics/{sample_type}/{sample}/metaeuk/{sample}_tax"
    conda:  "config/conda_yaml/metaeuk.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/metaeuk/{sample}",
        tmp_dir = "/home/kiledal/scratch_gdick1/metaEuk/{sample}",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        #tmp_dir = "/scratch/kiledal/metaEuk/{sample}",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/metaEuk/{sample_type}-{sample}.txt"
    log: "logs/metaEuk/{sample_type}-{sample}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        mem_mb = 160000, cpus=32, time_min=20000
    shell:
        """
        export MMSEQS_NUM_THREADS={resources.cpus}
        export TMPDIR={params.tmp_dir}
        
        # Different tmp dir if running on lab servers
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/metaEuk/{wildcards.sample}
        
        mkdir -p $TMPDIR

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}
        
        metaeuk \
            easy-predict \
            --threads {resources.cpus} \
            --split-memory-limit 50G \
            {input.assembly} \
            {params.unirefDB} \
            {params.out_prefix} \
            $TMPDIR
            2>&1 | tee -a {log} &&

        metaeuk taxtocontig \
            $TMPDIR/latest/contigs \
            {output.aa_seqs} \
            {output.header_map} \
            {params.unirefDB} \
            {output.tax_res} \
            $TMPDIR \
            --majority 0.5 \
            --tax-lineage 1 \
            --lca-mode 2

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r $TMPDIR 2>&1 | tee -a {log}
        #printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """

rule reads_unirefLCA_mmseqs_metaT:
    input:
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/noDedup_decon_rev_reads_fastp.fastq.gz"
    output:
        report = "data/omics/{sample_type}/{sample}/uniref_readmap_noDedup/{sample}_report"
    conda:  "config/conda_yaml/mmseqs.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/uniref_readmap_noDedup/{sample}",
        tmp_dir = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/tmp",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        #tmp_dir = "/tmp/kiledal/mmseqs2/{sample}",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        #tmp_dir = "/scratch/kiledal/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/reads_unirefLCA_mmseqs/{sample_type}-{sample}.txt"
    log: "logs/reads_unirefLCA_mmseqs/{sample_type}-{sample}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        mem_mb = 160000, cpus=32, time_min=20000
    shell:
        """
        export TMPDIR={params.tmp_dir}
        
        # Different tmp dir if running on lab servers
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/mmseqs2/{wildcards.sample}
        
        mkdir -p $TMPDIR
        mkdir $(basename {params.out_prefix})

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        # Previous copied reads to temp space, but no longer think this is nescessary. Should be deleted.
        #cp {input.fwd_reads} {params.tmp_fwd_reads}
        #cp {input.rev_reads} {params.tmp_rev_reads}

        #mmseqs touchdb {params.unirefDB} # loads  the database into memory, only use in large memory systems / nodes
        
        mmseqs \
            easy-taxonomy \
            {input.fwd_reads} {input.rev_reads} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            $TMPDIR \
            --lca-mode 3 \
            --orf-filter 1 \
            --orf-filter-s 3.5 \
            -s 4 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit 300G \
            2>&1 | tee -a {log}

            #{params.tmp_fwd_reads} {params.tmp_rev_reads} \
            # --db-load-mode 2 # this loads the database into memory, was an attempt to speed up on Great Lakes but limited memory & scratch space were limiting factors

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h $TMPDIR 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r {params.tmp_fwd_reads} {params.tmp_rev_reads} $TMPDIR 2>&1 | tee -a {log}
        rm -r $TMPDIR 2>&1 | tee -a {log}
        printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """

rule add_lineage_to_unirefLCAtax_metaT:
    input:
        db = ancient("data/reference/ncbi_tax"),
        report = "data/omics/{sample_type}/{sample}/uniref_readmap_noDedup/{sample}_report"
    output: 
        inspect_w_lineage = "data/omics/{sample_type}/{sample}/uniref_readmap_noDedup/{sample}_report_w_full_lineage",
        inspect_w_7_lev_lineage = "data/omics/{sample_type}/{sample}/uniref_readmap_noDedup/{sample}_report_w_standardized_lineage"
    conda: "config/conda_yaml/taxonkit.yaml"
    resources: cpus=1, mem_mb=10000, time_min=5440, mem_gb = 10
    log: "logs/unirefLCA_mmseqs_add_lineage/{sample_type}-{sample}.log"
    shell:
        """
        taxonkit lineage \
            {input.report} \
            --data-dir {input.db} \
            -i 5 \
            -o {output.inspect_w_lineage} 2>&1 | tee {log}

        taxonkit reformat \
            {output.inspect_w_lineage} \
            --data-dir {input.db} \
            -i 7 \
            -P \
            -o {output.inspect_w_7_lev_lineage} 2>&1 | tee -a {log}
        """

rule contig_unirefLCA_mmseqs_generic:
    input:
        contigs = "data/omics/{sample_type}/{sample}/assembly/{assembly}/contigs.renamed.fasta"
    output:
        report = "data/omics/{sample_type}/{sample}/assembly/{assembly}/{sample}_contig_report"
    conda:  "config/conda_yaml/mmseqs.yaml"
    params:
        #unirefDB = "/home/kiledal/scratch_gdick1/mmseqs_unirefdb/mmseqs2/uniref100",
        unirefDB = "data/reference/mmseqs2/uniref100",
        out_prefix = "data/omics/{sample_type}/{sample}/assembly/{assembly}/{sample}_contig",
        tmp_dir = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/tmp",
        # tmp_dir = "/dev/shm/akiledal/mmseqs2",
        # tmp_fwd_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/home/kiledal/scratch_gdick1/mmseqs/{sample}/{sample}__rev.fastq.gz",
        #tmp_dir = "/tmp/kiledal/mmseqs2/{sample}_contigs",
        #tmp_dir = "/scratch/kiledal/mmseqs2/{sample}_contigs",
        # tmp_fwd_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__fwd.fastq.gz",
        # tmp_rev_reads = "/tmp/kiledal/mmseqs2/{sample}/{sample}__rev.fastq.gz"
        #tmp_dir = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}",
        tmp_fwd_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__fwd.fastq.gz",
        tmp_rev_reads = "/old-geomicro/kiledal-extra/mmseqs_tmp/{sample}/{sample}__rev.fastq.gz"
    benchmark: "benchmarks/contig_unirefLCA_mmseqs_generic/{sample_type}-{sample}--{assembly}.txt"
    log: "logs/contig_unirefLCA_mmseqs_generic/{sample_type}-{sample}--{assembly}.log"
    resources:
        #mem_mb = 1450000, cpus=32, time_min=20000, partition = "largemem"
        mem_mb = 160000, cpus=32, time_min=7200
    shell:
        """
        export TMPDIR={params.tmp_dir}
        
        # Different tmp dir if running on lab servers
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/mmseqs2/{wildcards.sample}_contigs
        
        mkdir -p $TMPDIR

        # Log how much temp space is available, can cause job to fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee {log}
        df -h $ 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #mmseqs touchdb {params.unirefDB} # loads the database into memory, only use in large memory systems / nodes

        mmseqs \
            easy-taxonomy \
            {input.contigs} \
            {params.unirefDB} \
            ./{params.out_prefix} \
            $ \
            --lca-mode 3 \
            --orf-filter 1 \
            --orf-filter-s 3.5 \
            -s 4 \
            --tax-lineage 1 \
            --threads {resources.cpus} \
            --split-memory-limit $(echo "scale=0;({resources.mem_mb}*0.8)/1024" | bc)G \
            --db-load-mode 1 \
            2>&1 | tee -a {log}

            # --db-load-mode 2 # this loads the database into memory, was an attempt to speed up on Great Lakes but limited memory & scratch space were limiting factors

        # Prior to clearing temp files, log temp space to see if potential reason for fail
        printf "\nTemp directory storage available:\n" 2>&1 | tee -a {log}
        df -h $ 2>&1 | tee -a {log}
        printf "\n\n" 2>&1 | tee -a {log}

        #rm -r {params.tmp_fwd_reads} {params.tmp_rev_reads} $ 2>&1 | tee -a {log}
        rm -r $ 2>&1 | tee -a {log}
        printf "Done, and deleted temp dir" 2>&1 | tee -a {log}
        """

rule plass:
    input:
        fwd_reads = "data/omics/{sample_type}/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/{sample_type}/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        assembly = "data/omics/{sample_type}/{sample}/assembly/plass/assembly.faa",
    params:
        tmpdir = "/home/kiledal/scratch_gdick1/plass/{sample}"
    conda: "config/conda_yaml/plass.yaml"
    log: "logs/plass/{sample}-{sample_type}.log"
    benchmark: "benchmarks/plass/{sample}-{sample_type}.txt"
    resources: cpus=24, mem_mb=120000, time_min=2880 # standard assemblies
    shell:
        """
        export TMPDIR={params.tmpdir}
        [[ "${{HOSTNAME}}" == "cayman" || "${{HOSTNAME}}" == "vondamm" ]] && export TMPDIR=/scratch/$USER/plass/{wildcards.sample}
        mkdir -p $TMPDIR

        plass assemble \
            --threads {resources.cpus} \
            {input.fwd_reads} \
            {input.rev_reads} \
            {output.assembly} \
            $TMPDIR
        """

rule kofam_scan_plass:
    input:
        genes = "data/omics/metatranscriptomes/{sample}/assembly/plass/assembly.faa",
        profile = "data/reference/kegg/kofamscan/profiles",
        ko_list = "data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/omics/metatranscriptomes/{sample}/assembly/plass/{sample}_kofam_results.txt"
    conda: "config/conda_yaml/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofam_scan_plass/{sample}.txt"
    log: "logs/kofam_scan_plass/{sample}.log"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --format=detail-tsv \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.sample}_kofamscan \
            --ko-list {input.ko_list} {input.genes} | tee {log}
        """


rule map_reads_to_plass:
    input:
        plass = "data/omics/metatranscriptomes/{sample}/assembly/plass/assembly.faa",
        fwd_reads = "data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_rev_reads_fastp.fastq.gz"
    output:
        interleaved = temp("data/omics/metatranscriptomes/{sample}/reads/noDedup_decon_interleaved_reads_fastp.fastq.gz"),
        temp_bam = temp("data/omics/metatranscriptomes/{sample}/assembly/plass/mapped_temp.bam"),
        sam = temp("data/omics/metatranscriptomes/{sample}/assembly/plass/mapped.sam"),
        unsorted_bam = temp("data/omics/metatranscriptomes/{sample}/assembly/plass/mapped_unsorted.bam"),
        bam = "data/omics/metatranscriptomes/{sample}/assembly/plass/mapped.bam",
        diamond_db = temp("data/omics/metatranscriptomes/{sample}/assembly/plass/diamond_db")
    conda: "config/conda_yaml/diamond_read_map.yaml"
    benchmark: "benchmarks/map_reads_to_plass/{sample}.txt"
    log: "logs/map_reads_to_plass/{sample}.log"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 150000
    shell:
        """
        seqtk mergepe {input.fwd_reads} {input.rev_reads} > {output.interleaved}

        diamond makedb --in {input.plass} --db {output.diamond_db}

        # Step 2: Map metatranscriptomic reads to the PLASS assembly using DIAMOND (blastx mode) and output to SAM format
        diamond blastx -d {output.diamond_db} \
               -q {output.interleaved} \
               -o {output.sam} \
               --outfmt 101 \
               --max-target-seqs 1 \
               --evalue 1e-5 \
               --threads 8 \
               --sensitive

        # Convert sam to bam
        samtools view -bS {output.sam} > {output.temp_bam}
        
        # Filter mappings
        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 50 \
            --minId 80
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}       
        """

rule virsorter2:
    input:
        contigs = rules.rename_contigs.output.contigs
    output:
        virsorter_dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/virsorter2"),
        virsorter_done = touch("data/projects/{project}/{sample_type}/{sample}/bins/.virsorter2_done")
    params:
        out_prefix = "data/projects/{project}/{sample_type}/{sample}/bins/metadecoder/bins/{sample}"
    conda: "config/conda_yaml/virsorter.yaml"
    benchmark: "benchmarks/virsorter2/{sample_type}-{project}__{sample}.txt"
    log: "logs/virsorter2/{sample_type}-{project}__{sample}.log"
    resources: cpus=32, mem_mb=150000, time_min=10080 # standard samples
    priority: 3
    shell:
        """
        virsorter run \
            -w {output.virsorter_dir} \
            -i {input.contigs} \
            --db-dir data/reference/virsorter2 \
            --include-groups all \
            -j {resources.cpus} \
            all | tee {log}
        """

rule genomad:
    input:
        contigs = rules.rename_contigs.output.contigs
    output:
        genomad_dir = directory("data/projects/{project}/{sample_type}/{sample}/bins/genomad"),
        genomad_done = touch("data/projects/{project}/{sample_type}/{sample}/bins/.genomad_done")
    params:
        db_path = "data/reference/genomad/genomad_db",
        splits = "8",
        out_prefix = "data/projects/{project}/{sample_type}/{sample}/bins/genomad"
    conda: "config/conda_yaml/genomad.yaml"
    benchmark: "benchmarks/genomad/{sample_type}-{project}__{sample}.txt"
    log: "logs/genomad/{sample_type}-{project}__{sample}.log"
    resources: cpus=32, mem_mb=150000, time_min=10080 # standard samples
    priority: 3
    shell:
        """
        genomad end-to-end \
            --cleanup \
            --splits {params.splits} \
            {input.contigs} \
            {output.genomad_dir} \
            {params.db_path} | tee {log}
        """

rule gutsmash:
    input:
        "data/omics/{sample_type}/{sample}/bakta_assembly",
        genome = "data/omics/{sample_type}/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
    output:
        dir = directory("data/omics/{sample_type}/{sample}/gutsmash_assembly")
    params:
        db = "data/reference/antismash8_mibig4",
        genome = "data/omics/{sample_type}/{sample}/bakta_assembly/final.contigs.renamed.gbff"
    singularity: "docker://nmendozam/gutsmash"
    log: "logs/gutsmash/{sample_type}__{sample}.tsv"
    benchmark: "benchmarks/gutsmash/{sample_type}-{sample}.tsv"
    resources: cpus=1, mem_mb=100000, time_min=14400
    shell:
        """
        run_gutsmash.py \
            --minimal \
            --cb-knownclusters \
            --enable-genefunctions \
            --cpus {resources.cpus} \
            --output-dir {output.dir} \
            {params.genome} | tee {log}
        """



