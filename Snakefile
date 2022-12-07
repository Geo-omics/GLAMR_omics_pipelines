import os
import re
import snakemake.io
from glob import glob

configfile: "config.yaml"
report: "code/report/workflow.rst"

# Substitute bash $USER environment variable with actual user id, otherwise some steps fail
param_work_dir = config["work_dir"] #get working directory from config file
userid = env_var = os.environ['USER'] #get bash $USER variable
work_dir = param_work_dir.replace("$USER", userid) #sub $USER for actual username
current_dir = os.getcwd()

# Get import sample names
#metaG_samples = glob_wildcards("import/metagenomes/{sample}/").sample
#metaT_samples = glob_wildcards("import/metatranscriptomes/{sample}/").sample
#metabolome_samples = glob_wildcards("import/metabolomes/{sample}/").sample
#amplicon_samples = glob_wildcards("import/amplicons/{sample}/").sample

# Get sample names
#metaG_samples = glob_wildcards("data/omics/metagenomes/{sample}/reads").sample
jgi_samples = glob_wildcards("import/staging/jgi_2022/all_sample_filtered_reads/{sample}_interleaved.fastq.gz").sample
metaG_samples = glob_wildcards("data/projects/2022_geomicro_JGI_CSP/metagenomes/{sample}/").sample
#metaG_samples = glob_wildcards("data/projects/PRJNA464361/metagenomes/{sample}/").sample
#metaG_samples = ["E20212019","E20212012","E20212013","E20212010"]
read_download_samples = glob_wildcards("data/omics/metagenomes/{sample}/reads/accession").sample
metaT_samples = glob_wildcards("data/omics/metatranscriptomes/{sample}/").sample
metabolome_samples = glob_wildcards("data/omics/metabolomes/{sample}/").sample
amplicon_samples = glob_wildcards("data/omics/amplicons/{sample}/").sample

rule assemble:
    input: 
        #expand("data/omics/metagenomes/{sample}/assembly/metaspades/contigs.fasta",sample = metaG_samples),
        #expand("data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa",sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/assembly/metaspades_noNORM/contigs.fasta",sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.fa",sample = metaG_samples)

rule run_megahit:
    input:
        expand("data/omics/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.fa",sample = metaG_samples)

rule run_metaspades:
    input:
        expand("data/omics/metagenomes/{sample}/assembly/metaspades_noNORM/contigs.fasta",sample = metaG_samples)

rule run_prodigal:
    input: expand("data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa", sample = metaG_samples)

rule test:
    input: expand("data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",sample = metaG_samples)
    output: "test.out"


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
        snakemake run_drep_sep metaG_annotation run_humann_fastp run_sourmash data/sample_data/bracken_counts.tsv --rulegraph --dry-run | dot -Tpdf > {output.pdf}
        snakemake run_drep_sep metaG_annotation run_humann_fastp run_sourmash data/sample_data/bracken_counts.tsv --rulegraph --dry-run | dot -Tpng > {output.png}
        """

# rule import:
#     input:
#     output:
#     resources: cpus=1, mem_mb=8000, time_min=2880, mem_gb = 8
#     shell:
#         """

#         """

# rule gzip_fastx:

# rule unzip_fastx:


# rule zip_fastqs:
#     input:  SCRATCH + "/01RAW_fqs/{sample}"
#     output: temp(SCRATCH + "/02ZIPPED_fqs/{sample}")
#     params: outdir = SCRATCH + "/02ZIPPED_fqs/"
#     run:

#         if wildcards.sample.endswith('.fastq'):
#             shell("echo gzip {input}")
#             shell("echo mv {input}.gz {params.outdir}")
#         else:
#             shell("mv {input} {params.outdir}")


rule get_reads:
    input:
        acccession = "data/omics/metagenomes/{sample}/reads/accession"
    output:
        fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        touch = touch("data/omics/metagenomes/{sample}/reads/.reads_downloaded")
    params:
        read_dir = "data/omics/metagenomes/{sample}/reads/"
    conda: "config/conda_yaml/kingfisher.yaml"
    #benchmark: 
    #log:
    resources: time_min = 5000, heavy_network = 1
    shell:
        """
        export PATH=$PWD/code/kingfisher/bin:$PATH

        cd {params.read_dir}
        echo $(cat ./accession)

        kingfisher get \
            -r $(cat ./accession) \
            -m ena-ascp aws-http prefetch \
            --output-format-possibilities fastq.gz

        if (($(ls *_1.fastq.gz | wc -l)>1)); then
            echo "Error: Multiple files were downloaded."
            exit 1
        fi

        mv -f *_1.fastq.gz raw_fwd_reads.fastq.gz ||true
        mv -f *_2.fastq.gz raw_rev_reads.fastq.gz ||true

        # Download methods other than ascp can create .fastq files instead
        if (($(ls *_1.fastq | wc -l)>1)); then
            echo "Error: Multiple files were downloaded."
            exit 1
        fi

        mv -f *_1.fastq raw_fwd_reads.fastq ||true
        mv -f *_2.fastq raw_rev_reads.fastq ||true

        # If ASCP download fails, other methods can output .fastq, compress before finishing if that's the case
        if [ -e raw_fwd_reads.fastq ]; then gzip -1 raw_fwd_reads.fastq; fi
        if [ -e raw_rev_reads.fastq ]; then gzip -1 raw_rev_reads.fastq; fi
        """

rule run_get_reads:
    input:
        expand("data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz", sample = read_download_samples)

rule clumpify:
    input: 
        fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz"
    output: 
        done = touch("data/omics/metagenomes/{sample}/reads/done.touch")
    params:
        clumped_fwd_reads = "data/omics/metagenomes/{sample}/reads/clumped_raw_fwd_reads.fastq.gz",
        clumped_rev_reads = "data/omics/metagenomes/{sample}/reads/clumped_raw_rev_reads.fastq.gz",
    conda: "config/conda_yaml/main.yaml"
    benchmark:
        "benchmarks/clumpify/{sample}.txt"
    log: "logs/clumpify/{sample}_initial.log"
    resources: cpus=16, mem_mb = lambda wildcards, attempt: attempt * 175000,
        time_min=2880
        #partition = "largemem"
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
            2>&1 | tee {log}

        rm {input.fwd_reads} {input.rev_reads}
        mv {params.clumped_fwd_reads} {input.fwd_reads}
        mv {params.clumped_rev_reads} {input.rev_reads}
        """


rule deduplicate:
    input: 
        fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        clumpify = "data/omics/metagenomes/{sample}/reads/done.touch"
    output: 
        dedup_interleaved = temp("data/omics/metagenomes/{sample}/reads/dedup_interleaved.fastq.gz"),
        dedup_reads_fwd = "data/omics/metagenomes/{sample}/reads/dedup_reads_fwd.fastq.gz",
        dedup_reads_rev = "data/omics/metagenomes/{sample}/reads/dedup_reads_rev.fastq.gz"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/dedup/{sample}_dedup.log"
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


rule de_interleave:
    input: 
        reads = "import/staging/jgi_2022/all_sample_filtered_reads/{sample}_interleaved.fastq.gz"
    output: 
        reads_fwd = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        reads_rev = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/de_interleave/{sample}.log"
    # resources: cpus=36, time_min=2880,
    #     mem_mb = lambda wildcards, attempt: attempt * 170000,
    #     #partition = "largemem"
    resources: 
        #partition = "largemem",
        cpus = 24, 
        time_min = 7200,
        mem_mb = 120000
    shell:
        """
        #BBtools use more memory than given, amount given by 20% to stay within job specs.
        bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.6)/1" | bc)

        # Dedup only outputs interleaved files, this just converts back to paired
        reformat.sh in={input.reads} \
            out1={output.reads_fwd} \
            out2={output.reads_rev} 2>&1 | tee -a {log}
        """

rule JGI_SAMPLES_kraken_summarize_fastp:
    input:
        script = "code/merge_bracken.R",
        kraken_results = expand("data/omics/metagenomes/{sample}/kraken_fastp/{database}_{sample}_bracken.txt", database = ["refseq","gtdb"], sample = glob_wildcards("import/staging/jgi_2022/all_sample_filtered_reads/{sample}_interleaved.fastq.gz").sample),
        combined_tax_info = "data/reference/kraken_tax_info_merged.tsv"
    output:
        counts = "data/sample_data/JGI_bracken_counts.tsv",
        rel_abund = "data/sample_data/JGI_bracken_rel_abund.tsv"
    resources: cpus=1, mem_mb=5000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        ./{input.script} --taxonomy={input.combined_tax_info} --counts-out={output.counts} --rel-out={output.rel_abund}
        """

ruleorder: remove_contaminants_fastp > de_interleave

rule deinterleave_jgi:
    input:
        expand("data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz", sample = glob_wildcards("import/staging/jgi_2022/all_sample_filtered_reads/{sample}_interleaved.fastq.gz").sample)

rule fastp:
    input: 
        fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        clumped = "data/omics/metagenomes/{sample}/reads/done.touch"
    output: 
        tmp_fwd = temp("data/omics/metagenomes/{sample}/reads/temp_fastp_fwd_reads.fastq.gz"),
        tmp_rev = temp("data/omics/metagenomes/{sample}/reads/temp_fastp_rev_reads.fastq.gz"),
        fwd_reads = "data/omics/metagenomes/{sample}/reads/fastp_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/fastp_rev_reads.fastq.gz",
        html_dedup = "data/omics/metagenomes/{sample}/reads/qc/fastp_dedup.html",
        json_dedup = "data/omics/metagenomes/{sample}/reads/qc/fastp_dedup.json",
        html = "data/omics/metagenomes/{sample}/reads/qc/fastp.html",
        json = "data/omics/metagenomes/{sample}/reads/qc/fastp.json"
    conda: "config/conda_yaml/fastp.yaml"
    benchmark:
        "benchmarks/fastp/{sample}.txt"
    log: "logs/fastp/{sample}_fastp.log"
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

rule fastqc_fastp:
    input:
        fwd_reads = "data/omics/metagenomes/{sample}/reads/fastp_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/fastp_rev_reads.fastq.gz",
    output:
        #fwd_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/fwd_fastqc.html",
        #rev_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/rev_fastqc.html"
        touch("data/omics/metagenomes/{sample}/reads/fastqc_fastp/.done")
    conda:
          "config/conda_yaml/fastqc.yaml"
    resources: time_min = 7200, cpus = 24, mem_mb = 60000
    shell:
        """
        mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_fastp
        fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_fastp -t {resources.cpus} {input.fwd_reads}
        fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_fastp -t {resources.cpus} {input.rev_reads}
        """

rule fastqc_raw:
    input:
        fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
    output:
        #fwd_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/fwd_fastqc.html",
        #rev_report = "data/omics/metagenomes/{sample}/reads/fastqc_raw/rev_fastqc.html"
        touch("data/omics/metagenomes/{sample}/reads/fastqc_raw/.done")
    conda:
          "config/conda_yaml/fastqc.yaml"
    resources: time_min = 7200, cpus = 24, mem_mb = 60000
    shell:
        """
        mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_raw
        fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_raw -t {resources.cpus} {input.fwd_reads}
        fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_raw -t {resources.cpus} {input.rev_reads}
        """

rule multiqc:
    input: 
        "data/omics/metagenomes/{sample}/reads/fastqc_fastp/.done",
        "data/omics/metagenomes/{sample}/reads/fastqc_decontam/.done",
        "data/omics/metagenomes/{sample}/reads/fastqc_raw/.done",
        #"data/omics/metagenomes/{sample}/reads/fastqc_teal_decon/.done"
    output: 
        multiqc_dir = directory("data/omics/metagenomes/{sample}/reads/qc/multiqc")
    conda: "config/conda_yaml/multiqc.yaml"
    benchmark:
        "benchmarks/multiqc/{sample}.txt"
    log: "logs/multiqc/{sample}.log"
    resources: cpus=1, mem_mb = 20000, time_min=2880
    shell:
        """
        multiqc --interactive -d data/omics/metagenomes/{wildcards.sample}/reads/fastqc_* -o {output.multiqc_dir}
        """

rule run_fastqc_fastp:
    input: expand("data/omics/metagenomes/{sample}/reads/fastqc_fastp/.done", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/reads/qc/multiqc",sample = metaG_samples)



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
        human_genome = rules.get_contaminants.output.human_genome
    output:
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
        trimmed_fwd = "data/omics/metagenomes/{sample}/reads/trimmed_fwd_reads.fastq.gz",
        trimmed_rev = "data/omics/metagenomes/{sample}/reads/trimmed_rev_reads.fastq.gz",
        phix_rm_fwd = "data/omics/metagenomes/{sample}/reads/phix_fwd_reads.fastq.gz",
        phix_rm_rev = "data/omics/metagenomes/{sample}/reads/phix_rev_reads.fastq.gz",
        decon_fwd = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads.fastq.gz",
        decon_rev = "data/omics/metagenomes/{sample}/reads/decon_rev_reads.fastq.gz",
        cleaned_fwd = "data/omics/metagenomes/{sample}/reads/cleaned_fwd_reads.fastq.gz",
        cleaned_rev = "data/omics/metagenomes/{sample}/reads/cleaned_rev_reads.fastq.gz"
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/read_qc/{sample}.log"
    benchmark:
        "benchmarks/remove_contaminants/{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 120000, time_min = 2880
    shell:
        """        
        bbduk.sh -Xmx{resources.mem_mb}m \
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
        
        bbduk.sh -Xmx{resources.mem_mb}m \
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
        raw_reads_fwd = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        raw_reads_rev = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        deduped_reads_fwd = "data/omics/metagenomes/{sample}/reads/dedup_reads_fwd.fastq.gz",
        deduped_reads_rev = "data/omics/metagenomes/{sample}/reads/dedup_reads_rev.fastq.gz",
        qual_filt_and_trimmed_fwd = "data/omics/metagenomes/{sample}/reads/trimmed_fwd_reads.fastq.gz",
        qual_filt_and_trimmed_rev = "data/omics/metagenomes/{sample}/reads/trimmed_rev_reads.fastq.gz",
        decon_reads_fwd = "data/omics/metagenomes/{sample}/reads/cleaned_fwd_reads.fastq.gz",
        decon_reads_rev = "data/omics/metagenomes/{sample}/reads/cleaned_rev_reads.fastq.gz"
    output:
        "data/omics/metagenomes/{sample}/reads/{sample}_read_count.tsv"
    shell:
        """
        printf "read_state\tfwd_read_count\trev_read_count\n" > {output}
        printf "raw_reads\t$(($(zcat {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output}
        printf "deduped_reads\t$(($(zcat {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output}
        printf "filt_and_trimmed_reads\t$(($(zcat {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(zcat {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output}
        printf "decon_reads\t$(($(zcat {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
        """

# rule remove_spike_ins_fastp:
#     input:
#         spike_ins = rules.get_contaminants.output.spike_ins,
#         trimmed_fwd = rules.fastp.output.fwd_reads,
#         trimmed_rev = rules.fastp.output.rev_reads
#     output:
#         phix_rm_fwd = temp("data/omics/metagenomes/{sample}/reads/phix_fwd_reads_fastp.fastq.gz"),
#         phix_rm_rev = temp("data/omics/metagenomes/{sample}/reads/phix_rev_reads_fastp.fastq.gz")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/remove_spike_ins/{sample}_fastp.log"
#     params:
#         bbmap_index_path = "data/reference/contaminants"
#     benchmark:
#         "benchmarks/remove_spike_ins_fastp/{sample}.txt"
#     resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 150000, time_min = 2880
#     shell:
#         """        
#         bbmap_mem=$(echo "scale=-1; ({resources.mem_mb}*0.8)/1" | bc)
       
#         echo "Job memory= {resources.mem_mb}, bbmap allocated memory=$bbmap_mem because it is greedy"

#         echo "\n\n***Removing spike-in***\n\n" >> {log}
        
#         bbduk.sh -Xmx${{bbmap_mem}}m -eoom \
#             in1={input.trimmed_fwd} \
#             in2={input.trimmed_rev} \
#             outu1={output.phix_rm_fwd} \
#             outu2={output.phix_rm_rev} \
#             t={resources.cpus} k=31 hdist=1 \
#             ref={input.spike_ins} \
#             path={params.bbmap_index_path} \
#             1>>{log} 2>&1

#         echo "\n\n*** DONE ***\n\n" >> {log}
#         """

rule remove_contaminants_fastp:
    input:
        dedup_reads_fwd = rules.fastp.output.fwd_reads,
        dedup_reads_rev = rules.fastp.output.rev_reads,
        human_genome = rules.get_contaminants.output.human_genome,
        spike_ins = rules.get_contaminants.output.spike_ins,
        bbmap_index = "data/reference/contaminants/ref"
    output:
        phix_rm_fwd = temp("data/omics/metagenomes/{sample}/reads/phix_fwd_reads_fastp.fastq.gz"),
        phix_rm_rev = temp("data/omics/metagenomes/{sample}/reads/phix_rev_reads_fastp.fastq.gz"),
        decon_fwd = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        decon_rev = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    params:
        bbmap_index_path = "data/reference/contaminants"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/remove_contaminants_fastp/{sample}.log"
    benchmark:
        "benchmarks/remove_contaminants_fastp/{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 100000, time_min = 2880
    shell:
        """        
        bbduk.sh -Xmx{resources.mem_mb}m \
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
            in1={output.phix_rm_fwd} \
            in2={output.phix_rm_rev} \
            outu1={output.decon_fwd} \
            outu2={output.decon_rev} \
            t={resources.cpus} fast=t \
            ref={input.human_genome} \
            path={params.bbmap_index_path} \
            -Xmx{resources.mem_mb}m  \
            2>&1 | tee -a {log}
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

rule blast_Didymosphenia_geminata_chloroplast_16S:
    input:
        expand("data/omics/metagenomes/{sample}/BLAST/Didymosphenia_geminata_chloroplast_16S__{sample}.blastn", sample=metaG_samples)


rule fastqc_decontam:
    input:
        fwd_reads = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
        touch("data/omics/metagenomes/{sample}/reads/fastqc_decontam/.done")
    conda:
          "config/conda_yaml/fastqc.yaml"
    shell:
        """
        mkdir -p data/omics/metagenomes/{wildcards.sample}/reads/fastqc_decontam
        fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_decontam -t {resources.cpus} {input.fwd_reads}
        fastqc -o data/omics/metagenomes/{wildcards.sample}/reads/fastqc_decontam -t {resources.cpus} {input.rev_reads}
        """

rule run_fastqc_decontam:
    input: expand("data/omics/metagenomes/{sample}/reads/fastqc_decontam/.done", sample = metaG_samples)

rule bbnorm:
    input:
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        fwd_norm = temp("data/omics/metagenomes/{sample}/reads/bbnorm_fwd_reads.fastq.gz"),
        rev_norm = temp("data/omics/metagenomes/{sample}/reads/bbnorm_rev_reads.fastq.gz")
    params: "target=100 mindepth=2 bits=16 prefilter ecc=t"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/bbnorm/{sample}.log"
    benchmark:
        "benchmarks/bbnorm/{sample}.txt"
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
        raw_reads_fwd = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        raw_reads_rev = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        deduped_reads_fwd = "data/omics/metagenomes/{sample}/reads/fastp_fwd_reads.fastq.gz",
        deduped_reads_rev = "data/omics/metagenomes/{sample}/reads/fastp_rev_reads.fastq.gz",
        qual_filt_and_trimmed_fwd = rules.fastp.output.fwd_reads,
        qual_filt_and_trimmed_rev = rules.fastp.output.rev_reads,
        decon_reads_fwd = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        decon_reads_rev = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        #bbnorm_reads_fwd = rules.bbnorm.output.fwd_norm,
        #bbnorm_reads_rev = rules.bbnorm.output.rev_norm,
    output:
        "data/omics/metagenomes/{sample}/reads/{sample}_read_count_fastp.tsv"
    shell:
        """
        printf "read_state\tfwd_read_count\trev_read_count\n" > {output}
        printf "raw_reads\t$(($(zcat {input.raw_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.raw_reads_rev} | wc -l) / 4 ))\n" >> {output}
        printf "deduped_reads\t$(($(zcat {input.deduped_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.deduped_reads_rev} | wc -l) / 4 ))\n" >> {output}
        printf "filt_and_trimmed_reads\t$(($(zcat {input.qual_filt_and_trimmed_fwd} | wc -l) / 4 ))\t$(($(zcat {input.qual_filt_and_trimmed_rev} | wc -l) / 4 ))\n" >> {output}
        printf "decon_reads\t$(($(zcat {input.decon_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.decon_reads_rev} | wc -l) / 4 ))\n" >> {output}
        """
        #printf "bbnorm_reads\t$(($(zcat {input.bbnorm_reads_fwd} | wc -l) / 4 ))\t$(($(zcat {input.bbnorm_reads_rev} | wc -l) / 4 ))\n" >> {output}

rule run_count_reads:
    input: 
        expand("data/omics/metagenomes/{sample}/reads/{sample}_read_count.tsv", sample=metaG_samples),
        expand("data/omics/metagenomes/{sample}/reads/{sample}_read_count_fastp.tsv", sample=metaG_samples)

rule run_count_reads_fastp:
    input:
        expand("data/omics/metagenomes/{sample}/reads/{sample}_read_count_fastp.tsv", sample=metaG_samples)



rule assemble_metaspades:
    input:
        fwd_reads = rules.bbnorm.output.fwd_norm,
        rev_reads = rules.bbnorm.output.rev_norm
    output:
        assembly_dir = directory("data/omics/metagenomes/{sample}/assembly/metaspades"),
        contigs = "data/omics/metagenomes/{sample}/assembly/metaspades/contigs.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/metaspades/{sample}.log"
    benchmark: "benchmarks/metaspades/{sample}.txt"
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
        assembly_dir = directory("data/omics/metagenomes/{sample}/assembly/metaspades_noNORM"),
        contigs = "data/omics/metagenomes/{sample}/assembly/metaspades_noNORM/contigs.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/metaspades_noNORM/{sample}.log"
    benchmark: "benchmarks/metaspades_noNORM/{sample}.txt"
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

# rule assemble_megahit:
#     input:
#         fwd_reads = rules.bbnorm.output.fwd_norm,
#         rev_reads = rules.bbnorm.output.rev_norm
#     output:
#         touch("data/omics/metagenomes/{sample}/assembly/megahit/.done")
#     params:
#         assembly_dir = directory("data/omics/metagenomes/{sample}/assembly/megahit")
#     conda: "config/conda_yaml/main.yaml"
#     log: "logs/assembly/megahit/{sample}.log"
#     benchmark: "benchmarks/megahit/{sample}.txt"
#     resources: cpus = 16, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
#     shell:
#         """
#         rm -r {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
#         megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {params.assembly_dir} 2>&1 | tee {log}
#         """

# rule rename_megahit_contigs:
#     input: 
#         script = "code/rename_megahit_contigs.R",
#         assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit",
#         megahit_done = "data/omics/metagenomes/{sample}/assembly/megahit/.done"
#     output:
#         contigs = "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa",
#         done = touch("data/omics/metagenomes/{sample}/assembly/megahit/.contigs_renamed")
#     container: "docker://eandersk/r_microbiome"
#     resources: cpus = 1, time_min=200, mem_mb = 20000
#     shell:
#         """
#         pwd
#         ./{input.script} {output.contigs}
#         """

rule assemble_megahit:
    input:
        fwd_reads = rules.bbnorm.output.fwd_norm,
        rev_reads = rules.bbnorm.output.rev_norm
    output:
        contigs = "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa",
        #touch("data/omics/metagenomes/{sample}/assembly/megahit/.done")
    params:
        assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit/{sample}.log"
    benchmark: "benchmarks/megahit/{sample}.txt"
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
        contigs = "data/omics/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.fa",
        #touch("data/omics/metagenomes/{sample}/assembly/megahit/.done")
    params:
        assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit_noNORM"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit_noNORM/{sample}.log"
    benchmark: "benchmarks/megahit_noNORM/{sample}.txt"
    resources: cpus = 24, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        rm -r {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
        megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {params.assembly_dir} > {log}
        """

rule COassemble_megahit:
    input:
        fwd_reads = expand(rules.remove_contaminants_fastp.output.decon_fwd, sample = jgi_samples),
        rev_reads = expand(rules.remove_contaminants_fastp.output.decon_rev, sample = jgi_samples)
    output:
        concat_fwd = temp("tmp/fwd_concat.fastq"),
        concat_rev = temp("tmp/rev_concat.fastq"),
        contigs = "data/omics/metagenomes/jgi_coassembly/assembly/megahit/final.contigs.fa",
        #touch("data/omics/metagenomes/{sample}/assembly/megahit/.done")
    params:
        assembly_dir = "data/omics/metagenomes/jgi_coassembly/megahit"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit/coassembly.log"
    benchmark: "benchmarks/megahit/coassembly.txt"
    #resources: cpus = 92, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 150000
    resources: cpus = 32, time_min=20000, mem_mb = lambda wildcards, attempt: attempt * 1200000, partition = "largemem"
    shell:
        """
        
        zcat {input.fwd_reads} > {output.concat_fwd}
        zcat {input.rev_reads} > {output.concat_rev}

        rm -rf {params.assembly_dir} # for re-running, megahit doesn't overwrite automatically
        megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {output.concat_fwd} -2 {output.concat_rev} -o {params.assembly_dir} > {log}
        """


rule rename_megahit_contigs:
    input: 
        script = "code/rename_contigs.R",
        #assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit",
        contigs = "data/projects/{project}/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.fa",
        #assembly_done = "data/omics/metagenomes/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/projects/{project}/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        contig_info = "data/projects/{project}/metagenomes/{sample}/assembly/megahit_noNORM/contigs_info.tsv"
        #done = touch("data/omics/metagenomes/{sample}/assembly/megahit/.contigs_renamed")
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

rule rename_contigs:
    input: 
        script = "code/rename_contigs.R",
        #assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit",
        contigs = "data/omics/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.fa",
        #assembly_done = "data/omics/metagenomes/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/omics/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        contig_info = "data/omics/metagenomes/{sample}/assembly/megahit_noNORM/contigs_info.tsv"
        #done = touch("data/omics/metagenomes/{sample}/assembly/megahit/.contigs_renamed")
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

rule rename_metaspades_contigs:
    input: 
        script = "code/rename_contigs.R",
        contigs = "data/projects/{project}/metagenomes/{sample}/assembly/metaspades_noNORM/contigs.fasta",
        #assembly_done = "data/omics/metagenomes/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/projects/{project}/metagenomes/{sample}/assembly/metaspades_noNORM/contigs.renamed.fasta",
        contig_info = "data/projects/{project}/metagenomes/{sample}/assembly/metaspades_noNORM/contigs_info.tsv"
        #done = touch("data/omics/metagenomes/{sample}/assembly/megahit/.contigs_renamed")
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


rule merge_assemblies:
    input:
        metaspades_contigs = rules.assemble_metaspades.output.contigs,
        megahit_contigs = rules.rename_megahit_contigs.output.contigs,
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

rule run_quast:
    input: expand("data/omics/metagenomes/{sample}/assembly/quast", sample = metaG_samples)

rule prodigal:
    input:
        assembly = rules.rename_contigs.output.contigs
    output:
        proteins = "data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa",
        genes = "data/omics/metagenomes/{sample}/genes/{sample}_GENES.fna"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/prodigal/{sample}.log"
    benchmark: "benchmarks/prodigal/{sample}.txt"
    resources: cpus = 1, mem_mb = lambda wildcards, attempt: attempt * 16000, time_min = 2880
    shell:
        """
        prodigal -p meta -i {input.assembly} -a {output.proteins} -d {output.genes} 2>&1 | tee {log}
        """


rule calc_gene_abundance:
    input:
        genes = rules.prodigal.output.genes,
        proteins = rules.prodigal.output.proteins,
        fwd_reads = rules.remove_contaminants_fastp.output.decon_fwd,
        rev_reads = rules.remove_contaminants_fastp.output.decon_rev
    output:
        reads_vs_genes_rpkm = "data/omics/metagenomes/{sample}/genes/{sample}_READSvsGENES.rpkm",
        reads_vs_contigs_rpkm = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.rpkm",
        reads_vs_assembly_sam_gz = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam.gz"
    params:
        reads_vs_assembly_sam = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/calc_gene_abundance/{sample}.log"
    benchmark: "benchmarks/calc_gene_abundance/{sample}.txt"
    resources: cpus = 24, mem_mb = lambda wildcards, attempt: attempt * 64000, time_min = 2880
    shell:
        """
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {output.genes}) ref={output.genes} rpkm={output.reads_vs_genes_rpkm} 2>&1 | tee -a {log}
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.assembly}) ref={input.assembly} rpkm={output.reads_vs_contigs_rpkm} 32bit=t outm={params.reads_vs_assembly_sam} 2>&1 | tee -a {log}
        gzip {params.reads_vs_assembly_sam}
        """

rule download_uniref:
    output: 
        uniref100="data/reference/uniref/uniref100.fasta.gz"
    #conda: "config/conda_yaml/main.yaml"
    log: "logs/make_diamond_uniref_db/download_uniref.log"
    resources: cpus = 1, mem_mb=1000
    shell:
        """
        #cd data/reference/uniref
        #wget -N https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
        
        #mkdir data/reference/uniref
        cp data/reference/uniref100.fasta.gz data/reference/uniref/
        """

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
        gene_uniref_alignment = "data/omics/metagenomes/{sample}/{sample}_GENES.m8"
    params:
        "--top 0.5 --threads 10 --query-cover 50 --strand both -f 6 qseqid qlen sseqid slen qstart qend sstart send evalue pident mismatch qcovhsp scovhsp"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/align_to_uniref/{sample}_align_to_uniref.log"
    benchmark: "benchmarks/align_to_uniref/{sample}_align_to_uniref.txt"
    resources: cpus = 8, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        diamond blastx \
            -d {input.diamond_db} \
            -q {input.genes} \
            -o {output.gene_uniref_alignment} \
            {params} 2>&1 | tee {log}
        """

rule make_uniref_alignments:
    input: expand("data/omics/metagenomes/{sample}/{sample}_GENES.m8",sample = metaG_samples)

# Required input for the annotation script: 
# print "\t1. The gene or protein sequences:                      [sample]_GENES.fna\n";
# print "\t2. The alignment file of the sequences to UMRAD:       [sample]_GENES.m8 \n";
# print "\t3. The rpkm file for reads aligned to the genes:       [sample]_READSvsGENES.rpkm\n";
# print "\t4. The rpkm file for reads aligned to the contigs:     [sample]_READSvsCONTIGS.rpkm\n";
# print "\t5. The sample contigs sequences:                       [sample]_MCDD.fa\n\n";


#             /TAXONOMY\_DB.*$year\.txt/){	$intax	=$refdir.$file;}
# if($file =~ /UNIREF100\_INFO.*$year\.txt/){	$ininfo	=$refdir.$file;}
# if($file =~ /Function\_Names.*\.txt/){	


# Target rules for running kraken
rule calc_kraken_uniq:
    input: expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_out.txt", database = ["refseq","gtdb_r202"], sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_brackenMpa.txt", database = ["refseq","gtdb_r202"], sample = metaG_samples)

rule calc_kraken_uniq_gtdb:
    input: expand("data/omics/metagenomes/{sample}/kraken/{database}_{sample}_brackenMpa.txt", database = ["gtdb_r202"], sample = metaG_samples)

# Inspect Kraken2 databases, nescessary to produce nicely formatted taxonomy files for downstream analysis
rule kraken_inspect:
    input:
        db = ancient("/geomicro/data2/kiledal/references/kraken_databases/{database}")
    output: 
        inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/{database}/inspect.txt"
    conda: "config/conda_yaml/kraken.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        kraken2-inspect --db {input.db} > {output.inspect_file}
        """

rule add_lineage_to_inspect_gtdb:
    input:
        db = ancient("/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202"),
        inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect.txt"
    output: 
        inspect_w_lineage = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect_w_lineage.txt"
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
        db = ancient("/geomicro/data2/kiledal/references/kraken_databases/refseq"),
        inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect.txt"
    output: 
        inspect_w_lineage_unformatted = temp("/geomicro/data2/kiledal/references/kraken_databases/refseq/unformatted_inspect_w_lineage.txt"),
        inspect_w_lineage = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect_w_lineage.txt"
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
        # f_seq = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        # r_seq = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        db = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202",
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_report.txt",
        out = temp("data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_out.txt"),
        bracken = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_bracken.txt",
        bracken_report = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_brackenMpa.txt",
        unclass_f = temp("data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_unclassified_1.fasta"),
        unclass_r = temp("data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_unclassified_2.fasta"),
        bracken_input = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150,
        unclass_out = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_unclassified#.fasta"
    conda: "config/conda_yaml/kraken.yaml"
    log: "logs/kraken2_gtdb_w_uniq/{sample}.log"
    benchmark: "benchmarks/kraken2_gtdb_w_uniq/{sample}.txt"
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
        db = "/geomicro/data2/kiledal/references/kraken_databases/refseq",
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_report.txt",
        out = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_out.txt",
        bracken = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_bracken.txt",
        bracken_report = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_brackenMpa.txt",
        bracken_input = "data/omics/metagenomes/{sample}/kraken/refseq_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150
    conda: "config/conda_yaml/kraken.yaml"
    log: "logs/kraken2_refseq_w_uniq/{sample}.log"
    benchmark: "benchmarks/kraken2_refseq_w_uniq/{sample}.txt"
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
        db = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202"
    output:
        #temp(service("/dev/shm/gtdb_r202"))
        temp(directory("/dev/shm/gtdb_r202"))
    resources: cpus=1, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        mkdir {output}
        cp -r {input.db} /dev/shm/
        """

rule kraken2_gtdb_w_uniq_fastp:
    input:
        f_seq = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        # f_seq = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        # r_seq = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz",
        db = rules.kraken2_load_gtdb_DB.output,
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_report.txt",
        out = temp("data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_out.txt"),
        bracken = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_bracken.txt",
        bracken_report = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_brackenMpa.txt",
        unclass_f = temp("data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_unclassified_1.fasta"),
        unclass_r = temp("data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_unclassified_2.fasta"),
        bracken_input = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150,
        unclass_out = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_unclassified#.fasta"
    conda: "config/conda_yaml/kraken.yaml"
    benchmark: "benchmarks/kraken2_gtdb_w_uniq_fastp/{sample}.txt"
    resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
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
        db = "/geomicro/data2/kiledal/references/kraken_databases/refseq"
    output:
        #temp(service("/dev/shm/refseq"))
        temp(directory("/dev/shm/refseq"))
    resources: cpus=1, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
    shell:
        """
        mkdir {output}
        cp -r {input.db} /dev/shm/
        """

# Any reads not annotated with the GTDB database are then annotated with a RefSeq database
rule kraken2_refseq_w_uniq_fastp: ##Run kraken2
    input:
        f_seq = rules.kraken2_gtdb_w_uniq_fastp.output.unclass_f,
        r_seq = rules.kraken2_gtdb_w_uniq_fastp.output.unclass_r,
        db = rules.kraken2_load_refseq_DB.output,
        kreport2mpa = "code/kreport2mpa.py"
    output:
        report = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_report.txt",
        out = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_out.txt",
        bracken = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_bracken.txt",
        bracken_report = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_brackenReport.txt",
        bracken_mpa = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_brackenMpa.txt",
        bracken_input = "data/omics/metagenomes/{sample}/kraken_fastp/refseq_{sample}_for_bracken.txt"
    params:
        uniq_minimizer_threshold = 150
    benchmark: "benchmarks/kraken2_refseq_w_uniq_fastp/{sample}.txt"
    conda: "config/conda_yaml/kraken.yaml"
    resources: cpus=16, mem_mb=250000, time_min=1440, mem_gb = 250, partition = "largemem"
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

# Combine the kraken annotations and produce count table
rule kraken_summarize_fastp:
    input:
        script = "code/merge_bracken.R",
        kraken_results = expand("data/omics/metagenomes/{sample}/kraken_fastp/{database}_{sample}_bracken.txt", database = ["refseq","gtdb"], sample = metaG_samples),
        combined_tax_info = rules.kraken_database_tax_merge.output.combined_tax_info
    output:
        counts = "data/sample_data/bracken_counts.tsv",
        rel_abund = "data/sample_data/bracken_rel_abund.tsv"
    resources: cpus=1, mem_mb=5000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        ./{input.script} --taxonomy={input.combined_tax_info} --counts-out={output.counts} --rel-out={output.rel_abund}
        """

################


# Target rule to make all the metacodeR plots
rule plot_metacoders:
    input: expand("data/omics/metagenomes/{sample}/kraken/{sample}_kraken_metacodeR.pdf", sample=metaG_samples)


rule metacodeR:
    input:
        script = "code/plot_metacoder.R",
        abund = "data/sample_data/bracken_rel_abund.tsv"
        #metadata = "data/metadata.tsv"
    output: "data/omics/metagenomes/{sample}/kraken/{sample}_kraken_metacodeR.pdf"
    resources: cpus=1, mem_mb=8000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        {input.script} --abund={input.abund} --sample={wildcards.sample} --output={output}
        """


rule kofam_scan:
    input:
        genes = rules.prodigal.output.genes,
        profile = "data/reference/kegg/kofamscan/profiles",
        ko_list = "data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/omics/metagenomes/{sample}/kofam_scan/{sample}_kofam_results.txt"
    conda: "config/conda_yaml/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofamscan/{sample}.txt"
    log: "logs/kofamscan/{sample}.log"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --cpu={resources.cpus}  \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.sample}_kofamscan \
            --ko-list {input.ko_list} {input.genes}
        """
        

rule run_kofam_scan:
    input:
        expand("data/omics/metagenomes/{sample}/kofam_scan/{sample}_kofam_results.txt", sample = metaG_samples)


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
        done = touch("data/omics/metagenomes/{sample}/.annotation_done"),
    params:
        annotation_dir = directory("data/omics/metagenomes/{sample}/annotation")
    conda: "config/conda_yaml/main.yaml"
    resources: cpus=1, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 170000
    log: "logs/annotate_contigs/{sample}.log"
    benchmark: "benchmarks/annotate_contigs/{sample}.log"
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

#metaG_samples = "66737388b7e4c974ed413246b4408f1f"

rule metaG_annotation:
    input: 
        #expand("data/omics/metagenomes/{sample}/annotation", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/.annotation_done", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/reads/{sample}_read_count_fastp.tsv", sample = metaG_samples)
        #expand("data/omics/metagenomes/{sample}/bins/metabat2_bins", sample = metaG_samples),
        #expand("data/omics/metagenomes/{sample}/bins/gtdbtk", sample = metaG_samples),
        #expand("data/omics/metagenomes/{sample}/bins/checkm.txt", sample = metaG_samples)


# rule initial_binning:
#     input:
#     output:
#     conda: "config/conda_yaml/main.yaml"
#     resources: 
#     shell:
#         """
#        while read i; do echo "doing $i"; 
#             #metabat
#             awk -F'\t' '{print $1"\t"$2"\t"$6"\t"$6"\t"$6"\t"}' ${i}_READSvsCONTIGS.rpkm > ${i}_MB_abund.txt
#             metabat2 --maxP=97 --minS=93 --maxEdges=300 -m 1500 -i ${i}_MERGED_CONTIGS_COR.fasta -a ${i}_MB_abund.txt -o ${i}_MB_bins
#             metabat2 --maxP=99 --minS=97 --maxEdges=300 -m 1500 -i ${i}_MERGED_CONTIGS_COR.fasta -a ${i}_MB_abund.txt -o ${i}_MB_bins
#             #maxbin
#             awk -F'\t' '{print $1"\t"$6}' ${i}_READSvsCONTIGS.rpkm | grep "CLUSTER" > ${i}_MX_abund.txt
#             perl run_MaxBin.pl -thread 20 -contig ${i}_MERGED_CONTIGS_COR.fasta -abund ${i}_MX_abund.txt -out ${i}_MX_bins
#         done < samp_list.txt;
#         """



rule make_bins:
    input:
        #contigs = rules.MEC.output.corrected_assembly,
        contigs = rules.rename_megahit_contigs.output.contigs,
        fwd_reads = rules.remove_contaminants.output.decon_fwd,
        rev_reads = rules.remove_contaminants.output.decon_rev
    output:
        directory("data/omics/metagenomes/{sample}/bins/concoct_bins"),
        directory("data/omics/metagenomes/{sample}/bins/maxbin2_bins"),
        directory("data/omics/metagenomes/{sample}/bins/metabat2_bins"),
        directory("data/omics/metagenomes/{sample}/bins/work_files"),
        "data/omics/metagenomes/{sample}/bins/work_files/assembly.fa",
        fwd_reads = temp("data/omics/metagenomes/{sample}/reads/reads_1.fastq"),
        rev_reads = temp("data/omics/metagenomes/{sample}/reads/reads_2.fastq"),
        #out_dir = directory("data/omics/metagenomes/{sample}/bins")
    params:
        out_dir = "data/omics/metagenomes/{sample}/bins"
    conda: "config/conda_yaml/metawrap.yaml"
    resources: cpus=16, mem_mb=50000, time_min=2880, mem_gb = 50
    shell:
        """
        WORK_DIR=$PWD

        #data/omics/metagenomes/{wildcards.sample}/bins
        #mkdir -p data/omics/metagenomes/{wildcards.sample}/bins
        #cp data/qc_sequence_files/{wildcards.sample}[a-z]_R[12].fastq.gz data/bins/{wildcards.sample}/metaspades/
        gunzip -c {input.fwd_reads} > {output.fwd_reads}
        gunzip -c {input.rev_reads} > {output.rev_reads}
 
        #cd data/bins/{wildcards.sample}/metaspades
        #rename _R1.fastq _1.fastq *_R1.fastq
        #rename _R2.fastq _2.fastq *_R2.fastq

        #cd $WORK_DIR

        metawrap binning \
            -a {input.contigs} \
            -o {params.out_dir} \
            -t {resources.cpus} \
            -m {resources.mem_gb} \
            --metabat2 \
            --maxbin2 \
            --concoct \
            --universal \
            {output.fwd_reads} \
            {output.rev_reads}

        #rm -f data/bins/{wildcards.sample}/metaspades/*.fastq
        """


rule dastool:
    input:
        contigs = "data/omics/metagenomes/{sample}/bins/work_files/assembly.fa",
        #bin_folder = rules.make_bins.params.out_dir
        bins = "data/omics/metagenomes/{sample}/bins/metabat2_bins"
    params:
        bin_folder = rules.make_bins.params.out_dir
    output: 
        #summary = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_summary.txt",
        das_folder = directory("data/omics/metagenomes/{sample}/bins/DASTool"),
        das_bins_folder = directory("data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_bins")
    conda: "config/conda_yaml/das_tool.yaml"
    resources: cpus=8, mem_mb=50000, time_min=2880, mem_gb = 50
    shell:
        """
        touch {params.bin_folder}/ran_dastool.touch # for tracking that DAStool ran, even if unsuccessfully

        #Maxbin2
        Fasta_to_Contig2Bin.sh \
        -i {params.bin_folder}/maxbin2_bins \
        -e fa \
        > {params.bin_folder}/maxbin.scaffolds2bin.tsv

        #CONCOT
        Fasta_to_Contig2Bin.sh \
        -i {params.bin_folder}/concoct_bins \
        -e fa \
        > {params.bin_folder}/concoct.scaffolds2bin.tsv

        #Metabat2
        Fasta_to_Contig2Bin.sh \
        -i {params.bin_folder}/metabat2_bins \
        -e fa \
        > {params.bin_folder}/metabat2.scaffolds2bin.tsv

        DAS_Tool \
            -i {params.bin_folder}/maxbin.scaffolds2bin.tsv,{params.bin_folder}/concoct.scaffolds2bin.tsv,{params.bin_folder}/metabat2.scaffolds2bin.tsv \
            -l maxbin,concoct,metabat2 \
            -c {input.contigs} \
            -t {resources.cpus} \
            --write_bins \
            -o {output.das_folder}/

        # Rename bins to include sample name, make several downstream analyses easier
        cd {output.das_bins_folder}
        for f in *.fa; do mv -v -- "$f" "{wildcards.sample}_$f"; done
        """


rule checkm:
    input: rules.dastool.output.das_bins_folder
    output:
        dir = temp(directory("data/omics/metagenomes/{sample}/bins/checkm")),
        results = "data/omics/metagenomes/{sample}/bins/checkm.txt"
    conda: "config/conda_yaml/checkm.yaml"
    resources: cpus=8, mem_mb=80000, time_min=2880, mem_gb = 80
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {input} {output.dir}
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


rule gtdbtk:
    input:
        bins = rules.dastool.output.das_bins_folder,
        refs = "/geomicro/data2/kiledal/references/gtdbtk/release202"
    output: directory("data/omics/metagenomes/{sample}/bins/gtdbtk")
    conda: "config/conda_yaml/gtdbtk.yaml"
    resources: cpus=1, mem_mb=500000, time_min=2880, mem_gb = 500
    shell:
        """
        GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fa --genome_dir {input.bins} --out_dir {output} --cpus {resources.cpus}
        """

rule mag_coverage:
    input:
        fwd_reads = rules.remove_contaminants.output.cleaned_fwd,
        rev_reads = rules.remove_contaminants.output.cleaned_rev,
        bins = rules.dastool.output.das_bins_folder
    output: "data/omics/metagenomes/{sample}/bins/coverage.tsv"
    conda: "config/conda_yaml/coverm_env.yaml"
    resources: cpus=16, mem_mb=250000, time_min=2880
    shell:
        """
        coverm genome \
            -t {resources.cpus} \
            -m relative_abundance mean covered_bases variance length \
            --min-covered-fraction 0 \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            --genome-fasta-files {input.bins}/*.fa \
            -o {output}
        """

rule map_to_bins_for_reassembly:
    input:
        fwd_reads = rules.remove_contaminants.output.cleaned_fwd,
        rev_reads = rules.remove_contaminants.output.cleaned_rev,
        bins_folder = rules.dastool.output.das_bins_folder,
        bin = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_bins/{bin}.fa"
    output:
        #unfiltered_bam = "data/omics/metagenomes/{sample}/bins/reassembly/reads/{bin}_prefilt.bam",
        unfiltered_bam_dir = temp(directory("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}__bam_unfiltered")),
        filtered_bam = temp("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}.bam"),
        #filtered_bam_dir = directory("data/omics/metagenomes/{sample}/bins/reassembly/reads/{bin}__bam_filtered"),
        fwd_reads = temp("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}_R1.fastq.gz"),
        rev_reads = temp("data/omics/metagenomes/{sample}/bins/reassembly/mapped_reads/{bin}_R2.fastq.gz")
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=24, mem_mb=100000, time_min=2880
    shell:
        """
        coverm make \
            -t {resources.cpus} \
            -1 {input.fwd_reads} \
            -2 {input.rev_reads} \
            --reference {input.bin} \
            -o {output.unfiltered_bam_dir}

        coverm filter \
            -b {output.unfiltered_bam_dir}/*.bam \
            -o {output.filtered_bam} \
            --min-read-percent-identity 0.98 \
            --threads {resources.cpus}

        samtools fastq -@ {resources.cpus} \
            {output.filtered_bam} \
            -1 {output.fwd_reads} \
            -2 {output.rev_reads} \
            -0 /dev/null -s /dev/null -n
        """

rule reassemble_bin:
    input:
        fwd_reads = rules.map_to_bins_for_reassembly.output.fwd_reads,
        rev_reads = rules.map_to_bins_for_reassembly.output.rev_reads,
        bin = "data/omics/metagenomes/{sample}/bins/DASTool/_DASTool_bins/{bin}.fa"
    output:
        #assembly_dir = directory("data/omics/metagenomes/{sample}/bins/reassembly/{bin}"),
        contigs = "data/omics/metagenomes/{sample}/bins/reassembly/{bin}_reassembled_contigs.fasta"
    params: 
        assembly_dir = directory("data/omics/metagenomes/{sample}/bins/reassembly/{bin}")
    conda: "config/conda_yaml/main.yaml"
    log: "logs/bin_reassembly/{sample}/{bin}.log"
    benchmark: "logs/bin_reassembly/{sample}/{bin}.log"
    resources: cpus = 16, time_min=20160, mem_mb = 16000
        #mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        spades.py \
            -t {resources.cpus} \
            -m $(echo "scale=-1; ({resources.mem_mb}/1000)/1" | bc) \
            --careful \
		    --untrusted-contigs {input.bin} \
		    -1 {input.fwd_reads} \
		    -2 {input.rev_reads} \
		    -o {params.assembly_dir} > {log}

        mv {params.assembly_dir}/contigs.fasta {output.contigs}
        rm -r {params.assembly_dir}
        """

rule gather_bin:
    input: expand("data/omics/metagenomes/{sample}/bins/ran_dastool.touch",sample=metaG_samples)
    output: directory("data/omics/metagenome_bins")
    shell:
        """
        mkdir -p {output}
        
        ln data/omics/metagenomes/*/bins/DASTool/_DASTool_bins/*.fa {output}
        """

rule binning:
    input: 
        #"data/omics/metagenome_bins",
        #expand("data/omics/metagenomes/{sample}/bins/reassembly/{bin}_reassembled_contigs.fasta", sample = metaG_samples, bin = glob_wildcards("data/omics/metagenomes/{sample}/bins/reassembly/{bin}_reassembled_contigs.fasta").bin),
        "data/omics/metagenomes/0ae7941e0cc52de7e4913cf5defec020/bins/reassembly/0ae7941e0cc52de7e4913cf5defec020_concoct_bin.16_reassembled_contigs.fasta",
        expand("data/omics/metagenomes/{sample}/bins/coverage.tsv", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/bins/checkm.txt", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/bins/gtdbtk", sample = metaG_samples)


# rule gtdbtk_all:
#     input:
#         bins = rules.gather_bin.output,
#         refs = "/geomicro/data2/kiledal/references/gtdbtk/release202"
#     output: directory("data/gtdbtk")
#     conda: "code/gtdbtk.yaml"
#     resources: cpus=32, mem_mb=250000, time_min=2880, mem_gb = 250
#     shell:
#         """
#         GTDBTK_DATA_PATH={input.refs}

#         gtdbtk classify_wf --extension fa --genome_dir {input.bins} --out_dir {output} --cpus {resources.cpus} --pplacer_cpus 1
#         """

checkpoint drep:
    input: 
        bins = rules.gather_bin.output
    output:
        main_dir = directory("data/omics/metagenome_bins/derep"),
        MAGs= directory("data/omics/metagenome_bins/derep/dereplicated_genomes")
    conda: "config/conda_yaml/drep.yaml"
    resources: cpus=8, mem_mb=250000, time_min=2880,
    shell:
        """
        dRep dereplicate {output.main_dir} -g {input.bins}/*.fa
        """


rule humann:
    input:
        # Sample data #
        f_seq = rules.remove_contaminants.output.cleaned_fwd,
        r_seq = rules.remove_contaminants.output.cleaned_rev,
        bracken_mpa = "data/omics/metagenomes/{sample}/kraken/gtdb_{sample}_brackenMpa.txt",
        # Reference database #
        NUC_DB = "data/reference/humann/genome_reps_filt_annot.fna.gz",
        PROT_DB = "data/reference/humann/protein_database/uniref90_201901b.dmnd",
        NUC_fol = "data/reference/humann/",
        PROT_fol = "data/reference/humann/protein_database/"
    output:
        humann_output = directory("data/omics/metagenomes/{sample}/humann"),
        concat_unzipped_reads = temp("data/omics/metagenomes/{sample}/reads/for_humann.fastq")
    params:
        mem_use = "maximum"
    conda: "config/conda_yaml/humann.yaml"
    log: "logs/humann/{sample}.log"
    benchmark: "benchmarks/humann/{sample}.txt"
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

rule run_humann: 
    input: expand("data/omics/metagenomes/{sample}/humann", sample=metaG_samples)


rule humann_fastp:
    input:
        # Sample data #
        f_seq = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
        r_seq = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz",
        bracken_mpa = "data/omics/metagenomes/{sample}/kraken_fastp/gtdb_{sample}_brackenMpa.txt",
        # Reference database #
        NUC_DB = "/home/kiledal/scratch_gdick1/GVHD/data/reference/humann/genome_reps_filt_annot.fna.gz",
        PROT_DB = "/home/kiledal/scratch_gdick1/GVHD/data/reference/humann/protein_database/uniref90_201901b.dmnd",
        NUC_fol = "/home/kiledal/scratch_gdick1/GVHD/data/reference/humann/",
        PROT_fol = "/home/kiledal/scratch_gdick1/GVHD/data/reference/humann/protein_database/"
    output:
        humann_output = directory("data/omics/metagenomes/{sample}/humann_fastp"),
        concat_unzipped_reads = temp("data/omics/metagenomes/{sample}/reads/for_humann_fastp.fastq")
    params:
        mem_use = "maximum"
    conda: "config/conda_yaml/humann.yaml"
    log: "logs/humann/{sample}_fastp.log"
    benchmark: "benchmarks/humann/{sample}_fastp.txt"
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

rule run_humann_fastp: 
    input: expand("data/omics/metagenomes/{sample}/humann_fastp", sample=metaG_samples)


rule sourmash_sketch:
    input:
       f_seq = "data/omics/metagenomes/{sample}/reads/decon_fwd_reads_fastp.fastq.gz",
       r_seq = "data/omics/metagenomes/{sample}/reads/decon_rev_reads_fastp.fastq.gz"
    output:
       sig = "data/omics/metagenomes/{sample}/sourmash/{sample}.sig"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash_sketch/{sample}.log"
    benchmark: "benchmarks/sourmash_sketch/{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 20000
    shell:
        """
        sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge {wildcards.sample} -o {output.sig} {input.f_seq} {input.r_seq} 2>&1 | tee {log}
        """

rule sourmash_gather:
    input:
       sig = "data/omics/metagenomes/{sample}/sourmash/{sample}.sig",
       gtdb_refDB = "data/reference/sourmash/gtdb-rs207.dna.k31.zip",
       microcystis_refDB = "data/reference/sourmash/Microcystis_sigs.sig",
       taxDB = "data/reference/sourmash/gtdb_and_Microcystis_tax.db"
    output:
       reps = "data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.csv",
       tax = "data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.with-lineages.csv"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash/{sample}.log"
    benchmark: "benchmarks/sourmash/{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 20000
    shell:
        """
        sourmash gather {input.sig} {input.gtdb_refDB} {input.microcystis_refDB} -o {output.reps} 2>&1 | tee -a {log}

        sourmash tax annotate -g {output.reps} -t {input.taxDB} 2>&1 | tee -a {log}

        mv $(basename {output.tax}) data/omics/metagenomes/{wildcards.sample}/sourmash/
        """


rule run_sourmash: 
    input: expand("data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.csv", sample=metaG_samples)


rule sourmash_gather_Microcystis:
    input:
       sig = "data/omics/metagenomes/{sample}/sourmash/{sample}.sig",
       refDB = "data/reference/sourmash/Microcystis_sigs.sig"
    output:
       reps = "data/omics/metagenomes/{sample}/sourmash/{sample}_gather_Microcystis.csv",
       #tax = "data/omics/metagenomes/{sample}/sourmash/{sample}_gather_gtdbrs207_reps.with-lineages.csv"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash/{sample}.log"
    benchmark: "benchmarks/sourmash/{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 20000
    shell:
        """
        sourmash gather {input.sig} {input.refDB} -o {output.reps} 2>&1 | tee -a {log}

        #mv $(basename {output.reps}) data/omics/metagenomes/{wildcards.sample}/sourmash/
        """

rule run_sourmash_Microcystis: 
    input: expand("data/omics/metagenomes/{sample}/sourmash/{sample}_gather_Microcystis.csv", sample=metaG_samples)


# JGI seqs unprocessed

rule sourmash_gather_JGI:
    input:
       f_seq = "{sample}.fastq.gz",
       refDB = "data/reference/sourmash/gtdb-rs207.dna.k31.zip",
       taxDB = "data/reference/sourmash/gtdb-rs207.taxonomy.sqldb"
    output:
       sig = "{sample}.sig",
       reps = "{sample}_gather_gtdbrs207_reps.csv",
       tax = "{sample}_gather_gtdbrs207_reps.with-lineages.csv"
    conda: "config/conda_yaml/sourmash.yaml"
    log: "logs/sourmash_jgi/{sample}.log"
    benchmark: "benchmarks/sourmash_jgi/{sample}.txt"
    resources: cpus=1, time_min=4320, mem_mb = 20000
    shell:
        """
        sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge {wildcards.sample} -o {output.sig} {input.f_seq} 2>&1 | tee {log}

        sourmash gather {output.sig} {input.refDB} -o {output.reps} 2>&1 | tee -a {log}

        sourmash tax annotate -g {output.reps} -t {input.taxDB} 2>&1 | tee -a {log}

        mv $(basename {output.tax}) import/staging/jgi_2022/all_sample_filtered_reads/
        """


rule run_sourmash_jgi: 
    input: expand("import/staging/jgi_2022/all_sample_filtered_reads/{sample}_gather_gtdbrs207_reps.csv", sample=glob_wildcards("import/staging/jgi_2022/all_sample_filtered_reads/{sample}.fastq.gz").sample)





## New binning


rule link_reads_w_sample_names:
    input: "data/projects/{project}/metagenomes/{sample}/reads/decon_{dir}_reads_fastp.fastq.gz"
    output: "data/projects/{project}/metagenomes/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz"
    resources: cpus=1, mem_mb = 500
    shell:
        """
        ln {input} {output}
        """


rule map_to_contigs:
    input:
        expand("data/projects/{project}/metagenomes/{sample}/reads/fastp_decon/{sample}_{dir}.fastq.gz", sample = metaG_samples, dir = ["fwd", "rev"], project = "2022_geomicro_JGI_CSP"),
        contigs = rules.rename_megahit_contigs.output.contigs
    output: 
        bam_dir = directory(temp("data/projects/{project}/metagenomes/{sample}/bins/bam"))
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=16, mem_mb=50000, time_min=2880, mem_gb = 50
    shell:
        """
        coverm make -c data/projects/{wildcards.project}/metagenomes/*/reads/fastp_decon/*.fastq.gz \
            -r {input.contigs} \
            --discard-unmapped \
            -t {resources.cpus} \
            -o {output.bam_dir} 
        """

rule run_contig_map:
    input: expand("data/omics/metagenomes/{sample}/bins/bam", sample = metaG_samples)

rule contig_coverage:
    input:
        bam_dir = "data/projects/{project}/metagenomes/{sample}/bins/bam"
    output: 
        coverage = "data/projects/{project}/metagenomes/{sample}/bins/contig_coverage.tsv",
        coverage_metabat = "data/projects/{project}/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=8, mem_mb=100000, time_min=2880
    shell:
        """
        coverm contig \
            -b {input.bam_dir}/*.bam \
            -t {resources.cpus} \
            --output-file {output.coverage}

        coverm contig \
            -b {input.bam_dir}/*.bam \
            -t {resources.cpus} \
            --methods metabat \
            --output-file {output.coverage_metabat}
        """


rule index_contig_coverage:
    input:
        bam_dir = "data/projects/{project}/metagenomes/{sample}/bins/bam"
    output: 
        index_done = touch("data/projects/{project}/metagenomes/{sample}/bins/.bam_indexed")
    params:
        bam = "data/projects/{project}/metagenomes/{sample}/bins/bam/final.contigs.renamed.fa.decon_fwd_reads_fastp.fastq.gz.bam"
    conda: "config/conda_yaml/coverm.yaml"
    resources: cpus=4, mem_mb=20000, time_min=2880
    shell:
        """
        parallel -j {resources.cpus} samtools index -@ 1 ::: {input.bam_dir}/*.bam
        """

rule concoct:
    input:
        contigs = rules.rename_megahit_contigs.output.contigs,
        bam_index = rules.index_contig_coverage.output.index_done,
        bam_dir = "data/projects/{project}/metagenomes/{sample}/bins/bam"
    output:
        cut_contigs = "data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa",
        cut_contigs_bed = "data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/contigs_10K.bed",
        cut_coverage = "data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/cut_coverage_table.tsv"
    params:
        outdir = "data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/output",
        bam = "data/projects/{project}/metagenomes/{sample}/bins/bam/*.bam"
    benchmark: "benchmarks/concoct/{project}__{sample}.txt"
    conda: "config/conda_yaml/concoct.yaml"
    resources: cpus=16, mem_mb=100000, time_min=10080, mem_gb = 50
    shell:
        """
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {output.cut_contigs_bed} > {output.cut_contigs}
        
        concoct_coverage_table.py {output.cut_contigs_bed} {params.bam} > {output.cut_coverage}

        concoct --threads {resources.cpus} --composition_file {output.cut_contigs} --coverage_file {output.cut_coverage} -b {params.outdir}/
        
        merge_cutup_clustering.py {params.outdir}/clustering_gt1000.csv > {params.outdir}/clustering_merged.csv

        mkdir -p {params.outdir}/fasta_bins
        extract_fasta_bins.py {input.contigs} {params.outdir}/clustering_merged.csv --output_path {params.outdir}/fasta_bins
        """


rule metabat2:
    input:
        contigs = rules.rename_megahit_contigs.output.contigs,
        bam_index = rules.index_contig_coverage.output.index_done,
        bam_dir = "data/projects/{project}/metagenomes/{sample}/bins/bam",
        coverm_depth = "data/projects/{project}/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
    output:
        #depth = "data/omics/metagenomes/{sample}/bins/jgi_depth_summary.txt",
        done = touch("data/projects/{project}/metagenomes/{sample}/bins/METABAT2/.done")
    params:
        bin_name = directory("data/projects/{project}/metagenomes/{sample}/bins/METABAT2/metabat2")
    benchmark: "benchmarks/metabat2/{project}__{sample}.txt"
    singularity: "docker://metabat/metabat"
    resources: cpus=16, mem_mb=20000, time_min=2880, mem_gb = 50
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
        coverm_depth = "data/projects/{project}/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
    output:
        depths_file = "data/projects/{project}/metagenomes/{sample}/bins/maxbin/depths.txt"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb=50000, time_min=1000
    shell:
        """
        cd {current_dir}
        pwd

        ./{input.script} {input.coverm_depth}
        """

rule maxbin2:
    input:
        contigs = rules.rename_megahit_contigs.output.contigs,
        depth = rules.maxbin2_coverage.output.depths_file
    output:
        done = touch("data/projects/{project}/metagenomes/{sample}/bins/maxbin/.done")
    params:
        bin_dir = "data/projects/{project}/metagenomes/{sample}/bins/maxbin/maxbin"
    benchmark: "benchmarks/maxbin/{project}__{sample}.txt"
    conda: "config/conda_yaml/maxbin.yaml"
    resources: cpus=16, mem_mb=20000, time_min=10080, mem_gb = 50
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

#### For testing on only one sample ####
#metaG_samples = "c39841b318d0487eda9e0134e5c06381"
#metaG_samples = "coassembly"
########################################


rule calc_contig_coverage:
    input: 
        expand("data/projects/{project}/metagenomes/{sample}/bins/contig_coverage.tsv", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/semibin", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/METABAT2/.done/", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/maxbin/.done", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/VAMB", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/metadecoder/.done", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm.txt", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/das_tool/.done", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/.drep_done", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/.done_GTDB", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/.done_gunc", sample = metaG_samples, project = "2022_geomicro_JGI_CSP")

rule run_checkm_new_per_sample:
    input: expand("data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm.txt", sample =metaG_samples, project = "2022_geomicro_JGI_CSP")


rule gunc_and_drep:
    input: 
        expand("data/projects/{project}/metagenomes/{sample}/bins/.done_gunc", sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/.drep_done", sample = metaG_samples, project = "2022_geomicro_JGI_CSP")


rule bin_coassembly:
    input: 
        expand("data/projects/{project}/metagenomes/{sample}/bins/contig_coverage.tsv", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/assembly/megahit/final.contigs.renamed.fa", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/semibin", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/METABAT2/.done/", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/maxbin/.done", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/VAMB", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/metadecoder/.done", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm.txt", sample = "coassembly", project = "2022_geomicro_JGI_CSP"),
        #expand("data/omics/metagenomes/{sample}/bins/das_tool/.done", sample = "coassembly"),
        #expand("data/omics/metagenomes/{sample}/bins/.drep_done", sample = "coassembly")



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
        contigs = rules.rename_megahit_contigs.output.contigs,
        bam_dir = "data/projects/{project}/metagenomes/{sample}/bins/bam"
    output:
        out_dir = directory("data/projects/{project}/metagenomes/{sample}/bins/semibin"),
        done = touch("data/projects/{project}/metagenomes/{sample}/bins/semibin/.done")
    params:
    conda: "config/conda_yaml/semibin.yaml"
    benchmark: "benchmarks/semibin/{project}__{sample}.txt"
    log: "logs/semibin/{project}__{sample}.log"
    resources: cpus=16, mem_mb=170000, time_min=2880, mem_gb = 50
    shell:
        """
        WORK_DIR=$PWD

        SemiBin \
            single_easy_bin \
            --threads {resources.cpus} \
            --reference-db {input.ref_db} \
            -i {input.contigs} \
            -b {input.bam_dir}/*.bam \
            -o {output.out_dir} | tee {log}

            #  --environment mouse_gut ## can only be used for single sample binning
        """

rule run_semibin:
    input: expand("data/projects/{project}/metagenomes/{sample}/bins/semibin/.done", sample = metaG_samples, project = "2022_geomicro_JGI_CSP")


rule VAMB:
    input:
        contigs = rules.rename_megahit_contigs.output.contigs,
        coverm_depth = "data/projects/{project}/metagenomes/{sample}/bins/metabat_style_contig_coverage.tsv"
    output:
        outdir = directory("data/projects/{project}/metagenomes/{sample}/bins/VAMB")
    #conda: "config/conda_yaml/VAMB.yaml"
    conda: "VAMB"
    benchmark: "benchmarks/VAMB/{project}__{sample}.txt"
    log: "logs/VAMB/{project}__{sample}.log"
    resources: cpus=1, mem_mb=40000, time_min=1440, partition = "gpu", gpu = 1
    shell:
        """
        vamb -o _ --outdir {output.outdir} --fasta {input.contigs} --jgi {input.coverm_depth} --minfasta 200000 --cuda
        """

rule format_coverage_for_metadecoder:
    input:
        script = "code/make_metadecoder_coverage.R",
        coverage = "data/projects/{project}/metagenomes/{sample}/bins/contig_coverage.tsv",
        contigs = rules.rename_megahit_contigs.output.contigs
    output: "data/projects/{project}/metagenomes/{sample}/bins/metadecoder/coverage.tsv"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb = 50000, time_min=360
    shell:
        """
        pwd && cd {current_dir} && pwd

        ./{input.script} --coverage={input.coverage} --contigs={input.contigs} --out={output}
        """
        
rule metadecoder:
    input:
        contigs = rules.rename_megahit_contigs.output.contigs,
        coverage = "data/projects/{project}/metagenomes/{sample}/bins/metadecoder/coverage.tsv"
    output:
        touch("data/projects/{project}/metagenomes/{sample}/bins/metadecoder/.done"),
        seed = "data/projects/{project}/metagenomes/{sample}/bins/metadecoder/seed.txt",
        bins_dir = directory("data/projects/{project}/metagenomes/{sample}/bins/metadecoder/bins")
    params:
        out_prefix = "data/projects/{project}/metagenomes/{sample}/bins/metadecoder/bins/{sample}"
    #conda: "config/conda_yaml/VAMB.yaml"
    conda: "metadecoder"
    shadow: "minimal"
    benchmark: "benchmarks/metadecoder/{project}__{sample}.txt"
    log: "logs/metadecoder/{project}__{sample}.log"
    resources: cpus=1, mem_mb=150000, time_min=10080, partition = "gpu", gpu = 1
    shell:
        """
        mkdir -p {output.bins_dir}
        
        metadecoder seed --threads {resources.cpus} -f {input.contigs} -o {output.seed}

        metadecoder cluster -f {input.contigs} -c {input.coverage} -s {output.seed}  -o {params.out_prefix} | tee {log}
        """


rule standardize_bins:
    input:
        "data/projects/{project}/metagenomes/{sample}/bins/CONCOCT/cut_contigs_10K.fa",
        "data/projects/{project}/metagenomes/{sample}/assembly/megahit_noNORM/final.contigs.renamed.fa",
        "data/projects/{project}/metagenomes/{sample}/bins/semibin",
        "data/projects/{project}/metagenomes/{sample}/bins/METABAT2/.done/",
        "data/projects/{project}/metagenomes/{sample}/bins/maxbin/.done", 
        "data/projects/{project}/metagenomes/{sample}/bins/VAMB",
        "data/projects/{project}/metagenomes/{sample}/bins/metadecoder/.done",
        script = "code/standardize_bins.R"
    output: 
        contig_bin_mapping = "data/projects/{project}/metagenomes/{sample}/bins/contig_bins.rds",
        bins_linked = "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/.bins_linked"
    params: 
        sample = "{sample}",
        sample_dir = "data/projects/{project}/metagenomes/{sample}"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb = 50000, time_min=360
    shell:
        """
        pwd && cd cd {current_dir} && pwd

        ./{input.script} --sample_dir={params.sample_dir}
        """

rule checkm_new_per_sample:
    input: "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/.bins_linked"
    output:
        results = "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm.txt"
    params:
        in_dir = "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins",
        out_dir = "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm"
    conda: "config/conda_yaml/checkm.yaml"
    resources: cpus=16, mem_mb=80000, time_min=2880
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
        """


rule make_das_and_drep_inputs:
    input:
        "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/.bins_linked",
        "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm.txt",
        contig_bin_mapping = "data/projects/{project}/metagenomes/{sample}/bins/contig_bins.rds",
        script = "code/make_das_and_drep_inputs.R"
    output: 
        metabat2_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/metabat2_contigs.tsv",
        maxbin_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/maxbin_contigs.tsv",
        concoct_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/concoct_contigs.tsv",
        metadecoder_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/metadecoder_contigs.tsv",
        semibin_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/semibin_contigs.tsv",
        VAMB_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/VAMB_contigs.tsv",
        drep_bin_info = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/genome_info.csv",
        drep_bins_linked = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/.bins_linked"
    params: 
        sample = "{sample}"
    singularity: "docker://eandersk/r_microbiome"
    resources: cpus=1, mem_mb = 50000, time_min=360
    shell:
        """
        pwd && cd {current_dir} && pwd

        ./{input.script} --sample={params.sample}
        """


rule checkm_new:
    #input: "data/omics/metagenomes/metagenome_bins/raw_combined_bins"
    output:
        #dir = temp(directory("data/omics/metagenomes/metagenome_bins/raw_combined_bins")),
        results = "data/projects/{project}/metagenomes/metagenome_bins/raw_combined_bins/checkm.txt"
    params:
        in_dir = "data/projects/{project}/metagenomes/metagenome_bins/raw_combined_bins",
        out_dir = "data/projects/{project}/metagenomes/metagenome_bins/raw_combined_bins/checkm"
    conda: "config/conda_yaml/checkm.yaml"
    resources: cpus=24, mem_mb=120000, time_min=2880
    shell:
        """
        checkm lineage_wf --tab_table -f {output.results} -x fa -t {resources.cpus} {params.in_dir} {params.out_dir}
        """


rule dastool_new:
    input:
        contigs = rules.rename_megahit_contigs.output.contigs,
        checkm_res = "data/projects/{project}/metagenomes/{sample}/bins/all_raw_bins/checkm.txt",
        metabat2_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/metabat2_contigs.tsv",
        maxbin_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/maxbin_contigs.tsv",
        concoct_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/concoct_contigs.tsv",
        metadecoder_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/metadecoder_contigs.tsv",
        semibin_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/semibin_contigs.tsv",
        VAMB_contigs = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/VAMB_contigs.tsv"
    params:
        bin_folder = rules.make_bins.params.out_dir,
        das_prefix = "data/projects/{project}/metagenomes/{sample}/bins/das_tool/output/{sample}"
    output: 
        #summary = "data/projects/{project}/metagenomes/{sample}/bins/DASTool/_DASTool_summary.txt",
        das_done = touch("data/projects/{project}/metagenomes/{sample}/bins/das_tool/.done"),
        #das_bins_folder = directory("data/omics/metagenomes/{sample}/bins/das_tool/output/_DASTool_bins")
    conda: "config/conda_yaml/das_tool.yaml"
    benchmark: "benchmarks/dastool/{project}__{sample}.txt"
    log: "logs/dastool/{project}__{sample}.log"
    resources: cpus=8, mem_mb=50000, time_min=2880, mem_gb = 50
    shell:
        """
        mkdir -p $(dirname {params.das_prefix})
        
        DAS_Tool \
            -i {input.metabat2_contigs},{input.maxbin_contigs},{input.concoct_contigs},{input.metadecoder_contigs},{input.semibin_contigs},{input.VAMB_contigs} \
            -c {input.contigs} \
            -o {params.das_prefix} \
            -l metabat2,maxbin,concoct,metadecoder,semibin,VAMB \
            --threads {resources.cpus} \
            --write_bins \
             | tee {log}
        """


rule drep_new:
    input: 
        "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/.bins_linked"
    output:
        touch("data/projects/{project}/metagenomes/{sample}/bins/.drep_done")
    params:
        bins_linked = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/.bins_linked",
        input_bins = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/*.fa",
        genome_info = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/genome_info.csv", # Will need to be moved to inputs 
        main_dir = directory("data/projects/{project}/metagenomes/{sample}/bins/drep/"),
        MAGs= directory("data/projects/{project}/metagenomes/metagenome_bins/derep/dereplicated_genomes")
    conda: "config/conda_yaml/drep.yaml"
    benchmark: "benchmarks/drep/{project}__{sample}.txt"
    resources: cpus=8, mem_mb=150000, time_min=2880
    shell:
        """
        rm -rf {params.main_dir} # Clear any old drep output
        
        dRep dereplicate \
            {params.main_dir} \
            -p {resources.cpus} \
            --contamination 50 \
            --completeness 30 \
            --length 10000 \
            --genomeInfo {params.genome_info} \
            -g {params.input_bins}
        """

rule run_drep_sep:
    input:
        expand("data/projects/{project}/metagenomes/{sample}/bins/.drep_done",sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/das_tool/.done",sample = metaG_samples, project = "2022_geomicro_JGI_CSP"),
        expand("data/projects/{project}/metagenomes/{sample}/bins/.done_GTDB",sample = metaG_samples, project = "2022_geomicro_JGI_CSP")
        #expand("data/projects/{project}/metagenomes/{sample}/bins/.done_gunc",sample = metaG_samples)
        

rule GTDB:
    input:
        "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/.bins_linked",
        #"/home/kiledal/geomicro_home/references/.done_gtdb_refs_downloaded"
    params:
        input_bin_dir = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep",
        refs = "/home/kiledal/geomicro_home/references/gtdbtk/release207_V2",
        out_dir = "data/projects/{project}/metagenomes/{sample}/bins/GTDB"
    output:
        done = touch("data/projects/{project}/metagenomes/{sample}/bins/.done_GTDB")
    conda: "config/conda_yaml/gtdbtk.yaml"
    benchmark: "benchmarks/GTDB/{project}__{sample}.txt"
    log: "logs/GTDB/{project}__{sample}.log"
    resources: cpus=24, mem_mb=120000, time_min=2880
    shell:
        """
        export GTDBTK_DATA_PATH={params.refs}

        gtdbtk classify_wf --extension fa --genome_dir {params.input_bin_dir} --out_dir {params.out_dir} --cpus {resources.cpus}
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
        "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep/.bins_linked",
        ref_file = "data/reference/gunc_gtdb/gunc_db_gtdb95.dmnd",
        ref_dir = "data/reference/gunc_gtdb"
    output: 
        done = touch("data/projects/{project}/metagenomes/{sample}/bins/.done_gunc")
    params: 
        bin_dir = "data/projects/{project}/metagenomes/{sample}/bins/bins_for_drep",
        out_dir = "data/projects/{project}/metagenomes/{sample}/bins/gunc"
    resources: cpus = 24, mem_mb = 120000, time_min = 2880
    conda: "config/conda_yaml/gunc.yaml"
    benchmark: "benchmarks/gunc/{project}__{sample}.txt"
    log: "logs/gunc/{project}__{sample}.log"
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
        genes = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa",
        profile = "data/reference/kegg/kofamscan/profiles",
        ko_list = "data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/omics/metagenomes/coassembly/bins/kofamscan/{bin}_kofam_results.txt"
    conda: "config/conda_yaml/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofamscan/{bin}.txt"
    log: "logs/kofamscan/{bin}.log"
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

rule run_kofam_scan_bins:
    input:
        expand("data/omics/metagenomes/coassembly/bins/kofamscan/{bin}_kofam_results.txt", bin = glob_wildcards("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{BIN,[^/]+}.fa").BIN)

rule traitar:
    input:
        #bin_genes = expand("data/omics/metagenomes/metagenome_bins/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/metagenome_bins/{BIN,[^/]+}.fa").BIN),
        pfam_db = "data/reference/traitar_pfamDB",
        #bin_dir = "data/omics/metagenomes/metagenome_bins",
        sample_file = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar_sample_list.tsv"
    params:
        bin_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes",
        gene_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal",
        out_dir = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar"
    output: 
        #directory("data/omics/metagenomes/metagenome_bins/traitar")
        touch("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/traitar/.done")
    benchmark: "logs/traitar_coassembly/benchmark.txt"
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
        bin = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{bin}.fa"
    output: 
        proteins = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa",
        genes = "data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.gff"
    conda: "config/conda_yaml/main.yaml"
    resources: cpus = 1, mem_mb = 10000
    shell:
        """
        prodigal -p meta -i {input.bin} -a {output.proteins} -d {output.genes} #1>{log} 2>&1
        """

rule run_prodigal_mags_DREP:
    input: expand("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/prodigal/{bin}.faa", bin = glob_wildcards("data/omics/metagenomes/coassembly/bins/drep_ALL_SAMPLES/dereplicated_genomes/{BIN,[^/]+}.fa").BIN)
