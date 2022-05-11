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

# Get import sample names
#metaG_samples = glob_wildcards("import/metagenomes/{sample}/").sample
#metaT_samples = glob_wildcards("import/metatranscriptomes/{sample}/").sample
#metabolome_samples = glob_wildcards("import/metabolomes/{sample}/").sample
#amplicon_samples = glob_wildcards("import/amplicons/{sample}/").sample

# Get sample names
metaG_samples = glob_wildcards("data/omics/metagenomes/{sample}/reads").sample
metaT_samples = glob_wildcards("data/omics/metatranscriptomes/{sample}/").sample
metabolome_samples = glob_wildcards("data/omics/metabolomes/{sample}/").sample
amplicon_samples = glob_wildcards("data/omics/amplicons/{sample}/").sample

rule assemble:
    input: 
        expand("data/omics/metagenomes/{sample}/assembly/metaspades/contigs.fasta",sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa",sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/reads/dedup_reads.fastq",sample = metaG_samples)


rule prodigal:
    input: expand("data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa", sample = metaG_samples)

rule test:
    input: expand("data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",sample = metaG_samples)
    output: "test.out"

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


rule deduplicate:
    input: 
        fwd_reads = "data/omics/metagenomes/{sample}/reads/raw_fwd_reads.fastq.gz",
        rev_reads = "data/omics/metagenomes/{sample}/reads/raw_rev_reads.fastq.gz"
    output: 
        dedup_reads = "data/omics/metagenomes/{sample}/reads/dedup_reads.fastq"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/dedup/{sample}_dedup.log"
    resources: cpus=8, mem_mb = lambda wildcards, attempt: attempt * 100000,
        time_min=2880
    shell:
        """
        dedupe.sh \
            in1={input.fwd_reads} \
            in2={input.rev_reads} \
            out={output.dedup_reads} \
            t={resources.cpus} \
            1>{log} 2>&1
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
            1>{log} 2>&1
        """

rule remove_contaminants:
    input:
        dedup_reads = rules.deduplicate.output.dedup_reads,
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
    resources: cpus = 8, mem_mb = lambda wildcards, attempt: attempt * 60000, time_min = 2880
    shell:
        """        
        bbduk.sh -Xmx{resources.mem_mb}m \
            in={input.dedup_reads} \
            out1={output.trimmed_fwd} \
            out2={output.trimmed_rev} \
            t={resources.cpus} \
            minlen=50 \
            qtrim=rl \
            trimq=15 \
            ref={input.adapters} \
            path={params.bbmap_index_path} \
            ktrim=r k=23 mink=11 hdist=1 \
            1>{log} 2>&1
        
        echo "\n\n***doing spike-in removal***\n\n" >> {log}
        
        bbduk.sh -Xmx{resources.mem_mb}m \
            in1={output.trimmed_fwd} \
            in2={output.trimmed_rev} \
            outu1={output.phix_rm_fwd} \
            outu2={output.phix_rm_rev} \
            t={resources.cpus} k=31 hdist=1 \
            ref={input.spike_ins} \
            path={params.bbmap_index_path} \
            1>>{log} 2>&1
        
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
            1>>{log} 2>&1
        
        echo "\n\n***Running RemovePolyPairs.pl***\n\n" >> {log}
        perl code/RemovePolyPairs.pl {output.decon_fwd} {output.decon_rev} 50 {output.cleaned_fwd} {output.cleaned_rev} 1>>{log} 2>&1
        """

rule assemble_metaspades:
    input:
        fwd_reads = rules.remove_contaminants.output.cleaned_fwd,
        rev_reads = rules.remove_contaminants.output.cleaned_rev
    output:
        assembly_dir = directory("data/omics/metagenomes/{sample}/assembly/metaspades"),
        contigs = "data/omics/metagenomes/{sample}/assembly/metaspades/contigs.fasta"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/metaspades/{sample}.log"
    benchmark: "benchmarks/metaspades/{sample}.txt"
    resources: cpus = 16, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        metaspades.py -t {resources.cpus} -1 {input.fwd_reads} -2 {input.rev_reads} -o {output.assembly_dir} > {log}
        """

rule assemble_megahit:
    input:
        fwd_reads = rules.remove_contaminants.output.cleaned_fwd,
        rev_reads = rules.remove_contaminants.output.cleaned_rev
    output:
        assembly_dir = directory("data/omics/metagenomes/{sample}/assembly/megahit"),
        contigs = "data/omics/metagenomes/{sample}/assembly/megahit/final.contigs.fa"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/assembly/megahit/{sample}.log"
    benchmark: "benchmarks/megahit/{sample}.txt"
    resources: cpus = 16, time_min=7200, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        rm -r {output.assembly_dir} # for re-running, megahit doesn't overwrite automatically
        megahit -t {resources.cpus} --presets meta-sensitive -m 0.5 -1 {input.fwd_reads} -2 {input.rev_reads} -o {output.assembly_dir} > {log}
        """

rule merge_assemblies:
    input:
        metaspades_contigs = rules.assemble_metaspades.output.contigs,
        megahit_contigs = rules.assemble_megahit.output.contigs,
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
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.concat_contigs} out={output.dedup1} tuc mid=99 minscaf=200 rnc=f ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f sort=length absorbcontainment=t mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 1>{log} 2>&1
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup1} out={output.dedup2} tuc mid=99 minscaf=200 rnc=t ngn=t fo c pc=t fmj=t rc=t cc=t fcc=t mst=f sort=length absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 1>>{log} 2>&1
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup2} out={output.dedup3} tuc mid=99 minscaf=200 rnc=t ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f ordered=t absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t t={resources.cpus} 1>>{log} 2>&1
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup3} out={output.dedup4} tuc mid=99 minscaf=200 rnc=f ngn=f fo c pc=t fmj=t rc=t cc=t fcc=t mst=f ordered=t absorbcontainment=f mo=200 numaffixmaps=3 overwrite=t dot={output.dedup4_dot} t={resources.cpus} 1>>{log} 2>&1
        perl {input.merge_contigs_script} data/omics/metagenomes/{wildcards.sample}/assembly/{wildcards.sample} 99 1>>{log} 2>&1
        dedupe.sh -da -Xmx{resources.mem_mb}m -eoom in={output.dedup5} out={output.dedup6} t={resources.cpus} tuc mid=99 minscaf=200 overwrite=f 1>>{log} 2>&1
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
        bbmap.sh -da -Xmx${{use_mem}}m -eoom ref={input.assembly} path={params.assembly_dir} t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 idfilter=0.999 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} out={output.read_mapping} 1>{log} 2>&1
        samtools view -bShu {output.read_mapping} | samtools sort -m ${{mem_per_thread}}M -@ {resources.cpus} -o {output.read_mapping_sorted} 1>>{log} 2>&1
        samtools index {output.read_mapping_sorted} 1>>{log} 2>&1
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
            -o {output.corrected_assembly} 1>>{log} 2>&1
        """

rule predict_genes_and_calc_abundance:
    input:
        assembly = rules.MEC.output.corrected_assembly,
        fwd_reads = rules.remove_contaminants.output.cleaned_fwd,
        rev_reads = rules.remove_contaminants.output.cleaned_rev
    output:
        proteins = "data/omics/metagenomes/{sample}/proteins/{sample}_PROTEINS.faa",
        genes = "data/omics/metagenomes/{sample}/genes/{sample}_GENES.fna",
        reads_vs_genes_rpkm = "data/omics/metagenomes/{sample}/genes/{sample}_READSvsGENES.rpkm",
        reads_vs_contigs_rpkm = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.rpkm",
        reads_vs_assembly_sam_gz = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam.gz"
    params:
        reads_vs_assembly_sam = "data/omics/metagenomes/{sample}/assembly/{sample}_READSvsCONTIGS.sam"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/predict_genes_and_calc_abundance/{sample}_predict_genes_and_calc_abundance.log"
    benchmark: "benchmarks/predict_genes_and_calc_abundance/{sample}_predict_genes_and_calc_abundance.txt"
    resources: cpus = 4, mem_mb = lambda wildcards, attempt: attempt * 64000, time_min = 2880
    shell:
        """
        prodigal -p meta -i {input.assembly} -a {output.proteins} -d {output.genes} 1>{log} 2>&1
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {output.genes}) ref={output.genes} rpkm={output.reads_vs_genes_rpkm} 1>>{log} 2>&1
        bbmap.sh t={resources.cpus} ambig=random cigar=f maxindel=100 pairlen=600 minid=0.999 idtag=t printunmappedcount=t overwrite=t in1={input.fwd_reads} in2={input.rev_reads} path=$(dirname {input.assembly}) ref={input.assembly} rpkm={output.reads_vs_contigs_rpkm} 32bit=t outm={params.reads_vs_assembly_sam} 1>>{log} 2>&1
        gzip {params.reads_vs_assembly_sam}
        """

rule download_uniref:
    output: 
        uniref100="data/reference/uniref/uniref100.fasta.gz"
    conda: "config/conda_yaml/main.yaml"
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
        diamond makedb --threads {resources.cpus} --in {input.uniref100} -d data/reference/uniref/uniref100
        """

rule align_to_uniref:
    input:
        diamond_db = rules.make_diamond_uniref_db.output.uniref100_diamond,
        genes = rules.predict_genes_and_calc_abundance.output.genes
    output:
        gene_uniref_alignment = "data/omics/metagenomes/{sample}/{sample}_GENES.m8"
    conda: "config/conda_yaml/main.yaml"
    log: "logs/align_to_uniref/{sample}_align_to_uniref.log"
    benchmark: "benchmarks/align_to_uniref/{sample}_align_to_uniref.txt"
    resources: cpus = 8, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 100000
    shell:
        """
        diamond blastx -d {input.diamond_db} -q {input.genes} -o {output.gene_uniref_alignment}
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



rule annotate_contigs:
    input:
        genes = rules.predict_genes_and_calc_abundance.output.genes,
        gene_alignment = rules.align_to_uniref.output.gene_uniref_alignment,
        gene_rpkm = rules.predict_genes_and_calc_abundance.output.reads_vs_genes_rpkm,
        contig_rpkm = rules.predict_genes_and_calc_abundance.output.reads_vs_contigs_rpkm,
        contigs = rules.MEC.output.corrected_assembly,
        script = "code/Strain-Level_Metagenome_Analysis/AnnotateContigs.pl",
        UMRAD = "data/reference/UMRAD"
    output:
        annotation_dir = directory("data/omics/metagenomes/{sample}/annotation")
    conda: "config/conda_yaml/main.yaml"
    resources: cpus=1, time_min = 2880, mem_mb = lambda wildcards, attempt: attempt * 100000
    log: "logs/annotate_contigs/{sample}.log"
    benchmark: "benchmarks/annotate_contigs/{sample}.log"
    shell:
        """
        mkdir -p {output.annotation_dir}
        ln {input.genes} {input.gene_alignment} {input.gene_rpkm} {input.contig_rpkm} {input.contigs} {output.annotation_dir}/
        mv  {output.annotation_dir}/$(basename {input.contigs}) {output.annotation_dir}/{wildcards.sample}_MCDD.fa
        proj_dir=$(pwd)

        cp {input.script} {output.annotation_dir}/
        ln {input.UMRAD}/* {output.annotation_dir}/

        cd {output.annotation_dir}
        perl AnnotateContigs.pl -s {wildcards.sample} 
        
        #-d ${{proj_dir}}/{input.UMRAD}
        #> ${{proj_dir}}/{log}
        """

rule sample_annotation:
    input: 
        expand("data/omics/metagenomes/{sample}/annotation", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/bins/metabat2_bins", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/bins/gtdbtk", sample = metaG_samples),
        expand("data/omics/metagenomes/{sample}/bins/checkm.txt", sample = metaG_samples)


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
        contigs = rules.MEC.output.corrected_assembly,
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

rule coverm:
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

rule gather_bin:
    input: expand("data/omics/metagenomes/{sample}/bins/ran_dastool.touch",sample=metaG_samples)
    output: directory("data/omics/metagenome_bins")
    shell:
        """
        mkdir -p {output}
        
        ln data/omics/metagenomes/*/bins/DASTool/_DASTool_bins/*.fa {output}
        """

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