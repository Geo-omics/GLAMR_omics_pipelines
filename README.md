<h2>The Complete Metatranscriptome Pipeline</h2>
This is a direct read alignment pipeline. In the case of multiple conditions, it also includes optional differentially expressed gene (DEGs), differentally abundant function (DAFs), and differentially abundant organisms (DAOs) community analysis.

<h3>Purpose:</h3>
This pipeline was created so scientists could run the complete microbiome RNA-seq analysis in a single simplified process. I have a Universal Reference Database with genes from all sequenced organisms (Eukaryotes, Viruses, Archaea, Bacteria, Plasmids...) that have both functional (Kegg, COG, Pfam, GO, InterPro, MetaCyc, metal binding...etc) and phylogenetic annotations. After QC and alignment of your reads, various scripts are run to output:
<h5>1. Community Phylogenetic Tree(s) - which include quantitative and comparative (if specified) data
<h5>2. Community Functional Analysis - summary of reads, genes, top organism, and lowest common ancestor for each identified function
<h5>3. Gene Info Matrix - contains each matched gene from the universal database, with gene info, alignment score, various read counts (RPM, unique RPM, and RPKM), functional annotations and phylogeny. This info matrix can be used in R for comparative analysis.

<h3>Direct Alignment:</h3>
<h4>trimming -> read cleaning -> direct gene alignment -> differential analysis</h4>
<h4>OR</h4>
<h3>Metagenome Alignment:</h3>
You can first analyze your contig genes: https://github.com/TealFurnholm/Teals_Strain-Level_Metagenome_Pipeline/wiki/Contig-Analysis
<h4>trimming -> read cleaning -> align to annotated contigs -> differential analysis</h4>

<h3>Requirements:</h3>
Metatranscriptomics data is very large and requires substantial computing power. Hopefully you have server access and some familiarity working on a linux/unix system.
<br>* If you don't already have them on your system, install 
   <br>- trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
   <br>- bbtools: https://sourceforge.net/projects/bbmap/
   <br>- diamond: https://github.com/bbuchfink/diamond
   <br>- perl: https://www.perl.org/get.html
   <br>- R (if doing differential expression analysis): https://www.r-project.org/
   
<h2>How to Use:</h2>
https://github.com/TealFurnholm/Metatranscriptome/wiki



