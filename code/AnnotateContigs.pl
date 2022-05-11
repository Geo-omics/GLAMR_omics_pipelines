#use warnings;
use Statistics::Basic qw(:all nofill);
use Sort::Naturally 'nsort';


########################################################################
##########                      BEGIN                        ###########
########################################################################
$argv=join(" ", @ARGV);
if($argv =~ /\-s\s+(\S+)/){$samp=$1;} 
else{print "missing sample prefix ex: -s Toxin_Day2\nex:  Toxin_Day2.fastq\n"; die;}
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$month=$1; $year=$2;
$version=$1."_".$2;
$log=$samp."_AnalyzeContigs_".$version.".log";
open(LOG, ">", $log)||die;
print LOG "sample $samp version $version\n";



########################################################################
########## 	GIVE INSTRUCTIONS IF NEEDED	#######################
$argv = join(" ", @ARGV);
if($argv !~ /\w/ || $argv =~ /(^\-h|^\-help|\s\-h)/){
		print "\n\nPlease follow the instructions on GitHub:\n";
		print "\t1. construct UMRAD starting here -> https://github.com/TealFurnholm/Universal-Taxonomy-Database \n";
		print "\t2. run the metagenome through the pipeline, here -> https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis/wiki \n\n";
		print "Required sample data include:\n";
		print "(This script assumes you are running it where these sample data are.)\n";
		print "\t1. The gene or protein sequences:                      [sample]_GENES.fna\n";
		print "\t2. The alignment file of the sequences to UMRAD:       [sample]_GENES.m8 \n";
		print "\t3. The rpkm file for reads aligned to the genes:       [sample]_READSvsGENES.rpkm\n";
		print "\t4. The rpkm file for reads aligned to the contigs:     [sample]_READSvsCONTIGS.rpkm\n";
		print "\t5. The sample contigs sequences:                       [sample]_MCDD.fa\n\n";
		print "HOW TO USE ...\n";
		print "You must specify the following: \n";
		print "\t-d /path/to/UMRAD_reference/files/ \n";
		print "\t-s the_sample_prefix\n";
		print "\nEx: perl AnalyzeContigs.pl -s Control_Day_1 -d /reference_files/UMRAD/ \n\n";
		print "USER OPTIONS:\n";
		print "\t-c \t Minimum \% query coverage (default -c=70)\n";
		print "\t-m \t Minimum \% identity alignment (default -m=50)\n";
		print "\t-x \t Minimum number of genes to keep cellular organisms (default -x=3). Note: does not include viruses or plasmids\n";
		print "\t-k \t Top score margin, between 0 and 1. (default 0.9)\n";
		print "\t-n \t Minimum number of model genomes for functional models. Use 3 or greater (for averaging). \n";
		print "\t   \t More may result in lower rank/less specific model. (default 5)\n";
		print "\t\t Scores are \%identity x \%coverage (max possible score is 100x100=10000).\n";
		print "\t\t The -k margin keeps hits where: hit score >= top_score x [margin].\n";
		print "\t\t ex: 10000 x 0.75 = keep all hits with score >= 7500\n";
		print "\t\t For top score hits only = 1. To ignore hit scores = 0. (default = 0.75)\n\n";
		print "These options allow that the reference database may not have your sample's strains, just various relatives with mixed gene orthology\n";
		print "Lower stringency helps stop false negatives. Higher helps increase speed and reduce noise/false positives.\n";
		print "Many Excess/spurious hits and species are screened out later in the script.\n\n\n";
		die;
}
########################################################################
########################################################################



print "Checking user files \n\n";

########################################################################
##########		CHECK USER FILES		################
########################################################################
$dir = './';
opendir(DIR, $dir) or die "Could not open $dir\n";
@FILES = grep(/\./i, readdir DIR);
foreach my $file (@FILES){
	if($file =~ /$samp.*GENES\.m8/){ 		$ingvu=$file;}
	if($file =~ /$samp.*GENES.fna/){ 		$ingen=$file;}
	if($file =~ /$samp.*READSvsGENES.rpkm/){ 	$ingrc=$file;}
	if($file =~ /$samp.*READSvsCONTIGS.rpkm/){ 	$incov=$file;}
	if($file =~ /$samp.*MCDD.fa/){ 			$incon=$file;}
}
open(INGVU,$ingvu)||die "unable to open samp_GENES.m8 $ingvu:$!\n";
open(INGEN,$ingen)||die "unable to open samp_GENES.fna $ingen:$!\n";
open(INGRC,$ingrc)||die "unable to open samp_READSvsGENES.rpkm $ingrc:$!\n";
open(INCC, $incov)||die "unable to open samp_READSvsCONTIGS.rpkm $incov:$!\n";
open(INCON,$incon)||die "unable to open samp_MCDD.fa :$!\n";
print LOG "dir $dir\ningvu $ingvu\ningen $ingen\ningrc $ingrc\nincov $incov\nincon $incon\n";
########################################################################
########################################################################



print "Checking reference files \n \n";

########################################################################
##########		CHECK REFERENCE FILES		################
########################################################################
if($argv =~ /\s\-d\s+(\S+)/){ $refdir=$1;} else{ $refdir='./'; }
if($refdir !~ /\/$/){$refdir.='/';}
$refdir = '/geomicro/data22/teals_pipeline/BOSS/';
opendir(REFDIR, $refdir) or die "Could not open $dir\n";
@FILES = grep(/\./i, readdir REFDIR);
foreach my $file (@FILES){
	if($file =~ /TAXONOMY\_DB.*$year\.txt/){	$intax	=$refdir.$file;}
	if($file =~ /UNIREF100\_INFO.*$year\.txt/){	$ininfo	=$refdir.$file;}
	if($file =~ /Function\_Names.*\.txt/){		$infn	=$refdir.$file;}
}
open(INTAX, $intax)||die "unable to open intax $intax:$!\n";
if($ininfo=~/\.gz$/){open(INFO, "gunzip -c $ininfo |")||die "unable to open uniref info file:$!\n";}
else{open(INFO,$ininfo)||die "unable to open uniref info file:$!\n";}
open(INFN, $infn)||die "unable to open function names file:$!\n";
print LOG "refdir $refdir\nintax $intax\nininfo $ininfo\ninfn $infn\n";
########################################################################
########################################################################


print "Getting user settings \n \n";


########################################################################
##########		GET USER SETTINGS		################
########################################################################
if($argv =~/\s\-k\s+([\d\.]+)/){$top=$1;}else{$top = 0.9;}
if($argv =~/\s\-x\s+(\d+)/){$mingen=$1;} else{$mingen = 3;}
if($argv =~/\s\-m\s+(\d+)/){$minid=$1;}  else{$minid = 60;}
if($argv =~/\s\-c\s+(\d+)/){$mincov=$1}  else{$mincov = 60;}
if($argv =~/\s\-n\s+(\d+)/){$min_mod=$1} else{$min_mod = 5;}
$minsco=$minid*$mincov;
print LOG "top $top\nmingen $mingen\nminid $minid\nmincov $mincov\nmin_mod $min_mod\nminsco $minsco\n";
########################################################################
########################################################################





########################################################################
##########              GET FUNCTION NAMES              ################
########################################################################
	#INPUT FUNCTION NAMES
        open(INFN, $infn)||die;
        while(<INFN>){
                if($_ !~ /\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                (my $id, my $name)=split("\t",$_,-1);
                $id=~s/\s//g;
                $FUNC_NM{$id}=$name;
        }
        close(INFN);


	#GET CHEBI CPD NAMES
        print "RETRIEVE CPD NAMES\n";
        $start=147;
        $count=0;
        #while($count <=10){
        #	$count++;
        #	$file='https://ftp.ebi.ac.uk/pub/databases/chebi/archive/rel'.$start.'/Flat_file_tab_delimited/names.tsv.gz';
        #	qx{wget -O nx.tsv.gz $file};
        #	qx{cat nx.tsv.gz >> names.tsv.gz};
        #	qx{rm nx.tsv.gz};
        #	$start+=12;
        #	print "on count $count start $start\n";
        #}

	#INPUT NAMES
	print "INPUTTING NAMES\n";
	#qx{gunzip -f names.tsv.gz};
	open(INCPDNM, "names.tsv")||last;
	while(<INCPDNM>){
	        if($_!~/\w/){next;}
	        $_=uc($_);
	        $_=~s/[\r\n]+//;
	        @stuff=split("\t", $_);
	        $cpd="CHEBI:".$stuff[1];
	        $name=CleanNames($stuff[4]);
	        @PARTS = split("_", $name);
	        $mp=0; $mn='';
	        foreach $p (@PARTS){ $p=~/([A-Z]+)/; if(length($1)>$mp){ $mp=length($1); $mn=$1; }}
	        $CPD_NODD{$cpd}{$name}=@PARTS;
	        $CPD_NLEN{$cpd}{$name}=length($mn);
	}
	#qx{rm names.tsv*};

	#GET/CLEAN CPD NAMES -> MAKE HUMAN READABLE
	print "CLEAN CPD NAMES\n";
	foreach my $cpd (sort(keys %CPD_NODD)){
	        #FIRST CHECK FOR NAMES >= 7 CHARACTER STRINGS
	        #SORT BY INCREASING NON-WORDS
	        $mxo=10000; $mxl=0;
	        foreach my $name (sort{$CPD_NODD{$cpd}{$a}<=>$CPD_NODD{$cpd}{$b}} keys %{$CPD_NODD{$cpd}}){
	                $odd = $CPD_NODD{$cpd}{$name};
	                $len = $CPD_NLEN{$cpd}{$name};
	                if($len < 7){next;}
	                #NOW LOOK FOR NAME WITH LONGEST DESCRIPTOR
	                if($odd < $mxo){ $mxo=$odd; }
	                if($odd <= $mxo+1 ){
	                        if($len > $mxl){$mxl=$len; $CPD_NAMES{$cpd}=$name;}
	                }
	        }

	        #GET LONGEST NAME REGARDLESS OF ODD IF NO STRING LONGER THAN 7
	        if(!exists($CPD_NAMES{$cpd})){
	                foreach my $name (sort{$CPD_NLEN{$cpd}{$b}<=>$CPD_NLEN{$cpd}{$a}} keys %{$CPD_NLEN{$cpd}}){
	                        $CPD_NAMES{$cpd}=$name; last;
	                }
	        }

	        #FIX STUPID LONG NON-NAMES
	        if(length($CPD_NAMES{$cpd}) > 49){
	                $CPD_NAMES{$cpd} =~ s/\d+[A-Z]\_//g;
	                $CPD_NAMES{$cpd} =~ s/^[^A-Z]+|[^A-Z]+$//g;
	                @GN=();
	                @NM=split("_",$CPD_NAMES{$cpd});
	                while($NM[0]=~/./){
	                        $n=shift(@NM);
	                        if($n=~/^\d+$/){next;}
	                           if(grep /$n/, @NM){}
	                        elsif(grep /$n/, @GN){}
	                        else{push(@GN,$n);}
	                }
	                $CPD_NAMES{$cpd}=join("_",@GN);

	                while($CPD_NAMES{$cpd} =~ /\_.{1,2}\_/){$CPD_NAMES{$cpd} =~ s/\_.{1,2}\_/\_/; if(length($CPD_NAMES{$cpd})<7){last;}}
	                $CPD_NAMES{$cpd}="MOD-".$CPD_NAMES{$cpd};
	        }
	        $CPD_NAMES{$cpd}=~s/WURCS/COMPLEX_CARBOHYDRATE/;
		#print "cpd $cpd name $CPD_NAMES{$cpd}\n";
	}
	undef(%CPD_NODD);
	undef(%CPD_NLEN);
########################################################################
########################################################################
#$CPD_NAMES{$cpd}=$name;
#$FUNC_NM{$id}=$name;




########################################################################
##########	INPUT CONTIG AND GENE BASICS		################
########################################################################

## INPUT CONTIG SEQS
$/=">";
$time=localtime;
print "INPUT CONTIG SEQUENCES time $time $incon\n";
open(INCON, $incon)||die;
while(<INCON>){
		if($_ !~ /\w/){next;}
		$_=uc($_);
		@stuff=split("\n",$_);
		$contig=shift(@stuff);
		$seq=join("",@stuff);
		$seq=~s/[^A-Z]//g;
		$gc = $seq =~ tr/GC/gc/;
		$seq =uc($seq);
		$len = length($seq);
		$pgc = $gc*100/$len; 
		$pgc =~ s/(\.\d\d).*/$1/;
		$CON_PGC{$contig}=$pgc;
		$CON_LEN{$contig}=$len;
		$CON_SEQ{$contig}=$seq;
}
$/="\n";
close(INCON);

## INPUT CONTIG RPKM
$time=localtime;
print "INPUT CONTIG RPKMS time $time $incov\n";
open(INCC, $incov)||die;
while(<INCC>){
		if($_ !~ /\d/){next;}
		$_=uc($_);
		$_ =~ s/[\r\n]//;
		@stuff = split("\t", $_, -1);
		$contig=$stuff[0];
		if(!exists($CON_PGC{$contig})){next;}
		$depth=$stuff[3]; 
		$depth =~ s/(\.0*\d\d).*/$1/;
		$rpkm=$stuff[5]; 
		$rpkm =~ s/(\.0*\d\d).*/$1/;
		$CON_DEPTH{$contig}=$depth;
		$CON_RPKM{$contig}=$rpkm;
}
close(INCC);

## INPUT CONTIG GENES
$/=">";
$time=localtime;
print "INPUT CONTIG GENES time $time $ingen\n";
open(INGEN,$ingen)||die;
while(<INGEN>){
		if($_ !~ /\w/){next;}
		$_=uc($_);
		@stuff=split("\n",$_);
		$header = shift(@stuff);
		$header =~ /^((\S+)\_\d+)[\s\#]+(\d+)[\s\#]+(\d+)[\s\#]+.*PARTIAL\=(\d+)/i;
		$gene   =$1; 
		$contig =$2; 
		$start  =$3; 
		$len	=$4-$3; 
		$part   =$5;  

		if($gene!~/\w/ || $contig!~/\w/ || $start!~/\w/ || $len!~/\w/ || $part !~/\w/){die;}
		$seq=join("",@stuff);
		$seq=~s/[^A-Z]//g;
		my $count = $seq =~ tr/GC/gc/;
		$len=length($seq);
		$prgc=$count*100/$len; 
		$prgc =~ s/(\.\d\d).*/$1/;
		$CON_GENES{$contig}{$gene}=1;  
		$GENE_START{$gene}=$start;
		$GENE_PART{$gene}=$part;
		$GENE_PGC{$gene}=$prgc;
		$GENE_LEN{$gene}=$len;
}
$/="\n";
close(INGEN);

## INPUT GENE COVERAGE
$time=localtime;
print "INPUT GENE RPKMS time $time $ingrc\n";
open(INGRC,$ingrc)||die;
while(<INGRC>){
		if($_ !~ /\w/){next;}
		if($_ =~ /^\#/){next;}
		$_=uc($_);
		$_=~s/[\r\n]+//;
		@stuff=split("\t",$_,-1);
		$depth=$stuff[3]; 
		$depth =~ s/(\.0*\d\d).*/$1/;
		$rpkm=$stuff[5]; 
		$rpkm =~ s/(\.0*\d\d).*/$1/;
		$stuff[0] =~ s/\s.*//;
		$GENE_DEPTH{$stuff[0]}=$depth;
		$GENE_RPKM{$stuff[0]}=$rpkm; 
}
close(INGRC);
#####################################################
#####################################################
$genestart= keys %GENE_START;
$genepart=  keys %GENE_PART;
$genepgc=   keys %GENE_PGC;
$genelen=   keys %GENE_LEN;
$genedepth= keys %GENE_DEPTH;
$generpkm=  keys %GENE_RPKM;
$congen=    keys %CON_GENES;
$condepth=  keys %CON_DEPTH;
$conrpkm=   keys %CON_RPKM;
$conpgc=    keys %CON_PGC;
$conlen=    keys %CON_LEN;
print "genestart $genestart genepart $genepart genepgc $genepgc genelen $genelen genedepth $genedepth generpkm $generpkm\n";
print "congen $congen condepth $condepth conrpkm $conrpkm conpgc $conpgc conlen $conlen\n";
#####################################################
#####################################################





 

########################################################################
##########	INPUT ALIGNMENT AND DATABASE FILES			################
########################################################################

## INPUT TAXONOMY ##
$time=localtime;
print "INPUT TAXONOMY time $time $intax\n";
open(INTAX, $intax)||die;
while(<INTAX>){
	if($_ !~ /\w/){next;}
	$_ =~ s/[\r\n]//;
        $_=uc($_);
        @stuff = split("\t", $_, -1);
        $tid =  shift(@stuff);
        $lin = join(";",@stuff);
        if(!exists($LIN_TID{$lin})){$LIN_TID{$lin}=$tid;}
        $PHY{$tid}=$lin;
}


## INPUT DIAMOND HITS ##
$on=0;
$time=localtime;
print "input $ingvu time $time\n";
while(<INGVU>){
	if($_ !~ /\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_,-1);
	$gene = $stuff[0]; 
	$hit  = $stuff[2];
	$pid  = $stuff[9];
	if($pid < $minid){next;}
	$cov  = $stuff[11];
	if($cov < $mincov){next;}
	$sco  = $pid*$cov;
	$AG{$gene}=1;						#LIN SCREENING - track all sample genes
	if($sco > $TOP_HITS{$hit}){$TOP_HITS{$hit}=$sco;}	#LIN SCREENING - get best score for hit
	if($sco > $GENE_SCO{$gene}){
		$GENE_SCO{$gene}=$sco;  			#set new top score hit for the gene
		foreach my $h (keys %{$GENE_HITS{$gene}}){ 	#remove hits not in the $top X percentile
			if($GENE_HITS{$gene}{$h} < $GENE_SCO{$gene}*$top){ 
				delete($GENE_HITS{$gene}{$h});
				delete($HIT_GENES{$h}{$gene});
	}	}	}
	if($sco >= $GENE_SCO{$gene}*$top){ 
		$GENE_HITS{$gene}{$hit}=$sco; 
		$HIT_GENES{$hit}{$gene}=$sco;
}	}
$start_hits = keys %TOP_HITS;
undef(%GENE_SCO);


## INPUT HIT INFO ##
$on=0;
$time=localtime;
print "input $ininfo time $time\n";
open(INFO,"TEST_ANCON.txt")||die; ### !!!!
while(<INFO>){
	if($_ !~ /\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_,-1);
	$hit = $stuff[0];
	$name = $stuff[1];

	#taxa
	@tids = split(";",$stuff[5]); 
	%LINS=();
	foreach my $tid (@tids){
		if($PHY{$tid}!~/^(BACT|ARCH|MONA|EUKA|MICRO)/){next;}
		$lin=$PHY{$tid};
		$TOT_LIN_HITS{$lin}{$hit}=1; 			#LIN SCREENING - for genome % completeness
		if(!exists($TOP_HITS{$hit})){next;} 		#NOT A HIT, SO DONE
		$LINS{$lin}=1;					#FOR FUNCTIONS NEXT
		$LIN_HITS{$lin}{$hit}=$TOP_HITS{$hit};		#LIN SCREENING - for genome % completeness and sum sco
		$HIT_LINS{$hit}{$lin}=1;			#FOR GET UNIQS AND LATER LOOP CON-GENE-HIT-LINS
		foreach my $gene (keys %{$HIT_GENES{$hit}}){
			$LIN_GENES{$lin}{$gene}=1;		#LIN SCREENING - to account for good lin %AG genes
			$gene=~/(.*?)\_\d+$/;
			$contig=$1;
			$LIN_CONS{$lin}{$contig}=1;		#LATER OUTPUT PHYLOBIN CONTIGS
			$TOSS{$lin}=1;				#LIN SCREENING - for lin tracking
	}	}
	if(!exists($TOP_HITS{$hit})){next;}			#NOT A HIT, SO DONE
	$HIT_NAMES{$hit}=$name;

	#functions
	@IDS=();
	for my $i (7..21){
		if($stuff[$i]!~/\w/){next;}
		if($i==7 || $i==8 || $i==9){ $IDS[0]=$stuff[$i];}	#AA DOMAINS
		else{@IDS=split(";", $stuff[$i]);}
		foreach my $id (@IDS){
			if($i==10){$id="METAL-BINDING:".$id;}
			if($i==12){$id="LOCATION:".$id;}
			$HIT_FUNCS{$hit}{$i}{$id}=1; 			#FOR MAKING GENE FUNC OUTPUTS LATER
	}	}

	#compounds
	for my $i (22..24){
		if($stuff[$i]!~/\w/){next;}
		@IDS=split(";", $stuff[$i]);
		foreach my $id (@IDS){if($id=~/\w/){$HIT_CPDS{$hit}{$i}{$id}=1;}}	#FOR MAKING GENE CPDS LATER
	}
	if($on%100000==0){print "on $on hit $hit tophitsco $TOP_HITS{$hit}\n";} $on++;
	delete($TOP_HITS{$hit}); if(keys %TOP_HITS < 1){last;}
}
undef(%HIT_GENES);
undef(%TOP_HITS);
undef(%LINS);
undef(%PHY);

## GET UNIQUE GENES 
$time=localtime;
print "MAKE UNIQUE GENE LIST $time\n";
foreach my $gene (keys %GENE_HITS){
	$uniq=0;
	foreach my $hit (keys %{$GENE_HITS{$gene}}){
		#only 1 hit for the gene, only 1 lin for hit 
		if(keys %{$GENE_HITS{$gene}}==1 && keys %{$HIT_LINS{$hit}}==1){$uniq=1;}
		else{$uniq=0;}
	}
	$GENE_UNIQ{$gene}=$uniq;
}
########################################################################
########################################################################
$end_hits   = keys %HIT_GENES;
$sampl_hits = keys %HIT_LINS;
$func_hits  = keys %HIT_FUNCS;
$cpd_hits   = keys %HIT_CPDS;
$start_genes= keys %GENE_HITS;
$uniq_genes = keys %GENE_UNIQ;
$total_lins = keys %TOT_LIN_HITS;
$sampl_lins = keys %LIN_HITS;
$lingn	    = keys %LIN_GENES;
print "start_hits $start_hits end_hits $end_hits sampl_hits $sampl_hits func_hits $func_hits cpd_hits $cpd_hits \n";
print "start_genes $start_genes uniq_genes $uniq_genes\n";
print "total_lins $total_lins sampl_lins $sampl_lins lingn $lingn\n";
########################################################################
########################################################################






########################################################################
##########              SCREEN LINEAGES                 ################
########################################################################

$on=0;
$time=localtime;
print "GET METRIC 1: % GENOME COMPLETE time $time\n";
foreach my $lin (keys %LIN_HITS){ #every gene-hit lineage

        # GET PERCENT GENOME COMPLETE
        $tlg = keys %{$TOT_LIN_HITS{$lin}};	#all possible/sequenced lin genes in ref database
        $lg  = keys %{$LIN_HITS{$lin}};		#sequenced lin genes found in sample
        $per_gnome_comp = $lg*100/$tlg;
        $per_gnome_comp =~ s/(\.\d\d).*/$1/;
        $LIN_GNM_PGC{$lin}=$per_gnome_comp; 	#%complete of lineages in sample by reference

        # GET SUM ALN SCO
        $sco=0; $uniq=0;
        foreach my $hit (keys %{$LIN_HITS{$lin}}){
		#if the hit has only 1 lin, then uniq to lin
		if(keys %{$HIT_LINS{$hit}}==1){$uniq++;}
		else{$uniq+=0;}
                $sco+=$LIN_HITS{$lin}{$hit}; 	#top gene-hit score
        }
        $LIN_SUM_SCO{$lin}=$sco;
        $LIN_SUM_UNIQ{$lin}=$uniq;
}


#IN DECENDING QUALITY, ACCOUNT FOR THE SAMPLE GENES
%GOOD_LINS=();
$start_ag  =keys %AG;
$start_toss=keys %TOSS;
foreach my $lin (sort{$LIN_SUM_UNIQ{$b}<=>$LIN_SUM_UNIQ{$a} or $LIN_SUM_SCO{$b}<=>$LIN_SUM_SCO{$a} or $LIN_GNM_PGC{$b}<=>$LIN_GNM_PGC{$a}} keys %LIN_HITS){
	if(keys %AG < 1){last;}							#stop if no more unaccounted for sample genes left
	$count=0;
	foreach my $gene (keys %{$LIN_GENES{$lin}}){ if(exists($AG{$gene})){delete($AG{$gene}); $count++; }}
	if($count>0){ $GOOD_LINS{$lin}=1; delete($TOSS{$lin});}				#
}
$end_ag	 =keys %AG;
$end_toss=keys %TOSS;
$gkc=keys %GOOD_LINS;
print "start_toss $start_toss start_ag $start_ag end_toss $end_toss end ag $end_ag goodlins $gkc\n";

#OUTPUT LINEAGE TIP BINS
# get tips
foreach my $lin (%GOOD_LINS){ push(@ALL_LINS,$lin); }
@TIPS=();
$al=@ALL_LINS;
while($al>0){ 
	$lin=shift(@ALL_LINS);
	$al=@ALL_LINS;
	if(grep /$lin/, @ALL_LINS){next;}
	if(grep /$lin/, @TIPS){next;}
	push(@TIPS,$lin);
}

# get reference genomes
open(OUTTIPSTAT, ">", $samp."_TIP_STATS.txt")||die;
print OUTTIPSTAT "tip\tn50\tlin_gc\tlin_sco\tlin_uniq\tlin_pgc\n";
foreach my $lin (@TIPS){
	$lin_gc		=keys %{$LIN_GENES{$lin}};
	$lin_sco	=$LIN_SUM_SCO{$lin};
	$lin_uniq	=$LIN_SUM_UNIQ{$lin};
	$lin_pgc	=$LIN_GNM_PGC{$lin};
	if($lin_gc < $mingen ){next;}
	if($lin_gc < 100 && $lin=~/^[ABE]/){next;}
	if($lin_pgc < 10 && $lin=~/^[ABE]/){next;}

	$tip=$lin; $tip=~s/.*\;//g;
	$sum_con_len=0;
	%TIP_CONS=();
	open(OUTTIP, ">", "PHYLOBIN_".$tip.".fasta")||die;
	foreach my $contig (keys %{$LIN_CONS{$lin}}){
		print OUTTIP ">$contig\n$CON_SEQ{$contig}\n";
		$sum_con_len += length($CON_SEQ{$contig});
		$TIP_CONS{$contig}=$CON_LEN{$contig};
	}
	$min_n50=$sum_con_len/2;
	$nextlen=0;
	foreach my $contig (sort{$TIP_CONS{$b}<=>$TIP_CONS{$a}} keys %TIP_CONS){
		$len = $TIP_CONS{$contig};
		$nextlen += $len;
		if($nextlen > $min_n50){ $n50=$track_len; last;}
		$track_len = $len;
	}
	print OUTTIPSTAT "$tip\t$n50\t$lin_gc\t$lin_sco\t$lin_uniq\t$lin_pgc\n";
}
undef(%AG);
undef(@TIPS);
undef(%TOSS);
undef(%LIN_GENES);
undef(%LIN_CONS);
undef(%LIN_HITS);
undef(%LIN_SUM_SCO);
undef(%LIN_SUM_UNIQ);
undef(%LIN_GNM_PGC);
undef(%TOT_LIN_HITS);
########################################################################
########################################################################






########################################################################
##########	COMPILE GENE AND CONTIG STATS		################
########################################################################
open(OUTINFO, 	">", $samp."_contigs_".$version.".txt")||die "unable to open outinfo: $!\n";
open(OUTCYTO,    ">", $samp."_community_".$version.".txt")||die "unable to open outcyto: $!\n";
open(OUTFID, 	">", $samp."_functions_".$version.".txt")||die "unable to open outfid: $!\n";
open(OUTCPD, 	">", $samp."_compounds_".$version.".txt")||die "unable to open outcpd: $!\n";
open(OUTCINT, 	">", $samp."_cpd_int_".$version.".txt")||die "unable to open outcint: $!\n";


## LOOP THROUGH CONTIGS
$nov_con=0;
$nog_con=0;
foreach my $contig (sort(keys %CON_LEN)){
	
		#SET BASICS
		$conlen		=$CON_LEN{$contig};
		$conpgc		=$CON_PGC{$contig};
		$conrpkm	=$CON_RPKM{$contig};
		$condepth	=$CON_DEPTH{$contig};
		$congcnt	=keys %{$CON_GENES{$contig}};
		$gene		=0;
		$genepos	=0;
		$genepart	=0;
		$geneuniq	=0;
		$genelen	=0;
		$genepgc	=0;
		$generpkm	=0;
		$genedepth	=0;
		$genetids	=0;
		$genelca	="NO_GENE";
		$contids	=0;
		$conlca		="NO_GENE";
		@FUNCS=();
		for my $i (7..21){push(@FUNCS,0);} 	#set empty funcs
		$funcs		=join(";",@FUNCS);
		$reac=0; $prod=0; $tran=0;		#set emtpy cpds
		

		### GET CONTIG INFO ###
		#######################
		%CON_TIDS_G=(); %CON_TIDS_I=(); %CON_TIDS_B=();
		@CON_LINS_G=(); @CON_LINS_I=(); @CON_LINS_B=();
		foreach my $gene (keys %{$CON_GENES{$contig}}){ 
			foreach my $hit (keys %{$GENE_HITS{$gene}}){
				foreach my $lin (keys %{$HIT_LINS{$hit}}){
					$tid = $LIN_TID{$lin};
					if(!exists($GOOD_LINS{$lin}) || $lin!~/^[ABEM]/){next;}
					   if( $lin=~/UNCLASSIFIED.*UNCLASSIFIED.*UNCLASSIFIED/){ push(@CON_LINS_B,$lin); $CON_TIDS_B{$tid}=1;}
					elsif( $lin=~/^(MONA|BACT|ARCH|EUKA)/){ push(@CON_LINS_G,$lin); $CON_TIDS_G{$tid}=1;}
					elsif( $lin=~/^MICROBIOME/){ 		push(@CON_LINS_I,$lin); $CON_TIDS_I{$tid}=1;}
					else{					push(@CON_LINS_B,$lin); $CON_TIDS_B{$tid}=1;}
		}	}	}
		
		
		### MAKE CONTIG LCA AND CYTO ###
		################################
		@CON_TIDS=();
		#IF HAS 1. ONLY NOVEL GENES, 2. NO GENES, ELSE 3. GENES WITH SOME UMRAD HITS
		   if(keys %{$CON_GENES{$contig}}>0 && keys %CON_TIDS_G < 1 && keys %CON_TIDS_I < 1 && keys %CON_TIDS_B < 1){ 
				$conlca="NOVEL"; 
				if($conlen > 50000){
					$nov_con++; 
					$out_con = "LONG_NOVEL_CON_".$nov_con.".txt";
					open(OUTLC, ">", $out_con)||die;
					print OUTLC ">$contig\n$CON_SEQ{$contig}\n";
		}		}
		elsif(keys %{$CON_GENES{$contig}}<1){ 
				$conlca="N0_GENE";
				if($conlen > 50000){
					$nog_con++; 
					$out_con = "LONG_NO_GENE_CON_".$nog_con.".txt";
					open(OUTLC, ">", $out_con)||die;
					print OUTLC ">$contig\n$CON_SEQ{$contig}\n";
		}		}
		else{	   if($CON_LINS_G[0]=~/\w/){ 	$conlca=MakeLCA(@CON_LINS_G); foreach my $tid (keys %CON_TIDS_G){push(@CON_TIDS,$tid);}}
			elsif($CON_LINS_I[0]=~/\w/){    $conlca=MakeLCA(@CON_LINS_I); foreach my $tid (keys %CON_TIDS_I){push(@CON_TIDS,$tid);}}
			else{				$conlca=MakeLCA(@CON_LINS_B); foreach my $tid (keys %CON_TIDS_B){push(@CON_TIDS,$tid);}}
		}

		@CON_TIDS=nsort(@CON_TIDS);
		$contids=join(";",@CON_TIDS);
		@TAXON=split(";",$conlca);
		$conlca=~/^(.)/; $king=$1;
		for($i=0; $i<=$#TAXON; $i++){
		        $type = $king."_".$i;
		        if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
		        else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
			$TOT_PHYLA{$phyla}=1;
		        $CYTO_CON_LCA{$phyla}++;
		        $CYTO_CON_LEN{$phyla}+=$conlen;
		        $CYTO_CON_RPKM{$phyla}+=$conrpkm;
		}

				
		### GET GENE INFO ###
		#####################
		foreach my $gene (sort{$GENE_START{$a}<=>$GENE_START{$b}} keys %{$CON_GENES{$contig}}){ 
			%FUNC_LIST=(); 
			%CPD_LIST=();
			%GENE_TIDS_M=(); %GENE_TIDS_G=(); %GENE_TIDS_I=(); %GENE_TIDS_B=();
			@GENE_LINS_M=(); @GENE_LINS_G=(); @GENE_LINS_I=(); @GENE_LINS_B=();
			$maxsco=0;
			$besthit='';
			$name=''; 
			
			#BASIC
			$genepos	=$GENE_START{$gene};
			$genepart	=$GENE_PART{$gene};
			$geneuniq	=$GENE_UNIQ{$gene};
			$genepgc	=$GENE_PGC{$gene};
			$genelen	=$GENE_LEN{$gene};
			$generpkm	=$GENE_RPKM{$gene};
			$genedepth	=$GENE_DEPTH{$gene};
			$genelca	="NOVEL";

			### LOOP THROUGH UMRAD HITS
			foreach my $hit (keys %{$GENE_HITS{$gene}}){
				#best UMRAD-gene alignment score
				$sco = $GENE_HITS{$gene}{$hit};
				if($sco > $maxsco){$maxsco=$sco; $besthit=$hit; } 
				if(length($HIT_NAMES{$hit})>$name){$name=$HIT_NAMES{$hit};}
				
				#taxonomy
				foreach my $lin (keys %{$HIT_LINS{$hit}}){
					$tid = $LIN_TID{$lin};
                                        if(!exists($GOOD_LINS{$lin}) || $lin!~/^[ABEM]/){next;}
					   if( $lin=~/UNCLASSIFIED.*UNCLASSIFIED.*UNCLASSIFIED/){push(@GENE_LINS_B,$lin); $GENE_TIDS_B{$tid}=1; }
					elsif( $lin=~/^MONA/){ 			push(@GENE_LINS_M,$lin); $GENE_TIDS_M{$tid}=1; }
					elsif( $lin=~/^(BACT|ARCH|EUKA)/){ 	push(@GENE_LINS_G,$lin); $GENE_TIDS_G{$tid}=1; }
					elsif( $lin=~/^MICROBIOME/){ 		push(@GENE_LINS_I,$lin); $GENE_TIDS_I{$tid}=1; }
					elsif( $lin=~/^[ABEM]/){		push(@GENE_LINS_B,$lin); $GENE_TIDS_B{$tid}=1; }
					else{}
				}
				
				#functions
				foreach my $i (keys %{$HIT_FUNCS{$hit}}){
					foreach my $id ( keys %{$HIT_FUNCS{$hit}{$i}} ){
						if($id=~/\w/){$FUNC_LIST{$i}{$id}=$maxsco;}
				}	}	
				
				#compounds
				foreach my $i (keys %{$HIT_CPDS{$hit}}){
					foreach my $id (keys %{$HIT_CPDS{$hit}{$i}}){
						if($id=~/\w/){$CPD_LIST{$i}{$id}=$maxsco;}
			}	}	}


			##############################################
			#####	 	MAKE ALL LINS CYTO	######
			@GEN_TIDS=();
			@GENE_LINS=();
			   if($GENE_LINS_M[0]=~/^MONA/){		@GENE_LINS=@GENE_LINS_M; $genelca=MakeLCA(@GENE_LINS_M); foreach my $tid (keys %GENE_TIDS_M){push(@GEN_TIDS,$tid);}}
			elsif($GENE_LINS_I[0]=~/^MICRO/){		@GENE_LINS=@GENE_LINS_I; $genelca=MakeLCA(@GENE_LINS_I); foreach my $tid (keys %GENE_TIDS_I){push(@GEN_TIDS,$tid);}}
			elsif($GENE_LINS_G[0]=~/^(ARCH|BACT|EUKA)/){	@GENE_LINS=@GENE_LINS_G; $genelca=MakeLCA(@GENE_LINS_G); foreach my $tid (keys %GENE_TIDS_G){push(@GEN_TIDS,$tid);}}
			elsif($GENE_LINS_B[0]=~/^[ABEM]/){		@GENE_LINS=@GENE_LINS_B; $genelca=MakeLCA(@GENE_LINS_B); foreach my $tid (keys %GENE_TIDS_B){push(@GEN_TIDS,$tid);}}
			else{ $genelca="NOVEL"; }
			@GEN_TIDS=nsort(@GEN_TIDS);
			$genetids = join(";", @GEN_TIDS);
			if($GEN_TIDS[0] !~ /\d/){$genetids=0;}
			%seen=(); @GENE_LINS = grep { !$seen{$_}++ } @GENE_LINS;
			foreach my $lin (@GENE_LINS){
				@TAXON=split(";",$lin);
				$lin=~/^(.)/; $king=$1;
				for($i=0; $i<=$#TAXON; $i++){
					$type = $king."_".$i;
					if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
					else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
					$TOT_PHYLA{$phyla}=1;
					$CYTO_TOT_GEN{$phyla}++;
					$CYTO_GEN_LEN{$phyla}+=$genelen;
					$CYTO_GEN_RPK{$phyla}+=$generpkm;
					$CYTO_GEN_CON_RPK{$phyla}+=$conrpkm;
					$CYTO_GEN_DEP{$phyla}+=$genedepth;
					$CYTO_GEN_UNI{$phyla}+=$geneuniq;
					if($genepart=~/^0+$/){$CYTO_GEN_PAR{$phyla}++;}
			}	}
			#### MAKE GENE LCA CYTO	
			@TAXON=split(";",$genelca);
			$genelca=~/^(.)/; $king=$1;
			for($i=0; $i<=$#TAXON; $i++){
			        $type = $king."_".$i;
			        if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
			        else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
			        $CYTO_GEN_LCA{$phyla}++;
				$TOT_PHYLA{$phyla}=1;
			}
			#############################################
			#############################################
				
		
		
			### FUNCTIONS
			#############################################
			@FOUT=();
			%RXNS=();
			%TCDB=();
			for my $i (7..21){
				$j=$i-7;
				@FUNCS=(); 
				foreach my $id (keys %{$FUNC_LIST{$i}}){
					if($id!~/\w/){next;}
					push(@FUNCS,$id); #TRACK GENE FUNCS i 7-21
					$FUNC_GENES{$id}++; #track all funcs for later output
					$FUNC_SCOS{$id}+=$FUNC_LIST{$i}{$id};
					$FUNC_LCA{$id}{$genelca}++;
					$FUNC_RPKM{$id}+=$generpkm;
					if($i>18){$RXNS{$id}=1; }
					if($i==11){						
						$id=~s/.*?(TCDB\:[A-Z\d\.]+).*?/$1/;
						$TCDB{$id}=1;}
					if($i=~/^(7|8|9)$/){last;}
				}
				@FUNCS=nsort(@FUNCS);
				$func=join(";",@FUNCS);
				if($func!~/\w/){$func='';}
				$FOUT[$j]=$func;
			}
			$funcs=join("\t",@FOUT);
			#############################################
			#############################################
				
				

			### COMPOUNDS
			#############################################
			@REACT=(); @PRODK=(); @TRANS=();
			@TAXON=split(";",$genelca);
			for my $i (22..24){
				if($i==22){$cpd_type="REACTANT"; }
				if($i==23){$cpd_type="PRODUCT"; }
				if($i==24){$cpd_type="TRANSPORT";}
				foreach my $id (keys %{$CPD_LIST{$i}}){ 
					#general cpd summary
					$CPD_GENES{$id}{$cpd_type}++;
					$CPD_SCO{$id}{$cpd_type}+=$CPD_LIST{$i}{$id};
					$CPD_LCA{$id}{$cpd_type}{$genelca}++;
					$CPD_RPKM{$id}{$cpd_type}+=$conrpkm;
					#gene cpd info
					if($i==22){push(@REACT,$id);}
					if($i==23){push(@PRODK,$id);}
					if($i==24){push(@TRANS,$id);}
					#cpd-rxn interactions
					$root=$id."_".$CPD_NAMES{$id};
					if($i==22 && $id =~ /CHEBI/){
						foreach my $rxn (keys %RXNS){
							for($j=1; $j<=$#TAXON; $j++){
				    			        $type = $king."_".$j;
					                        $source=$id."_R_".$rxn."_".$TAXON[$j-1];
				            			$target=$id."_R_".$rxn."_".$TAXON[$j];
				            			$phyla = "$type\t$source\t$target";
				            			$CPD_INT_CYTO{$id}{R}{$phyla}{$gene}=$conrpkm;
								if($i==1){
									$type = $king."_0";
									$phyla = "$type\t$root\t$source";
									$CPD_INT_CYTO{$id}{R}{$phyla}{$gene}=$conrpkm;
								}
				        }	}	}
					if($i==23 && $id =~ /CHEBI/){
						foreach my $rxn (keys %RXNS){
							for($i=1; $i<=$#TAXON; $i++){
			                       		        $type = $king."_".$i;
					                        $source=$id."_P_".$rxn."_".$TAXON[$i-1];
			                            		$target=$id."_P_".$rxn."_".$TAXON[$i];
			                              		$phyla = "$type\t$source\t$target";
			                              		$CPD_INT_CYTO{$id}{P}{$phyla}{$gene}=$conrpkm;
								if($i==1){
									$type = $king."_0";
									$phyla = "$type\t$root\t$source";
									$CPD_INT_CYTO{$id}{P}{$phyla}{$gene}=$conrpkm;
								}
					}	}	}
					if($i==24 && $id =~ /CHEBI/){
						foreach my $tcdb (keys %TCDB){
							for($i=1; $i<=$#TAXON; $i++){
			                       		        $type = $king."_".$i;
					                        $source=$id."_T_".$tcdb."_".$TAXON[$i-1];
			                             		$target=$id."_T_".$tcdb."_".$TAXON[$i];
			                              		$phyla = "$type\t$source\t$target";
			                              		$CPD_INT_CYTO{$id}{T}{$phyla}{$gene}=$conrpkm;
								if($i==1){
									$type = $king."_0";
									$phyla = "$type\t$root\t$source";
									$CPD_INT_CYTO{$id}{T}{$phyla}{$gene}=$conrpkm;
								}
			}	}	}	}	}
			@REACT=nsort(@REACT);
			@PRODK=nsort(@PRODK);
			@TRANS=nsort(@TRANS);
			$reac=join(";",@REACT);
			$prod=join(";",@PRODK);
			$tran=join(";",@TRANS);
			#############################################
			#############################################
		
			print OUTINFO "$contig\t$conlen\t$congcnt\t$conpgc\t$conrpkm\t$condepth\t";
			print OUTINFO "$gene\t$genepos\t$genepart\t$geneuniq\t$genelen\t$genepgc\t$generpkm\t$genedepth\t";
			print OUTINFO "$genetids\t$genelca\t$contids\t$conlca\t";
			print OUTINFO "$name\t$maxsco\t$besthit\t$funcs\t$reac\t$prod\t$tran\n";
		}	
}


#OUTPUT CYTOSCAPE
print OUTCYTO "type\tSource\tTarget\tsum_con_len\tsum_gene_len\tsum_tot_gene\tcomplete_genes\tuniq_genes\tcon_rpkm\tgene_rpkm\tgene-con_rpkm\tgene_depth\tgene_lca\tcon_lca\n";
foreach my $phyla (%TOT_PHYLA){
	if($phyla!~/[A-Z]/){next;}
	$con_len    =$CYTO_CON_LEN{$phyla};		if($con_len !~/\w/){$con_len =0;}
	$gen_len    =$CYTO_GEN_LEN{$phyla};		if($gen_len !~/\w/){$gen_len =0;}
	$tot_gen    =$CYTO_TOT_GEN{$phyla};		if($tot_gen !~/\w/){$tot_gen =0;}
	$part_gen   =$CYTO_GEN_PAR{$phyla};		if($part_gen !~/\w/){$part_gen =0;}	
	$uniq_gen   =$CYTO_GEN_UNI{$phyla};		if($uniq_gen !~/\w/){$uniq_gen =0;}
	$con_rpkm   =$CYTO_CON_RPKM{$phyla};		if($con_rpkm !~/\w/){$con_rpkm =0;} 		$con_rpkm=~s/(\.\d\d\d).*/$1/;
	$gen_rpkm   =$CYTO_GEN_RPK{$phyla};		if($gen_rpkm !~/\w/){$gen_rpkm =0;}		$gen_rpkm=~s/(\.\d\d\d).*/$1/;
	$gen_con_rpkm=$CYTO_GEN_CON_RPK{$phyla};	if($gen_con_rpkm !~/\w/){$gen_con_rpkm =0;}	$gen_con_rpkm=~s/(\.\d\d\d).*/$1/;
	$gen_dept   =$CYTO_GEN_DEP{$phyla};		if($gen_dept !~/\w/){$gen_dept =0;}		$gen_dept=~s/(\.\d\d\d).*/$1/;
	$gen_lca    =$CYTO_GEN_LCA{$phyla};		if($gen_lca !~/\w/){$gen_lca =0;}
	$con_lca    =$CYTO_CON_LCA{$phyla};		if($con_lca !~/\w/){$con_lca =0;}
	print OUTCYTO "$phyla\t$con_len\t$gen_len\t$tot_gen\t$part_gen\t$uniq_gen\t$con_rpkm\t$gen_rpkm\t$gen_con_rpkm\t$gen_dept\t$gen_lca\t$con_lca\n";
}


#OUTPUT FUNCTION TABLE
print OUTFID "fid\tfname\tFUNC_GENES\tFUNC_SCOS\tFUNC_RPKM\tfidlca\n";
foreach my $fid (sort(keys %FUNC_GENES)){
	@FIDLINS=();
	foreach my $glca (keys %{$FUNC_LCA{$fid}}){if($glca=~/^(ARCH|BACT|EUKA|MONA)/){push(@FIDLINS,$glca);}}
	$fidlca=MakeLCA(@FIDLINS);
	$fname = $FUNC_NM{$fid};
	print OUTFID "$fid\t$fname\t$FUNC_GENES{$fid}\t$FUNC_SCOS{$fid}\t$FUNC_RPKM{$fid}\t$fidlca\n";
}


#OUTPUT COMPOUND TABLE
print OUTCPD "cpd\tcpd_name\ttype\tCPD_GENES\tCPD_SCO\tCPD_RPKM\tcpdlca\n";
foreach my $cpd (sort(keys %CPD_GENES)){
	foreach my $type (sort(keys %{$CPD_GENES{$cpd}})){
		$cpd_name =$CPD_NAMES{$cpd};
		$cpd_gc	  =$CPD_GENES{$cpd}{$type};
		$cpd_sco  =$CPD_SCO{$cpd}{$type};
		$cpd_rpkm =$CPD_RPKM{$cpd}{$type};
		@CPDLINS=();
		foreach my $glca (keys %{$CPD_LCA{$cpd}{$type}}){push(@CPDLINS,$glca);}
		$cpdlca=MakeLCA(@CPDLINS);
		print OUTCPD "$cpd\t$cpd_name\t$type\t$cpd_gc\t$cpd_sco\t$cpd_rpkm\t$cpdlca\n";
}	}


#OUTPUT CPD INTERACTION
print OUTCINT "cpd\tprot\ttype\tsource\ttarget\tgc\tsum_rpkm\n";
foreach my $id (sort(keys %CPD_INT_CYTO)){
	if(!exists($CPD_INT_CYTO{$id}{T})){next;} #no transporters, no care
	foreach my $dir (sort(keys %{$CPD_INT_CYTO{$id}})){
		foreach my $phyla (sort(keys %{$CPD_INT_CYTO{$id}{$dir}})){
			$gc=keys %{$CPD_INT_CYTO{$id}{$dir}{$phyla}};
			$sum_rpkm=0;
			foreach my $gene (keys %{$CPD_INT_CYTO{$id}{$dir}{$phyla}}){ $sum_rpkm+=$CPD_INT_CYTO{$id}{$dir}{$phyla}{$gene}; }
			print OUTCINT "$id\t$dir\t$phyla\t$gc\t$sum_rpkm\n";
}	}	}









###########################
### 	SUBROUTINES 	###
###########################

sub MakeLCA{
        my @ARR;
        %seen=();
        @array1 = @_;
        @array1 = grep { !$seen{$_}++ } @array1;
        #get the kingdoms, JIC lineage is NCA
        %LET=();
        foreach my $lin (@array1){
                if($lin !~ /^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/i){next;}
                $lin =~ /^(\w)/; $LET{$1}=1;
        }
        @LET=();
        foreach my $let (sort(keys %LET)){push(@LET,$let);}
        $let=join("",@LET);
        $LCA='';
        $len1 = @array1;
        if($len1 == 1){$LCA = $array1[0]; }
        elsif($len1 > 1){
                $first = $array1[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @array1);
                        $len2 = @matched;
                        if($len2 == $len1){push(@ARR, $alevel);}
                        else{last;}
                }
                $len3 = @ARR;
                if($len3 > 1){$LCA = join(";", @ARR);}
                elsif($len3==1){$LCA = $ARR[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }
        #add kingdoms to NCA
        if($LCA eq "NCA"){ $LCA.="-".$let; }
        return($LCA);
}
		

sub CleanNames{
        $name = $_[0];

        #remove pointless ambiguators
        $name =~ s/(CANDIDATUS|CANDIDATUAS|CANDIDATE|VOUCHERED|UNDESCRIBED|UNSCREENED|UNKNOWN|UNCULTIVATED|UNCULTURED)\s*/UNCHARACTERIZED\_/g;
        $name =~ s/(UNIDENTIFIED|UNCLASSIFIED|CONTAMINATION|SCREENED|UNASSIGNED|PUTATIVE|HYPOTHETICAL)\s/UNCHARACTERIZED\_/g;
        $name =~ s/(UNCHARACTERIZED\_)+/UNCHARACTERIZED\_/g;
        $name =~ s/\-*LIKE\s*/\_/g;

        #remove junk punctuation/standardize
        $name =~ s/\s+/_/g;
        $name =~ s/[^\w]+/_/g;
        $name =~ s/\_+/\_/g;
        $name =~ s/(^\_+|\_+$)//g;

        return($name);
}
