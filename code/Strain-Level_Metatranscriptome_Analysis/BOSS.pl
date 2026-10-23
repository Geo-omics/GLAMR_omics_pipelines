use Cwd;
my $dir = getcwd;
$dir = $dir."/";
print "dir is $dir\n";
$dir2=$dir."BOSS/";
print "dir2 is $dir2\n";

#GET PBS stuff
$pbs1 = "\#PBS -l nodes=1:ppn=cores";
$pbs2 = "\#PBS -l pmem=corememGB";
$pbs3 = "\#PBS -l walltime=runtime";
$pbs4 = "\#PBS -l qos=flux";
$pbs5 = "\#PBS -A account";
$pbs6 = "\#PBS -q fluxtype";
$pbs7 = "\#PBS -N input";
$pbs8 = "\#PBS -j oe";
$pbs9 = "\#PBS -M email";
$pbs10 = "\#PBS -m abe";
$pbs11 = "\#PBS -V";
$pbs12 = "\#PBS -o ".$dir2."run_logs/sample_type.outfile\n";
$pbs13 = "if [ -n \"\$PBS_O_WORKDIR\" ]; then cd \$PBS_O_WORKDIR; fi";

#extra for muscato:
	$pbs14 = "\#PBS -l ddisk=40GB";
	$pbs15 = "export LC_ALL=C";
	$pbs16 = "export GOPATH=\${HOME}/go";
	$pbs17 = "export GOBIN=\${GOPATH}/bin";
	$pbs18 = "export PATH=\${GOBIN}:\${PATH}";
 	$pbs19 = "export GOGC=20";
 	
$unixlines = "mkdir ".$dir2."run_logs\n";


#GET USER PARAMS
print "Muscato will process all the fastq files listed in design.matrix that are located in the FASTQ folder.\n";
print "If you plan on trimming, make sure your adapter sequences are in the adapters.fa file.\n";
print "Are you ready to proceed? (y/n)\n"; 
$go = <STDIN>;
if($go =~ /n/){die;}

#GET USER EMAIL
$email='';
$try=0;
while($email !~ /\@/ || $email !~ /\./ || $email !~ /\w/){
	print "Please enter a valid email.\n";
	$email = <STDIN>;
	$email =~ s/\s+//g;
	$try++;
	if($try == 5){print "There seems to be a problem with entering your email. Please contact tealfurn@umich.edu for assistance\n"; die;}
}
$pbs9 =~ s/email/$email/;

#GET USER FLUX ACCT
$account='';
$try=0;
while($account !~ /flux/){
	print "Please enter you flux user account (eg. trfurn_fluxm).\n";
	$account = <STDIN>;
	$account =~ s/\s+//g;
	$try++;
	if($try == 5){print "There seems to be a problem with entering your flux account. Please contact tealfurn@umich.edu for assistance\n"; die;}
}
$pbs5 =~ s/account/$account/;

#GET USER FLUX TYPE
$fluxtype='';
$try=0;
while($fluxtype !~ /flux/){
	print "Please enter the flux type (eg. fluxm or fluxod).\n";
	$fluxtype = <STDIN>;
	$fluxtype =~ s/\s+//g;
	$try++;
	if($try == 5){print "There seems to be a problem with entering your flux type. Please contact tealfurn@umich.edu for assistance\n"; die;}
}
$pbs6 =~ s/fluxtype/$fluxtype/;

#GET READ LENGTH AND MIN LENGTH CUTOFF
print "What is your read length?\n";
$rlen = <STDIN>;
$rlen =~ s/\D//g;
$try=0;
while($rlen !~ /\d{2,10000}/ ){ print "Read length must be 25 or greater. Please enter your read length\n";
	$rlen = <STDIN>;
	$rlen =~ s/\D+//g;
	$try++;
	if($try == 5){print "There seems to be a problem with entering your read length. Please contact tealfurn@umich.edu for assistance\n"; die;}
}
	
print "What minimum read length cutoff?\n";
$minrl = <STDIN>;
$minrl =~ s/\D+//g;
$try=0;
while($minrl !~ /\d{2,10000}/ || $minrl > $rlen){ print "Read length is $rlen. What minimum read length cutoff, less than or equal to $rlen?\n";
	$minrl = <STDIN>;
	$minrl =~ s/\D+//g;
	$try++;
	if($try == 5){print "There seems to be a problem with entering a minimum length. Please contact tealfurn@umich.edu for assistance\n"; die;}
}
if($minrl < 25){print "Muscato not designed for reads less than 25nt\n"; die;}

#EITHER USER DOES WHOLE PIPELINE OR CHOOSES/MODS PIECES
	print "Do you want run the full default pipeline? (y/n)\n";
	$doall = <STDIN>;
	$doall =~ s/.*(y|n).*/\1/i;
	while($doall !~ /(y|n)/i){ print "Do you want run the full default pipeline? (y/n)\n"; 
		$doall = <STDIN>;
		$doall =~ s/.*(y|n).*/\1/i;}
		
if($doall =~ /y/){
	$dotrim = "y";
	$doclean = "y";
	$doribo = "y";
	$doalign = "y";
	$docount = "y";
	$doRDEG = "y";
	$pmatch = 0.95;
	$unixlines .= "mkdir ".$dir."ALIGNED\n";
	$unixlines .= "mkdir ".$dir."NOTALIGNED\n";
	$maxm = 1000;
}
else{
	#GET USER JOBS
	print "Do you want to trim your reads? (y/n)\n";
	$dotrim = <STDIN>;
	$dotrim =~ s/.*(y|n).*/\1/i;
	while($dotrim !~ /(y|n)/i){ print "Do you want to trim your reads? (y/n)\n"; 
		$dotrim = <STDIN>;
		$dotrim =~ s/.*(y|n).*/\1/i;}
		
	print "Do you want to remove polyA/T and low-complexity reads? (y/n)\n";
	$doclean = <STDIN>;
	$doclean =~ s/.*(y|n).*/\1/i;
	while($doclean !~ /(y|n)/i){ print "Do you want to remove polyA/T and low-complexity reads? (y/n)\n"; 
		$doclean = <STDIN>;
		$doclean =~ s/.*(y|n).*/\1/i;}	
		
	print "Do you want to remove ribosomal contamination? (y/n)\n";
	$doribo = <STDIN>;
	$doribo =~ s/.*(y|n).*/\1/i;
	while($doribo !~ /(y|n)/i){ print "Do you want to remove ribosomal contamination? (y/n)\n";
	        $doribo = <STDIN>;
	        $doribo =~ s/.*(y|n).*/\1/i;}

	print "Do you want to align your reads? (y/n)\n";
	$doalign = <STDIN>;
	$doalign =~ s/.*(y|n).*/\1/i;
	while($doalign !~ /(y|n)/i){ print "Do you want to align your reads? (y/n)\n";
		$doalign = <STDIN>;
		$doalign =~ s/.*(y|n).*/\1/i;}
		
	if($doribo =~ /y/i || $doalign =~ /y/i){
		
		print "What percent identity do you want? (enter value between 0.3 and 1)\n";
		$pmatch = <STDIN>;
		$pmatch =~ s/\s+//g;
		while($pmatch !~ /\d/ || $pmatch < 0.3 || $pmatch > 1){ print "What percent identity do you want? (enter value between 0.3 and 1)\n";
			$pmatch = <STDIN>;
			$pmatch =~ s/\s+//g;}
			
		print "What is the maximum number of matches to find? (1000 for faster, 10000 for more accuracy)\n";
		$maxm = <STDIN>;
		$maxm =~ s/\s+//g;
		while($maxm !~ /\d/ || $maxm < 10 || $maxm > 50000){ print "please pick a number between 10 and 50,000\n";
			$maxm = <STDIN>;
			$maxm =~ s/\s+//g;}
		$unixlines .= "mkdir ".$dir."ALIGNED\n";
		$unixlines .= "mkdir ".$dir."NOTALIGNED\n";
	}
			
	print "Do you want to convert read matches to gene counts? (y/n)\n";
	$docount = <STDIN>;
	$docount =~ s/.*(y|n).*/\1/i;
	while($docount !~ /(y|n)/i){ print "Do you want to align your reads? (y/n)\n";
		$doalign = <STDIN>;
		$doalign =~ s/.*(y|n).*/\1/i;}
		
	print "Do you want to run differential gene expression analysis? (y/n)\n";
	$doRDEG = <STDIN>;
	$doRDEG =~ s/.*(y|n).*/\1/i;
	while($doRDEG !~ /(y|n)/i){ print "Do you want to run differential gene expression analysis? (y/n)\n";
		$doRDEG = <STDIN>;
		$doRDEG =~ s/.*(y|n).*/\1/i;}
}



#GET FASTQ FILES
$filedir = $dir.'FASTQ/';
opendir(DIR, $filedir) || die "Can't opendir $filedir: $!";
my @files= readdir DIR;
foreach my $file (@files) { 
	$file = $filedir.$file;
	if(-f $file && $file =~ /\.(fastq|fq)$/i ){ 
	push(@FQS, $file);}
}


#GET DESIGN 
$design = $dir.'design.matrix';
open(DESIGN, $design)||die "unable to open $design:$!\n";
$i=0;
while(<DESIGN>){
	$_ =~ s/[\n\r]+//;
	if($_ !~ /\w/){next;}
	@stuff = split("\t", $_);
	if($i==0){for $j (1..$#stuff){$CONDS{$j}=$stuff[$j];}}
	else{ 
		@fc = grep (/\Q$stuff[0]\E/, @FQS);
		push(@goodfiles, $fc[0]);
		$FILE2SAMP{$fc[0]}=$stuff[0];
		$flc='';
		$flc = @fc;
		if($flc != 1){print "sample $stuff[0] in your design.matrix is either missing or there are two files with $stuff[0] in the name.\nFix the fastq and/or design.matrix and retype perl BOSS.pl into the bash prompt.\n"; die;}
		for $j (1..$#stuff){
			if($stuff[$j]==1){$FILECOND{$CONDS{$j}}.="$fc[0]\n";}
		}
	}
	$i++;
}
@FQS = @goodfiles;


##TRIMMOMATIC##
if($dotrim =~ /y/i){ 
	#TRIM LOW QUALITY BP
	print "Please enter a minimum Phred score (10-35)\n";     
	$prs = <STDIN>;
    $prs =~ s/\D+//g;
    while($prs < 10 || $prs > 35){ print "Please enter a Phred score between 10 and 35\n";
        $prs = <STDIN>;
        $prs =~ s/\D+//g;}

	#CROP FIRST BP OF READ
	print "How many bp would you like to trim from the start of the sequence? (0 to 10)\n";
    	$hedc = <STDIN>;
    	$hedc =~ s/\D+//g;
    	while($hedc < 0 || $hedc > 10){ print "Please enter a trim number between 0 to 10\n";
        	$hedc = <STDIN>;
        	$hedc =~ s/\D+//g;}
		if($minrl-$hedc < 25){print "Muscato not designed for reads less than 25nt\n"; die;}
	
	#CHECK COMMAND
	$unixlines .= "mkdir ".$dir."CLEANED\n";  ### goes in BIGBOSS.pbs
	$command = "java -jar \$TRIMM_JAR/trimmomatic-0.36.jar SE -phred33 samp_fastq ".$dir."CLEANED/sample.trim ILLUMINACLIP:".$dir."adapters.fa:2:21:10 LEADING:prs TRAILING:prs SLIDINGWINDOW:4:prs HEADCROP:hedc MINLEN:minrl";
	$command =~ s/minrl/$minrl/g;
	$command =~ s/hedc/$hedc/g;
	$command =~ s/prs/$prs/g;
	print "current trim program settings are:\n$command\nContinue? (y/n)\n";
	$incmd = <STDIN>;
	if($incmd !~ /y/i){
		print "paste your command or type \"n\"\n";
		$incmd = <STDIN>;
		if($incmd !~ /trimmomatic/i){
			print "Copy (highlight and CTRL+C) and modify the above command, then re-run 'perl BOSS.pl' (do not change input/output files!).\n";
			print "If you have paired end reads - merge using JGI bbmerge before this step and make sure only things you want aligned are in the FASTQ folder.\n";
			print "If not using ScriptSeq epidemiology adapters, put fasta-formatted adapters into adapters.fa in the BOSS folder.\n";
			print "Visit: http://www.usadellab.org/cms/?page=trimmomatic for more assistance\n";
			die;}
		else{$command = <STDIN>;}
	}		
	
	#MAKE PBS TRIM FILES
	foreach(@FQS){
		$type = "trim";
		$sample = $FILE2SAMP{$_};
		$jobname = $dir2.$sample."_trim.pbs";
		$jobid = $sample."_trim";
		$mycommand= $command;
		$mycommand =~ s/samp_fastq/$_/g;
		$mycommand =~ s/sample/$sample/g;
		
		$mypbs1 = $pbs1;
		$mypbs2 = $pbs2;
		$mypbs3 = $pbs3;
		$mypbs7 = $pbs7;
		$mypbs12 = $pbs12;
		
		$mypbs1 =~ s/cores/1/;
		$mypbs2 =~ s/coremem/4/;
		$mypbs3 =~ s/runtime/8:00:00/;
		$mypbs7 =~ s/input/$jobid/;		
		$mypbs12 =~ s/sample/$sample/;
        $mypbs12 =~ s/type/trim/;

		open(TRIMJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
		print TRIMJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$mypbs12\n$pbs13\n\n$mycommand\n";
		print TRIMJOB "qsub ".$dir2.$sample."_check.pbs\n";
	}
}




##CLEAN READS##
if($doclean =~ /y/i){
	if($dotrim !~ /y/i){$unixlines .= "mkdir ".$dir."CLEANED\n";}
	$cleanpl = $dir2."RemovePoly.pl";
	open(INCLEAN, $cleanpl)||die "unable to open $cleanpl:$!\n";
	@inclean = <INCLEAN>;
	$inclean = join('', @inclean);
	foreach(@FQS){
		$type = "clean";
		$sample = $FILE2SAMP{$_};
		$jobname = $dir2.$sample."_clean.pbs";
 		$jobid = $sample."_clean";
 		$command = "perl ".$dir2."Clean_".$sample.".pl";
		if($dotrim =~ /y/i){$file = $dir."CLEANED/$sample\.trim";}
		else{$file = $_;}
		
		#Make Clean_sample.pl
		$myclean = $inclean;
		$myclean =~ s/samp_fastq/$file/g;
		$myclean =~ s/sample/$sample/g;
		$outpl = $dir2."Clean_$sample\.pl";  
		open(OUTCLEAN, ">", $outpl)||die "unable to open $outpl:$!\n";
		print OUTCLEAN "$myclean\n";
		
		#Make Sample_clean.pbs	
		$mypbs1 = $pbs1;
		$mypbs2 = $pbs2;
		$mypbs3 = $pbs3;
		$mypbs7 = $pbs7;
    	$mypbs12 =$pbs12;

		$mypbs1 =~ s/cores/1/;
		$mypbs2 =~ s/coremem/1/;
		$mypbs3 =~ s/runtime/10:00:00/;

		$mypbs7 =~ s/input/$jobid/;
		$mypbs12 =~ s/sample/$sample/;
    	$mypbs12 =~ s/type/clean/;
    	
		open(CLEANJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
		print CLEANJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$mypbs12\n$pbs13\n\n$command\n";
		print CLEANJOB "qsub ".$dir2.$sample."_check.pbs\n";
	}
}





#ADJUST WINDOW SETTINGS FOR USER READ LENGTH
if($doribo =~ /y/i || $doalign =~ /y/i){
	if($rlen<=50){
		$wl = 15;
		$sub = ($rlen - $wl - 7)/3;
		$one = int($sub);
		$two = int($sub*2);
		$three = int($sub*3);
		$winnums = "$one,$two,$three";
		$command1 =~ s/winnum/$winnums/;}
	elsif($rlen<=75){
		$wl = 20;
		$sub = ($rlen - $wl - 7)/3;
		$one = int($sub);
		$two = int($sub*2);
		$three = int($sub*3);
		$winnums = "5,$one,$two,$three";
		$command1 =~ s/winnum/$winnums/;}	
	elsif($rlen<=200){
		$wl = 20;
		$sub = ($rlen - $wl - 7)/3;
		$one = int($sub);
		$two = int($sub*2);
		$three = int($sub*3);
		$winnums = "5,$one,$two,$three";
		$command1 =~ s/winnum/$winnums/;}	
	else{
		$wl = 25;
		$sub = ($rlen - $wl - 10)/3;
		$one = int($sub);
		$two = int($sub*2);
		$three = int($sub*3);
		$winnums = "10,$one,$two,$three";
		$command1 =~ s/winnum/$winnums/;}
}
	
	

##REMOVE RIBOSOMAL CONTAMINATION##
if($doribo =~ /y/i){
	
	#CREATE ALIGNMENT COMMAND
	$command1 = "muscato --ReadFileName=samplein --ResultsFileName=".$dir."ALIGNED/sample_RRNA_matches.txt --GeneFileName=".$dir2."RRNADB2R.txt.sz --GeneIdFileName=".$dir2."RRNADB2R_ids.txt.sz --Windows=winnum --WindowWidth=winlen --BloomSize=400000000 --NumHash=20 --PMatch=pcid --MinDinuc=5 --MinReadLength=minrl --MaxMatches=1000 --MaxConfirmProcs=5 --MaxReadLength=rlen --MatchMode=best --MMTol=1 --TempDir=tmp/sample";
	$command1 =~ s/winlen/$wl/;
	$command1 =~ s/pcid/$pmatch/;
	$command1 =~ s/minrl/$minrl/;	
	$command1 =~ s/rlen/$rlen/;

	foreach(@FQS){ 
		$sample = $FILE2SAMP{$_};
		$jobname = $dir2.$sample."_ribo.pbs";
		$jobid = $sample."_ribo";
		
		if($doclean =~ /y/i ){$file = $dir."CLEANED/".$sample."_clean.fastq";}
		elsif($dotrim =~ /y/i){$file = $dir."CLEANED/".$sample.".trim";}
		else{$file = $_;}
		
		$mycommand1 = $command1;
		$mycommand1 =~ s/samplein/$file/g;
		$mycommand1 =~ s/sample/$sample/g;
		$mycommand1 =~ s/maxm/$maxm/g;
		
		$mypbs1 = $pbs1;
		$mypbs2 = $pbs2;
		$mypbs3 = $pbs3;
		$mypbs7 = $pbs7;
        $mypbs12 =$pbs12;
 
		$mypbs1 =~ s/cores/6/;
		$mypbs2 =~ s/coremem/1/;
		$mypbs3 =~ s/runtime/40:00:00/;
		$mypbs7 =~ s/input/$jobid/;
		$mypbs12 =~ s/sample/$sample/;
        $mypbs12 =~ s/type/ribo/;

		open(ALIGNJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
		print ALIGNJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$pbs14\n$mypbs12\n$pbs13\n$pbs15\n$pbs16\n$pbs17\n$pbs18\n$pbs19\n\n$mycommand1\n";
		print ALIGNJOB "qsub ".$dir2.$sample."_check.pbs\n";
        close(ALIGNJOB);	
	}
}




#MAKE ALIGNMENT PBS FILES
if($doalign =~ /y/i){
	
	#CREATE ALIGNMENT COMMAND
	$command1 = "muscato --ReadFileName=samplein --ResultsFileName=".$dir."ALIGNED/sample_matches.txt --GeneFileName=".$dir2."BIGDBSEQS.txt.sz --GeneIdFileName=".$dir2."BIGDBSEQS_ids.txt.sz --Windows=winnum --WindowWidth=winlen --BloomSize=400000000 --NumHash=20 --PMatch=pcid --MinDinuc=5 --MinReadLength=minrl --MaxMatches=maxm --MaxConfirmProcs=5 --MaxReadLength=rlen --MatchMode=best --MMTol=1 --TempDir=tmp/sample";

	#MAKE PBS
	$command1 =~ s/winlen/$wl/;
	$command1 =~ s/pcid/$pmatch/;
	$command1 =~ s/minrl/$minrl/;	
	$command1 =~ s/rlen/$rlen/;

	foreach(@FQS){ 
		$sample = $FILE2SAMP{$_};
		$jobname = $dir2.$sample."_align.pbs";
		$jobid = $sample."_align";
		
		if($doribo =~ /y/i){$file = $dir."ALIGNED/sample_RRNA_matches.nonmatch.txt.fastq";}
		elsif($doclean =~ /y/i ){$file = $dir."CLEANED/".$sample."_clean.fastq";}
		elsif($dotrim =~ /y/i){$file = $dir."CLEANED/".$sample.".trim";}
		else{$file = $_;}
		
		$mycommand1 = $command1;		
		$mycommand1 =~ s/samplein/$file/g;
		$mycommand1 =~ s/sample/$sample/g;
		$mycommand1 =~ s/maxm/$maxm/g;
		
		$mypbs1 = $pbs1;
		$mypbs2 = $pbs2;
		$mypbs3 = $pbs3;
		$mypbs7 = $pbs7;
		$mypbs12 =$pbs12;

		$mypbs1 =~ s/cores/6/;
		$mypbs2 =~ s/coremem/1/;
		$mypbs3 =~ s/runtime/80:00:00/;
		$mypbs7 =~ s/input/$jobid/;		
        $mypbs12 =~ s/sample/$sample/;
        $mypbs12 =~ s/type/align/;

		open(ALIGNJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
		print ALIGNJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$pbs14\n$mypbs12\n$pbs13\n$pbs15\n$pbs16\n$pbs17\n$pbs18\n$pbs19\n\n$mycommand1\n";
		print ALIGNJOB "mv ".$dir."ALIGNED/*".$sample."_matches.nonmatch.txt.fastq ".$dir."NOTALIGNED/\n";
		print ALIGNJOB "awk \'{print \$5}\' ".$dir."ALIGNED/".$sample."_matches.txt > ".$dir."ALIGNED/".$sample."_Genes.txt\n";
    	print ALIGNJOB "sed -i \'s/\\\(;\\\|_r\\\)\$//g\' ".$dir."ALIGNED/".$sample."_Genes.txt\n";
    	print ALIGNJOB "sed -i \'s/;/\\n/g\' ".$dir."ALIGNED/".$sample."_Genes.txt\n";
    	print ALIGNJOB "sort -u -o ".$dir."ALIGNED/".$sample."_Genes.txt ".$dir."ALIGNED/".$sample."_Genes.txt\n";
		print ALIGNJOB "qsub ".$dir2.$sample."_check.pbs\n";
    	close(ALIGNJOB);
    }
}	


#GET COUNTS
if($docount =~ /y/i){
	
	$countpl = $dir2."GetCounts.pl";
	open(INCOUNT, $countpl)||die "unable to open $cleanpl:$!\n";
	@incount = <INCOUNT>;
	$incount = join('', @incount);
	
	foreach(@FQS){
		$sample = $FILE2SAMP{$_};
		$file = $dir."ALIGNED/".$sample."_matches.txt";
		$jobname = $dir2.$sample."_count.pbs";
 		$jobid 	= $sample."_count";
 		$command = "perl ".$dir2."Count_".$sample.".pl";
 		
		#Make Count_sample.pl
		$mycount = $incount;
		$mycount =~ s/sample/$sample/g;
		$mycount =~ s/dir1/$dir/g;
		$mycount =~ s/dir2/$dir2/g;
		$outpl = $dir2."Count_".$sample.".pl";  
		open(OUTCOUNT, ">", $outpl)||die "unable to open $outpl:$!\n";
		print OUTCOUNT "$mycount\n";
				
		#Make Sample_count.pbs	
		$mypbs1 = $pbs1;
		$mypbs2 = $pbs2;
		$mypbs3 = $pbs3;
		$mypbs5 = $pbs5; 
		$mypbs6 = $pbs6; 
		$mypbs7 = $pbs7;
        $mypbs12 =$pbs12;

		
		if(-e $file){
			my $fsize = (-s $file) / (1024 * 1024 * 1024);
			$corenum = int($fsize/25*1.75 + 0.5);
		}
		else{$corenum =3;}	
		
 		$mypbs1 =~ s/cores/$corenum:largemem/;
		$mypbs2 =~ s/coremem/25/;
		$mypbs3 =~ s/runtime/20:00:00/;
		$mypbs5 =~ s/flux.*/fluxod/;
		$mypbs6 =~ s/flux.*/fluxod/;		
		$mypbs7 =~ s/input/$jobid/;
        $mypbs12 =~ s/sample/$sample/;
        $mypbs12 =~ s/type/count/;
        
		open(COUNTJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
		print COUNTJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$mypbs5\n$mypbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$mypbs12\n$pbs13\n\n$command\n";
		print COUNTJOB "qsub ".$dir2.$sample."_check.pbs\n";
	}	
}





#GET COUNTS
if($doRDEG =~ /y/i){
	
		#GET CONTRASTs
		$contrast = $dir.'contrast.matrix';
		open(CONTS, $contrast)||die "unable to open $contrast:$!\n";
		$i=0;
		while(<CONTS>){
		        $_ =~ s/[\n\r]+//;
		        if($_ !~ /\w/){next;}
		        @stuff = split("\t", $_);
		        if($i==0){for $j (1..$#stuff){$CONTRASTS{$j}=$stuff[$j];}}
		        else{
		                for $j (1..$#stuff){
		                        if($stuff[$j]<0){$NEGS{$CONTRASTS{$j}}{$stuff[0]}=1;}
		                        if($stuff[$j]>0){$POSI{$CONTRASTS{$j}}{$stuff[0]}=1;}
		                }
		        }
		        $i++;
		}
		foreach my $contrast (keys %NEGS){
		        my @negs;
		        my @posi;
		        foreach my $condition (keys %{$NEGS{$contrast}}){ @n1=split('\n', $FILECOND{$condition}); push(@negs, $condition); push(@negs,@n1);}
		        foreach my $condition (keys %{$POSI{$contrast}}){ @p1=split('\n', $FILECOND{$condition}); push(@posi, $condition); push(@posi,@p1);}

		        print "Your contrast $contrast compares \n";
		        for my $i (0..$#posi){print "$posi[$i]\n";}
		        print "vs\n";
		        for my $i (0..$#negs){print "$negs[$i]\n";}
		}
	


	#MAKE INFO PBS
	$unixlines .= "mkdir ".$dir."STATS\n";
	$command = "cat ".$dir."ALIGNED/*_Genes.txt | sort -u > ".$dir."ALIGNED/AllGenes.txt \nperl -e 'open(A, \"".$dir."ALIGNED/AllGenes.txt\"); while(<A>){/.+/ && \$k{\$&}++}; while(<>){/^[^\\t]+/ && do{print if defined(\$k{\$&})}}\' ".$dir2."BIGDB.txt >> ".$dir."STATS/GeneInfo.txt\ncut -d\$\'\\t\' -f1,3,4,7,9-15 ".$dir."STATS/GeneInfo.txt > ".$dir."STATS/info.matrix\n";
	$mypbs1 = $pbs1;
	$mypbs2 = $pbs2;
	$mypbs3 = $pbs3;
	$mypbs7 = $pbs7;
	$mypbs12 =$pbs12;
	$mypbs12 =~ s/sample/info/;
	$mypbs12 =~ s/type/info/;
	$mypbs1 =~ s/cores/20:largemem/;
	$mypbs2 =~ s/coremem/25/;
	$mypbs3 =~ s/runtime/6:00:00/;
	$jobname = $dir2."Get_info.pbs";
	$jobid = "Get_info";
	$mypbs7 =~ s/input/$jobid/;
	$command =~ s/samplein/$file/g;
	$command =~ s/sample/$sample/g;
	open(INFOJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
	print INFOJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$mypbs12\n$pbs13\n\n$command\n";
	print INFOJOB "echo \"Completed successfully\"\n";
	print INFOJOB "qsub ".$dir2.$sample."_check.pbs\n";
	
	
	#MAKE limma.R PBS
	$command = "R GetDEGs.R\n";
	$mypbs1 = $pbs1;
	$mypbs2 = $pbs2;
	$mypbs3 = $pbs3;
	$mypbs7 = $pbs7;
	$mypbs12 =$pbs12;
	$mypbs12 =~ s/sample/degs/;
	$mypbs12 =~ s/type/degs/;
	$mypbs1 =~ s/cores/2:largemem/;
	$mypbs2 =~ s/coremem/25/;
	$mypbs3 =~ s/runtime/48:00:00/;
	$jobname = $dir2."limmaR.pbs";
	$jobid = "limmaR";
	$mypbs7 =~ s/input/$jobid/;
	open(RJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
	print RJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$mypbs12\n$pbs13\n\n$command\n";
	print RJOB "qsub ".$dir2.$sample."_check.pbs\n";
}




#MAKE BIGBOSS.pbs
print "making BIGBOSS\n";
$jobname = $dir2."BIGBOSS.pbs";
$mypbs1 = $pbs1;
$mypbs2 = $pbs2;
$mypbs3 = $pbs3;
$mypbs7 = $pbs7;
$mypbs12 =$pbs12;
$mypbs12 =~ s/sample/boss/;
$mypbs12 =~ s/type/boss/;
$mypbs1 =~ s/cores/1/;
$mypbs2 =~ s/coremem/1/;
$mypbs3 =~ s/runtime/1:00:00/;
$mypbs7 =~ s/input/BIGBOSS/;
open(BOSSJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
print BOSSJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$mypbs12\n$pbs13\n\n$unixlines\n";

for my $i (0..$#FQS){
	$sample = $FILE2SAMP{$FQS[$i]};
	if($dotrim =~ /y/i){$type ="trim"; }
	elsif($doclean =~ /y/i){$type ="clean";}
	elsif($doribo =~ /y/i){$type ="ribo";}
	elsif($doalign =~ /y/i){$type ="align";}
	else{$type ="info";}
	print BOSSJOB "JOB0=\$(qsub -V ".$dir2.$sample."_".$type.".pbs)\n";
}


#MAKE JOB CHECKLIST FOR CHECKJOB.PL
$j2 = $dir2."joblist.txt";
open(JOBLST, ">", $j2)||die "unable to open $j2:$!\n";
if($dotrim =~ /y/i){print JOBLST "trim\n";}
if($doclean =~ /y/i){print JOBLST "clean\n";}
if($doribo =~ /y/i){print JOBLST "ribo\n";}
if($doalign =~ /y/i){print JOBLST "align\n";}
if($docount =~ /y/i){print JOBLST "count\n";}
if($doRDEG =~ /y/i){print JOBLST "info\ndegs\n";}
print JOBLST "done\n";


#MAKE SAMPLE_CHECK.PBS
print "making check\n";
for my $i (0..$#FQS){
	$sample = $FILE2SAMP{$FQS[$i]};
	$jobname = $dir2.$sample."_check.pbs";
	$input = $sample."_check";
	
	$mypbs1 = $pbs1;
	$mypbs2 = $pbs2;
	$mypbs3 = $pbs3;
	$mypbs7 = $pbs7;
    $mypbs12 =$pbs12;
  
	$mypbs1 =~ s/cores/1/;
	$mypbs2 =~ s/coremem/1/;
	$mypbs3 =~ s/runtime/1:00:00/;
	$mypbs7 =~ s/input/$input/;	
	$mypbs12 =~ s/sample/$sample/;
    $mypbs12 =~ s/type/check/;

	if($dotrim =~ /y/i){$type ="trim";}
	elsif($doclean =~ /y/i){$type ="clean";}
	elsif($doribo =~ /y/i){$type ="ribo";}
	elsif($doalign =~ /y/i){$type ="align";}
	elsif($docount =~ /y/i){$type ="count";}
	else{$type ="degs";}
	
	open(CHECKJOB, ">", $jobname)||die "unable to open $jobname:$!\n";
	print CHECKJOB "$mypbs1\n$mypbs2\n$mypbs3\n$pbs4\n$pbs5\n$pbs6\n$mypbs7\n$pbs8\n$pbs9\n$pbs10\n$pbs11\n$pbs12\n$pbs13\n\n";
	print CHECKJOB "perl ".$dir2."checkjob.pl $sample $type\n"; 	
	print CHECKJOB "bash ".$dir2.$sample."_check.sh\n";
}



