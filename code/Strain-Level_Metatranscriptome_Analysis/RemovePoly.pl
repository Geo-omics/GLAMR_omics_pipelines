use warnings;

$input1 = $ARGV[0];
$minlen = $ARGV[1]; 
$minlen =~ s/\D+//g;
$output1 = $ARGV[2];
while($ARGV[0]!~/\.(fastq|fq)/ || $minlen !~ /^\d+$/ || $output1!~/\.(fastq|fq)/){ 
	print "please run the perl script with input fastq, minimum read length, and output fastq. Can be compressed.\nEx: perl RemovePoly.pl in.fastq.gz 100 out.fastq.gz\n";
}

if($input1 =~ /\.gz$/i){open(INPUT1, "gunzip -c $input1 |") or die "gunzip $input1: $!";}
else{open(INPUT1, $input1)||die;}

if($output1 =~ /\.gz$/i){open(OUTPUT1,'>>:gzip', $output1)||die;}
else{open(OUTPUT1, ">", $output1)||die; }


$on=0;
$badline=0;
$start = localtime;
while( my $line = <INPUT1> . <INPUT1> . <INPUT1> . <INPUT1> ) {

	$totreads++;
   	$line = uc($line);
	@stuff = split("\n", $line);
	if($stuff[1] !~ /^[A-Z]+$/){print "non-nucleotide detected on $on at $stuff[1]\n"; $badline++; next;}
	if($badline > 3){print "File parsing issue, make sure all lines in file match fastq format and pairs are matched\n"; die;}


	#REMOVE AMBIGS
	$stuff[1] =~ s/[^ATGC]/N/g;
	if($stuff[1] =~ /^N|N$/){
		@seq1 = split("", $stuff[1]);
		@sco1 = split("", $stuff[3]);
		while($seq1[0] eq "N"){ shift(@seq1); shift(@sco1); }
		@seq1 = reverse(@seq1);
		while($seq1[0] eq "N"){ shift(@seq1); shift(@sco1); }
		@seq1 = reverse(@seq1);
		$stuff[1] = join("", @seq1);
		$stuff[3] = join("", @sco1);
	}


	#REMOVE LOW ENTROPY READS
	%DI = ();
	@SEQ = split("", $stuff[1]);
	for my $x (0..$#SEQ-1){
	        $y = $x+1;
	        $di = join("", @SEQ[$x..$y]);
	        $DI{$di}=1;
	}
	$kc1 = keys %DI;

	
	#REMOVE POLY AND SHORT REPEATS
	$seq = $stuff[1];
	$sco = $stuff[3];
	while($seq =~ /^.?A{20,}|T{20,}|G{20,}|C{20,}|(.{2,6})\1{5,}/){
        if($seq =~ /^(A{20,}|T{20,}|G{20,}|C{20,})/){ $beg = @+[0]; }
        elsif( $seq =~ /^(.{2,6})\1{5,}/){ $beg = @+[0]; }
        elsif($seq =~ /^.?(A{20,}|T{20,}|G{20,}|C{20,})/){ $beg = @+[0]; }
        elsif( $seq =~ /^.?(.{2,6})\1{5,}/){ $beg = @+[0]; }
        else{last;}
        $beg = @+[0];
        @seq = split("", $seq);
        @seq = @seq[$beg..$#seq];
        @sco = split("", $sco);
        @sco = @sco[$beg..$#sco];
        $seq = join("", @seq);
        $sco = join("", @sco);
	}

	while($seq =~ /A{20,}|T{20,}|G{20,}|C{20,}|(.{2,6})\1{5,}.?$/){
        if($seq =~ /(A{20,}|T{20,}|G{20,}|C{20,})$/){ $beg = @+[0]; }
        elsif( $seq =~ /(.{2,6})\1{5,}$/){ $beg = @+[0]; }
        elsif($seq =~ /(A{20,}|T{20,}|G{20,}|C{20,}).?$/){ $beg = @+[0]; }
        elsif( $seq =~ /(.{2,6})\1{5,}.?$/){ $beg = @+[0]; }
        else{last;}
        $end = @-[0];
        $end--;
        @seq = split("", $seq);
        @seq = @seq[0..$end];
        @sco = split("", $sco);
        @sco = @sco[0..$end];
        $seq = join("", @seq);
        $sco = join("", @sco);
	}
	$stuff[1] = $seq;
	$stuff[3] = $sco;

	#OUTPUT READS
	my $NFs = $stuff[1] =~ tr/Nn/NN/; 
	   if($kc1 < 7){$lowent++; next;} #both low entropy
	elsif($NFs > length($stuff[1])*0.05){$ambig++; next;} #both too many ambigs
	elsif(length($stuff[1]) < $minlen){$short++; next;} #both too short
	else{$good++;}

	$SEQS{$stuff[1]}=$stuff[0];
	$SCOS{$stuff[1]}=$stuff[3];
	$CNT{$stuff[1]}++;


	if($on%1000000==0){$time = localtime; print "on $on time $time tot $totreads bgood $good lowent $lowent ambig $ambig short $short\n";} 
	$on++;
}
$end = localtime;
print "on $on start $start end $end minlen $minlen\n";
print "tot $totreads bothgood $bothgood forgood $forgood revgood $revgood bothlowent $lowent bothambig $ambig bothshort $short mixed $mixed\n";


foreach my $seq1 (keys %SEQS){
		$nm1=$SEQS{$seq1}."_".$CNT{$seq1};
		print OUTPUT1 "$nm1\n$seq1\n+\n$SCOS{$seq1}\n";
}
