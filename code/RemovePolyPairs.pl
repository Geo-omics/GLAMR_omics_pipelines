#use warnings;

$input1 = $ARGV[0];
$input2 = $ARGV[1];
$minlen = $ARGV[2];
$output1 = $ARGV[3];
$output2 = $ARGV[4];

if($input1 !~ /(\_1|fwd)/i){ print "Please specify forward read file with either _1_ or _fwd in the file name\n"; }
if($input2 !~ /(\_2|rev)/i){ print "Please specify reverse read file with either _2_ or _rev in the file name\n"; }
if($output1 !~ /\w/ || $output2 !~ /\w/ || $output1 eq $output2){ print "Please specify forward and reverse output file names\n"; }

if($input1 =~ /\.gz$/i){
        open(INPUT1, "gunzip -c $input1 |") or die "gunzip $input1: $!";
        open(INPUT2, "gunzip -c $input2 |") or die "gunzip $input2: $!";
}
else{   open(INPUT1, $input1)||die; open(INPUT2, $input2)||die; }
if($output1 =~ /\.gz$/i){
        open(OUTPUT1,'>:gzip', $output1)||die;
        open(OUTPUT2,'>:gzip', $output2)||die;
}
else{   open(OUTPUT1, ">", $output1)||die;
        open(OUTPUT2, ">", $output2)||die; }
$log=$output1;
$log=~s/\_*(fwd|for|1\_).*?$/_RPP.log/;
open(LOG, ">>", $log)||die;

$on=0;
$badline=0;
$start = localtime;
while( my $line = <INPUT1> . <INPUT1> . <INPUT1> . <INPUT1> . <INPUT2> . <INPUT2> . <INPUT2> . <INPUT2> ) {
        if($badline > 3){print "File parsing issue\n"; die;}
        $totreads++;
        $line = uc($line);
        @stuff = split("\n", $line);
        if($stuff[1] !~ /^[A-Z]+$/ || $stuff[5] !~ /^[A-Z]+$/){print "non-nucleotide detected on $on at $stuff[1] \n or \n $stuff[5]\n"; $badline++; next;}
        $stuff[4] =~ /^(\S+)(\s|\_).*?/;
        $read=$1;
        if($stuff[0] !~ /$read/ || $read !~ /\w/){ print "paired name issues $stuff[0] ne $stuff[4] read $read badline $badline\n"; $badline++; next;}
        if(exists($SEQS{$stuff[1]}{$stuff[5]})){$dup++; next;}
        $save1 = $stuff[1]; $save3 = $stuff[3];
        $save5 = $stuff[5]; $save7 = $stuff[7];
        $len1 = length($stuff[1]);
        $len2 = length($stuff[5]);
        if($len1 < $minlen && $len2 < $minlen){$short++; next;} #both too short

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
        $stuff[5] =~ s/[^ATGC]/N/g;
        if($stuff[5] =~ /^N|N$/){
                @seq2 = split("", $stuff[5]);
                @sco2 = split("", $stuff[7]);
                while($seq2[0] eq "N"){ shift(@seq2); shift(@sco2); }
                @seq2 = reverse(@seq2);
                while($seq2[0] eq "N"){ shift(@seq2); shift(@sco2); }
                @seq2 = reverse(@seq2);
                $stuff[5]=join("", @seq2);
                $stuff[7] = join("", @sco2);
        }
        if(exists($SEQS{$stuff[1]}{$stuff[5]})){$dup++; next;}
        my $NFs = $stuff[1] =~ tr/Nn/NN/;
        my $NRs = $stuff[5] =~ tr/Nn/NN/;
        if($NFs > $len1*0.05 && $NRs > $len2*0.05){ $ambig++; next;} #both too many ambigs


        #REMOVE POLY AND SHORT REPEATS
        $seq = $stuff[1];
        $sco = $stuff[3];
        while($seq =~ /^.?(A{30,}|T{30,}|G{30,}|C{30,}|(.{2})\2{14,}|(.{3})\3{9,}|(.{4})\4{6,}|(.{5})\5{5,}|(.{6})\6{4,})/){
                $beg = @+[0];
                @seq = split("", $seq);
                @seq = @seq[$beg..$#seq];
                @sco = split("", $sco);
                @sco = @sco[$beg..$#sco];
                $seq = join("", @seq);
                $sco = join("", @sco);
        }
        while($seq =~ /(A{30,}|T{30,}|G{30,}|C{30,}|(.{2})\2{14,}|(.{3})\3{9,}|(.{4})\4{6,}|(.{5})\5{5,}|(.{6})\6{4,}).?$/){
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

        #REMOVE POLY REVERSE
        $seq = $stuff[5];
        $sco = $stuff[7];
        while($seq =~ /^.?(A{30,}|T{30,}|G{30,}|C{30,}|(.{2})\2{14,}|(.{3})\3{9,}|(.{4})\4{6,}|(.{5})\5{5,}|(.{6})\6{4,})/){
                $beg = @+[0];
                @seq = split("", $seq);
                @seq = @seq[$beg..$#seq];
                @sco = split("", $sco);
                @sco = @sco[$beg..$#sco];
                $seq = join("", @seq);
                $sco = join("", @sco);
        }
        while($seq =~ /(A{30,}|T{30,}|G{30,}|C{30,}|(.{2})\2{14,}|(.{3})\3{9,}|(.{4})\4{6,}|(.{5})\5{5,}|(.{6})\6{4,}).?$/){
                $end = @-[0];
                $end--;
                @seq = split("", $seq);
                @seq = @seq[0..$end];
                @sco = split("", $sco);
                @sco = @sco[0..$end];
                $seq = join("", @seq);
                $sco = join("", @sco);
        }
        $stuff[5] = $seq;
        $stuff[7] = $sco;
        $len1 = length($stuff[1]);
        $len2 = length($stuff[5]);
        if($len1 < $minlen && $len2 < $minlen){$short++; next;} #both too short


        #REMOVE LOW ENTROPY READS
        %DI = ();
        @SEQ = split("", $stuff[1]);
        for my $x (0..$#SEQ-1){
                $y = $x+1;
                $di = join("", @SEQ[$x..$y]);
                $DI{$di}=1;
        }
        $kc1 = keys %DI;
        %DI2 = ();
        @SEQ = split("", $stuff[5]);
        for my $x (0..$#SEQ-1){
                $y = $x+1;
                $di = join("", @SEQ[$x..$y]);
                $DI2{$di}=1;
        }
        $kc2 = keys %DI2;
        if($kc1 < 7 && $kc2 < 7){$lowent++; next;} #both low entropy

        #OUTPUT READS
        if($kc1 < 7 || $NFs > $len1*0.05 || $len1 < $minlen){ #forward read bad
           if($kc2 < 7 || $NRs > $len2*0.05 || $len2 < $minlen){ $mixed++; next;} #both bad
           $revgood++; $stuff[1]=$save1; $stuff[3]=$save3;
        }
        if($kc2 < 7 || $NRs > $len2*0.05 || $len2 < $minlen){ #reverse read bad
           if($kc1 < 7 || $NFs > $len1*0.05 || $len1 < $minlen){ $mixed++; next;} #both bad
           $forgood++; $stuff[5]=$save5; $stuff[7]=$save7;
        }
        else{$bothgood++;}


        $SEQS{$stuff[1]}{$stuff[5]}=$read;
        $SCOS{$stuff[1]}=$stuff[3];
        $SCOS{$stuff[5]}=$stuff[7];


        if($on%1000000==0){$time = localtime;
                print "on $on time $time tot $totreads dup $dup bothgood $bothgood forgood $forgood revgood $revgood bothlowent $lowent bothshort $short mixed $mixed $read\n";
                print LOG "on $on time $time tot $totreads dup $dup bothgood $bothgood forgood $forgood revgood $revgood bothlowent $lowent bothshort $short mixed $mixed $read\n";
        }
        $on++;
}


foreach my $seq1 (sort(keys %SEQS)){
        if($seq1 !~ /^[A-Z]+$/){$blr++; next;}
        foreach my $seq2 (sort(keys %{$SEQS{$seq1}})){
                if($seq2 !~ /^[A-Z]+$/){$brr++; next;}
                $name = $SEQS{$seq1}{$seq2};
                if($name !~ /\w/){$nn++; next; }
                $nm1=$SEQS{$seq1}{$seq2}."_1";
                $nm2=$SEQS{$seq1}{$seq2}."_2";
                $totout++;
                print OUTPUT1 "$nm1\n$seq1\n+\n$SCOS{$seq1}\n";
                print OUTPUT2 "$nm2\n$seq2\n+\n$SCOS{$seq2}\n";
        }
}

$end = localtime;
print LOG "$input1 $output1 on $on start $start end $end minlen $minlen tot $totreads totout $totout ";
print LOG "noname $nn dup $dup bothgood $bothgood forgood $forgood revgood $revgood bothlowent $lowent bothambig $ambig bothshort $short mixed $mixed blr $blr brr $brr\n";
