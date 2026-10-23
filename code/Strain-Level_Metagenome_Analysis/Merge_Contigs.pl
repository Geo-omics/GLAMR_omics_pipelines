#use warnings;

#BEFORE YOU BEGIN
#USE bbtools dedupe to remove containments and duplicates and orient the sequences into overlaps
   #   comics dedupe in=multi_assembly_contigs.fasta out=dedup.fasta tuc mid=99 e=2 minscaf=200 fo mo=200 c pc=t cc=t fcc=t fmj=t absorbcontainment=t numaffixmaps=3 overwrite=t
#USE bbtools again on the dedupe output to generate the .dot file - use more  permissive edit to allow for indels, will screen later
   #   comics dedupe in=dedup.fasta out=dedup2.fasta tuc e=10 rnc=t ngn=t fo c pc=f absorbcontainment=f mo=200 numaffixmaps=3 dot=graph.dot overwrite=t
#the .dot file outputs based on input contigs and will not keep contig names if you cluster (bug in dedupe)


$sample = $ARGV[0];
$sample =~ s/\s+//g;

$mid=$ARGV[1];
$mid =~ s/\D+//g;
if($mid !~/\d\d/){ print "need to add %identity: perl Merge_ContigsY.pl sample 95 \n"; die;}
$userpid = $mid/100;
$percid = 1-$userpid; #edit distance


$incon  = $sample."_dedup4.fa";
$indot  = $sample."_graph4.dot";
$prod   = $sample."_prod_X.fa";
$output = $sample."_MERGED_CONTIGS.fa";
$log    = $sample."_MERGED_CONTIGS.log";
$starttime = localtime;
open(INCON, $incon)||die "unable to open $incon: $!\n";
open(INDOT, $indot)||die "unable to open $indot: $!\n";
open(OUTPUT, ">", $output)||die;
open(LOG, ">", $log)||die;


#FIRST: GET OVERLAP INFO AND SEQ IDS
$on=0;
$max=0;
$totseqlen=0;
$/=">";
print "INPUT CONTIGS $starttime\n";
while(<INCON>){
        if($_ !~ /\w/){next;}
        $on++;
        $_ = uc($_);
        @stuff = split("\n",$_);
        $header = shift(@stuff);
        $header =~ /(CLUSTER\_(\d+)\,CONTIG\_\d+)/;
        $contig = $1;
        $cluster = $2;

        $seq = join("", @stuff);
        $seq =~ s/[^A-Z]//g;
        $CON2SEQ{$contig}=$seq;
        $ALLSEQS{$seq}=$cluster;
        $NAME{$seq}=$contig;
        $CLUC{$cluster}=1;

        $totseqlen+=length($seq);
        if(length($seq)>$max){$max=length($seq);}
}
$sc = keys %ALLSEQS;
$ckc = keys %CLUC;
$avglen = $totseqlen/$on;
print     "contigs:\t$on\ndistinct:\t$sc\nlongest:\t$max\naverage:\t$avglen\ntotallen:\t$totseqlen\nclusters:\t$ckc\n";
print LOG "contigs:\t$on\ndistinct:\t$sc\nlongest:\t$max\naverage:\t$avglen\ntotallen:\t$totseqlen\nclusters:\t$ckc\n";




#GET OVERLAPS
$/="\n";
$time = localtime;
print "INPUT OVERLAP INFO $time\n";
open(INDOT, $indot)||die;
while(<INDOT>){
        if($_ !~ /\-\>/){next;}
        $_ = uc($_);
        $_ =~ /.*(CLUSTER\_(\d+)\S+)\s.*(CLUSTER\S+)\s.*\=\"(.+)\"/i;
        $contig1 = $1;
        $cluster = $2;
        $contig2 = $3;
        @stuff = split(",",$4);
        $seq1 = $CON2SEQ{$contig1};
        $seq2 = $CON2SEQ{$contig2};
        $s1len = length($seq1)-1;
        $s2len = length($seq2)-1;
        $olen = $stuff[1];
        $edit = $stuff[2];

        #screen out bad overlaps
        $maxbad = $olen*$percid;
        if($edit > $maxbad){next;}
        if($stuff[4] !=0 && $stuff[6] !=0 && $stuff[5] !=0 && $stuff[7] !=0 ){ $NO{$cluster}=1; next; }
        if($stuff[4] >= $stuff[5] || $stuff[6] >= $stuff[7]){ $NO{$cluster}=1; next; } #reverse_dir
        if($stuff[4] == 0 && $stuff[5] == $s1len){ $NO{$cluster}=1; next; } #containment
        if($stuff[6] == 0 && $stuff[7] == $s2len){ $NO{$cluster}=1; next; } #containment

        #orient so it is always head (end of first seq) tail (start of second seq) overlap position
        if($stuff[4] < $stuff[6] && $stuff[5] < $stuff[7]){
                $CLUSTERS{$cluster}{$seq2}{$seq1}=$stuff[6]."|".$stuff[7]."|".$s2len."|".$stuff[4]."|".$stuff[5]."|".$s1len."|".$edit."|".$olen;
                print "cluster $cluster seq2>seq1 $CLUSTERS{$cluster}{$seq2}{$seq1}\n";
                $KEEP{$cluster}{$seq1}++; $KEEP{$cluster}{$seq2}+=0;
                $ML{$seq1}=1; $ML{$seq2}=1;
        }
        elsif($stuff[4] > $stuff[6] && $stuff[5] > $stuff[7]){
                $CLUSTERS{$cluster}{$seq1}{$seq2}=$stuff[4]."|".$stuff[5]."|".$s1len."|".$stuff[6]."|".$stuff[7]."|".$s2len."|".$edit."|".$olen;
                print "cluster $cluster seq1>seq2 $CLUSTERS{$cluster}{$seq1}{$seq2}\n";
                $KEEP{$cluster}{$seq2}++; $KEEP{$cluster}{$seq1}+=0;
                $ML{$seq1}=1; $ML{$seq2}=1;
        }
        else{ $NO{$cluster}=1; next; }
}
$ckc = keys %CLUSTERS;
$nkc = keys %NO;
$skc = keys %ML;
print LOG "not overlapped clusters $nkc ol clusters $ckc olseqs $skc\n";
print     "not overlapped clusters $nkc ol clusters $ckc olseqs $skc\n";
undef(%CON2SEQ);
undef(%NO);


##############################################
#######    OUTPUT SEQS W/O OVERLAPS    #######
##############################################
$cntout=0;
$time = localtime;
$sakc = keys %ALLSEQS;
print "OUTPUT NON-OVERLAPPING CONTIGS\n";
foreach my $seq (sort(keys %ALLSEQS)){
        if(!exists($ML{$seq})){
                $cntout++;
                $name = $NAME{$seq};
                $name=~s/[^\d\_\|]+//g;
                $name=~s/^\_|\_$//g;
                print OUTPUT ">CLUSTER:$name\n$seq\n";
                delete($CON2SEQ{$NAME{$seq}});
                delete($ALLSEQS{$seq});
                delete($NAME{$seq});
        }
}
$akc = keys %ALLSEQS;
undef(%ML);
print LOG "total contigs $sakc, output $cntout unoverlaped contigs, there are $akc contigs remaining to be merged\n";
print     "total contigs $sakc, output $cntout unoverlaped contigs, there are $akc contigs remaining to be merged\n";
##############################################
##############################################
##############################################






##############################################
#########    MERGEOVERLAP CLUSTERS   #########
##############################################
print "looping through clusters\n";
$on=0;
$totclust = keys %CLUSTERS;
$time = localtime;
foreach my $cluster (sort{$CLUSTERS{$a} <=> $CLUSTERS{$b}} keys %CLUSTERS){

        #### ~ SCREEN %KEEP FOR HEADS ~ ####
        $min=100000000; %HEADS=();
        foreach my $seq (sort{$KEEP{$cluster}{$a} <=> $KEEP{$cluster}{$b}} keys %{$KEEP{$cluster}}){
                if($KEEP{$cluster}{$seq} < $min ){ $min = $KEEP{$cluster}{$seq}; }
                if($KEEP{$cluster}{$seq}==$min){push(@HEADS,$seq);}
                else{last;}
        }
        $sthedkc = @HEADS;
        $CLHC{$sthedkc}++;
        print "cluster $cluster heads $sthedkc\n";

        #### ~ EXTEND FROM HEADS ~ ####
        while($HEADS[0]=~/\w/){
                $head=shift(@HEADS);
                $name=$NAME{$head};
                $name=~s/[^\d\_\|ab]+//g;
                $name=~s/^\_|\_$//g;
                @HEAD=split("",$head);
                print "cluster $cluster head $NAME{$head}\n";
                foreach my $tail (keys %{$CLUSTERS{$cluster}{$head}}){
                        $name2=$NAME{$tail};
                        $name2=~s/[^\d\_\|ab]+//g;
                        $name2=~s/^\_|\_$//g;
                        @TAIL=split("",$tail);

                        $regex=$CLUSTERS{$cluster}{$head}{$tail};
                        $regex =~ /^(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)/;
                        $headstol=$1; $headeol=$2;
                        $tailstol=$4; $taileol=$5;

                        if($headstol>0){$header=join("",@HEAD[0..$headstol-1]);}  else{$header='';}
                        if($taileol < $#TAIL){$tailer=join("",@TAIL[$taileol+1..$#TAIL]);}  else{$tailer='';}
                        $mid1=join("",@HEAD[$headstol..$headeol]);
                        $mid2=join("",@TAIL[$tailstol..$taileol]);


                        @NEW=(); @NNS=();
                        if($mid1 eq $mid2){
                                $ns1 = $header.$mid1.$tailer; push(@NEW, $ns1);
                                $nn1 = $name."|".$name2; push(@NNS, $nn1);
                        }
                        else{   $ns1 = $header.$mid1.$tailer; push(@NEW, $ns1);
                                $nn1 = $name."|".$name2."a";  push(@NNS, $nn1);
                                $ns2 = $header.$mid2.$tailer; push(@NEW, $ns2);
                                $nn2 = $name."|".$name2."b";  push(@NNS, $nn2);
                        }

                        print "cluster $cluster nn1 $nn1 nn2 $nn2 regex $CLUSTERS{$cluster}{$head}{$tail}\n";

                        #check if more tails to add to current tail
                        if(exists($CLUSTERS{$cluster}{$tail})){
                                foreach my $t2 (keys %{$CLUSTERS{$cluster}{$tail}}){
                                        for my $x (0..$#NEW){
                                                $nx=$NEW[$x];
                                                $nn=$NNS[$x];
                                                push(@HEADS,$nx);  #new seq is a head
                                                $NAME{$nx}=$nn;
                                                #make new head/tail in %CLUSTERS
                                                #since new seq ends with current tail, set regex to current tail and next tail
                                                $CLUSTERS{$cluster}{$nx}{$t2}=$CLUSTERS{$cluster}{$tail}{$t2};
                                        }
                                }
                        }
                        else{   #nothing more to add to these sequences, they are done
                                for my $x (0..$#NEW){
                                        $nx=$NEW[$x]; $nn=$NNS[$x];
                                        print "done $nn\n";
                                        $DONE{$nx}=$nn;
                                }
                        }
                        delete($CLUSTERS{$cluster}{$head}{$tail}); #get rid of done seq in %CLUSTERS
                }
                delete($CLUSTERS{$cluster}{$head});
        }
}



$dkc = keys %DONE;
foreach my $seq (keys %DONE){
        $name = $DONE{$seq};
        $name=~s/[^\d\_\|]+//g;
        $name=~s/^\_|\_$//g;
        $PROD_FIX{$name}{$seq}=$DONE{$seq};
}

$on=0;
$tbf = keys %PROD_FIX;
foreach my $ns (keys %PROD_FIX){
        $on++;
        $kc = keys %{$PROD_FIX{$ns}};
        print "on $on of $tbf kc $kc\n";
        if($kc < 2){next;}
        %CONF=();
        %COMP=();
        %NLEN=();
        foreach my $seq (keys %{$PROD_FIX{$ns}}){
                my ($conf, $compl, $nlen) = do_prodigal($seq);
                $CONF{$seq}=$conf;
                $COMP{$seq}=$compl;
                $NLEN{$seq}=$nlen;
        }
        $mxconf=0; $mxcomp=0; $mxnlen=0;
        foreach my $seq (sort{$COMP{$b}<=>$COMP{$a} || $NLEN{$b}<=>$NLEN{$a} || $CONF{$b}<=>$CONF{$a}} keys %COMP){
                print "$ns name $DONE{$seq} conf $CONF{$seq} comp $COMP{$seq} nlen $NLEN{$seq}\n";
                if($COMP{$seq}>$mxcomp && $NLEN{$seq}>$mxnlen){ $mxconf=$CONF{$seq}; $mxcomp=$COMP{$seq}; $mxnlen=$NLEN{$seq};}
                if($NLEN{$seq}>=$mxnlen && $COMP{$seq}>=$mxcomp && $CONF{$seq}>=$mxconf){next;}
                else{delete($DONE{$seq});}
        }
}

$dkc_red = keys %DONE;
foreach my $seq (keys %DONE){
        $name = $DONE{$seq};
        @mj = ( $name =~ /[ab]/g );
        $mj = @mj;
        $MULTIMID{$mj}++;
        @ol = ( $name =~ /\|/g );
        $ol=@ol;
        $SOLC{$ol}++;

        print OUTPUT ">CLUSTER:$name\n$seq\n";
        print LOG "type 2\tCLUSTER:$name\n";
}

foreach my $x (sort(keys %CLHC)){print "cluster hc $x\t$CLHC{$x}\n";}
foreach my $x (sort(keys %MULTIMID)){print "MULTIMID $x\t$MULTIMID{$x}\n";}
foreach my $x (sort(keys %SOLC)){print "multiol $x\t$SOLC{$x}\n";}

foreach my $x (sort(keys %CLHC)){       print LOG "cluster hc $x\t$CLHC{$x}\n";}
foreach my $x (sort(keys %MULTIMID)){   print LOG "MULTIMID $x\t$MULTIMID{$x}\n";}
foreach my $x (sort(keys %SOLC)){       print LOG "multiol $x\t$SOLC{$x}\n";}
print "total clusters to be overlapped: $totclust\n";
print "initial donemerge seqs: $dkc\n";
print "total mergealts: $tbf\n";
print "final donemerge seqs: $dkc_red\n";



##### ~~ SUBROUTINES ~~ #####
#############################

sub do_prodigal{
        my $nsX = $_[0];
        open(PROD, ">", $prod)||die;
        print PROD ">nsX\n$nsX\n";
        $compl=0; $conf=0; $nlen =0;
        $Q = qx(comics prodigal -c -q -p meta -f gff -i $prod);
        @Q = split("\n", $Q);
        foreach my $l (@Q){
                if($l =~ /^\#/){next;}
                $l =~ /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+).*partial\=(\d+).*conf\=(\d+)/i;
                $len = $3-$2;
                if($4 !~ /1/){$compl++;}
                $conf+=$5;
                $nlen += $len;
        }
        return($conf, $compl, $nlen);
}





