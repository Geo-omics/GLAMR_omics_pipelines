#use warnings;
#THIS SCRIPT IS FOR USE WITH DIAMOND ALIGNED NGS SEQUENCING READS DESCRIBED HERE:


#FILES LIST
$indiam = $ARGV[0];
$samp = $ARGV[1];


$infn    = '/geomicro/data22/teals_pipeline/Function_Names.txt';
$intax   = '/geomicro/data22/teals_pipeline/TAXONOMY_DB_2020.txt';
$ininfo  = $samp.'_URDB.txt.gz';

open(INFN, $infn)||die;
open(INDI, $indiam)||die;
open(INTAX, $intax)||die;
if( -e $ininfo && -s $ininfo){ open(ININFO, "gunzip -c $ininfo |")||die; }
else{   print "First Time, looping through the URDB\n";
        $ininfo = '/geomicro/data22/teals_pipeline/URDB_MAY_2019.txt.gz';
        $outinfo = $samp."_URDB.txt.gz";
        open(ININFO, "gunzip -c $ininfo |") or die "problem open gunzip $ininfo: $!";
        open(OUTINFO, '>:gzip', $outinfo)||die;
}

#OUTPUT FILES
$outlog  = $samp.".log";
$outcyto = $samp.".cyto";
$outcent = $samp."_info.matrix";
$outfunc = $samp."_funcs.txt";
open(LOG, ">>", $outlog)||die;
open(OUTCENT,">",$outcent)||die;
open(OUTCYTO,">",$outcyto)||die;
open(OUTFUNC,">",$outfunc)||die;


#INPUT TAXONOMY
print "INPUT TAXONOMY\n";
while(<INTAX>){
        $_=uc($_);
        if($_ !~ /\w/){next;}
        @stuff = split("\t", $_,-1);
        $stuff[$#stuff]=~ s/[\n\r\;]+//g;
        $tid = shift(@stuff);
        $PHY{$tid}=join(";",@stuff);
	if($tid=~/00000$/){ print "tid $tid lin $PHY{$tid}\n";}
}

print "INPUT FUNCTION NAMES\n";
while(<INFN>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $id, my $name)=split("\t",$_,-1);
        $id=~s/\s//g;
        $FUNC2NAME{$id}=$name;
}


#DIAMOND MATCHES
print "INPUT DIAMOND MATCHES\n";
while(<INDI>){
    	$_=uc($_);
    	$_=~s/[\r\n\>\@]//g;
    	@stuff=split("\t",$_, -1);
	$read = $stuff[0];
    	$hit = $stuff[2];
    	$hit =~ s/\_R//;
    	$HIT_LEN{$hit}=$stuff[3];
    	$sco = $stuff[9]*$stuff[11];
	$READ_HITS{$read}{$hit}=$sco;
	$GENE_HITS{$hit}{$read}=$sco;
	if($on%1000000==0){$gc = keys %GENE_HITS; $rc = keys %READ_HITS; print "on $on INDIAM rc $rc gc $gc hit $hit read $read\n";}
        $on++;
}
$gc = keys %GENE_HITS; $rc = keys %READ_HITS; 
print LOG "on $on INDIAM rc $rc gc $gc hit $hit read $read\n";

#GET UNIQUE READ COUNTS PER HIT
foreach my $read (keys %READ_HITS){ if(keys %{$READ_HITS{$read}}==1){ foreach my $hit (keys %{$READ_HITS{$read}}){ $HIT_URC{$hit}++; } } }
$totreads= keys %READ_HITS;
$ugc = keys %HIT_URC;
print LOG "total reads $totreads genes w/uniq reads $ugc\n";

#INPUT GENE INFORMATION
while(<ININFO>){
        if($_ !~ /\w/){next;}
	$_ =~ s/[\n\r]+//;
	$_ = uc($_);
        @stuff = split("\t", $_, -1);
	if(exists($HIT_LEN{$stuff[0]})){
		$gene = $stuff[0];
		$out = join("\t", @stuff[0,1,3,9..13,18]);
		$OINF{$out}=1;

	    if($ininfo=~/$samp/){
		#GET TIDS PER GENE and GOOD TIDS
		if( $stuff[1] =~ /HUMAN\_/ ){$HIT_TID{$gene}{9606}=1;}
		else{   @TIDS = split(";", $stuff[2]);
                        foreach my $tid (@TIDS){
                                if($PHY{$tid} =~ /PRIMATE/){
                                        delete($HIT_TID{$gene});
                                        $HIT_TID{$gene}{9606}=1;
                                        if($HIT_URC{$gene}>0){$TID_URC{9606}+=$HIT_URC{$gene};}
                                        else{$TID_URC{9606}+=0;}
                                        last;
                                }
                                if($PHY{$tid} !~ /\w/){next;}
                                $HIT_TID{$gene}{$tid}=1;
                                if($HIT_URC{$gene}>0){$TID_URC{$tid}+=$HIT_URC{$gene};}
                                else{$TID_URC{$tid}+=0;}
                        }
                }

                #GET FUNCS
                @FUNC = split(";", $stuff[3]);
                foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^K/){     $FUNCS{K}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[4]);
                foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^P/){     $FUNCS{P}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[5]);
                foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^[KC]/){  $FUNCS{C}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[6]);
                foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^G/){     $FUNCS{G}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[7]);
                foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^I/){     $FUNCS{I}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[8]);
                foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /RXN/){    $FUNCS{R}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}

	    }
	    else{
		$out = join("\t", @stuff[0,1,3,9..13,18]);
		$OINF{$out}=1;
                #GET TIDS PER GENE and GOOD TIDS
		if( $stuff[1] =~ /HUMAN\_/ ){$HIT_TID{$gene}{9606}=1;}
                else{ 	@TIDS = split(";", $stuff[3]);
                	foreach my $tid (@TIDS){ 
				if($PHY{$tid} =~ /PRIMATE/){
					delete($HIT_TID{$gene}); 
					$HIT_TID{$gene}{9606}=1;
					if($HIT_URC{$gene}>0){$TID_URC{9606}+=$HIT_URC{$gene};}
					else{$TID_URC{9606}+=0;} 
					last;
				} 
				if($PHY{$tid} !~ /\w/){next;}
				$HIT_TID{$gene}{$tid}=1;
				if($HIT_URC{$gene}>0){$TID_URC{$tid}+=$HIT_URC{$gene};}
				else{$TID_URC{$tid}+=0;}
			}
		}

                #GET FUNCS
                @FUNC = split(";", $stuff[9]);  
		foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^K/){	$FUNCS{K}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[10]); 
		foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^P/){	$FUNCS{P}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[11]); 
		foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^[KC]/){	$FUNCS{C}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[12]); 
		foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^G/){	$FUNCS{G}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[13]); 
		foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /^I/){	$FUNCS{I}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
                @FUNC = split(";", $stuff[18]); 
		foreach my $id (@FUNC){ if($id =~ /[\w\d]/ && $id =~ /RXN/){	$FUNCS{R}{$gene}{$id}++; if(!exists($FUNC2NAME{$id}) || $FUNC2NAME{$id} !~ /\w/){$FUNC2NAME{$id}=$stuff[1];}}}
	    }

                #GET NAME
                $G2NAME{$gene}=$stuff[1];
        }
}
$tgc = keys %GENE_HITS;
$gkc = keys %HIT_TID;
$tkc = keys %TID_URC;
print LOG "totgene hits $tgc gene hits w/tids $gkc tot tids $tkc\n"; 
foreach my $out (sort(keys %OINF)){print OUTINFO "$out\n";}
undef(%OINF);



print "PROCESS GENE OUTPUT\n";
print OUTCENT "URDB_gene\tname\tgene_len\tsum_score\ttot_rpm\tuniq_rpm\trpkm\tkegg\tpfam\tcog\tgo\tinterpro\tmetacyc\tlca\n";
foreach my $gene (sort(keys %HIT_TID)){

	$name = $G2NAME{$gene};
	$tr = keys %{$GENE_HITS{$gene}};
	$trpm = $tr*1000000/$totreads; $trpm =~ s/(\.\d\d).*/$1/;
	$rpkm = $trpm*1000/$HIT_LEN{$gene}; $rpkm =~ s/(\.\d\d).*/$1/;
	$urpm = $HIT_URC{$gene}/$totreads*1000000; $urpm =~ s/(\.\d\d).*/$1/;
	$scor = 0; foreach my $read (keys %{$GENE_HITS{$gene}}){ $scor+=$GENE_HITS{$gene}{$read}; } $scor =~ s/(\.\d\d).*/$1/;
	$glen = $HIT_LEN{$gene};

	#remove less likely tids
	$max=0;
	@LINS=();
	foreach my $tid (sort{ $TID_URC{$b} <=> $TID_URC{$a} } keys %{$HIT_TID{$gene}}){
		if($TID_URC{$tid}>$max){$max=$TID_URC{$tid};}
		if($TID_URC{$tid}<$max*0.9 && $PHY{$tid}!~/^MONA/){next;}
		push(@LINS,$PHY{$tid});
	}

	#make LCA and CYTO
	my %seen = (); @LINS = grep { ! $seen{ $_ }++ } @LINS;
	$lca = MakeLCA(@LINS);
	$lca =~ /^(.)/; $king=$1;
	@TAXON = split(";",$lca);
	for($i=0; $i<=$#TAXON; $i++){
		$type = $king."_".$i;
		if($i == 0){ $phyla = "$type\tROOT\t$TAXON[$i]"; }
		else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]"; }
		$CYTO_GC{$phyla}++;
		$CYTO_RPKM{$phyla}+=$rpkm;
		foreach my $read (keys %{$GENE_HITS{$gene}}){
			$CYTO_TRPM{$phyla}{$read}=1;
			if(keys %{$READ_HITS{$read}}==1){$CYTO_URPM{$phyla}{$read}=1;}
		}
	}

	$kid=''; $max=0; foreach my $id (sort{$FUNCS{K}{$gene}{$b} <=> $FUNCS{K}{$gene}{$a} } keys %{$FUNCS{K}{$gene}}){ 
		$cnt = $FUNCS{K}{$gene}{$id}; if($cnt>$max){$max=$cnt;} 
		if($cnt >= $max*0.9){ $kid.=$id.";"; $FUNC_SCOR{$id}{$gene}=$scor; 
			foreach my $lin (@LINS){ $FUNC_LINS{$id}{$lin}+=$trpm;}
			foreach my $read (keys %{$GENE_HITS{$gene}}){ $FUNCS_TRC{$id}{$read}=1; } 
		} 
	}
        $pid=''; $max=0; foreach my $id (sort{$FUNCS{P}{$gene}{$b} <=> $FUNCS{P}{$gene}{$a} } keys %{$FUNCS{P}{$gene}}){
                $cnt = $FUNCS{P}{$gene}{$id}; if($cnt>$max){$max=$cnt;} if($cnt >= $max*0.9){ $pid.=$id.";"; $FUNC_SCOR{$id}{$gene}=$scor; foreach my $lin (@LINS){ $FUNC_LINS{$id}{$lin}+=$trpm;} foreach my $read (keys %{$GENE_HITS{$gene}}){ $FUNCS_TRC{$id}{$read}=1; } } }
        $cid=''; $max=0; foreach my $id (sort{$FUNCS{C}{$gene}{$b} <=> $FUNCS{C}{$gene}{$a} } keys %{$FUNCS{C}{$gene}}){
                $cnt = $FUNCS{C}{$gene}{$id}; if($cnt>$max){$max=$cnt;} if($cnt >= $max*0.9){ $cid.=$id.";"; $FUNC_SCOR{$id}{$gene}=$scor; foreach my $lin (@LINS){ $FUNC_LINS{$id}{$lin}+=$trpm;} foreach my $read (keys %{$GENE_HITS{$gene}}){ $FUNCS_TRC{$id}{$read}=1; }} }
        $gid=''; $max=0; foreach my $id (sort{$FUNCS{G}{$gene}{$b} <=> $FUNCS{G}{$gene}{$a} } keys %{$FUNCS{G}{$gene}}){
                $cnt = $FUNCS{G}{$gene}{$id}; if($cnt>$max){$max=$cnt;} if($cnt >= $max*0.9){ $gid.=$id.";"; $FUNC_SCOR{$id}{$gene}=$scor; foreach my $lin (@LINS){ $FUNC_LINS{$id}{$lin}+=$trpm;} foreach my $read (keys %{$GENE_HITS{$gene}}){ $FUNCS_TRC{$id}{$read}=1; }} }
        $iid=''; $max=0; foreach my $id (sort{$FUNCS{I}{$gene}{$b} <=> $FUNCS{I}{$gene}{$a} } keys %{$FUNCS{I}{$gene}}){
                $cnt = $FUNCS{I}{$gene}{$id}; if($cnt>$max){$max=$cnt;} if($cnt >= $max*0.9){ $iid.=$id.";"; $FUNC_SCOR{$id}{$gene}=$scor; foreach my $lin (@LINS){ $FUNC_LINS{$id}{$lin}+=$trpm;} foreach my $read (keys %{$GENE_HITS{$gene}}){ $FUNCS_TRC{$id}{$read}=1; }} }
        $rid=''; $max=0; foreach my $id (sort{$FUNCS{R}{$gene}{$b} <=> $FUNCS{R}{$gene}{$a} } keys %{$FUNCS{R}{$gene}}){
                $cnt = $FUNCS{R}{$gene}{$id}; if($cnt>$max){$max=$cnt;} 
		if($cnt >= $max*0.9){ $rid.=$id.";"; $FUNC_SCOR{$id}{$gene}=$scor; 
			foreach my $lin (@LINS){ $FUNC_LINS{$id}{$lin}+=$trpm;} 
			foreach my $read (keys %{$GENE_HITS{$gene}}){ $FUNCS_TRC{$id}{$read}=1; }
		}}

	$nums = $glen."\t".$scor."\t".$trpm."\t".$urpm."\t".$rpkm;
	$funcs = $kid."\t".$pid."\t".$cid."\t".$gid."\t".$iid."\t".$rid;
	$funcs =~ s/\;\t/\t/g;
	$funcs =~ s/\;$//;
	print OUTCENT "$gene\t$name\t$nums\t$funcs\t$lca\n";
	
	delete($GENE_HITS{$gene});
        delete($G2NAME{$gene});
	delete($HIT_LEN{$gene});
	delete($HIT_TID{$gene});
	delete($HIT_URC{$gene});
}


print "PROCESS CYTOSCAPE PHYLOGENY\n";
print OUTCYTO "type\tSource\tTarget\tgene_count\ttotal_rpm\tuniq_rpm\trpkm\n";
foreach my $phyla (sort{$CYTO_RPKM{$b} <=> $CYTO_RPKM{$a}} keys %CYTO_RPKM){
	$tc = keys %{$CYTO_TRPM{$phyla}};
	$trpm = $tc*1000000/$totreads;
	$trpm =~ s/(\.\d\d).*/$1/;
	$uc = keys %{$CYTO_URPM{$phyla}};
	$urpm = $uc*1000000/$totreads;
	$urpm =~ s/(\.\d\d).*/$1/;
	$CYTO_RPKM{$phyla}=~ s/(\.\d\d).*/$1/;
	print OUTCYTO "$phyla\t$CYTO_GC{$phyla}\t$trpm\t$urpm\t$CYTO_RPKM{$phyla}\n";
}
undef(%CYTO_GC);
undef(%CYTO_TRPM);
undef(%CYTO_URPM);
undef(%CYTO_RPKM);





print "PROCESS FUNCTIONS\n";
print OUTFUNC "function\tname\tgenes\ttotal_rpm\tuniq_rpm\tsum_score\ttop_rpm\ttop_lin\tlca\n";
foreach my $id (sort(keys %FUNC_SCOR)){

	$name = $FUNC2NAME{$id};

	$gc = keys %{$FUNC_SCOR{$id}};
	$rc = keys %{$FUNCS_TRC{$id}};
	$trpm = $rc*1000000/$totreads; $trpm =~ s/(\.\d\d).*/$1/;

	$scor = 0; foreach my $gene (keys %{$FUNC_SCOR{$id}}){ $scor+=$FUNC_SCOR{$id}{$gene}; } $scor =~ s/(\.\d\d).*/$1/;
	$uc = 0; foreach my $read (keys %{$FUNCS_TRC{$id}}){ if(keys %{$READ_HITS{$read}}==1){ $uc++;}}
	$urpm = $uc*1000000/$totreads; $urpm =~ s/(\.\d\d).*/$1/;	

	@LINS=(); $toplin=''; $toprpm=0;
	foreach my $lin (sort{$FUNC_LINS{$id}{$b} <=> $FUNC_LINS{$id}{$a}} keys %{$FUNC_LINS{$id}}){ 
		push(@LINS, $lin); if($toplin eq ''){$toplin=$lin; $toprpm=$FUNC_LINS{$id}{$lin}; $toprpm =~ s/(\.\d\d).*/$1/;}
	}

	$lca = MakeLCA(@LINS); 

	print OUTFUNC "$id\t$name\t$gc\t$trpm\t$urpm\t$scor\t$toprpm\t$toplin\t$lca\n";
}


sub MakeLCA{
        my @TAXON;
        %seen=();
        @PHYL = @_;
        @PHYL = grep { !$seen{$_}++ } @PHYL;
        if(grep {/^MONA/} @PHYL){
                @TMP=();
                foreach my $lin (@PHYL){
                        if($lin =~ /^MONA/){ push(@TMP,$lin); }
                }
                 @PHYL=@TMP;
        }
        $len1 = @PHYL;
        if($len1 == 1){$LCA = $PHYL[0]; }
        elsif($len1 >1){
                $first = $PHYL[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @PHYL);
                        $len2 = @matched;
                        if($len2 == $len1){push(@TAXON, $alevel);}
                        else{last;}
                }
                $len3 = @TAXON;
                if($len3 > 1){$LCA = join(";", @TAXON);}
                elsif($len3==1){$LCA = $TAXON[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }
        return($LCA);
}


