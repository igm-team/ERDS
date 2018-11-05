#!/usr/bin/perl
#Author: Mingfu Zhu, mingfu.zhu@duke.edu
use List::Util qw[min max];

##input and directory setup
##################################################
our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
my $parameter_sample=$ARGV[0];
my $functions=$ARGV[1];
require "$parameter_sample";
require "$functions";
my $int_prob=0.99;

$cnv_dir="$cnv_dir"."_"."$win_size_dup";
$gcbin_dir="$cnv_dir/gc_bin";
system("mkdir $cnv_dir") if (!-e $cnv_dir);
system("mkdir -p $gcbin_dir") if (!-e $gcbin_dir);
my $lam_all="$cnv_dir/lam.txt";
my $gc_weight_file="$cnv_dir/gc_weight.txt";
my $mq_weight_file="$cnv_dir/mq_weight.txt";

my $win_count=$win_size_dup/$win_size_del;

#Chromosome setup
my @aut_chr;
for(my $j=0; $j<=$#normal_chr; $j++){
	if($normal_chr[$j]=~m/^c*h*r*\d+$/){
		push(@aut_chr,$normal_chr[$j]);
	}
}
my %chr_hash;
for(my $j=0; $j<=$#normal_chr; $j++){
    $chr_hash{$normal_chr[$j]}=$j;
}

my (@chr_length, @chr_start, @row_len, @rrow_len);
my (@list,$in);
my $chr_max;
my $chr_max_length=0;
open(FAI,"$ref_fai_file") or die "Can't open $ref_fai_file";
while($in=<FAI>){
	$in=~/^(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)/;
	if(exists($chr_hash{$1})){
		$chr_length[$chr_hash{$1}]=$2;
		$chr_start[$chr_hash{$1}]=$3;
		$row_len[$chr_hash{$1}]=$4;
		$rrow_len[$chr_hash{$1}]=$5;
		if($chr_max_length<$2){
			$chr_max_length=$2;
			$chr_max=$1;
		}
	}
}
close(FAI) or die "Can't close $ref_fai_file";

##rc
##################################################
my ($seq, $sub_seq, $gc_count, $n_count, @gcb);
my ($rc, $map, $start, $end, $length);

for(my $j=0; $j<=$#normal_chr; $j++){
	open(REF,"$ref_file") or die "Can't open $ref_file";	
	seek(REF, $chr_start[$j],0);
	read(REF, $seq, $chr_length[$j]*$rrow_len[$j]/$row_len[$j]);
	$seq=~ s/\n//g;
	$seq=substr($seq,0,$chr_length[$j]);
	close(REF) or die "Can't close $ref_file";
	my (@start, @end, @rc, @map);

	my $chr_rc="$cnv_dir/$normal_chr[$j].rd";
	my $rc_file="$rc_dir/$normal_chr[$j].rc";
	open(RC, "$rc_file") or die "Can't open $rc_file";
	while($in=<RC>){
		$in=~m/(\d+)\t(\d+)\t(\S+)\t(\d+)\n/;
		push(@start,$1);
		push(@end,$2);
		push(@rc,$3);
		push(@map,$4);
	}
	close(RC) or die "Can't close $rc_file file";
	open(OUT, ">$chr_rc") or die "Can't open $chr_rc";
	for(my $k=0; $k<=$#start; $k+=$win_count){
		$start=$start[$k];
		$end=$end[min($k+$win_count-1,$#start)];
		$length=$end-$start+1;
		$rc=0;
		$map=0;
		my (@rc_tmp,@map_tmp);
		for(my $l=0; $l<=$win_count-1; $l++){
			push(@rc_tmp,$rc[$k+$l]);
			push(@map_tmp,$map[$k+$l]);
		}
		$rc=midhalf(@rc_tmp)/$win_count;
		$map=midhalf(@map_tmp)/$win_count;
		$gc_count=0;
		$n_count=100;
		$sub_seq=substr($seq,$start-1,$length);
		$gc_count=($sub_seq =~ tr/C|G|c|g//);
		$n_count=($sub_seq =~ tr/N|n//);
		$gc_count=int($gc_count/$length*100+0.5);
		$n_count=int($n_count/$length*100+0.5);		

		if($j<=$#aut_chr and $n_count==0 and $map>=$mq_thh){
			push(@gcb,[$gc_count,int($rc+0.5)]);
		}
		printf OUT ("$start\t$end\t%.2f\t$gc_count\t$n_count\t%d\n",$rc,$map);				
	}
	close(OUT) or die "Can't close $chr_rc";
}

#######################GCbin
my $num_gcb=$#gcb;
my @GCB;
my $fh;
for(my $i=0; $i<=100; $i++){
	my $gcbin_file="$gcbin_dir/$i.gcbin";
	push(@GCB,"GCB$i");
	open($GCB[$i],">$gcbin_file");
}

for(my $k=0; $k<=$num_gcb; $k++){
	$fh=$GCB[$gcb[$k][0]];
	printf $fh "%d\n",int($gcb[$k][1]+0.5);
}
for(my $i=0; $i<=100; $i++){
	close($GCB[$i]);
}

#######################gender
my @X_files=<$cnv_dir/*X.rd>;
my @Y_files=<$cnv_dir/*Y.rd>;
my $X_file=$X_files[0];
my $Y_file=$Y_files[0];
if(-e $X_file and -e $Y_file){
	my @X_rc;
	open(X,"$X_file") or die "Can't open $X_file";
	while ($in=<X>) {
		$in=~s/[\r\n]//g;
		@list=split /\s/, $in;
		push(@X_rc,$list[2]);
	}
	close(X) or die "Can't close $X_file";
	
	my @Y_rc;
	open(Y,"$Y_file");
	while ($in=<Y>) {
		$in=~s/[\r\n]//g;
		@list=split /\s/, $in;
		push(@Y_rc,$list[2]);
	}
	close(Y) or die "Can't close $Y_file";

	if(mean(@X_rc)>6*mean(@Y_rc)){
		if(-e $Y_file){
			system("rm $Y_file");
		}
	}
	undef @X_rc;
	undef @Y_rc;
}

##GC estimation
my $phmm_dir="$script_dir/phmm";
open(B,">$lam_all")||die"cannot open $lam_all\n";
for(my $i=0; $i<=100; $i++){
	my $gcbin_file="$gcbin_dir/$i.gcbin";
	my $lam="$gcbin_dir/$i.lam";
	if (-z $gcbin_file){
		print B "0\t1\n";
	}
	else{
		system("$phmm_dir/em $gcbin_file > $lam");
		open(A,"$lam")||die"cannot open $lam\n";
		$in=<A>;
		$in=~s/[\r\n]//g;
		@list=split /\s/, $in;
		printf B "%.2f\t%.2f\n",$list[0],$list[1];
		close(A)||die"cannot close $lam\n";	
	}
}
close(B)||die"cannot close $lam_all\n";

print_weight($gc_weight_file,100);
print_weight($mq_weight_file,$max_mq);

my $chr=$chr_max;
my $hmm_file="$cnv_dir/$chr.hmm.txt";
print_hmm();
system("perl $script_dir/hmm_dup.pl $parameter_sample $functions $chr");

my $gc_col=3;
my $mq_col=5;
my (@count_gc,@count_mq);

for(my $k=0; $k<=$num_state+1; $k++){
	for(my $j=0; $j<=100; $j++){
		$count_gc[$j][$k]=0;
	}
	for(my $j=0; $j<=$max_mq; $j++){
		$count_mq[$j][$k]=0;
	}
}

my $chr_vit="$cnv_dir/$chr.vit";
open(A,"$chr_vit") or die "Can't open $chr_vit";
my $chr_rd="$cnv_dir/$chr.rd";
open(B,"$chr_rd") or die "Can't open $chr_rd";
while ($in=<A>) {
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	$vit=$list[0]-1;
	$in=<B>;
	@list=split /\s/, $in;
	$gc=int($list[$gc_col]);
	$mq=int($list[$mq_col]);
	$count_gc[$gc][$vit]++;
	$count_mq[$mq][$vit]++;
}
close(A) or die "Can't open $chr_vit";
close(B) or die "Can't open $chr_rd";

open(A,">$gc_weight_file")||die"cannot open $gc_weight_file\n";
for(my $j=0; $j<=100; $j++){
	my $row_total=0;
	my $del_total=0;
	my $dup_total=0;
	my $del_weight=1;
	my $dup_weight=1;
	for(my $k=0; $k<=$num_state; $k++){
		$row_total+=$count_gc[$j][$k];
	}
	for(my $k=0; $k<=1; $k++){
		$del_total+=$count_gc[$j][$k];	
	}
	for(my $k=3; $k<=$num_state; $k++){
		$dup_total+=$count_gc[$j][$k];	
	}
	if($row_total>0){
		$del_weight=$exp_pct/($del_total/$row_total) if($del_total>0);
		$dup_weight=$exp_pct/($dup_total/$row_total) if($dup_total>0);
	}
	printf A "%.3f\t%.3f\n",$del_weight,$dup_weight;
}
close(A)||die"cannot close $gc_weight_file\n";

open(A,">$mq_weight_file")||die"cannot open $mq_weight_file\n";
for(my $j=0; $j<=$max_mq; $j++){
	my $row_total=0;
	my $del_total=0;
	my $dup_total=0;
	my $del_weight=1;
	my $dup_weight=1;
	for(my $k=0; $k<=$num_state; $k++){
		$row_total+=$count_mq[$j][$k];
	}
	for(my $k=0; $k<=1; $k++){
		$del_total+=$count_mq[$j][$k];	
	}
	for(my $k=3; $k<=$num_state; $k++){
		$dup_total+=$count_mq[$j][$k];	
	}
	if($row_total>0){
		$del_weight=$exp_pct/($del_total/$row_total) if($del_total>0);
		$dup_weight=$exp_pct/($dup_total/$row_total) if($dup_total>0);
	}
	printf A "%.3f\t%.3f\n",$del_weight,$dup_weight;
}
close(A)||die"cannot close $mq_weight_file\n";

##functions
sub midhalf{
	my @array=@_;
	my $count = $#array+1;
	my $sum;
	my @sarray = sort { $a <=> $b } @array; 
	my $upper=int($count*0.75+0.5);
	my $bottom=int($count*0.25+0.5);
	for(my $i=$bottom; $i<=$upper; $i++){
		$sum+=$sarray[$i];
	}	
	return $sum*$count/($upper-$bottom+1);
}

sub print_weight{
	my $weight_file=$_[0];
	my $total_line=$_[1];
	open(W,">$weight_file") or die "Can't open $weight_file";
	for(my $k=0; $k<=$total_line; $k++){
		printf W ("1\t1\n",1);
	}
	close(W) or die "Can't close $weight_file";
}

sub print_hmm{
	open(HMM,">$hmm_file") or die "Can't open $hmm_file";
	printf HMM ("N=\t%d\nA:\n",$num_state+1);
	my @p;
	my $pa=$int_prob;
	my $pb=(1-$pa)/($num_state);
	for(my $i=0; $i<=$num_state; $i++){
		for(my $k=0; $k<=$num_state; $k++){
			$p[$i][$k]=$pb;
		}
		$p[$i][$i]=$pa;
	}
	for(my $i=0; $i<=$num_state; $i++){
		for(my $k=0; $k<=$num_state; $k++){
			printf HMM ("\%.4f\t",$p[$i][$k]);
		}
		print HMM "\n";
	}

	print HMM "pi:\n";
	my $int_p=1/($num_state+1);
	for(my $k=0; $k<=$num_state; $k++){
		printf HMM ("\%.4f\t",$int_p);
	}
	print HMM "\n";
	close(HMM) or die "Can't close $hmm_file";
}
