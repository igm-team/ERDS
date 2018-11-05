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
open(FAI,"$ref_fai_file") or die "Can't open $ref_fai_file";
while($in=<FAI>){
	$in=~/^(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)/;
	if(exists($chr_hash{$1})){
		$chr_length[$chr_hash{$1}]=$2;
		$chr_start[$chr_hash{$1}]=$3;
		$row_len[$chr_hash{$1}]=$4;
		$rrow_len[$chr_hash{$1}]=$5;
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

	my $chr_rc="$cnv_dir/$normal_chr[$j].rd";
	my $rc_file="$rc_dir/$normal_chr[$j].rc";
	open(RC, "$rc_file") or die "Can't open $rc_file";
	open(OUT, ">$chr_rc") or die "Can't open $chr_rc";

	while($in=<RC>){
		$in=~m/(\d+)\t(\d+)\t(\S+)\t(\d+)\n/;
		$start=$1;
		$end=$2;
		$rc=$3;
		$map=$4;
		$length=$end-$start+1;
		$gc_count=0;
		$n_count=100;
		if($end<$chr_length[$j]){
			$sub_seq=substr($seq,$start-1,$length);
			$gc_count=($sub_seq =~ tr/C|G|c|g//);
			$n_count=($sub_seq =~ tr/N|n//);
			$gc_count=int($gc_count/$length*100+0.5);
			$n_count=int($n_count/$length*100+0.5);		
		}
		if($j<=$#aut_chr and $n_count==0 and $map>=$mq_thh){
			push(@gcb,[$gc_count,$rc]);
		}
		print OUT "$start\t$end\t$rc\t$gc_count\t$n_count\t$map\n";				
	}
	close(RC) or die "Can't close $rc_file file";
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
	printf $fh "%d\n", int($gcb[$k][1]+0.5);
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
