#!/usr/bin/perl -w
#Author: Mingfu Zhu, mingfu.zhu@duke.edu
use List::Util qw[min max];
use strict;

our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
my $parameter_sample=$ARGV[0];
my $functions=$ARGV[1];
require "$parameter_sample";
require "$functions";
my $type="DUP";
my $chr=$ARGV[2];
my $int_prob=0.99;
$num_state++;

##input and directory setup
##################################################
$cnv_dir="$cnv_dir"."_"."$win_size_dup";
$gcbin_dir="$cnv_dir/gc_bin";
my $lam_all="$cnv_dir/lam.txt";
my $gc_weight_file="$cnv_dir/gc_weight.txt";
my $mq_weight_file="$cnv_dir/mq_weight.txt";
my $phmm_dir="$script_dir/phmm";
my $hmm_dir="$script_dir/hmm";
my $chr_rd="$cnv_dir/$chr.rd";
my $chr_sd="$sd_dir/$chr.sd";
my $chr_dat="$cnv_dir/$chr.dat";
my $chr_vit="$cnv_dir/$chr.vit";
if(!-s $chr_rd){
	print "no chromosome $chr_rd\n"  if ($chr!~/Y/);
	exit;
}

my (@list, $in);
my ($rd, $rd_gc, $gc_count, $gc_weight, $n_count);
my ($map_q, $mq_weight, $weight, $start, $end, $length, $cn, $logp_sum, $ratio);
my @het_count;
my $locus;

#gender setup
my $gender;
my %gender;
my @aut_chr;
my @X_files=<$cnv_dir/*X.rd>;
my @Y_files=<$cnv_dir/*Y.rd>;
if($#X_files==0 and $#Y_files==0){
	$gender="M";
	for(my $j=0; $j<=$#chr; $j++){
		if($chr[$j]=~m/^c*h*r*\d+$/){
			push(@aut_chr,$chr[$j]);
			$gender{$chr[$j]}=2;
		}
		elsif($chr[$j]=~m/^c*h*r*X|Y$/){
			$gender{$chr[$j]}=1;
		}
		else{
			$gender{$chr[$j]}=0;
		}
	}
}
else{
	$gender="F";
	for(my $j=0; $j<=$#chr; $j++){
		if($chr[$j]=~m/^c*h*r*\d+|X$/){
			push(@aut_chr,$chr[$j]);
			$gender{$chr[$j]}=2;
		}
		else{
			$gender{$chr[$j]}=0;
		}
	}
}

my (@lam, @factor);
open(L,$lam_all) or die "Can't open $lam_all";
while($in=<L>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	push(@lam,$list[0]);
	push(@factor,$list[1]);
}
close(L) or die "Can't close $lam_all";

####### weight
my @gc_weight;
my @mq_weight;
open(W,$gc_weight_file) or die "Can't open gc_weight_file";
while($in=<W>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	push(@gc_weight,$list[0]) if ($type eq "DEL");
	push(@gc_weight,$list[1]) if ($type eq "DUP");
}
close(W) or die "Can't close gc_weight_file";

open(W,$mq_weight_file) or die "Can't open mq_weight_file";
while($in=<W>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	push(@mq_weight,$list[0]) if ($type eq "DEL");
	push(@mq_weight,$list[1]) if ($type eq "DUP");
}
close(W) or die "Can't close mq_weight_file";

#end_hash
my %end_hash;
open(SD,$chr_sd) or die "Can't open $chr_sd";
while($in=<SD>){
	if($in=~ m/^\w+:\d+-(\d+)/){
		$end=$1;
		$end_hash{$end}=1;	
	}
}
close(SD) or die "Can't close $chr_sd";

##data
open(RD,$chr_rd) or die "Can't open $chr_rd";
open(DAT,">$chr_dat") or die "Can't open $chr_dat";
while($in=<RD>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	$start=$list[0];
	$end=$list[1];
	$rd=$list[2];
	$gc_count=$list[3];
	$n_count=$list[4];
	$map_q=$list[5];
	$weight=($gc_weight[$gc_count]*$mq_weight[$map_q])**(1/2);
	$rd=$rd/$factor[$gc_count]*$weight;
	$rd_gc=$lam[$gc_count]/($factor[$gc_count]*2)*$weight;		

	my $end_rounded=int(($end-1)/$sd_win_size+1)*$sd_win_size;
	if(exists($end_hash{$end_rounded})){
		$rd=-1;
		$rd_gc=-1;
	}
	elsif($n_count>10){
		$rd=-1;
		$rd_gc=-1;
	}
	elsif($n_count<10){
		$rd=$rd*100/(100-$n_count);
	}
	my $het_count=0;
	printf DAT ("$start\t$end\t%.2f\t%.2f\t$het_count\n",$rd,$rd_gc);
}
close(RD) or die "Can't close $chr_rd";
close(DAT) or die "Can't close $chr_dat";


##HMM file
###########################################################################
my $hmm_file="$cnv_dir/$chr.hmm.txt";
print_hmm();

###Call HMM##
if($gender{$chr}==2){
	system("$phmm_dir/phmm $hmm_file $chr_dat>$chr_vit");
}
else{
	system("$hmm_dir/hmm $hmm_file $chr_dat>$chr_vit");
}
system("rm $hmm_file");


##functions
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

