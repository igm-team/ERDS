#!/usr/bin/perl
#Author: Mingfu Zhu, mingfu.zhu@duke.edu
use List::Util qw[min max];
use strict;

our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
my $parameter_sample=$ARGV[0];
my $functions=$ARGV[1];
require "$parameter_sample";
require "$functions";
my $chr=$ARGV[2];
$num_state++;

##input and directory setup
##################################################
system("cp $cnv_dir/$chr.sd.rd \"$cnv_dir\"\"_\"\"$win_size_dup\"/") if (-e "$cnv_dir/$chr.sd.rd");
$cnv_dir="$cnv_dir"."_"."$win_size_dup";
my $phmm_dir="$script_dir/phmm";
my $hmm_dir="$script_dir/hmm";
my $chr_rd="$cnv_dir/$chr.rd";
my $chr_dat="$cnv_dir/$chr.dat";
my $chr_vit="$cnv_dir/$chr.vit";
my $chr_evt="$cnv_dir/$chr.evt";
my $chr_sd_rd="$cnv_dir/$chr.sd.rd";
if(!-s $chr_vit){
	print "no chromosome $chr_vit\n"  if ($chr!~/Y/);
	exit;
}

my $chr_evt_dup="$cnv_dir/$chr.dup.evt";

my (@list, $in);
my ($start, $end, $length, $cn, $logp_sum, $ratio, $rd, $rd_gc);

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
my $normal_cn=$gender{$chr};

my $lam="$cnv_dir/lam.txt";
my @lam;
open(L,$lam) or die "Can't open $lam";
while($in=<L>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	push(@lam,$list[0]/2);
}
close(L) or die "Can't close $lam";

my (@raw_rd,@raw_rd_gc,@ref_cn);
open(RD,"$chr_rd") or die "Can't open $chr_rd";
while ($in=<RD>) {
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	my $raw_rd=$list[2];
	my $gc=$list[3];
	push(@raw_rd,$raw_rd);
	push(@raw_rd_gc,$lam[$gc]);
	push(@ref_cn,$normal_cn);
}
close(RD) or die "Can't close $chr_rd";

##events
open(EVT,">$chr_evt") or die "Can't open $chr_evt";
my (@start, @end, @rd, @rd_gc, @state, @logp);
open(DAT,"$chr_dat") or die "Can't open $chr_dat";
open(VIT,"$chr_vit") or die "Can't open $chr_vit";
while ($in=<DAT>) {
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	$rd=$list[2];
	$rd_gc=$list[3];
	push(@start,$list[0]);
	push(@end,$list[1]);
	push(@rd,$list[2]);
	push(@rd_gc,$list[3]);

	$in=<VIT>;
	$in=~s/[\r\n]//g;
	$in=~/^(\d+)\t(\S+)/;
	my $q=$1-1;
	push(@logp,$2);
	if($q==$num_state){
		push(@state,$normal_cn);
	}
	elsif($q==$num_state-1){
		$ratio=1;
		if($rd_gc>0){
			$ratio=int($rd/$rd_gc+0.5);
		}
		push(@state,max($num_state,$ratio));
	}
	else{
		push(@state,$q);
	}
}
close(VIT) or die "Can't close $chr_vit";			
close(DAT) or die "Can't close $chr_dat";

my $sd_region=0;
for(my $i=0;$i<=$#state;$i++){
	while($state[$i]==$normal_cn and $i<$#state){
		$i++;
	}
	$start=$start[$i];
	$cn=$state[$i];
	$logp_sum=$logp[$i];
	while($state[$i+1]==$cn and $i<$#state){
		$logp_sum+=$logp[$i+1];
		$i++;
	}
	$end=$end[$i];
	$length=$end-$start+1;
	if($length>0){
		printf EVT ("$chr\t$start\t$end\t$length\t$normal_cn\t$cn\t%.2f\t$sd_region\n",$logp_sum);
	}
}

##SD
my (@start_sd, @end_sd, @ref_cn_sd, @cn_sd, @logp_sd, @type_sd);
open(SDRD,$chr_sd_rd) or die "cannot open $chr_sd_rd";
$sd_region=1;
while($in=<SDRD>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	my ($chr_sd, $start_sd, $end_sd, $ref_cn_sd, $rd_sd, $rd_gc_sd, $cn_sd, $logp_sd, $type_sd)=@list; 
	$rd_sd/=($sd_win_size/$win_size_dup);
	$rd_gc_sd/=($sd_win_size/$win_size_dup);
	for(my $j=int($start_sd/$win_size_dup);$j<=int(($end_sd-1)/$win_size_dup);$j++){
		$raw_rd[$j]=$rd_sd;
		$raw_rd_gc[$j]=$rd_gc_sd;
		$ref_cn[$j]=$ref_cn_sd;
	}
	
	if($logp_sd>1){
		push (@start_sd, $start_sd);
		push (@end_sd, $end_sd);
		push (@ref_cn_sd, $ref_cn_sd);
		push (@cn_sd, $cn_sd);
		push (@logp_sd, $logp_sd);
		push (@type_sd, $type_sd);	
	}
}
close(SDRD) or die "Can't close $chr_sd_rd";

for(my $k=0;$k<=$#start_sd;$k++){
	my $tmp=$k;
	my $logp_sum=$logp_sd[$k];
	while($ref_cn_sd[$k+1]==$ref_cn_sd[$tmp] and $start_sd[$k+1]-$end_sd[$k]<=1 and $k<$#start_sd){	
		if($cn_sd[$k+1]==$cn_sd[$tmp]){
			$k++;
			$logp_sum+=$logp_sd[$k];
		}
		else{
			last;
		}
	}
	if($type_sd[$tmp] ne "NORM"){
		printf EVT ("$chr\t$start_sd[$tmp]\t$end_sd[$k]\t%d\t$ref_cn_sd[$tmp]\t$cn_sd[$tmp]\t$logp_sum\t$sd_region\n",$end_sd[$k]-$start_sd[$tmp]+1);
	}
}
close(EVT) or die "Can't close $chr_evt";
sortfile($chr_evt);
mergefile($chr_evt);

####################merge events
open(EVTDUP,">$chr_evt_dup") or die "Can't open $chr_evt_dup";
my ($left, $right, $prec, @region_cn);
my (@start_evt, @end_evt, @len_evt, @cn_evt, @ref_cn_evt, @logp_sum, @type_evt, @sd_region);
open(EVT,"$chr_evt") or die "Can't open $chr_evt";
while($in=<EVT>){
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	my $start=$list[1];
	my $end=$list[2];
	my $length=$list[3];
	my $ref_cn_evt=$list[4];
	my $cn_evt=$list[5];
	my $logp_sum=$list[6];
	my $sd_region=$list[7];
	push (@start_evt, $start);
	push (@end_evt, $end);
	push (@len_evt, $length);
	push (@cn_evt, $cn_evt);
	push (@ref_cn_evt, $ref_cn_evt);
	push (@logp_sum, $logp_sum);
	push (@sd_region, $sd_region);
	if($cn_evt==$ref_cn_evt or ($cn_evt>$ref_cn_evt and $logp_sum<2)){
		push (@type_evt, "NORM");
	}
	elsif($cn_evt<$ref_cn_evt){
		push (@type_evt, "DEL");
	}
	elsif($cn_evt>$ref_cn_evt){
		push (@type_evt, "DUP");
	}
}
close(EVT) or die "Can't close $chr_evt";
#system("rm $chr_evt");
my $evt_num=$#start_evt;

my @gap;
$gap[0]=1000000;
for(my $j=1; $j<=$evt_num; $j++){
	if($type_evt[$j-1] eq "DEL" or $type_evt[$j] eq "DEL"){
		$gap[$j]=1000000;
	}
	else{
		$gap[$j]=$start_evt[$j]-$end_evt[$j-1]-1;
	}
}
$gap[$evt_num+1]=1000000;

for(my $j=1; $j<=$evt_num; $j++){
	if($gap[$j]<=$gap_events){
		$len_evt[$j-1]=$end_evt[$j]-$start_evt[$j-1]+1;
		$logp_sum[$j-1]+=$logp_sum[$j];
		$sd_region[$j-1]=max($sd_region[$j-1], $sd_region[$j]);
		splice (@start_evt, $j,1);
		splice (@end_evt, $j-1,1);
		splice (@len_evt, $j,1);
		splice (@logp_sum, $j,1);
		splice (@gap, $j,1);
		splice (@sd_region, $j,1);
		$evt_num--;
		$j-- if($j>=1);
	}
}

for(my $j=0; $j<=$evt_num; $j++){
	my $start=$start_evt[$j];
	my $end=$end_evt[$j];
	my $length=$end-$start+1;
	my $logp_sum=$logp_sum[$j];
	my $nor_logp=$logp_sum/$length*$win_size_dup;
	my $prec=0;
	if($length>=$win_size_dup){
		@region_cn=region_cn($start,$end);
		next if($nor_logp<2 or ($region_cn[1]-$region_cn[0]<1 and $region_cn[0]<=2) or ($region_cn[1]-$region_cn[0]<2 and $region_cn[0]>2));
		print EVTDUP "$chr\t$start\t$end\t$length\tDUP\t$logp_sum\t$prec\t$region_cn[0]\t$region_cn[1]\n";
	}
}
close(EVTDUP) or die "Can't close $chr_evt_dup";
sortfile($chr_evt_dup);


#function
sub region_cn{
	my ($start,$end)=@_;
	my $length=$end-$start+1;
	my $exp_rd=0;
	my $region_cn=0;
	my $ref_cn=0;
	my $s=int(($start-1)/$win_size_dup);
	my $e=int(($end-1)/$win_size_dup);	
	if($s==$e){
		$exp_rd=$raw_rd_gc[$s];
		$region_cn=$raw_rd[$s];
		$ref_cn=$ref_cn[$s];
	}
	else{
		$exp_rd+=($win_size_dup*($s+1)-$start+1)/$win_size_dup*$raw_rd_gc[$s];
		$region_cn+=($win_size_dup*($s+1)-$start+1)/$win_size_dup*$raw_rd[$s];
		$ref_cn+=($win_size_dup*($s+1)-$start+1)/$win_size_dup*$ref_cn[$s];

		$exp_rd+=($end-$win_size_dup*$e)/$win_size_dup*$raw_rd_gc[$e];
		$region_cn+=($end-$win_size_dup*$e)/$win_size_dup*$raw_rd[$e];
		$ref_cn+=($end-$win_size_dup*$e)/$win_size_dup*$ref_cn[$e];
		
		for(my $i=$s+1; $i<=$e-1; $i++){
			$exp_rd+=$raw_rd_gc[$i];
			$region_cn+=$raw_rd[$i];
			$ref_cn+=$ref_cn[$i];
		}
		$ref_cn=int($ref_cn*$win_size_dup/$length+0.5);
	}
	return ($ref_cn,0) if($region_cn==0 and $exp_rd==0);
	return ($ref_cn,$normal_cn) if($region_cn>0 and $exp_rd==0);
	$region_cn=int($region_cn/$exp_rd+0.5);
	return ($ref_cn, $region_cn);
}
