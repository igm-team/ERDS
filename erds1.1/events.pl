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
my $dist_thh=($ins_pem-2*$read_length)/2;
my $inser_half=int($ins_pem/2+0.5);

##input and directory setup
##################################################
my $phmm_dir="$script_dir/phmm";
my $hmm_dir="$script_dir/hmm";
my $chr_pem_bam="$rc_dir/$chr.pem.bam";
my $chr_rd="$cnv_dir/$chr.rd";
my $chr_dat="$cnv_dir/$chr.dat";
my $chr_vit="$cnv_dir/$chr.vit";
my $chr_evt="$cnv_dir/$chr.evt";
my $chr_sd_rd="$cnv_dir/$chr.sd.rd";
if(!-s $chr_vit){
	print "no chromosome $chr_vit\n" if ($chr!~/Y/);
	exit;
}

my $chr_evt_del="$cnv_dir/$chr.del.evt";
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
my (@start, @end, @rd, @rd_gc, @hets, @state, @logp);
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
	push(@hets,$list[4]);

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
	$rd_sd/=($sd_win_size/$win_size_del);
	$rd_gc_sd/=($sd_win_size/$win_size_del);
	for(my $j=int($start_sd/$win_size_del);$j<=int(($end_sd-1)/$win_size_del);$j++){
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

####################merge events and del check
open(EVTDEL,">$chr_evt_del") or die "Can't open $chr_evt_del";
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
	if($cn_evt<$ref_cn_evt){
		push (@start_evt, $start);
		push (@end_evt, $end);
		push (@len_evt, $length);
		push (@cn_evt, $cn_evt);
		push (@ref_cn_evt, $ref_cn_evt);
		push (@logp_sum, $logp_sum);
		push (@type_evt, "DEL");
		push (@sd_region, $sd_region);
	}
}
close(EVT) or die "Can't close $chr_evt";
system("rm $chr_evt");
my $evt_num=$#start_evt;

my @gap;
$gap[0]=1000000;
for(my $j=1; $j<=$evt_num; $j++){
	if($type_evt[$j-1] eq $type_evt[$j]){
		$gap[$j]=$start_evt[$j]-$end_evt[$j-1]-1;
	}
	else{
		$gap[$j]=1000000;
	}
}
$gap[$evt_num+1]=1000000;

for(my $j=1; $j<=$evt_num-1; $j++){
	if($gap[$j]<=$gap_events){
		$len_evt[$j-1]=$end_evt[$j]-$start_evt[$j-1]+1;
		$logp_sum[$j-1]+=$logp_sum[$j];
		$sd_region[$j-1]=max($sd_region[$j-1], $sd_region[$j]);
		splice (@start_evt, $j,1);
		splice (@end_evt, $j-1,1);
		splice (@len_evt, $j,1);
		splice (@logp_sum, $j,1);
		splice (@type_evt, $j,1);
		splice (@gap, $j,1);
		splice (@sd_region, $j,1);
		$evt_num--;
		$j-- if($j>=1);
	}
}

for(my $j=0; $j<=$evt_num; $j++){
	next if($type_evt[$j] eq "NORM");
	my $start=$start_evt[$j];
	my $end=$end_evt[$j];
	my $length=$end-$start+1;
	my $logp_sum=$logp_sum[$j];
	if($type_evt[$j] eq "DEL"){
		my $temp_pem_pairs=$pem_pairs;
		if($length<=$small_pem){
			$temp_pem_pairs++;
		}
		my @pem=pem_check($start,$end,$temp_pem_pairs,0);
		if($pem[0]>0){
			$left=$pem[1];
			$right=$pem[2];
			$prec=$pem[3];
			$length=$right-$left+1;
			@region_cn=region_cn($left,$right);
			next if($region_cn[0]<=$region_cn[1]);
			print EVTDEL "$chr\t$left\t$right\t$length\tDEL\t$logp_sum\t$prec\t$region_cn[0]\t$region_cn[1]\n";
		}
		elsif($length>=$large_pem){
			@region_cn=region_cn($start,$end);
			next if($region_cn[0]<=$region_cn[1]);
			$prec=0;
			print EVTDEL "$chr\t$start\t$end\t$length\tDEL\t$logp_sum\t$prec\t$region_cn[0]\t$region_cn[1]\n";
		}
	}
}

##low rd scan
if ($gender{$chr}==1){
	$low_scan_ratio--;
}
for(my $i=0; $i<=$#start; $i++){
	if($rd_gc[$i]>=0 and $raw_rd[$i]<=$low_scan_ratio*$raw_rd_gc[$i] and $state[$i]==$normal_cn and $hets[$i]<=1){
		$start=$start[$i];
		while(($rd_gc[$i+1]>=0 and $raw_rd[$i+1]<=$low_scan_ratio*$raw_rd_gc[$i+1] and $state[$i+1]==$normal_cn and $hets[$i+1]<=1) or ($rd_gc[$i+2]>=0 and $raw_rd[$i+2]<=$low_scan_ratio*$raw_rd_gc[$i+2] and $state[$i+2]==$normal_cn and $hets[$i+2]<=1) and $i<$#start-1){
			$i++;
		}
		$end=$end[$i];
		#require length>=2*$win_size_del
		next if ($end-$start+1<2*$win_size_del);
		my @pem=pem_check($start,$end,$low_pem_pairs,1);
		if($pem[0]>0){
			$left=$pem[1];
			$right=$pem[2];
			$prec=$pem[3];
			$length=$right-$left+1;
			@region_cn=region_cn($left,$right);
			next if($region_cn[0]<=$region_cn[1]);
			print EVTDEL "$chr\t$left\t$right\t$length\tDEL\tNA\t$prec\t$region_cn[0]\t$region_cn[1]\n";
		}
	}
}
close(EVTDEL) or die "Can't close $chr_evt_del";
sortfile($chr_evt_del);

#function
sub pem_check{
	my (@left, @right, @rate, $left, $right, @left_clip, @right_clip, $left_clip, $right_clip);
	my ($start,$end,$count_thh,$MQ)=@_;
	my $length=$end-$start+1;
	if($length>1000){
		$rcp_rate=0.4;
	}
	if($length>5000){
		$rcp_rate=0.5;
	}
	if($length>10000){
		$rcp_rate=0.8;
	}
	my $flank_rate=1/$rcp_rate-1;

	my $start_flank=max(1,$start-int($length*$flank_rate));
	my $end_flank=$end+int($length*$flank_rate);
	my $region="$chr:$start_flank-$end_flank";
	my $line_count;
	open(ST, "$samtools view $chr_pem_bam $region | ")||die"$!\n";
	while(my $line=<ST>) {
		$line_count++;
		if($line_count>=$count_thh*2){
			last;
		}
	}
	close(ST);
	if($line_count<$count_thh*2){
		return 0;
	}
	
	my $ins_cut=$ins_pem+max($min_det_size,$length*$rcp_rate);
	my $middle=int(($start+$end)/2+0.5);
	my (%ID,%paired_ID);
	my $idx=0;
	$prec=0;

	my (@new_left, @new_right);
	open(ST, "$samtools view $chr_pem_bam $region | ")||die"$!\n";
	while(my $line=<ST>) {
		$line=~s/[\r\n]//g;
		@list=split /\s/, $line;
		my $flag=$list[1];
		my $pos=$list[3];
		my $mq=$list[4];
		my $CIAGR=$list[11];
		if($CIAGR=~/MD:Z:(\S+)/){
			$CIAGR=$1;
		}		
		my $pos2=$list[7];
		my $ins_read=abs($list[8]);
		my $read_len=$read_length;
		my $ID=$list[0];
		if($ID=~/(\S+)\/\d+/){
			$ID=$1;
		}
		next if($mq<$MQ);
		my $M=0;
		while($CIAGR=~/(\d+)[MI]/g){
			$M+=$1;
		}
		if($pos<=$middle and $CIAGR=~/(\d+)S$/ and $1>=$clip_len){
			push(@left_clip,$pos+$M);
		}
		if($pos>=$middle and $CIAGR=~/^(\d+)S/ and $1>=$clip_len){
			push(@right_clip,$pos-1);
		}
		next if(length($CIAGR)>$CIAGR_length or $M<$aligned_length_thh);
		my @b=bin($flag);
		#opposite direction and point to the inner
		if($ins_read>$ins_cut and $ins_read<2*$large_pem and $b[4]!=$b[5]){
			if(!exists($ID{$ID})){
				$ID{$ID}=1;
				if($pos2>$pos){
					$left=$pos+$M-$read_len+$inser_half;
					$right=$pos2+$read_len-$inser_half;
				}
				else{
					$left=$pos2+$inser_half;
					$right=$pos+$M-$inser_half;
				}				

				my $overlap=min($end,$right)-max($start,$left)+1;
				if($overlap/$length>=$rcp_rate and $overlap/($right-$left+1)>=$rcp_rate){
					push(@left, $left);
					push(@right, $right);
					$paired_ID{$ID}=$idx;
					$idx++;
				}
			}
			elsif(exists($paired_ID{$ID})){
				push(@new_left,$left[$paired_ID{$ID}]);
				push(@new_right,$right[$paired_ID{$ID}]);				
			}
		}
	}
	close(ST);
	$left=0;
	
	if($#new_left>=$count_thh-1){
		my ($left_prec, $right_prec);
		($left,$right)=cluster($start,$end,$count_thh,$dist_thh,@new_left,@new_right);
		#$left==0 if no signature from cluster.
		if($left>0){
			#search nearby softclip
			$left_clip=clip_mode($left,$right,$min_det_size-$read_length,$clip_spt_pem,@left_clip) if($#left_clip>=$clip_spt_pem-1);		
			$right_clip=clip_mode($right,$left,$min_det_size-$read_length,$clip_spt_pem,@right_clip) if($#right_clip>=$clip_spt_pem-1);		
			if($left_clip>0 or $right_clip>0){
				if($left_clip>0 and $right_clip>0){
					$left=$left_clip;
					$right=$right_clip;	
					$prec=1;
				}
				elsif($left_clip>0){
					$left=$left_clip;
				}
				elsif($right_clip>0){
					$right=$right_clip;
				}
		
			}
			return ($#left+1, $left, $right, $prec);
		}
	}
	
	#only softclip signiture
	elsif($length>=$softclip_check_thh){
		#higher softclip_spt
		$left_clip=clip_mode($start,$end,$win_size_del,$clip_spt_npem,@left_clip) if($#left_clip>=$clip_spt_npem-1);		
		$right_clip=clip_mode($end,$start,$win_size_del,$clip_spt_npem,@right_clip) if($#right_clip>=$clip_spt_npem-1);		
		if($left_clip>0 and $right_clip>0){
			$left=$left_clip;
			$right=$right_clip;	
			return (1111, $left, $right, $prec);
		}	
	}
	return 0;
}

sub clip_mode{
	my ($center,$bound_mate,$distance,$clip_spt,@clip)=@_;
	my @new_clip;
	for(my $i=0; $i<=$#clip; $i++){
		if(abs($clip[$i]-$center)<$distance){
			push(@new_clip,$clip[$i]);
		}
	}
	my %mode;
	for(my $i=0; $i<=$#new_clip; $i++){
		if(!exists($mode{$new_clip[$i]})){
			$mode{$new_clip[$i]}=1;
		}
		else{
			$mode{$new_clip[$i]}++;
		}
	}
	my $max=0;
	my $overlap_max=0;
	my $clip;
	foreach my $key (keys %mode) {
		my $overlap=abs($key-$bound_mate);
		if ($max*$overlap_max<$mode{$key}*$overlap){
			$max=$mode{$key};
			$overlap_max=$overlap;
			$clip=$key;
		}
	}
	return $clip if($max>=$clip_spt);
	return 0 if($max<$clip_spt);
}

sub region_cn{
	my ($start,$end)=@_;
	my $length=$end-$start+1;
	my $exp_rd=0;
	my $region_cn=0;
	my $ref_cn=0;
	my @rd_windows;
	my $s=int(($start-1)/$win_size_del);
	my $e=int(($end-1)/$win_size_del);
	if($s==$e){
		$exp_rd=$raw_rd_gc[$s];
		$ref_cn=$ref_cn[$s];
	}
	else{
		$exp_rd+=($win_size_del*($s+1)-$start+1)/$win_size_del*$raw_rd_gc[$s];
		$ref_cn+=($win_size_del*($s+1)-$start+1)/$win_size_del*$ref_cn[$s];
		$exp_rd+=($end-$win_size_del*$e)/$win_size_del*$raw_rd_gc[$e];
		$ref_cn+=($end-$win_size_del*$e)/$win_size_del*$ref_cn[$e];
		for(my $i=$s+1; $i<=$e-1; $i++){
			$exp_rd+=$raw_rd_gc[$i];
			$ref_cn+=$ref_cn[$i];
		}
		$ref_cn=int($ref_cn*$win_size_del/$length+0.5);
	}
	if($length>=$large_pem){
		push(@rd_windows,($win_size_del*($s+1)-$start+1)/$win_size_del*$raw_rd[$s]);
		push(@rd_windows,($end-$win_size_del*$e)/$win_size_del*$raw_rd[$e]);	
		for(my $i=$s+1; $i<=$e-1; $i++){
			push(@rd_windows,$raw_rd[$i]);
		}
	}
	else{
		@rd_windows=rd_bam_region_bp($chr,$start,$end);
	}
	$region_cn=midhalf(@rd_windows);
	return ($ref_cn,0) if($region_cn==0 and $exp_rd==0);
	return ($ref_cn,$normal_cn) if($region_cn>0 and $exp_rd==0);
	$region_cn=int($region_cn/$exp_rd+0.5);
	return ($ref_cn, $region_cn);
}

sub rd_bam_region_bp{
	my ($chr,$start,$end)=@_;
	my $length=$end-$start+1;
	my @rd;
	open(ST, "$samtools view $bam_file $chr:$start-$end| ")||die"$!\n";
	while(my $row=<ST>) {
		if($row=~m/^\S+\t\d+\t\S+\t(\d+)\t\d+\t\w+\t\S+\t\S+\t\S+\t(\S+)\t/){
			my $read_start=$1;
			my $rd_length=length($2);
			my $read_end=$read_start+$rd_length-1;
			for(my $i=max($read_start,$start); $i<=min($read_end,$end); $i++){
				$rd[$i-$start]++;
			}			
		}
	}
	close(ST) ||die"$!\n";
	for(my $i=0; $i<=$#rd; $i++){
		$rd[$i]/=$win_size_del;
	}
	return @rd;
}

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