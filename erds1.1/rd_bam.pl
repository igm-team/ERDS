#!/usr/bin/perl
use strict;
use List::Util qw[min max];

our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
my $parameter_sample=$ARGV[0];
my $functions=$ARGV[1];
require "$parameter_sample";
require "$functions";

my $chr=$ARGV[2];
my $chr_rc="$rc_dir/$chr.rc";
my $chr_sd_single="$sd_dir/$chr.single.sd";
my $chr_sd_rc="$rc_dir/$chr.sd.rc";
my $chr_pem_sam="$rc_dir/$chr.pem.sam";
my $chr_pem_bam="$rc_dir/$chr.pem.bam";

my ($row, $pos, $CIAGR, $alg_length, $ins_read, $mq, $bq);
my (@rc, @mq, @rd_bp);
my ($rc_avg, $mq_avg);
my (@list, $in);

#read bam file
if($chr=~m/^c*h*r*X|c*h*r*Y|c*h*r*\d+/){
	system("$samtools view -H $bam_file >$chr_pem_sam");
	open(PSAM,">>$chr_pem_sam") or die "Can't open $chr_pem_sam";
	read_bam($bam_file);
	close(PSAM) or die "Can't close $chr_pem_sam";
	system("$samtools view -Sb $chr_pem_sam >$chr_pem_bam");
	system("$samtools index $chr_pem_bam");
	system("rm $chr_pem_sam");
}
else{
	read_bam_simple($bam_file);
}

#print out
my ($start, $end);
open(RC,">$chr_rc") or die "Can't open $chr_rc";
for(my $i=0; $i<=$#rc; $i++){
	$mq_avg=60;
	if(defined($mq[$i])){
		$mq_avg=$mq[$i]/$rc[$i];
	}
	my @rd_window;
	$start=$i*$win_size_del+1;
	$end=($i+1)*$win_size_del;
	for(my $i=$start; $i<=$end; $i++){
		push(@rd_window, $rd_bp[$i]);
	}	
	$rc_avg=midhalf(@rd_window)/$win_size_del;
	printf RC ("$start\t$end\t%.2f\t%d\n",$rc_avg,int($mq_avg+0.5));	
}
close(RC) or die "Can't close $chr_rc";

open(SDS,"$chr_sd_single") or die "Can't open $chr_sd_single";
open(SDRC,">$chr_sd_rc") or die "Can't open $chr_rc";
while(<SDS>){
	if(~ m/^\w+:(\d+)-(\d+)\n/){
		$start=$1;
		$end=$2;
	}
	my @rd_window;
	for(my $i=$start; $i<=$end; $i++){
		push(@rd_window, $rd_bp[$i]);
	}	
	$rc_avg=midhalf(@rd_window)/($end-$start+1);
	printf SDRC ("$start\t$end\t%.2f\n",$rc_avg);	
}
close(SDS) or die "Can't close $chr_sd_single";
close(SDRC) or die "Can't close $chr_rc";

##functions
sub read_bam{
	$bam_file=$_[0];
	my $ins_cut=$ins_pem+$min_det_size;
	my (@row,@pos,@pem);
	open(ST, "$samtools view $bam_file $chr| ")||die"$!\n";
	while(my $in=<ST>) {
		$in=~s/[\r\n]//g;
		@list=split /\s/, $in;
		$pos=$list[3];
		$mq=$list[4];
		$CIAGR=$list[5];
		$alg_length=0;
		while($CIAGR=~/(\d+)[M]/g){
			$alg_length+=$1;
		}
		next if($alg_length<$aligned_length_thh or $CIAGR eq "*" or length($CIAGR)>$CIAGR_length);
		$ins_read=abs($list[8]);
		$bq=$list[10];
		$list[5]="1M";

		if($ins_read>$ins_cut and $ins_read<2*$large_pem){
			$row="$list[0]";
			for(my $j=1; $j<=8; $j++){
				$row=$row."\t$list[$j]";
			}
			$row=$row."\tA\tH\tMD:Z:$CIAGR";
			push(@row,$row);
			push(@pos,$pos);
			push(@pem,1);
		}
		elsif($in!~/XC:i:/ and $CIAGR=~/S/){
			$row="$list[0]";
			for(my $j=1; $j<=8; $j++){
				$row=$row."\t$list[$j]";
			}
			$row=$row."\tA\tH\tMD:Z:$CIAGR";
			push(@row,$row);
			push(@pos,$pos);
			push(@pem,0);
		}

		$alg_length=0;
		while($CIAGR=~/(\d+)[MI]/g){
			$alg_length+=$1;
		}
		rc($pos,$alg_length,$mq);		
	}
	close(ST) ||die"$!\n";
	for(my $j=0; $j<=$#row; $j++){
		if($pem[$j]==1 or $pos[$j]-$pos[$j-1]<$read_length or $pos[$j+1]-$pos[$j]<$read_length){
			print PSAM "$row[$j]\n";
		}
	}
}

sub read_bam_simple{
	$bam_file=$_[0];
	my $ins_cut=$ins_pem+$min_det_size;
	open(ST, "$samtools view $bam_file $chr| ")||die"$!\n";
	while(my $row=<ST>) {
		if($row=~m/^\S+\t\d+\t\S+\t(\d+)\t(\d+)\t(\w+)\t/){
			$pos=$1;
			$mq=$2;
			$CIAGR=$3;
			$alg_length=0;
			while($CIAGR=~/(\d+)[M]/g){
				$alg_length+=$1;
			}
			if($alg_length>=$aligned_length_thh){
				$alg_length=0;
				while($CIAGR=~/(\d+)[MI]/g){
					$alg_length+=$1;
				}
				rc($pos,$alg_length,$mq);		
			}
		}
	}
	close(ST) ||die"$!\n";
}

sub rc{
	my ($start, $alg_length, $mq)=@_;
	my $end=$start+$alg_length-1;
	
	for(my $i=$start; $i<=$end; $i++){
		$rd_bp[$i]++;
	}	
	
	my $s=int(($start-1)/$win_size_del);
	my $e=int(($end-1)/$win_size_del);	
	if($s==$e){
		$rc[$s]+=$alg_length;
		$mq[$s]+=$mq*$alg_length;
	}
	else{
		$rc[$s]+=($win_size_del*($s+1)-$start+1);
		$rc[$e]+=($end-$win_size_del*$e);
		$mq[$s]+=$mq*($win_size_del*($s+1)-$start+1);
		$mq[$e]+=$mq*($end-$win_size_del*$e);
		for(my $i=$s+1; $i<=$e-1; $i++){
			$rc[$i]+=$win_size_del;
			$mq[$i]+=$mq*$win_size_del;
		}				
	}
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
