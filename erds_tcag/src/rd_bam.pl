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
my $SAMTOOLS_EXTRA_FLAGS = "-F 256";

my ($row, $pos, $CIAGR, $alg_length, $ins_read, $mq, $bq);
my (@rc, @mq, @rd_bp);
my ($rc_avg, $mq_avg);
my (@list, $in);

my $time;

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


sub read_bam_simple{
	my $bfile_to_read=$_[0];
	my $ins_cut=$ins_pem+$min_det_size;
	open(ST, "$samtools view $SAMTOOLS_EXTRA_FLAGS '$bfile_to_read' $chr| ")||die"$!\n";
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

my $bam_file_to_read = $bam_file eq "ERDS_SEPARATE_BAMS" ? "$chr.bam" : $bam_file;

#read bam file
if($chr=~m/^c*h*r*X|c*h*r*Y|c*h*r*\d+/){
	system("$samtools view $SAMTOOLS_EXTRA_FLAGS -H '$bam_file_to_read' > '$chr_pem_sam'");
	
	$time = localtime;

	print "[$time] Starting $ERDS_dir/read_bam script with parameters '$bam_file_to_read' $ins_pem $min_det_size '$chr' $aligned_length_thh $CIAGR_length $large_pem $win_size_del $read_length '$chr_rc' '$chr_sd_rc' '$chr_pem_sam'\n";
	system("$ERDS_dir/read_bam '$bam_file_to_read' $ins_pem $min_det_size '$chr' $aligned_length_thh $CIAGR_length $large_pem $win_size_del $read_length '$chr_rc' '$chr_sd_rc' '$chr_pem_sam' > '$sample_dir/debug/read_bam_chr$chr.txt' 2> '$sample_dir/debug/read_bam_chr$chr.err'");
	print "read_bam script complete for chr$chr\n";

	system("$samtools view $SAMTOOLS_EXTRA_FLAGS -Sb '$chr_pem_sam' > '$chr_pem_bam'");
	system("$samtools index '$chr_pem_bam'");
	system("rm '$chr_pem_sam'");
} else {
	print "Doing read_bam_simple\n";
	read_bam_simple($bam_file_to_read);

	#print out
	my ($start, $end);
	open(RC,">$chr_rc") or die "Can't open $chr_rc";
	for(my $i=0; $i<=$#rc; $i++){
		$mq_avg=60;
		if(defined($mq[$i])){
			$mq_avg=$mq[$i]/$rc[$i];
		}

		print "mq_avg for $i is $mq_avg\n";

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
	`touch $chr_sd_rc`;
}


