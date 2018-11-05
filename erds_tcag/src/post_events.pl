#Authors: Mingfu Zhu & David B. Goldstein
#!/usr/bin/perl -w
use strict;
use List::Util qw[min max];

##input and directory setup
##################################################
our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
my $parameter_sample=$ARGV[0];
my $functions=$ARGV[1];
require "$parameter_sample";
require "$functions";

my $del="$sample_dir/$samplename.del.events";

##events
my ($chr, $start, $end, $length, $ref_cn, $cn);
my (@list, $in);
my $overlap_rate_thh1=0.5;
my $overlap_rate_thh2=0.9;
my $del_thh=100;

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
		}
		if($chr[$j]!~m/^c*h*r*X|Y$/){
			$gender{$chr[$j]}=2;
		}
		else{
			$gender{$chr[$j]}=1;
		}
	}
}
else{
	$gender="F";
	for(my $j=0; $j<=$#chr; $j++){
		if($chr[$j]=~m/^c*h*r*\d+|X$/){
			push(@aut_chr,$chr[$j]);
		}
		if($chr[$j]!~m/^c*h*r*Y$/){
			$gender{$chr[$j]}=2;
		}
		else{
			$gender{$chr[$j]}=0;
		}
	}
}

open(EVTDEL,">$del") or die "Can't open $del";
for(my $j=0; $j<=$#normal_chr; $j++){
	$chr=$normal_chr[$j];
	my $chr_evt_del="$cnv_dir/$chr.del.evt";
	next if (!-e $chr_evt_del);

	my ($log_sum,@chr,@start,@end,@log_sum,@rows);
	open(CDEL,"$chr_evt_del") or die "Can't open $chr_evt_del";
	while ($in=<CDEL>) {
		$in=~s/[\r\n]//g;
		@list=split /\s/, $in;
		$length=$list[3];
		$ref_cn=$list[7];
		$cn=$list[8];
		next if($length<$del_thh or ($ref_cn>$gender{$chr} and $cn>=$ref_cn));

		$start=$list[1];
		$end=$list[2];
		$log_sum=$list[5];
		push(@start,$start);
		push(@end,$end);
		push(@log_sum,$log_sum);		
		push(@rows,"$in\n");
	}
	close(CDEL) or die "Can't close $chr_evt_del";
	
	my $num_rows=$#start;
	for(my $j=0; $j<=$num_rows-1; $j++){
		my $overlap=min($end[$j], $end[$j+1])-max($start[$j], $start[$j+1]);
		if($overlap>0){
			my $lengtha=$end[$j]-$start[$j]+1;
			my $lengthb=$end[$j+1]-$start[$j+1]+1;
			if($overlap/$lengtha>=$overlap_rate_thh1 and $overlap/$lengthb>=$overlap_rate_thh1){
				if($log_sum[$j]>$log_sum[$j+1]){
					$rows[$j+1]="";
				}
				elsif($log_sum[$j]<=$log_sum[$j+1]){
					$rows[$j]="";
				}
				else{
					if($lengtha<$lengthb){
						$rows[$j+1]="";
					}
					else{
						$rows[$j]="";
					}				
				}
			}
			elsif($overlap/$lengtha>=$overlap_rate_thh2){
				$rows[$j]="";
			}
			elsif($overlap/$lengthb>=$overlap_rate_thh2){
				$rows[$j+1]="";
			}
		}
	}
	
	for(my $j=0; $j<=$#start; $j++){
		print EVTDEL "$rows[$j]";
	}
}
close(EVTDEL) or die "Can't close $del";
sortfile($del);
