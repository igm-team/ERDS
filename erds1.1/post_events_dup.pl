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

my $dup="$sample_dir/$samplename.dup.events";
system("rm $dup") if(-e $dup);
$cnv_dir="$cnv_dir"."_"."$win_size_dup";

##events
my (@list,$in);
open(EVTDUP,">$dup") or die "Can't open $dup";
for(my $j=0; $j<=$#normal_chr; $j++){
	my $chr=$normal_chr[$j];
	my $chr_evt_dup="$cnv_dir/$chr.dup.evt";
	next if(!-e $chr_evt_dup);
	open(CDUP,$chr_evt_dup) or die "Can't open $chr_evt_dup";
	while($in=<CDUP>){
		print EVTDUP $in;
	}
	close(CDUP) or die "Can't close $chr_evt_dup";
}
close(EVTDUP) or die "Can't close $dup";
sortfile($dup);
