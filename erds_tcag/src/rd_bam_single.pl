#!/usr/bin/perl -w
use strict;

our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
my $parameter_sample=$ARGV[0];
my $functions=$ARGV[1];
my $num_processors = $ARGV[2];
require "$parameter_sample";
require "$functions";

sub partition_array {
	my $num_processors = shift;
	my @array = reverse @_;

	my @partitions;

	for my $i (0..$num_processors - 1) {
		my @this_partition;

		for (my $j = $i; $j <= $#array; $j += $num_processors) {
			push @this_partition, $array[$j];
		}

		push @partitions, [ @this_partition ];
	}

	return @partitions;
}
my @partitions = partition_array($num_processors, @chr);
my @child_pids;
for my $i (0..$#partitions) {
	if (my $pid = fork) {
		push @child_pids, $pid;
	} else {
		foreach my $chr (@{$partitions[$i]}) {
			my $time = localtime;
			print "[$time] PID $$ is performing computation for $chr\n";
			system("perl $script_dir/rd_bam.pl $parameter_sample $functions '$chr' > '$sample_dir/debug/rd_bam_chr$chr.txt' 2> '$sample_dir/debug/rd_bam_chr$chr.err'");
			print "Error code for $chr is $?\n" if $? != 0;
		}
		exit;
	}
}

for my $i (0..$#partitions) {
	waitpid($child_pids[$i], 0);
}
