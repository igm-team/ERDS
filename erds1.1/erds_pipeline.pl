#!/usr/bin/perl
#Author: Mingfu Zhu, mingfu.zhu@duke.edu
#Version 1.1
#http://www.duke.edu/~mz34/erds.htm

use strict;
use Cwd qw(abs_path);
use List::Util qw[min max];
use Getopt::Long qw(:config bundling);

usage() if($#ARGV<0 or $ARGV[0] eq "h" or $ARGV[0] eq "help");
my $erds_version="ERDS1.1";

###########setup arguments
#line a, copy to every script
our ($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
our (@chr, @normal_chr);
#order indentical to line a.
my @para_name=qw(samplename bam_file variant_file aligned_length_thh ERDS_dir mq_thh exp_pct num_state ref_file ref_fai_file win_size_del win_size_dup sample_dir script_dir sd_dir log_dir rc_dir cnv_dir gcbin_dir samtools max_mq read_length min_det_size sd_win_size ins_pem gap_events low_scan_ratio pem_pairs low_pem_pairs rcp_rate clip_len clip_spt_pem clip_spt_npem softclip_check_thh large_pem small_pem CIAGR_length);

my $begin_step=0;
my $end_step=10;
my $cluster=0;
my $parameter;
my $sd="empty";
my $num_stdd=4;
my @para;

print "\n############################$erds_version############################\n";
my $path = abs_path($0);
if($path=~/([\S\/]*)\/erds_pipeline\.pl$/){
	$ERDS_dir=$1;
}
print "Your input command line is\nperl $path ";

my @log;
for(my $j=0; $j<=$#ARGV; $j++){
    print "$ARGV[$j] ";
	push(@log, "$ARGV[$j] ");
}
print "\n";

GetOptions(
	#You must define the following parameters.
	'o=s' => \$sample_dir, #-o: followed by the path of output directory.(-o output_dir)
	'b=s' => \$bam_file, #-b: followed by specified bam_file. (-b bam_file)
	'v=s' => \$variant_file, #-v: followed by specified variant_file file called by SNV calling tools. (-v variant_file)
	'r=s' => \$ref_file, #-r: followed by specified reference file. (-r reference_fasta_file)

	#You can optionally define the following parameters.
	'cluster|c' => \$cluster, #--cluster or -c: if indicated, run cluster. By default single cpu. (-c)
	'samtools=s' => \$samtools, #--samtools: followed by path_to_samtools. By default the attached version 1.12. (--samtools path_to_samtools)
	'sd=s' => \$sd, #--sd: followed by the SD version of b36, b37 or empty. By default empty. (--sd b36|b37|empty)
	'exp_pct|p=s' => \$exp_pct, #--exp_pct or -p: followed by expected percentage of deletions or duplications. By default 0.0025. (-p value)
	'n=s' => \$samplename, #--name: followed by the string of the sample name. (--name sample_name)
	
	#The following parameters are for advanced users only. You can optionally define them only if you are sure.
	'parameter=s' => \$parameter, #--parameter: followed by specified parameter file. (--parameter parameter_file)
	'del=s' => \$win_size_del, #--del: followed by specified win_size_del. By default 200. (--del win_size_del)
	'dup=s' => \$win_size_dup, #--dup: followed by specified win_size_dup. By default 1000. (--dup win_size_dup)
	'large=s' => \$large_pem, #--large: followed by specified length (in bp). Deletions larger than this value will escape from PEM check. By default 50000.(--large integer)
	'pem=s' => \$pem_pairs, #--pem: followed by specified value for PEM support required. By default 1. (--pem integer)
	'small=s' => \$small_pem, #--small: followed by specified length (in bp). Deletions smaller than this value will require more PEM support (\$pem+1). By default 1000. (--small integer)
	'scan=s' => \$low_scan_ratio, #--scan: followed by specified value for low RD ratio scan. By default 1.4. (--scan value)
	'pairs=s' => \$low_pem_pairs, #--pairs: followed by specified value for PEM support required for low_pem_pairs. By default 3. (--pairs integer)
	'clip=s' => \$clip_spt_pem, #--clip: followed by specified value for clip_spt_pem. By default 2. (--clip integer)
	'soft=s' => \$softclip_check_thh, #--soft: followed by specified length for softclip_check. By default 1000. (--soft integer)
	'num_stdd=s' => \$num_stdd, #--num_stdd: followed by specified number. Mapping distance larger than mean library size plus \$num_stdd*\$std_dev is regarded as discordant. By default 4. (--num_stdd integer)
	'begin=s' => \$begin_step, #--begin: followed by the step to run. Use for control the program during debuging or rerun. By default 0. (--begin integer) 
	'end=s' => \$end_step, #--end: followed by the step to stop. By default 10. (--end integer)
);

die "Cannot find bam_file $bam_file\n" if (!-e $bam_file);
die "Cannot find variant_file $variant_file\n" if (!-e $variant_file);
die "Cannot find ref_file $ref_file\n" if (!-e $ref_file);
$ref_fai_file="$ref_file.fai";
die "Cannot find ref_fai_file $ref_fai_file\n" if (!-e $ref_fai_file);

if(!defined($parameter)){
	$parameter="$ERDS_dir/parameter.txt";
}
my $functions="$ERDS_dir/functions.pl";
require "$functions";
system("mkdir -p $sample_dir") if(!-e $sample_dir);
$script_dir="$sample_dir/script";
system("rm -r $script_dir") if (-e $script_dir);
system("cp -r $ERDS_dir $script_dir");
$sd_dir="$script_dir/database/SD_list_"."$sd";

if(!defined($samtools)){
	$samtools="$script_dir/software/samtools";
}
if(!defined($samplename)){
	$samplename=get_samplename();
}

my $sample_log="$sample_dir/sample.log";
open(LOG,">$sample_log")||die"cannot open $sample_log\n";
print LOG "Your input command line is\nperl $path ";
for(my $j=0; $j<=$#log; $j++){
    print LOG "$log[$j]";
}
print LOG "\n";
close(LOG)||die"cannot close $sample_log\n";

my @para=($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
read_parameter($parameter);

my @chr_length;
open(FAI,"$ref_fai_file") or die "Can't open $ref_fai_file";
while(my $in=<FAI>){
	$in=~m/^(\S+)\t/;
	push(@chr,$1);
	if($in=~m/^(c*h*r*X|c*h*r*Y|c*h*r*\d+)\t(\d+)\t/){
		push(@normal_chr,$1);
		push(@chr_length,$2);	
	}
}
close(FAI) or die "Can't close $ref_fai_file";

my $parameter_sample="$sample_dir/parameter-$samplename.pl";
my $num_ram=19;

test_bam();
@para=($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length);
open(PP,">$parameter_sample") or die "Can't open $parameter_sample";
print PP "our \$parameter=\"$parameter_sample\";\n";
print PP "our \$cluster=$cluster;\n";
for(my $j=0; $j<=$#para_name; $j++){
	if($para[$j]!~/^\d+$/ or $para[$j]=~/^0\d+/){
		print PP "our \$$para_name[$j]=\"$para[$j]\";\n";
	}
	else{
		print PP "our \$$para_name[$j]=$para[$j];\n";
	}
}
printf PP "our \@chr=qw(%s);\nour \@normal_chr=qw(%s);\n", join(" ",@chr), join(" ",@normal_chr);
close(PP);

require "$parameter_sample" or die "$parameter_sample does not exist!";
if(-e $log_dir){
	system("rm -r $log_dir");
}
system("mkdir -p $log_dir");
if(!-e $rc_dir){
	mkdir $rc_dir;
}
if(!-e $gcbin_dir){ 
   system("mkdir -p $gcbin_dir");
}
$functions="$script_dir/functions.pl";

my %chr;
for(my $j=0; $j<=$#chr; $j++){
    $chr{$chr[$j]}=$j;
}
sd_scan();

#####cluster
if($cluster){
	my $chr_range=$#chr+1;
	my $normal_chr_range=$#normal_chr+1;
	#rd_bam, step 0
	script_gen($samplename,"$ERDS_dir/rd_bam_cluster.txt","$script_dir/rd_bam_cluster.pl","rb_"."$samplename","1-$chr_range","");
	system("qsub $script_dir/rd_bam_cluster.pl $parameter_sample $functions") if($begin_step<=0 and $end_step>=0);
	
	#prehmm, step 1
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/pre_hmm_cluster.pl","pre_"."$samplename", "", "rb_"."$samplename");
	system("qsub $script_dir/pre_hmm_cluster.pl \"pre_hmm.pl\" $parameter_sample $functions") if($begin_step<=1 and $end_step>=1);
	
	#em, step 2
	script_gen($samplename,"$ERDS_dir/em_cluster.txt","$script_dir/em_cluster.pl","em_"."$samplename","","pre_"."$samplename");
	system("qsub $script_dir/em_cluster.pl $parameter_sample $functions") if($begin_step<=2 and $end_step>=2);
	
	#hmm, step 3
	script_gen($samplename,"$ERDS_dir/hmm_cluster.txt","$script_dir/hmm_cluster.pl","hmm_"."$samplename","1-$normal_chr_range","em_"."$samplename");
	system("qsub $script_dir/hmm_cluster.pl $parameter_sample $functions") if($begin_step<=3 and $end_step>=3);
	
	#events, step 4
	script_gen($samplename,"$ERDS_dir/events_cluster.txt","$script_dir/events_cluster.pl","evt_"."$samplename","1-$normal_chr_range","hmm_"."$samplename");
	system("qsub $script_dir/events_cluster.pl $parameter_sample $functions") if($begin_step<=4 and $end_step>=4);
	
	#post_events, step 5
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/post_events_cluster.pl","post_"."$samplename", "", "evt_"."$samplename");
	system("qsub $script_dir/post_events_cluster.pl \"post_events.pl\" $parameter_sample $functions") if($begin_step<=5 and $end_step>=5);
	
	#prehmm_dup, step 6
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/pre_hmm_dup_cluster.pl","pred_"."$samplename", "", "hmm_"."$samplename");
	system("qsub $script_dir/pre_hmm_dup_cluster.pl \"pre_hmm_dup.pl\" $parameter_sample $functions") if($begin_step<=6 and $end_step>=6);
	
	#hmm_dup, step 7
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/hmm_dup_cluster.pl","hmmd_"."$samplename","1-$normal_chr_range","pred_"."$samplename");
	system("qsub $script_dir/hmm_dup_cluster.pl \"hmm_dup.pl\" $parameter_sample $functions") if($begin_step<=7 and $end_step>=7);
	
	#events_dup, step 8
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/events_dup_cluster.pl","evtd_"."$samplename","1-$normal_chr_range","hmmd_"."$samplename");
	system("qsub $script_dir/events_dup_cluster.pl \"events_dup.pl\" $parameter_sample $functions") if($begin_step<=8 and $end_step>=8);
	
	#post_events_dup, step 9
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/post_events_dup_cluster.pl","postd_"."$samplename", "", "evtd_"."$samplename");
	system("qsub $script_dir/post_events_dup_cluster.pl \"post_events_dup.pl\" $parameter_sample $functions") if($begin_step<=9 and $end_step>=9);
	
	#merge, step 10
	script_gen($samplename,"$ERDS_dir/cluster.txt","$script_dir/merge_cluster.pl","merge_"."$samplename", "", "post_"."$samplename,postd_"."$samplename");
	system("qsub $script_dir/merge_cluster.pl \"merge.pl\" $parameter_sample $functions") if($begin_step<=10 and $end_step>=10);
}

#####single CPU
else{
	#rd_bam, step 0
	printtime("ERDS stared step 0 in $samplename...");
	system("perl $script_dir/rd_bam_single.pl $parameter_sample $functions") if($begin_step<=0 and $end_step>=0);
	printtime("ERDS finished step 0 in $samplename...");
	print ;
	
	#prehmm, step 1
	printtime("ERDS stared step 1 in $samplename...");	
	system("perl $script_dir/pre_hmm.pl $parameter_sample $functions") if($begin_step<=1 and $end_step>=1);
	printtime("ERDS finished step 1 in $samplename...");
	
	#em, step 2
	printtime("ERDS stared step 2 in $samplename...");
	system("perl $script_dir/em_single.pl $parameter_sample $functions") if($begin_step<=2 and $end_step>=2);
	printtime("ERDS finished step 2 in $samplename...");
	
	#hmm, step 3
	printtime("ERDS stared step 3 in $samplename...");
	system("perl $script_dir/single.pl hmm.pl $parameter_sample $functions") if($begin_step<=3 and $end_step>=3);
	printtime("ERDS finished step 3 in $samplename...");
	
	#events, step 4
	printtime("ERDS stared step 4 in $samplename...");
	system("perl $script_dir/single.pl events.pl $parameter_sample $functions") if($begin_step<=4 and $end_step>=4);
	printtime("ERDS finished step 4 in $samplename...");
	
	#post_events, step 5
	printtime("ERDS stared step 5 in $samplename...");
	system("perl $script_dir/post_events.pl $parameter_sample $functions") if($begin_step<=5 and $end_step>=5);
	printtime("ERDS finished step 5 in $samplename...");
	
	#pre_hmm_dup, step 6
	printtime("ERDS stared step 6 in $samplename...");
	system("perl $script_dir/pre_hmm_dup.pl $parameter_sample $functions") if($begin_step<=6 and $end_step>=6);
	printtime("ERDS finished step 6 in $samplename...");
	
	#hmm_dup, step 7
	printtime("ERDS stared step 7 in $samplename...");
	system("perl $script_dir/single.pl hmm_dup.pl $parameter_sample $functions") if($begin_step<=7 and $end_step>=7);
	printtime("ERDS finished step 7 in $samplename...");
	
	#events_dup, step 8
	printtime("ERDS stared step 8 in $samplename...");
	system("perl $script_dir/single.pl events_dup.pl $parameter_sample $functions") if($begin_step<=8 and $end_step>=8);
	printtime("ERDS finished step 8 in $samplename...");
	
	#post_events_dup, step 9
	printtime("ERDS stared step 9 in $samplename...");
	system("perl $script_dir/post_events_dup.pl $parameter_sample $functions") if($begin_step<=9 and $end_step>=9);
	printtime("ERDS finished step 9 in $samplename...");
	
	#merge, step 10
	printtime("ERDS stared step 10 in $samplename...");
	system("perl $script_dir/merge.pl $parameter_sample $functions") if($begin_step<=10 and $end_step>=10);
	printtime("ERDS finished all steps in $samplename. The results are at $sample_dir.");
}
print "\n############################$erds_version############################\n\n";

#functions
sub get_samplename{
	open(ST, "$samtools view -H $bam_file|"); 
	while(my $in=<ST>) {
		if($in=~/\sSM:(\S+)\s/){
			$samplename=$1;
			last;
		}
	}
	close(ST);
	return $samplename;
}

sub read_parameter{
	my $parameter=$_[0];
	my %para;
	open(P,"$parameter") or die "Cannot open $parameter\n";
	while (my $in=<P>) {
		$in=~s/\s|,|;//g;
		if($in=~m/^(\w+):\s*(\S+)/){
			$para{$1}=$2;
		}
	}
	for(my $j=0; $j<=$#para_name; $j++){
		if(exists($para{$para_name[$j]}) and !defined($para[$j])){
			$para[$j]=$para{$para_name[$j]};
		}
	}
	close(P);
	#in an order identical to line a
	($samplename, $bam_file, $variant_file, $aligned_length_thh, $ERDS_dir, $mq_thh, $exp_pct, $num_state, $ref_file, $ref_fai_file, $win_size_del, $win_size_dup, $sample_dir, $script_dir, $sd_dir, $log_dir, $rc_dir, $cnv_dir, $gcbin_dir, $samtools, $max_mq, $read_length, $min_det_size, $sd_win_size, $ins_pem, $gap_events, $low_scan_ratio, $pem_pairs, $low_pem_pairs, $rcp_rate, $clip_len, $clip_spt_pem, $clip_spt_npem, $softclip_check_thh, $large_pem, $small_pem, $CIAGR_length)=@para;
}

sub test_bam{
	my (@mq, @CIAGR, @ins, @len, @smq, @sins, @slen, @bq);
	for(my $j=0; $j<=$num_ram; $j++){
		my $ran=int(rand($#normal_chr+1));
		my $chr_ran=$normal_chr[$ran];
		my $region_length=100000;
		my $start_ran=int(rand($chr_length[$ran]-$region_length));
		my $end_ran=$start_ran+$region_length-1;		
		open(ST, "$samtools view -q 10 $bam_file $chr_ran:$start_ran-$end_ran| awk '{print \$2 \"\t\" \$5 \"\t\" \$6 \"\t\" \$9 \"\t\" length(\$10) \"\t\" \$11'} |"); 
		while(<ST>) {
			if(m/^(\d+)\t(\d+)\t(\S+)\t(\S+)\t(\d+)\t(\S+)\n/){
				my $flag=$1;
				my @b=bin($flag);
				my $mq=$2;
				my $CIAGR=$3;
				my $insert=$4;
				my $len=$5;
				my $bq=$6;
				if($b[0]==1 and $b[1]==1){
					push(@mq, $mq);
					push(@ins, abs($insert));
					push(@len, $len);			
				}
				if($CIAGR=~/^(\d+)S/){
					push(@CIAGR, $CIAGR);
					push(@bq, $bq);
				}
			}
		}
		close(ST);
	}
	
	@smq = sort {$a <=> $b} @mq;
	@slen = sort {$a <=> $b} @len;
	if($#smq>1000){
		$max_mq=$smq[$#smq];
		$ins_pem=int(med(@ins)+0.5);
		$read_length=$slen[int($#slen*0.9)];
		my $ins_cut=$num_stdd*std_dev(@ins);
		my @new_ins;
		for(my $j=0; $j<=$#ins; $j++){
			if(abs($ins[$j]-$ins_pem)<=$ins_cut){
				push(@new_ins, $ins[$j]);		
			}
		}
		$min_det_size=int($num_stdd*std_dev(@new_ins)+0.5);
	}
	else{
		die "no enough reads\n";
	}
}

sub sd_scan{
	my @chr_sd_single;
	for(my $i=0; $i<=$#chr; $i++){
		my $chr=$chr[$i];
		my $chr_sd;

		if($chr=~m/(c*h*r*)(M\w?)/){
			my $short_chr="MT";
			$chr_sd="$sd_dir/$short_chr.sd";
			system("mv $chr_sd $sd_dir/$chr.sd") if $short_chr ne $chr;
		}
		elsif($chr=~m/(?:chr|ch|c)*(\w+\.*\w*)/){
			my $short_chr=$1;
			$chr_sd="$sd_dir/$short_chr.sd";
			if ($short_chr ne $chr and $sd ne "empty"){
				system("mv $chr_sd $sd_dir/$chr.sd");
			}
		}

		$chr_sd="$sd_dir/$chr.sd";
		if($sd eq "empty"){
			open(SD,">$chr_sd");
			close(SD);
		}
		open(SD,"$chr_sd") or die "Can't open $chr_sd";
		while(my $in=<SD>){
			$in=~s/[\r\n]//g;
			my @list=split /\s/, $in;
			for(my $k=0; $k<=$#list; $k++){
				if($list[$k]=~ m/^(\w+):(\d+)-(\d+)/){
					my $chr_sd=$1;
					my $start_sd=int($2/$win_size_del+0.5)*$win_size_del+1;
					my $end_sd=int($3/$win_size_del+0.5)*$win_size_del;
					my $locus="$chr_sd:$start_sd-$end_sd";
					$chr_sd_single[$chr{$1}].="$locus\n";
				}
			}
		}
		close(SD) or die "Can't close $chr_sd";
	}
	for(my $i=0; $i<=$#chr; $i++){
		my $chr=$chr[$i];
		my $chr_sd_single="$sd_dir/$chr.single.sd";
		open(SDS,">$chr_sd_single") or die "Can't open $chr_sd_single";
		if(defined($chr_sd_single[$i])){
			print SDS $chr_sd_single[$i];
		}
		close(SDS) or die "Can't close $chr_sd_single";
	}
}

sub script_gen{
	my $sample_name=$_[0];
	my $blue_print=$_[1];
	my $script=$_[2];
	my $job_id=$_[3];
	my $job_range=$_[4];
	my $hold_id=$_[5];
	open(BP,$blue_print) or die "Can't open $blue_print";
	open(SC,">$script") or die "Can't open $script";
	print SC "\#\$\t-N\t$job_id\n";
	print SC "\#\$\t-o\t$log_dir/$job_id.log\n";
	print SC "\#\$\t-e\t$log_dir/$job_id.err\n";
	if($job_range ne ""){
		print SC "\#\$\t-t\t$job_range\n";	
	}	
	if($hold_id ne ""){
		print SC "\#\$\t-hold_jid\t$hold_id\n";
	}
	while(<BP>){
		print SC $_;
	}
	close(BP) or die "Can't close $blue_print";
	close(SC) or die "Can't close $script";
}

sub usage{
    print "\n############################\nUsage:\n############################\n";
    print STDERR <<"	_EOT_";
	
perl \$ERDS_dir/erds_pipeline.pl -b \$bam_file -v \$snps_indels_file -o \$output_dir -r \$ref_file OPTIONS
Below is an example when I ran ERDS in sample na12878chgv.
perl /home/mz34/erds1.1/erds_pipeline.pl -b /nfs/seqsata06/ALIGNMENT/BUILD37/GENOME/na12878chgv/combined/na12878chgv_final.bam -v /nfs/seqsata06/ALIGNMENT/BUILD37/GENOME/na12878chgv/combined/na12878chgv.recalSNPs.vcf -o /nfs/seqsata07/Mingfu/na12878chgv -r /nfs/goldsteindata/refDB/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa

	#You must define the following parameters.
	-o: followed by the path of output directory.(-o output_dir)
	-b: followed by specified bam_file. (-b bam_file)
	-v: followed by specified variant_file file called by SNV calling tools. (-v variant_file)
	-r: followed by specified reference file. (-r reference_fasta_file)

	#You can optionally define the following parameters.
	--cluster or -c: if indicated, run cluster. By default single cpu. (-c)
	--samtools: followed by the specified samtools. By default the attached version 1.12. (--samtools path_to_samtools)
	--sd: followed by the SD version of b36, b37 or empty. By default empty.(--sd b36|b37|empty)
	--exp_pct or -p: followed by expected percentage of deletions or duplications. By default 0.0025.(-p value)
	--name: followed by the string of the sample name.(-name sample_name)

	#The following parameters are for advanced users only. You can optionally define them only if you are sure.
	--parameter: followed by specified parameter file. (--parameter parameter_file)
	--del: followed by specified win_size_del. By default 200. (--del win_size_del)
	--dup: followed by specified win_size_dup. By default 1000. (--dup win_size_dup)
	--large: followed by specified length (in bp). Deletions larger than this value will escape from PEM check. By default 50000.(--large integer)
	--pem: followed by specified value for PEM support required. By default 1. (--pem integer)
	--small: followed by specified length (in bp). Deletions smaller than this value will require more PEM support (\$pem+1). By default 1000. (--small integer)
	--scan: followed by specified value for low RD ratio scan. By default 1.4. (--scan value)
	--pairs: followed by specified value for PEM support required for low_pem_pairs. By default 3. (--pairs integer)
	--clip: followed by specified value for clip_spt_pem. By default 2. (--clip integer)
	--soft: followed by specified length for softclip_check. By default 1000. (--soft integer)
	--num_stdd: followed by specified number. Mapping distance larger than mean library size plus \$num_stdd*\$std_dev is regarded as discordant. By default 4. (--num_stdd integer)
	--begin: followed by the step to run. Use for control the program during debuging or rerun. By default 0. (--begin integer) 
	--end: followed by the step to stop. By default 10. (--end integer)

For more information, please visit
http://www.duke.edu/~mz34/erds.htm
	
	_EOT_
	exit(1);
}