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

my $events="$sample_dir/$samplename.events";
my $del="$sample_dir/$samplename.del.events";
my $dup="$sample_dir/$samplename.dup.events";
system("cp $dup $events");
system("cat $del>>$events");

sortfile($events);

##VCF
my $vcf_file="$sample_dir/$samplename.erds.vcf";
my $date=date();
my @list=split /\//, $ref_file;
my $ref_version=$list[$#list];
@list=split /\./, $ref_version;
$ref_version=$list[0];
open(VCF,">$vcf_file") or die "Can't open $vcf_file";
print VCF "##fileformat=VCFv4.1\n";
print VCF "##fileDate=$date\n";
print VCF "##reference=$ref_version\n";
print VCF "##assembly=$ref_file\n";
print VCF "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n";
print VCF "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
print VCF "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
print VCF "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print VCF "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print VCF "##ALT=<ID=DUP,Description=\"Duplication\">\n";
print VCF "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">\n";
print VCF "##FORMAT=<ID=REFCN,Number=1,Type=Integer,Description=\"Copy number genotype for reference genome\">\n";
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$samplename\n";

open(EVT,"$events") or die "Can't open $events";
while (my $in=<EVT>) {
	$in=~s/[\r\n]//g;
	@list=split /\s/, $in;
	my $chr=$list[0];
	my $start=$list[1];
	my $END=$list[2];
	my $SVLEN=$list[3];
	my $SVTYPE="DUP";
	if($list[4] eq "DEL"){
		$SVTYPE="DEL";
		$SVLEN=-$SVLEN;
	}
	my $prec=$list[6];
	if($prec==1){
		$prec="PRECISE";
	}
	else{
		$prec="IMPRECISE";
	}
	my $REFCN=$list[7];
	my $CN=$list[8];
	print VCF "$chr\t$start\t.\t.\t<$SVTYPE>\t.\tPASS\t$prec;SVTYPE=$SVTYPE;END=$END;SVLEN=$SVLEN\tREFCN:CN\t$REFCN:$CN\n";
}
close(EVT) or die "Can't close $events";
close(VCF) or die "Can't close $vcf_file";

