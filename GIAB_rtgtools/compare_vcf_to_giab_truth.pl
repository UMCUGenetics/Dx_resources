#!/usr/bin/perl -w
use strict;
use Getopt::Long;


my ($nist,$sample,$bed,$high_conf,$vcf,$threads,$vcf2,$sample2) = ("v2.19","U175754CFGIAB12878", '/hpc/diaggen/software/production/Dx_tracks/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank.bed', '', '',6,'', '');
my ($true_path,$truevcf);
my $ref = '/hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2/Homo_sapiens.GRCh37.GATK.illumina.SDF';
my $GATK = "/hpc/local/CentOS7/cog/software/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /hpc/diaggen/data/databases/ref_genomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta ";


my $opt_help = "This tools compares a vcf to the GIAB consensus truth set and calculates with RTGtools sensitivity and precision

		Usage:  compare_vcf_to_giab_truth.pl 	
			-vcf [v] vcf file to analyse [$vcf],
			-nist [n] NIST consensus version (v2.19|v3.3.2) [$nist],
			-sample [s] Samplename in test vcf [$sample],
			-vcf2 [v2] vcffile to use as truthset for duplo measurements; overlap (true positives) will be compared to nist,
			-sample2 [s2] samplename in vcf2, 
			-bed bed [b] file to analyse, ie exome [$bed],
			-high_conf [c] bed file with high confident regions [$high_conf],
			-t threads to use during analyses [$threads],
			-h help\n";

my $result = GetOptions("vcf|v=s" => \$vcf,
			"nist|n=s" => \$nist,
			"sample|s=s" => \$sample,
			"vcf2|v2=s" => \$vcf2,
			"sample2|s2=s" => \$sample2,
			"bed|b=s" => \$bed,
			"high_conf|c=s" => \$high_conf,
			"t=i" => \$threads,
			"help|?") or die $opt_help;

die $opt_help unless ($vcf);

if ($vcf2 and !$sample2) { die "need sample2!\n"};

if ($nist eq "v2.19") {
    $true_path ='/hpc/diaggen/data/databases/GIAB/NIST_v2.19/';
    $truevcf = 'GIAB12878_nist2.19_truth.vcf.gz';
    $high_conf = '/hpc/diaggen/data/databases/GIAB/NIST_v2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_noMT.bed';
}elsif ($nist eq "v3.3.2") {
    $true_path ='/hpc/diaggen/data/databases/GIAB/NIST_v3.3.2/';
    $truevcf = 'nistregions_HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz';
    $high_conf = '/hpc/diaggen/data/databases/GIAB/NIST_v3.3.2/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed';
}


if ($vcf2) {
    system "java -Xmx4g -jar $GATK -T SelectVariants -se $sample --excludeNonVariants --removeUnusedAlternates -V $vcf -o $sample\_$vcf";
    system "java -Xmx4g -jar $GATK -T SelectVariants -se $sample2 --excludeNonVariants --removeUnusedAlternates -V $vcf2 -o $sample2\_$vcf2";
    #run rtgtools in duplomode with vcf and vcf2
    my $outputfoldername = "RTG_$nist\_$vcf\_$vcf2";
	    my $command = "java -Xmx4G -jar  /hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2/RTG.jar vcfeval -t $ref -T $threads --baseline=$sample2\_$vcf2 --calls=$sample\_$vcf -o $outputfoldername --all-records";
            system "$command" ;


    #set vcf file to TP vcf from rtgeval
    my $vcf_inclpath = "$outputfoldername/tp.vcf.gz";
    $vcf = "giab_duplo_true_pos.vcf.gz";
    #splitout snps
    system "java -Xmx4g -jar $GATK -T SelectVariants -L $high_conf --excludeNonVariants --removeUnusedAlternates -V $vcf_inclpath -selectType SNP -o SNP_nistregions_$vcf";
    #splitout indels
    system "java -Xmx4G -jar $GATK -T SelectVariants -L $high_conf --excludeNonVariants --removeUnusedAlternates -V $vcf_inclpath -selectType INDEL -o INDEL_nistregions_$vcf";
}else{
    #splitout snps
    print("java -Xmx4g -jar $GATK -T SelectVariants -se $sample -L $high_conf --excludeNonVariants --removeUnusedAlternates -V $vcf -selectType SNP -o SNP_nistregions_$vcf");
    system "java -Xmx4g -jar $GATK -T SelectVariants -se $sample -L $high_conf --excludeNonVariants --removeUnusedAlternates -V $vcf -selectType SNP -o SNP_nistregions_$vcf";
    #splitout indels
    print("java -Xmx4G -jar $GATK -T SelectVariants -se $sample -L $high_conf --excludeNonVariants --removeUnusedAlternates -V $vcf -selectType INDEL -o INDEL_nistregions_$vcf");
    system "java -Xmx4G -jar $GATK -T SelectVariants -se $sample -L $high_conf --excludeNonVariants --removeUnusedAlternates -V $vcf -selectType INDEL -o INDEL_nistregions_$vcf";
}

#Run RTGTOOLS
my @passedonly = ('no','yes');
my @exomeonly = ('no','yes');
my $outfile = "rtg_results_$vcf\_$nist.txt";

foreach my $passedonly (@passedonly) {
    foreach my $exomeonly (@exomeonly) {
	foreach my $vcf (<*_nistregions*.gz>) {
	    next if -d $vcf;
	    my $outputfoldername = "RTG_exome-$exomeonly\_passed-$passedonly\_$nist\_$vcf";
	    my $type = $1 if $vcf=~ /(SNP|INDEL)_/;
	    
	    system "echo $type $outputfoldername >>$outfile ";
	    my $command = "java -Xmx4G -jar  /hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2/RTG.jar vcfeval -t $ref -T $threads --baseline=$true_path$type\_$truevcf --calls=$vcf -o $outputfoldername";
	    
	    $command .= " --all-records " if ($passedonly eq "no");
	    $command .= " --bed-regions=$bed " if ($exomeonly eq "yes");
            system "$command >> $outfile" ;
	}
    }
}

# parse rtgeval results into table
open IN, $outfile;
open OUT2, ">parsed_$outfile";
my @header = ('Type', 'Sample','Threshold','True-pos','False-pos','False-neg','Precision','Sensitivity','F-measure'); 
print OUT2 join("\t",@header),"\n";
while (my $line=<IN>) {
    chomp($line);
    my @OUT;
    my @line = split(/\s+/,$line);
    next if !$line[1];
    if ($line[1] =~ /SNP|INDEL/) {
	#print OUT2 join("\t",@line),"\t";
        print OUT2 join("\t",@line);
    }
    elsif ($line[1] eq 'None') { #get values without threshold
	print OUT2 join("\t",@line),"\n";
    }
}
