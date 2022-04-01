#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopt ('s:',\%opts);
my $sample=$opts{"s"};
my $path="/hpf/projects/sbowdin/andyding/yiming";

unless ($sample){
 die "Error: sample ID missing, please specify sample ID using -s [sample].\n";
}

my %tool;
open file1, "<$path/wgs_pipeline_config.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 $tool{$split1[0]}=$split1[1];
}
close file1;

my $invcf="";
for (my $i=0; $i<=39; $i++){
 if (length($i)==1){
  $i="0$i";
 }
 $invcf="$invcf INPUT=$path/intermediate/$sample/$sample.chunk$i.filtered.g.vcf ";
}

# MergeVCFs
system ("java -Xms2000m -Xmx2500m -jar $tool{'picard'} MergeVcfs $invcf OUTPUT=$path/intermediate/$sample/$sample.filtered.g.vcf");

# Reblock
system ("$tool{'gatk'} --java-options \"-Xms3000m -Xmx3000m\" ReblockGVCF -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -V $path/intermediate/$sample/$sample.filtered.g.vcf -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 -O $path/intermediate/$sample/$sample.g.vcf");

# ValidateVCF
system ("$tool{'gatk'} --java-options \"-Xms32000m -Xmx32000m\" ValidateVariants -V $path/intermediate/$sample/$sample.g.vcf -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -L $tool{'reference'}/wgs_calling_regions.hg38.interval_list -gvcf --validation-type-to-exclude ALLELES --dbsnp $tool{'reference'}/Homo_sapiens_assembly38.dbsnp138.vcf");

# CollectVariantCallingMetrics
system ("java -Xms2000m -Xmx2500m -jar $tool{'picard'} CollectVariantCallingMetrics INPUT=$path/intermediate/$sample/$sample.g.vcf OUTPUT=$path/intermediate/$sample/$sample.g DBSNP=$tool{'reference'}/Homo_sapiens_assembly38.dbsnp138.vcf SEQUENCE_DICTIONARY=$tool{'reference'}/Homo_sapiens_assembly38_masked.dict TARGET_INTERVALS=$tool{'reference'}/wgs_evaluation_regions.hg38.interval_list GVCF_INPUT=true");

# Clean-up
system ("mv $path/intermediate/$sample/$sample.g.vcf $path/single_sample_wgs_final/$sample");
system ("$tool{'bgzip'} -f $path/single_sample_wgs_final/$sample/$sample.g.vcf");
system ("$tool{'tabix'} -p vcf $path/single_sample_wgs_final/$sample/$sample.g.vcf.gz");

for (my $i=0; $i<=39; $i++){
 if (length($i)==1){
  $i="0$i";
 }
 system ("rm $path/intermediate/$sample/$sample.chunk$i.filtered.g.vcf $path/intermediate/$sample/$sample.chunk$i.filtered.g.vcf.idx");
}
system ("rm $path/intermediate/$sample/$sample.filtered.g.vcf $path/intermediate/$sample/$sample.filtered.g.vcf.idx");

exit 2;
