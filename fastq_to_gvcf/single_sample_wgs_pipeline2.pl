#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopt ('s:c:',\%opts);
my $sample=$opts{"s"};
my $chunk=$opts{"c"};
my $path="/hpf/projects/sbowdin/andyding/yiming";

unless ($sample){
 die "Error: sample ID missing, please specify sample ID using -s [sample].\n";
}
unless ($chunk){
 die "Error: chunk missing, please specify using -c [chunk] (from 1-40).\n";
}
$chunk--;
if (length($chunk)==1){
 $chunk="0$chunk";
}

my %tool;
open file1, "<$path/wgs_pipeline_config.txt";
while (<file1>){
 chomp;
 my @split1=split /\t/,$_;
 $tool{$split1[0]}=$split1[1];
}
close file1;

system ("mkdir -p $path/intermediate/$sample");
system ("mkdir -p $path/single_sample_wgs_final/$sample");

# HaplotypeCaller_GATK4_VCF
system ("$tool{'gatk'} --java-options \"-Xmx14000m -Xms14000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\" HaplotypeCaller -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -I $path/single_sample_wgs_final/$sample/$sample.aligned.duplicate_marked.sorted.bam -O $path/intermediate/$sample/$sample.chunk$chunk.init.g.vcf -contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation --dragen-mode --dragstr-params-path $path/intermediate/$sample/$sample.CalibrateDragstrModel -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF -L $path/hg38_scattered_intervals_40/00$chunk-scattered.interval_list");

# DragenHardFilterVcf
system ("$tool{'gatk'} --java-options \"-Xms2000m -Xmx2500m\" VariantFiltration -V $path/intermediate/$sample/$sample.chunk$chunk.init.g.vcf --filter-expression \"QUAL < 10.4139\" --filter-name \"DRAGENHardQUAL\" -O $path/intermediate/$sample/$sample.chunk$chunk.filtered.g.vcf");

# Clean-up
system ("rm $path/intermediate/$sample/$sample.chunk$chunk.init.g.vcf $path/intermediate/$sample/$sample.chunk$chunk.init.g.vcf.idx");

exit 2;
