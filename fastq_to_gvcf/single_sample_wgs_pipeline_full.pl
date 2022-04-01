#!/usr/bin/perl -w
use strict;
use Getopt::Std;
my %opts;
getopt ('s:',\%opts);
my $sample=$opts{"s"};
my $path="/hpf/projects/sbowdin/andyding/yiming";
my $temp_path="/localhd/$ENV{'PBS_JOBID'}";

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

system ("mkdir -p $path/intermediate/$sample");
system ("mkdir -p $path/single_sample_wgs_final/$sample");

# Determine the number of read groups
system ("ls $path/fastq/$sample/*.gz > $path/intermediate/$sample/$sample.fastq_gz_files.txt");
open file1, "<$path/intermediate/$sample/$sample.fastq_gz_files.txt";
my (%rg,%fq);
while (<file1>){
 chomp;
 my @split1=split /\_/,$_;
 $fq{$_}=1;
 foreach my $split1 (@split1){
  if ($split1=~/^readgroup/){
   $rg{$split1}=1;
  }
 }
}
close file1;

my $markduplicate_input="";
foreach my $rg (keys %rg){
 # FASTQ to uBAM
 my $fq1="$path/fastq/$sample/$sample\_$rg\_R1.fastq.gz";
 my $fq2="$path/fastq/$sample/$sample\_$rg\_R2.fastq.gz";
 if (exists $fq{$fq2}){
  system ("java -Xmx32g -jar $tool{'picard'} FastqToSam FASTQ=$fq1 FASTQ2=$fq2 OUTPUT=$temp_path/$sample.$rg.unmapped.bam READ_GROUP_NAME=$sample SAMPLE_NAME=$sample LIBRARY_NAME=$sample\_$rg PLATFORM=Illumina SEQUENCING_CENTER=TCAG");
 }
 else{
  system ("java -Xmx32g -jar $tool{'picard'} FastqToSam FASTQ=$fq1 OUTPUT=$temp_path/$sample.$rg.unmapped.bam READ_GROUP_NAME=$sample SAMPLE_NAME=$sample LIBRARY_NAME=$sample\_$rg PLATFORM=Illumina SEQUENCING_CENTER=TCAG");
 }

 # CollectQualityYieldMetrics
 system ("java -Xms2000m -Xmx3000m -jar $tool{'picard'} CollectQualityYieldMetrics INPUT=$temp_path/$sample.$rg.unmapped.bam OQ=TRUE OUTPUT=$path/intermediate/$sample/$sample.$rg.unmapped.quality_yield_metrics");

 # DragmapAlignment SamToFastqAndDragmapAndMba
 system ("$tool{'dragen-os'} -b $temp_path/$sample.$rg.unmapped.bam -r $tool{'reference'} --interleaved=1 | $tool{'samtools'} view -h -O BAM - > $temp_path/$sample.$rg.aligned.bam");
 system ("java -Dsamjdk.compression_level=2 -Xms1000m -Xmx1000m -jar $tool{'picard'} MergeBamAlignment VALIDATION_STRINGENCY=SILENT EXPECTED_ORIENTATIONS=FR ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_REMOVE=RG ATTRIBUTES_TO_REMOVE=NM ATTRIBUTES_TO_REMOVE=MD ALIGNED_BAM=$temp_path/$sample.$rg.aligned.bam UNMAPPED_BAM=$temp_path/$sample.$rg.unmapped.bam REFERENCE_SEQUENCE=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta PAIRED_RUN=true SORT_ORDER=\"unsorted\" IS_BISULFITE_SEQUENCE=false ALIGNED_READS_ONLY=false CLIP_ADAPTERS=false MAX_RECORDS_IN_RAM=2000000 ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant PROGRAM_RECORD_ID=\"dragen-os\" PROGRAM_GROUP_VERSION=\"dragen-os\" PROGRAM_GROUP_COMMAND_LINE=\"dragen-os -b input_bam -r dragen_reference --interleaved=1\" PROGRAM_GROUP_NAME=\"dragen-os\" UNMAPPED_READ_STRATEGY=COPY_TO_TAG ALIGNER_PROPER_PAIR_FLAGS=true UNMAP_CONTAMINANT_READS=true ADD_PG_TAG_TO_READS=false OUTPUT=$temp_path/$sample.$rg.MergeBamAlignment.bam");

 # CollectUnsortedReadgroupBamQualityMetrics
 system ("java -Xms5000m -Xmx6500m -jar $tool{'picard'} CollectMultipleMetrics INPUT=$temp_path/$sample.$rg.MergeBamAlignment.bam OUTPUT=$path/intermediate/$sample/$sample.$rg ASSUME_SORTED=true PROGRAM=null PROGRAM=CollectBaseDistributionByCycle PROGRAM=CollectInsertSizeMetrics PROGRAM=MeanQualityByCycle PROGRAM=QualityScoreDistribution METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=ALL_READS");
 $markduplicate_input="$markduplicate_input INPUT=$temp_path/$sample.$rg.MergeBamAlignment.bam";
}
print "MarkDuplicates input $markduplicate_input\n";

# MarkDuplicates
system ("java -Dsamjdk.compression_level=2 -Xms32g -jar $tool{'picard'} MarkDuplicates $markduplicate_input OUTPUT=$temp_path/$sample.MarkDuplicates.bam METRICS_FILE=$path/intermediate/$sample/$sample.duplicate_metrics VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=\"queryname\" CLEAR_DT=\"false\" ADD_PG_TAG_TO_READS=false");

# Clean-up, remove read group BAM files
foreach my $rg (keys %rg){
 system ("rm $temp_path/$sample.$rg.MergeBamAlignment.bam $temp_path/$sample.$rg.unmapped.bam $temp_path/$sample.$rg.aligned.bam");
}

# SortSampleBam
system ("java -Dsamjdk.compression_level=2 -Xms32g -Xmx32g -jar $tool{'picard'} SortSam INPUT=$temp_path/$sample.MarkDuplicates.bam OUTPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam SORT_ORDER=\"coordinate\" CREATE_INDEX=true CREATE_MD5_FILE=true MAX_RECORDS_IN_RAM=300000");

# Clean-up, remove MarkDuplicates BAM file
system ("rm $temp_path/$sample.MarkDuplicates.bam");

# CrossCheckFingerprints
system ("java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx3000m -jar $tool{'picard'} CrosscheckFingerprints OUTPUT=$path/intermediate/$sample/$sample.CrosscheckFingerprints HAPLOTYPE_MAP=$tool{'reference'}/Homo_sapiens_assembly38.haplotype_database.txt EXPECT_ALL_GROUPS_TO_MATCH=true INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam LOD_THRESHOLD=-20.0 CROSSCHECK_BY=\"READGROUP\"");

# CollectReadgroupBamQualityMetrics
system ("java -Xms5000m -Xmx6500m -jar $tool{'picard'} CollectMultipleMetrics INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam REFERENCE_SEQUENCE=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta OUTPUT=$path/intermediate/$sample/$sample ASSUME_SORTED=true PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=\"CollectGcBiasMetrics\" METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=READ_GROUP");

# CollectAggregationMetrics
system ("java -Xms5000m -Xmx6500m -jar $tool{'picard'} CollectMultipleMetrics INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam REFERENCE_SEQUENCE=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta OUTPUT=$path/intermediate/$sample/$sample ASSUME_SORTED=true PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=CollectSequencingArtifactMetrics PROGRAM=QualityScoreDistribution PROGRAM=\"CollectGcBiasMetrics\" METRIC_ACCUMULATION_LEVEL=null METRIC_ACCUMULATION_LEVEL=SAMPLE METRIC_ACCUMULATION_LEVEL=LIBRARY");

# CalculateReadGroupChecksum
system ("java -Xms1000m -Xmx3500m -jar $tool{'picard'} CalculateReadGroupChecksum INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam OUTPUT=$path/intermediate/$sample/$sample.CalculateReadGroupChecksum");

# ConvertToCram
system ("$tool{'samtools'} view -C -T $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta $temp_path/$sample.aligned.duplicate_marked.sorted.bam | tee $temp_path/$sample.aligned.duplicate_marked.sorted.cram | md5sum | awk '{print \$1}' > $path/intermediate/$sample/$sample.aligned.duplicate_marked.sorted.cram.md5");
system ("$tool{'samtools'} index $temp_path/$sample.aligned.duplicate_marked.sorted.cram");

# ValidateSamFile
system ("java -Xms32000m -Xmx32000m -jar $tool{'picard'} ValidateSamFile INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.cram OUTPUT=$path/intermediate/$sample/$sample.aligned.duplicate_marked.sorted.cram.ValidateSamFile REFERENCE_SEQUENCE=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta MAX_OUTPUT=1000000000 IGNORE=\"MISSING_TAG_NM\" MODE=VERBOSE IS_BISULFITE_SEQUENCED=false");

# CollectWgsMetrics
system ("java -Xms2000m -Xmx2500m -jar $tool{'picard'} CollectWgsMetrics INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta INCLUDE_BQ_HISTOGRAM=true INTERVALS=$tool{'reference'}/wgs_coverage_regions.hg38.interval_list OUTPUT=$path/intermediate/$sample/$sample.aligned.duplicate_marked.sorted.bam.CollectWgsMetrics USE_FAST_ALGORITHM=true READ_LENGTH=250");

# CollectRawWgsMetrics
system ("java -Xms16000m -jar $tool{'picard'} CollectRawWgsMetrics INPUT=$temp_path/$sample.aligned.duplicate_marked.sorted.bam VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta INCLUDE_BQ_HISTOGRAM=true INTERVALS=$tool{'reference'}/wgs_coverage_regions.hg38.interval_list OUTPUT=$path/intermediate/$sample/$sample.aligned.duplicate_marked.sorted.bam.CollectRawWgsMetrics USE_FAST_ALGORITHM=true READ_LENGTH=250");

# CalibrateDragstrModel
system ("$tool{'gatk'} --java-options \"-Xmx16000m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Dsamjdk.reference_fasta=$tool{'reference'}/Homo_sapiens_assembly38_masked.fasta\" CalibrateDragstrModel -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -I $temp_path/$sample.aligned.duplicate_marked.sorted.bam -str $tool{'reference'}/Homo_sapiens_assembly38.str -O $path/intermediate/$sample/$sample.CalibrateDragstrModel");

# HaplotypeCaller_GATK4_VCF
system ("$tool{'gatk'} --java-options \"-Xmx32000m -Xms32000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\" HaplotypeCaller -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -I $temp_path/$sample.aligned.duplicate_marked.sorted.bam -L $tool{'reference'}/wgs_calling_regions.hg38.interval_list -O $temp_path/$sample.init.g.vcf -contamination 0 -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation --dragen-mode --dragstr-params-path $path/intermediate/$sample/$sample.CalibrateDragstrModel -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 -ERC GVCF");

# DragenHardFilterVcf
system ("$tool{'gatk'} --java-options \"-Xms2000m -Xmx2500m\" VariantFiltration -V $temp_path/$sample.init.g.vcf --filter-expression \"QUAL < 10.4139\" --filter-name \"DRAGENHardQUAL\" -O $temp_path/$sample.filtered.g.vcf");

# Reblock
system ("$tool{'gatk'} --java-options \"-Xms3000m -Xmx3000m\" ReblockGVCF -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -V $temp_path/$sample.filtered.g.vcf -do-qual-approx --floor-blocks -GQB 20 -GQB 30 -GQB 40 -O $path/intermediate/$sample/$sample.g.vcf");

# ValidateVCF
system ("$tool{'gatk'} --java-options \"-Xms32000m -Xmx32000m\" ValidateVariants -V $path/intermediate/$sample/$sample.g.vcf -R $tool{'reference'}/Homo_sapiens_assembly38_masked.fasta -L $tool{'reference'}/wgs_calling_regions.hg38.interval_list -gvcf --validation-type-to-exclude ALLELES --dbsnp $tool{'reference'}/Homo_sapiens_assembly38.dbsnp138.vcf");

# CollectVariantCallingMetrics
system ("java -Xms2000m -Xmx2500m -jar $tool{'picard'} CollectVariantCallingMetrics INPUT=$path/intermediate/$sample/$sample.g.vcf OUTPUT=$path/intermediate/$sample/$sample.g DBSNP=$tool{'reference'}/Homo_sapiens_assembly38.dbsnp138.vcf SEQUENCE_DICTIONARY=$tool{'reference'}/Homo_sapiens_assembly38_masked.dict TARGET_INTERVALS=$tool{'reference'}/wgs_evaluation_regions.hg38.interval_list GVCF_INPUT=true");

# Clean-up
system ("mv $temp_path/$sample.aligned.duplicate_marked.sorted.cram $path/single_sample_wgs_final/$sample");
system ("mv $temp_path/$sample.aligned.duplicate_marked.sorted.cram.crai $path/single_sample_wgs_final/$sample");
system ("mv $path/intermediate/$sample/$sample.aligned.duplicate_marked.sorted.cram.md5 $path/single_sample_wgs_final/$sample");
system ("mv $temp_path/$sample.aligned.duplicate_marked.sorted.bam $path/single_sample_wgs_final/$sample");
system ("mv $temp_path/$sample.aligned.duplicate_marked.sorted.bai $path/single_sample_wgs_final/$sample");
system ("mv $temp_path/$sample.aligned.duplicate_marked.sorted.bam.md5 $path/single_sample_wgs_final/$sample");
system ("mv $path/intermediate/$sample/$sample.g.vcf $path/single_sample_wgs_final/$sample");
system ("$tool{'bgzip'} -f $path/single_sample_wgs_final/$sample/$sample.g.vcf");
system ("$tool{'tabix'} -p vcf $path/single_sample_wgs_final/$sample/$sample.g.vcf.gz");

system ("rm $temp_path/$sample.init.g.vcf $temp_path/$sample.init.g.vcf.idx $temp_path/$sample.filtered.g.vcf $temp_path/$sample.filtered.g.vcf.idx $path/intermediate/$sample/$sample.fastq_gz_files.txt");

exit 2;
