data = /path/to/data/directory
results=/path/to/results/directory
tissue= # blood or liver
dnarna= # DNA or RNA
sample= # sample name
bowtie=/share/apps/genomics/bowtie2-2.4.1/bowtie2
samtools=/share/apps/genomics/samtools-1.9/bin/samtools
picard=/share/apps/genomics/picard-2.26.9/picard.jar

#### AAV2 expression ####

# Make output directories
mkdir -p ${results}/${tissue}/${sample}/${dnarna}/bowtie/ ${results}/${tissue}/${sample}/${dnarna}/bam

# Align human-filtered fastqs to AAV2 genome
${bowtie}-build ${data}/genomes/NC_001401.fasta ${data}/genomes/NC_001401

(${bowtie} -x ${data}/genomes/NC_001401 --very-sensitive \
-1 ${data}/fastq/${tissue}/${dnarna}/${sample}_${dnarna}_filtered_1.fastq -2 ${data}/fastq/${tissue}/${dnarna}/${sample}_${dnarna}_filtered_2.fastq \
-S ${results}/${tissue}/${sample}/${dnarna}/bowtie/${sample}_aav2.sam) 2> ${results}/${tissue}/${sample}/${dnarna}/bowtie/${sample}_aav2_log.txt

# Convert to bam
${samtools} view -bS ${results}/${tissue}/${sample}/${dnarna}/bowtie/${sample}_aav2.sam > ${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2.bam

# Sort
${samtools} sort ${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2.bam > ${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_sorted.bam                                                                                                 

# Index
${samtools} index  ${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_sorted.bam                                                                                                         

# Remove duplicates
java -Xmx4g -jar ${picard} MarkDuplicates \
I=${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_sorted.bam \
O=${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_sorted_dedup.bam \
M=${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_picard.dup.metrics REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT 

# Index
java -Xmx4g -jar ${picard}  BuildBamIndex I=${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_sorted_dedup.bam  VALIDATION_STRINGENCY=SILENT

# Make depth files for import into R
${samtools} depth ${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_aav2_sorted_dedup.bam \
> ${results}/${tissue}/${sample}/${dnarna}/bam/${sample}_${dnarna}_aav2_sorted_dedup.depth

#### Transcript assembly ####
 mkdir ${results}/star/${tissue}/${sample}/${dnarna}/ ${results}/stringtie/${tissue}/${sample}/${dnarna}/

${star} --genomeDir ${genome_dir} \
--runThreadN 6 \
--runMode alignReads \
--readFilesIn ${results}/${tissue}/${sample}/${dnarna}/aav_mapped_fastq/${tissue}_${sample}_${dnarna}_1.fastq \
${results}/${tissue}/${sample}/${dnarna}/aav_mapped_fastq/${tissue}_${sample}_${dnarna}_2.fastq \
--outFileNamePrefix ${results}/star/${tissue}/${sample}/${dnarna}/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard

${stringtie} -o ${results}/stringtie/${tissue}/${sample}/${dnarna}/${tissue}_${sample}_${dnarna}_aav2.gtf \
-s 1 ${results}/star/${tissue}/${sample}/${dnarna}/Aligned.sortedByCoord.out.bam

