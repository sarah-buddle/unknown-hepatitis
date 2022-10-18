#### Nanopore integration ####

for sample in sample1 sample2 sample4 sample5; do

# See if any of the insertions in human genome map to AAV2

${minimap2} -d ${aav2_genome}.mmi ${aav2_genome}.fasta

grep 'SVTYPE=INS' ${results}/epi2me_sv/${sample}/${sample}_sv.vcf | grep 'SVLEN=[1-5][0-9][0-9][0-9]' \
> ${results}/epi2me_sv/${sample}/${sample}_insertions.vcf

awk '{print ">" $2 "\n" $5}' ${results}/epi2me_sv/${sample}/${sample}_insertions.vcf \
> ${results}/epi2me_sv/${sample}/${sample}_insertions.fasta

${minimap2} -ax map-ont --secondary-no ${aav2_genome}.mmi ${results}/epi2me_sv/${sample}/${sample}_insertions.fasta \
> ${results}/epi2me_sv/${sample}/${sample}_insertions.sam

samtools fasta ${results}/epi2me_sv/${sample}/${sample}_insertions.sam -F 4 > ${results}/epi2me_sv/${sample}/${sample}_aav2.fasta

# See if any of the insertions in human genome map to HHV6

${minimap2} -d ${hhv6_genome}.mmi ${hhv6_genome}.fasta

${minimap2} -ax map-ont --secondary=no ${hhv6_genome}.mmi ${results}/epi2me_sv/${sample}/${sample}_insertions.fasta \
> ${results}/epi2me_sv/${sample}/${sample}_hhv6_insertions.sam

samtools fastq ${results}/epi2me_sv/${sample}/${sample}_hhv6_insertions.sam -F 4 \
> ${results}/epi2me_sv/${sample}/${sample}_hhv6_insertions.fastq

# See if any of the reads WIMP classified as AAV2 have human elements 
grep -i 'adeno-associated virus - 2' ${results}/wimp/${sample}/${sample}_wimp.csv | awk -F, '{print $2}' \
> ${results}/wimp/${sample}/aav2_reads.txt

fgrep -F -f ${results}/wimp/${sample}/aav2_reads.txt ${data}/${sample}/${sample}_all.fastq -A 3 | 
grep -v -- "^--$" \
> ${results}/wimp/${sample}/aav2_reads.fastq

${minimap2} -ax map-ont --secondary=no ${human_genome}.mmi ${results}/wimp/${sample}/aav2_reads.fastq \
> ${results}/wimp/${sample}/aav2_reads_human.sam

samtools fastq ${results}/wimp/${sample}/aav2_reads_human.sam -F 4 \
> ${results}/wimp/${sample}/aav2_reads_human.fastq

# BLAST to human genome
${blastn}/makeblastdb -in ${human_genome}.fa -parse_seqids -title "human_genome" -dbtype nucl

${blastn}/blastn -db ${human_genome}.fa -query ${results}/wimp/${sample}/aav2_reads.fastq \
-out ${results}/wimp/${sample}/aav2_reads_human_blast.txt -outfmt 6

# See if any of the unclassfied reads map to the AAV2 genome
grep ',Unclassified,' ${results}/wimp/${sample}/${sample}_wimp.csv | awk -F, '{print $2}' \
> ${results}/wimp/${sample}/unclassified_reads.txt

fgrep -F -f ${results}/wimp/${sample}/unclassified_reads.txt ${data}/${sample}/${sample}_all.fastq -A 3 | 
grep -v -- "^--$" \
> ${results}/wimp/${sample}/unclassified_reads.fastq

${minimap2} -ax map-ont --secondary=no ${aav2_genome}.mmi ${results}/wimp/${sample}/unclassified_reads.fastq \
> ${results}/wimp/${sample}/unclassified_reads_aav2.sam

samtools fastq -F 4 ${results}/wimp/${sample}/unclassified_reads_aav2.sam \
> ${results}/wimp/${sample}/unclassified_reads_aav2.fastq

# BLAST unclassified reads to the AAV2 genome
${blastn}/makeblastdb -in ${aav2_genome}.fasta -parse_seqids -title "aav2_genome" -dbtype nucl

${blastn}/blastn -db ${aav2_genome}.fasta -query ${results}/wimp/${sample}/unclassified_reads.fastq \
-out ${results}/wimp/${sample}/unclassified_reads_aav2_blast.txt -outfmt 6

# See if any HHV6 reads map to human genome
grep -i 'Human betaherpesvirus 6B' ${results}/wimp/${sample}/${sample}_wimp.csv | awk -F, '{print $2}' \
> ${results}/wimp/${sample}/hhv6_reads.txt

fgrep -F -f ${results}/wimp/${sample}/hhv6_reads.txt ${data}/${sample}/${sample}_all.fastq -A 3 | 
grep -v -- "^--$" \
> ${results}/wimp/${sample}/hhv6_reads.fastq

${minimap2} -ax map-ont --secondary=no ${human_genome}.mmi ${results}/wimp/${sample}/hhv6_reads.fastq \
> ${results}/wimp/${sample}/hhv6_reads_human.sam

samtools fastq ${results}/wimp/${sample}/hhv6_reads_human.sam -F 4 \
> ${results}/wimp/${sample}/hhv6_reads_human.fastq

# See if any unclassified reads map to HHV6
${minimap2} -ax map-ont --secondary=no ${hhv6_genome}.mmi ${results}/wimp/${sample}/unclassified_reads.fastq \
> ${results}/wimp/${sample}/unclassified_reads_hhv6.sam

samtools fastq -F 4 ${results}/wimp/${sample}/unclassified_reads_hhv6.sam \
> ${results}/wimp/${sample}/unclassified_reads_hhv6.fastq

done




