############ Analyssi of Nanopore data ############

### Add paths to data and software
data=
btru=
minimap2=
aav2_rc=
human_genome=
results=
samtools=

sample=samplename

mkdir -p ${results}/${sample}/aav2_reads

############ Map raw reads to AAV2 genome ############

### Map to reverse complement of AAV2 genome
${minimap2} -ax map-ont -t 4 ${aav2_rc}.mmi \
${btru}/Nanopore/${sample}/fastq_pass/* \
${data}/${sample}/fastq_pass/* \
${samtools} view -bh - > ${results}/${sample}/${sample}_aav2_rc_2.bam

### Extract aligned reads and convert to fastq
${samtools} view -bF 4 -h ${results}/${sample}/${sample}_aav2_rc.bam |
${samtools} fastq - \
> ${results}/${sample}/${sample}_aav2_rc.fastq

############ Map reads that mapped to AAV2 to human genome ############

### Map AAV2 reads to humangenome
${minimap2} -ax map-ont ${human_genome}.mmi \
${results}/${sample}/${sample}_aav2_rc.fastq |
${samtools} view -bh - > ${results}/${sample}/${sample}_aav2_human.bam

### Extract aligned reads and convert to fastq
${samtools} view -bF 4 -h ${results}/${sample}/${sample}_aav2_human.bam |
${samtools} fastq - \
> ${results}/${sample}/${sample}_aav2_human.fastq

### Extract unaligned reads and convert to fastq
${samtools} view -bf 4 -h ${results}/${sample}/${sample}_aav2_human.bam |
${samtools} fastq - \
> ${results}/${sample}/${sample}_aav2_nohuman.fastq

############ Run Kraken on the raw reads ############
${kraken}/kraken2 --threads 4 --db ${kraken_db} \
${data}/${sample}/${sample}.fastq.gz \
--report ${results}/${sample}/${sample}_kraken.report \
> ${results}/${sample}/${sample}_kraken.txt

############ Analysis of the AAV2 reads ############

### Get the kraken AAV2 read IDs and remove any human that are that length
grep -P "\t10804\t" ${results}/${sample}/${sample}_aav2_kraken.txt |
grep -vP "\t9606\t" \
> ${results}/${sample}/${sample}_aav2_kraken_only.txt

awk '{print $2}' ${results}/${sample}/${sample}_aav2_kraken_only.txt \
> ${results}/${sample}/${sample}_aav2_kraken_reads.txt

### Get the minimap read IDs
awk 'NR % 4 == 1' ${results}/${sample}/${sample}_aav2_nohuman.fastq |
sed 's/@//' > ${results}/${sample}/${sample}_aav2_minimap_nohuman_reads.txt

### Get reads classified as AAV2 by both minimap and kraken
grep -F -f ${results}/${sample}/${sample}_aav2_minimap_nohuman_reads.txt \
${results}/${sample}/${sample}_aav2_kraken_reads.txt | sort -u \
> ${results}/${sample}/${sample}_aav2_reads_both.txt

grep -F -f ${results}/${sample}/${sample}_aav2_reads_both.txt -A 3 \
${results}/${sample}/${sample}_aav2_nohuman.fastq |
sed 's/^--$//' | awk 'NF' \
> ${results}/${sample}/${sample}_aav2_both.fastq

rm ${results}/${sample}/${sample}_aav2_both.fastq
while read line; do
  grep "$line" -A 3 -m 1 ${results}/${sample}/${sample}_aav2_nohuman.fastq \
  >> ${results}/${sample}/${sample}_aav2_both.fastq
done < ${results}/${sample}/${sample}_aav2_reads_both.txt

### Split resulting file into individual reads
mkdir aav2_reads
split -l 4 --numeric-suffixes --suffix-length=3 ${results}/${sample}/${sample}_aav2_both.fastq \
${results}/${sample}/aav2_reads/AAV2_read_

### Convert to fasta
for file in ${results}/${sample}/aav2_reads/*; do

    cat ${file} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${file}.fasta

done

### Move fastqs to a separate file
mkdir ${results}/${sample}/aav2_reads/fastq
find ${results}/${sample}/aav2_reads/ -type f ! -iname "*.fasta" -maxdepth 1 -exec mv -t ${results}/${sample}/aav2_reads/fastq {} \;

### Get length of reads
for file in ${results}/${sample}/aav2_reads/AAV2_read*.fasta; do

    grep -v '>' ${file} | wc -m >> ${results}/${sample}/aav2_reads/length.txt

done

### Calculate total bases and N50
sort -nr ${results}/${sample}/aav2_reads/length.txt > ${results}/${sample}/aav2_reads/length_sorted.txt

total_length=$(awk '{SUM+=$1}END{print SUM}' ${results}/${sample}/aav2_reads/length_sorted.txt)

echo Total bases $total_length > ${results}/${sample}/aav2_reads/stats.txt

cumulative_length=0

N50=0

while read length; do
  cumulative_length=$((cumulative_length + length))

  if [ $cumulative_length -ge $((total_length / 2)) ]; then
    N50=$length
    break
  fi
done < ${results}/${sample}/aav2_reads/length_sorted.txt

echo N50 $N50 >> ${results}/${sample}/aav2_reads/stats.txt

### Make dot plots (local machine)
mkdir plots
for file in *.fasta; do
  redotable --window 20 ${file} \
  ~/mnt/BTRU-scratch/sarah/data/hepatitis_preprocessed/genomes/NC_001401.fasta \
  plots/${file}_redotable.png
done
