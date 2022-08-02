 #### Assemble the reads with no matches to our database used in metaMix ####

 data=/path/to/data/directory
 results/path/to/results/directory
 dnarna= # either DNA or RNA
 tissue= # either liver or blood
 sample= # sample ID

# Make directories
mkdir -p ${results}/${tissue}/${sample}/${dnarna}/${options}spades_blast_removed/output \
${results}/${tissue}/${sample}/${dnarna}/blast_removed/
 
 # # Make a list of read IDs found by BLAST
awk '{print $1}' ${data}/blastn/${tissue}/${dnarna}/${sample}/nucleotide/${sample}_${dnarna}_megaBLAST_norRNA.tab \
> ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_genes.txt

# Select any fasta entries that don't match those gene names and extract those entries from the fasta files
cat ${data}/blastn/${tissue}/${dnarna}/${sample}/blastn_rRNA/NM_paired_${sample}_${dnarna}_filtered_rRNA.fasta \
| paste - - | grep -v -F -f ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_genes.txt \
> ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_removed.txt

awk  'BEGIN {FS=OFS=" "} {print $1}'  ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_removed.txt | \
sed 's/>//' | sed 's/\/1//' | sed 's/\/2//' | sort | uniq -d > ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids.txt

sed 's/$/\/1/' ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids.txt \
> ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids_1.txt

sed 's/$/\/2/' ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids.txt \
> ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids_2.txt

seqtk subseq ${data}/qc_fastq/${tissue}/${dnarna}/${sample}/${sample}_${dnarna}_1.fastq \
${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids_1.txt \
> ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_removed_1.fastq

seqtk subseq ${data}/qc_fastq/${tissue}/${dnarna}/${sample}/${sample}_${dnarna}_2.fastq \
${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_ids_2.txt \
> ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_removed_2.fastq


# Sort data ready for assembly
seqkit sort ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_removed_1.fastq > ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_sorted_1.fastq
seqkit sort ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_blast_removed_2.fastq > ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_sorted_2.fastq


# Assemble using spades
spades.py --${options} \
-o ${results}/${tissue}/${sample}/${dnarna}/${options}spades_blast_removed/output \
-1 ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_sorted_1.fastq \
-2 ${results}/${tissue}/${sample}/${dnarna}/blast_removed/${sample}_sorted_2.fastq