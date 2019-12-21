# Genome Assembly

## Using all reads and *SPAdes*

- Reads in *rstudio.iu.edu.tr* server, folder `/home/bioinfo/reads/4-F20-96_S2_L001_R?_001.fastq.gz`
- quality control using *fastqc*
  - put the results in `public_html` folder, so you can see them on-line
- clean and trim using some tool like *trimmomatic*
  - Suggest some alternative
- quality control again using *fastqc*. Is it better?
- assemble using *SPAdes*
- visualize using *Bandage*
- Calculate N50

###### Quality control

```bash
fastqc /home/bioinfo/reads/4-F20-96_S2_L001_R* -o outputfastqc/
```

```bash
cp *html /home/eskutkaan/public_html/
```

Copied to the public_html, we can view through the link;

```bash
rstudio.iu.edu.tr/~eskutkaan/
```

```bash
TrimmomaticPE /home/eskutkaan/finalhomework/fastqs/4-F20-96_S2_L001_R1_001.fastq.gz /home/eskutkaan/finalhomework/fastqs/4-F20-96_S2_L001_R2_001.fastq.gz trimmed_1.fastq trimmed_2.fastq ILLUMINACLIP: /home/eskutkaan/TruSeq3-PE.fa:3:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

You can view our trimmed fastqs' quality files, named as trimmed_?_fastqc.html

The trimmed fastqs are better than the older ones as you can look for quality reports.

After looking for some alternatives i just found that the trimmomatic is the best for now.

```bash
spades -1 4-F20-96_S2_L001_R1_001.fastq.gz -2 4-F20-96_S2_L001_R2_001.fastq.gz -o .
```

We can use the command above to assemble with SPAdes.

## Using *phrap*

- Convert all reads from *fastq* format to *fasta* and *qual*. It can be done using Python
- assemble using *phrap*
- visualize using *consed*
- calculate N50

We can convert them to fasta files with the pyfastx.

```bash
pyfastx fq2fa 4-F20-96_S2_L001_R1_001.fastq.gz -o .
pyfastx fq2fa 4-F20-96_S2_L001_R2_001.fastq.gz -o .
```

```python
from Bio import SeqIO
SeqIO.convert("4-F20-96_S2_L001_R1_001.fastq", "fastq", "4-F20-96_S2_L001_R1_001.fastq.fa.qual", "qual")
SeqIO.convert("4-F20-96_S2_L001_R2_001.fastq", "fastq", "4-F20-96_S2_L001_R2_001.fastq.fa.qual", "qual")
```

```bash
phrap 4-F20-96_S2_L001_R1_001.fastq.fa
phrap 4-F20-96_S2_L001_R2_001.fastq.fa
```

```bash
pyfastx info 4-F20-96_S2_L001_R1_001.fastq.fa
```

```bash
Sequence counts: 125564
Total bases: 27633180
GC content: 51.02%
A counts: 6794775
C counts: 7061046
G counts: 7034780
N counts: 6755
T counts: 6735824
Mean length: 220.07
Median length: 250.00
Max length: 251
Min length: 35
N50, L50: 250, 55058
length >= 1000: 0
```

```bash
pyfastx info 4-F20-96_S2_L001_R2_001.fastq.fa
```

```
Sequence counts: 125564
Total bases: 27727284
GC content: 51.08%
A counts: 6712685
C counts: 7111917
G counts: 7047574
N counts: 6699
T counts: 6848409
Mean length: 220.82
Median length: 250.00
Max length: 251
Min length: 35
N50, L50: 251, 55234
length >= 1000: 0
```

## Using reference sequence

- Terje suggested that this primer could be similar to NC_0021750. Test this hypothesis by aligning all reads to this sequence.
  - Use *bowtie*, *bowtie2*, *bwa mem* and *bwa aln*
  - all these produce SAM files
  - are all results the same? How can you compare them?
  - Which alignment is “the best”
- prepare *fastq* files containing only the reads that align to the plasmid
- assemble these reads using *phrap*

```bash
efetch -db nuccore -format fasta -id NC_025175 > NC_025175.fna
```

Firstly, we will get our reference fasta to map our reads with the command above.

```bash
bwa index NC_025175.fna
```

```bash
bwa mem NC_025175.fna 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > mem.sam
```

We can compare the results from these algorithms with the looking for sam file in text mode with the command;

```bash
less mem.sam
```

They look the same except the minimal structure differences.

I think bwa mem is the best alignment, because the other sams were smaller than the mem one. And didn't look so obvious to me.

```bash
samtools view mem.sam -b -o mem.bam
```

With the command above we can convert our sam file to bam file.

```bash
samtools sort mem.bam -o sorted_mem.bam | samtools index sorted_mem.bam
```

With the command above we can sort and index our bam.

Position Flag 4 indicates an unmapped read in samtools view options. Hence,

```bash
samtools view -F4 -h sorted_mem.bam > mem.mapped.sam
samtools view -f4 -h sorted_mem.bam > mem.unmapped.sam
```

We are going to use mapped sam file to produce fastq.

```bash
samtools view mem.mapped.sam -b -o mem.mapped.bam | samtools fastq mem.mapped.bam -1 mapped1.fastq | samtools fastq mem.mapped.bam -2 mapped2.fastq
```

We will convert our fastqs to fasta and qual.

```python
from Bio import SeqIO
SeqIO.convert("mapped1.fastq", "fastq", "mapped1.fasta", "fasta")
SeqIO.convert("mapped1.fastq", "fastq", "mapped1.qual", "qual")
SeqIO.convert("mapped2.fastq", "fastq", "mapped2.fasta", "fasta")
SeqIO.convert("mapped2.fastq", "fastq", "mapped2.qual", "qual")
```

```bash
mv mapped1.qual mapped1.fasta.qual
mv mapped2.qual mapped2.fasta.qual
```

Assemble these reads using *phrap*

```bash
phrap mapped1.fasta
phrap mapped2.fasta
```

## Analysis of alignment

- what is the coverage of each nucleotide of the reference plasmid?
  - are there any missing regions? What genes are there?
  - are there any regions with atypical coverage? What genes are there?
- Which reads are
  - only in the plasmid
  - only in the chromosome
  - in both plasmid and chromosome

## using Mind-the-gap

[https://gatb.inria.fr/software/mind-the-gap/](https://gatb.inria.fr/software/mind-the-gap/)

# Multiple sequence alignment

- Microbiota genes
- clustal
- t-coffee

```bash
epos -input microbiota_genes/accessions.txt -db nuccore  
```

```bash
efetch -db nuccore -format fasta > phylotypes.fna
```

```bash
clustalw -infile=phylotypes.fna -output=DIFFERENT_OUTPUTS -outfile=tree
```

now we have .aln .dnd .nxs files for phylotypes.fna, .aln for alignments but not fasta, .dnd for text tree format, .nxs for NEXUS format.

```bash
clustalx
```

clustalx for graphical interface to interpret our alignments.

```bash
t_coffee phylotypes.fna
```

It gives phylotypes.html file and we will move it to public_html folder to view.

```bash
rstudio.iu.edu.tr/~eskutkaan/reyhanay
```

We can reach the html file through the link.

```bash
jalview
```

We can use jalview to make a tree, too.
