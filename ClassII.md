# Class II

## *in silico* Sorting

In the [BICA](https://github.com/lintian0616/bica) project, [Xenome](https://academic.oup.com/bioinformatics/article/28/12/i172/269972/Xenome-a-tool-for-classifying-reads-from-xenograft) was used to separate human cancer cell reads from the mouse stroma cell reads.

### Build Xenome Reference

Before sorting the reads, we need to build the xenome reference index from the **graft** (human) and **host** (mouse) reference sequences.

```
cd /project/zhang/Hai_training
mkdir XenomeRef
```

The index building will take more than several hours, and you need to submit the batch job script.

```
#PBS -l nodes=1:ppn=8
#PBS -l vmem=80gb
#PBS -l mem=80gb
#PBS -l walltime=10:00:00
#PBS -m ae
#PBS -M haiw@bcm.edu
#PBS -o /project/zhang/pbs_output
#PBS -e /project/zhang/pbs_output
#PBS -S /bin/bash
#PBS -j oe
#PBS -q medium
#########################
module load anaconda/2.5.0
module load bedtools/2.17.0
module load bedops/2.4.2
module load bismark/0.15.0
module load boost/1.60.0
module load bowtie/1.1.0
module load bowtie2/2.1.0
module load bwa/0.7.13
module load bseqc/master
module load bsmap/2.90
module load circos/0.64
module load cufflinks/2.2.1
module load fastqc/0.11.2
module load libgtextutils/0.7 fastx-toolkit/0.0.14
module load htseq/2.7
module load macs/2.1.0
module load moabs/1.3.2
module load java/1.8.0u71 picard/2.1.0
module load R/3.2.3
module load rsem/1.2.28
module load samtools/0.1.19
module load sratoolkit/2.3.5-2
module load star/2.5.2b
module load tophat/2.0.10
module load homer/4.8

cd /project/zhang/Hai_training/XenomeRef
genomeDir=/project/zhang/Hai_training/Genome_training

/home/ltian/xenome-1.0.1-r/xenome index -M 64 -T 8 -P kmer25 --tmp-dir ./ -v -H $genomeDir/mm.fa -G $genomeDir/hs.fa
```

### Sorting the Reads

Let's test one sample (e.g., **Orth_rep1**) first. Let's use interactive job for the test.

```
qsub -I -l vmem=32gb

cd /project/zhang/Hai_training/Orth_rep1

XenomeIdx=/project/zhang/Hai_training/XenomeRef/kmer25
/home/ltian/xenome-1.0.1-r/xenome classify -P $XenomeIdx --pairs -i *_1.fastq.gz -i *_2.fastq.gz --graft-name human --host-name mouse
```

You can write a **for** loop to go through all the samples.

```
#PBS -l nodes=1:ppn=8
#PBS -l vmem=80gb
#PBS -l mem=80gb
#PBS -l walltime=10:00:00
#PBS -m ae
#PBS -M haiw@bcm.edu
#PBS -o /project/zhang/pbs_output
#PBS -e /project/zhang/pbs_output
#PBS -S /bin/bash
#PBS -j oe
#PBS -q medium
#########################
module load anaconda/2.5.0
module load bedtools/2.17.0
module load bedops/2.4.2
module load bismark/0.15.0
module load boost/1.60.0
module load bowtie/1.1.0
module load bowtie2/2.1.0
module load bwa/0.7.13
module load bseqc/master
module load bsmap/2.90
module load circos/0.64
module load cufflinks/2.2.1
module load fastqc/0.11.2
module load libgtextutils/0.7 fastx-toolkit/0.0.14
module load htseq/2.7
module load macs/2.1.0
module load moabs/1.3.2
module load java/1.8.0u71 picard/2.1.0
module load R/3.2.3
module load rsem/1.2.28
module load samtools/0.1.19
module load sratoolkit/2.3.5-2
module load star/2.5.2b
module load tophat/2.0.10
module load homer/4.8

XenomeIdx=/project/zhang/Hai_training/XenomeRef/kmer25
cd /project/zhang/Hai_training/samples

for f in ./*;
  do
    cd "$f" && echo Entering into $f and begin analysis
    /home/ltian/xenome-1.0.1-r/xenome classify -P $XenomeIdx --pairs -i *_1.fastq.gz -i *_2.fastq.gz --graft-name ${PWD##*/}_human --host-name ${PWD##*/}_mouse
    cd ..
  done;
```

We only use human specific reads or mouse specific reads for downstream analysis. One problem of **Xenome 1.0** is that the `@` and `+` symbols in fastq file are missing.

```
Orth_rep1.15.1 15 length=76
CTATTATCATTCCCACTTTACAATGAGAAAATTGAGTTTCTTTCCACTCAACCACAACCCCCAGATGAGGAAGGAG
Orth_rep1.15.1 15 length=76
AAAAAEEEEAEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEE6EEEAEAAEEEEEEEEEEEEEEEEEEE
```

You can use the custom python code `XenomeCorrect.py`, first upload that python script file to the home directory.

```
scp XenomeCorrect.py bcm:~
```

You can use the command below to fix this problem and gzip the file simutaneously. This may require a little bit higher memory, and you can run `qsub -I -l vmem=8gb -d $PWD`.

```
python ~/XenomeCorrect.py Orth_rep1_human_1.fastq | gzip > Orth_rep1_human_c_1.fastq.gz
```

You can glimpse the file by the command:

```
zcat Orth_rep1_human_c_1.fastq.gz | head
```

## Read Mapping and Quantification

Since we are using **STAR -> RSEM** pipeline, we will skip the big **sam** or **bam** files and obtain the counting matrix table directly.

```
#PBS -l nodes=1:ppn=8
#PBS -l vmem=64gb
#PBS -l mem=64gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M haiw@bcm.edu
#PBS -o /project/zhang/pbs_output
#PBS -e /project/zhang/pbs_output
#PBS -S /bin/bash
#PBS -j oe
#PBS -q short
#########################
module load anaconda/2.5.0
module load bedtools/2.17.0
module load bedops/2.4.2
module load bismark/0.15.0
module load boost/1.60.0
module load bowtie/1.1.0
module load bowtie2/2.1.0
module load bwa/0.7.13
module load bseqc/master
module load bsmap/2.90
module load circos/0.64
module load cufflinks/2.2.1
module load fastqc/0.11.2
module load libgtextutils/0.7 fastx-toolkit/0.0.14
module load htseq/2.7
module load macs/2.1.0
module load moabs/1.3.2
module load java/1.8.0u71 picard/2.1.0
module load R/3.2.3
module load rsem/1.2.28
module load samtools/0.1.19
module load sratoolkit/2.3.5-2
module load star/2.5.2b
module load tophat/2.0.10
module load homer/4.8

cd /project/zhang/Hai_training/BICA_sample/Orth_rep1
RSEM_Idx=/project/zhang/Hai_training/Genome_training/hs_RSEM/hs

rsem-calculate-expression --star -p 8 --star-gzipped-read-file --no-bam-output --calc-ci --estimate-rspd --quiet --phred33-quals --paired-end *_human_c_1.fastq.gz *_human_c_2.fastq.gz $RSEM_Idx ${PWD##*/}_human_RSEM
```

If you want to output **bam** file for genome viewer (e.g., [IGV](http://software.broadinstitute.org/software/igv/)), you just need to replace `--no-bam-output` with `--output-genome-bam`.

You can also run a for-loop to analyze the samples continuously.

```
#PBS -l nodes=1:ppn=8
#PBS -l vmem=64gb
#PBS -l mem=64gb
#PBS -l walltime=3:00:00
#PBS -m ae
#PBS -M haiw@bcm.edu
#PBS -o /project/zhang/pbs_output
#PBS -e /project/zhang/pbs_output
#PBS -S /bin/bash
#PBS -j oe
#PBS -q short
#########################
module load anaconda/2.5.0
module load bedtools/2.17.0
module load bedops/2.4.2
module load bismark/0.15.0
module load boost/1.60.0
module load bowtie/1.1.0
module load bowtie2/2.1.0
module load bwa/0.7.13
module load bseqc/master
module load bsmap/2.90
module load circos/0.64
module load cufflinks/2.2.1
module load fastqc/0.11.2
module load libgtextutils/0.7 fastx-toolkit/0.0.14
module load htseq/2.7
module load macs/2.1.0
module load moabs/1.3.2
module load java/1.8.0u71 picard/2.1.0
module load R/3.2.3
module load rsem/1.2.28
module load samtools/0.1.19
module load sratoolkit/2.3.5-2
module load star/2.5.2b
module load tophat/2.0.10
module load homer/4.8

cd /project/zhang/Hai_training/BICA_sample
RSEM_Idx=/project/zhang/Hai_training/Genome_training/hs_RSEM/hs

for f in ./*;
  do
    cd "$f" && echo Entering into $f and begin analysis
    python ~/XenomeCorrect.py ${PWD##*/}_human_1.fastq | gzip > ${PWD##*/}_human_c_1.fastq.gz
    python ~/XenomeCorrect.py ${PWD##*/}_human_2.fastq | gzip > ${PWD##*/}_human_c_2.fastq.gz
    rsem-calculate-expression --star -p 8 --star-gzipped-read-file --no-bam-output --calc-ci --estimate-rspd --quiet --phred33-quals --paired-end ${PWD##*/}_human_c_1.fastq.gz ${PWD##*/}_human_c_2.fastq.gz $RSEM_Idx ${PWD##*/}_human_RSEM
    cd ..
  done;
```

## Data Normalization

We can use `rsem-generate-data-matrix` to merge individual RSEM samples. Basically, there are three counting outputs: **expected_count**, **TPM**, **FPKM**, **posterior_mean_count**, **pme_TPM**, **pme_FPKM**, the output table will merge the values of **expected_count**.

* If you want to compare the same gene across samples, you need to do some type of normalization, which I suggest the median normalization proposed by `DESeq2`. You can use the merged output table;
* If you want to compare genes within a sample, you do not need normalization. You can use **TPM** (Transcripts Per Kilobase Million) or **FPKM** (Fragments Per Kilobase Million);
* The posterior mean estimates (**posterior\_mean_count**, **pme_TPM**, **pme_FPKM**) are calculated using the Collapsed Gibbs sampler implemented as part of RSEM.

```
cd /project/zhang/Hai_training/BICA_sample

## Human Cancer Cell Gene Expression
rsem-generate-data-matrix */*human_RSEM.genes.results > BICA_human.count.txt

## Mouse Stromal Cell Gene Expression
rsem-generate-data-matrix */*mouse_RSEM.genes.results > BICA_mouse.count.txt
```

After obtaining the merged matrix table, you can transfer the merged table to Cyverse istance. Here is the example code.

```
scp BICA_human.count.txt lintian0616@128.196.64.106:~
```

We then use `DEseq2` to normalize the data in [R programming](https://www.r-project.org/).

```{r}
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)
human_count <- read.table("BICA_human.count.txt", header=TRUE, sep="\t", row.names=1)
human_count.ann <- data.frame(type=c("BICA", "IVBL", "2D", "2D", "Orth"), row.names=colnames(human_count))
human_count.dds <- DESeqDataSetFromMatrix(round(data.matrix(human_count)), colData=human_count.ann, design=~ type)
human_count.dds <- DESeq(human_count.dds)
sizeFactors(human_count.dds) ## This is the read coverage/sequence depth for each sample
human_count_rld <- rlog(human_count.dds) ## We use regularized log transformation (rld)
human_count_rld_gene <- assay(human_count_rld)
```

Then, we need to convert the **ENSEMBL IDs** to **Gene Symbols**, which can be done using `biomaRt` package.

```{r}
library(biomaRt)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
rownames(human_count_rld_gene) <- sub("\\..+", "", rownames(human_count_rld_gene))
map <- getBM(mart=ensembl, attributes=c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=rownames(human_count_rld_gene))
map <- map[map$hgnc_symbol!="", ]
genenames <- unique(map$hgnc_symbol[map$ensembl_gene_id %in% rownames(human_count_rld_gene)])

human_count_gene <- matrix(NA, nrow=length(genenames), ncol=ncol(human_count_rld_gene))
rownames(human_count_gene) <- genenames
colnames(human_count_gene) <- rownames(human_count.ann)

for(g in genenames) {
  p <- as.character(map$ensembl_gene_id[map$hgnc_symbol==g])
  if(length(p)==1) {
    human_count_gene[g, ] <- human_count_rld_gene[p, ]
  }
  else{
    human_count_gene[g, ] <- apply(human_count_rld_gene[p, ], 2, max)
  }
}
```

## Download Data from Basespace

Go to the illumina account and log in. Click the project folder, and download using the application.

You will see 8 files for each pair-end sample. You can merge the 4 channels using the command below. For example, if the data are stored in `/Volumes/LaCie/Basespace/`.

```
cd /Volumes/LaCie/Basespace/BICA_04252016/BICA-29983961/
## Go to the 1st sample:
cd 042216_01-35383443
cat *L001_R1* *L002_R1* *L003_R1* *L004_R1* > BICA1_S1_R1.fastq.gz
cat *L001_R2* *L002_R2* *L003_R2* *L004_R2* > BICA1_S1_R2.fastq.gz
```