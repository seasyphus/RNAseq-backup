# Class I

## ssh

```
ssh ltian@chemo.dldcc.bcm.edu
```

To analyze the samples, you can use either interactive mode or batch job mode.

### Interactive mode

`qsub -I -l vmem=16gb -d $PWD`

### Batch mode

Here is the head of the job script (**.sh**).

```
#PBS -l nodes=1:ppn=8
#PBS -l vmem=120gb
#PBS -l mem=120gb
#PBS -l walltime=60:00:00
#PBS -m ae
#PBS -M ltian@bcm.edu
#PBS -o /project/zhang/pbs_output
#PBS -e /project/zhang/pbs_output
#PBS -S /bin/bash
#PBS -j oe
#PBS -q long
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
```

### Load modules

You can check the available modules by `module avail`. To load the module, use `module load XXX`.

## Data download

### Genome and Annotation Download

You can download the reference from [GENCODE](https://www.gencodegenes.org/).

```
cd /project/zhang
mkdir Genome_Hai
## Mouse (Release M12, GRCm38.p5)
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/GRCm38.primary_assembly.genome.fa.gz

## Human (Release 25, GRCh38.p7)
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz

gunzip *

mv gencode.vM12.annotation.gtf mm.gtf
mv GRCm38.primary_assembly.genome.fa mm.fa
mv gencode.v25.annotation.gtf hs.gtf
mv GRCh38.primary_assembly.genome.fa hs.fa
```

### Download the RNA-seq Data

The data for BICA project are deposited to [GSE84114](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84114)

The raw data are stored in the [ftp](ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP078/SRP078006). Let's use the first sample (**SRR3747264**, **MCF-7, rep1**) as example.

```
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP078/SRP078006/SRR3747264/SRR3747264.sra

## Change sample name
mv SRR3747264.sra MCF7_rep1.sra

## Use interactive mode
qsub -I -l vmem=16gb
```

You need to load the module called **sratoolkit/2.3.5-2** by typing `module load sratoolkit/2.3.5-2`. But usually, this should have been appended to the end of `~/.bashrc`. Check it by `less ~/.bashrc`.

Then, you can convert the `.sra` to `fastq` by the command below.

```
fastq-dump -I --split-3 MCF7_rep1.sra
```

It will be better to compress the file before mapping/quantification as this saves lots of space.

```
gzip *fastq
```

It is helpful to check the read quality by using `fastqc`.

```
fastqc *fastq.gz
```

Once this quality control is finished, you can download the html reports to your local computer and open the reports using web browsers (Chrome, Firefox). To download the reports, open the terminal window in you local computer and run the commands below.

```
scp ltian@chemo.dldcc.bcm.edu:/project/zhang/Hai_training/sampleDirectory/*html .
```

## Build Reference

We will use the STAR-RSEM pipeline. This is ideal for the RNA-seq novice as **1)** it is pretty simple (one step, from fastq.gz to gene count table); **2)** it can omit the large intermediate files (.bam, .sam). This may require large memory, and it is better to use batch mode.

One very important parameter is `--star-sjdboverhang readLength-1`. If the read length is 75, you should use `--star-sjdboverhang 74`.

You can create an empty text file by typing:

```
nano RSEM_Ref_Build.sh
```

Paste the following codes to nano. 

```
#PBS -l nodes=1:ppn=8
#PBS -l vmem=100gb
#PBS -l mem=100gb
#PBS -l walltime=12:00:00
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

mkdir /project/zhang/Genome_Hai/mm_RSEM
mkdir /project/zhang/Genome_Hai/hs_RSEM
GenomesRef=/project/zhang/Genome_Hai

## Mouse
cd /project/zhang/Genome_Hai/mm_RSEM
rsem-prepare-reference -p 8 --gtf /$GenomesRef/mm.gtf --star --star-sjdboverhang 74 $GenomesRef/mm.fa mm

## Human
cd /project/zhang/Genome_Hai/hs_RSEM
rsem-prepare-reference -p 8 --gtf /$GenomesRef/hs.gtf --star --star-sjdboverhang 74 $GenomesRef/hs.fa hs
```

You can submit the batch script by typing:

```
qsub RSEM_Ref_Build.sh
```

You can check the status of the job by `qstat`. You can delete the job by `qdel jobID`.