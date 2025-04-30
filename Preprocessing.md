# Raw reads preprocessing

For this step we used 2 TB of memmory SSD, 64 GB DDR4, and an Intel 9 twelve-core processor.  

## Step 1: Setup directories and organize input data
Once installed the bioinformatic tools, downloaded and verificated all the samples put all of them in a folder called **data**. 

Change the format of sample names to eICh24_X_1.fq.gz and eICh24_X_2.fq.gz 
```bash
# Reads _1
for file in *1.fq.gz; do
  newname=$(echo "$file" | sed -E 's/^(eICh24_[0-9]+)_.+(_1\.fq\.gz)$/\1\2/')
  mv "$file" "$newname"
done

# Reads _2
for file in *2.fq.gz; do
  newname=$(echo "$file" | sed -E 's/^(eICh24_[0-9]+)_.+(_2\.fq\.gz)$/\1\2/')
  mv "$file" "$newname"
done
```

Create the folder tree and move the "data" folder.

```bash
.
├── proyect/
│    ├── README.md
│    ├── INSTALL.md
│    ├── bioinfo-tools.yml
│    ├── checksums.md5
│    ├── adapters/
│    ├── data/
│    │   ├── eICh24_1_1.fq.gz
│    │   └── eICh24_1_2.fq.gz
│    │   ├── eICh24_2_1.fq.gz
│    │   └── eICh24_2_2.fq.gz
│    ├── results/
│    │   └── fastqc/
│    │   │    ├── fastqc_one/
│    │   │    └── fastqc_final/
│    │   └── trimmomatic/
│    │   │    ├── paired/
│    │   │    └── unpaired/
│    │   └── trimgalore/
│    │   └── bowtie2/
│    │   │    ├── index/
│    │   │    ├── paired/
│    │   │    └── unpaired/
│    │   └── prinseq/
│    │   │    ├── filtered/
│    │   │    └── discarded/

```

## Step2: First quality control

Activate the environment
```bash
conda activate bioinfo-tools
```

**Replace user with the username, e.g. jdsl2009**

Run fastqc to the samples
```bash
# Move to the folder fastqc_one
cd /home/user/proyect/results/fastqc_one

fastqc /home/user/proyect/data/*.gz -o /home/user/proyect/results/fastqc/fastqc_one
```

# Step 3: Remove adapters and artifacts
Download the adapters folder and remove the adapters with the files.fa

```bash
for i in {1..48}; do
    line=$(printf "eICh24_%d" "$i")

    # Adapters pathway
    adapter_file="/home/user/proyect/adapters/${line}.fa"
    
    # Run Trimmomatic with the adapter file corresponding to each sample
    trimmomatic PE \
    /home/user/proyect/data/"$line"_1.fq.gz \
    /home/user/proyect/data/"$line"_2.fq.gz \
    /home/user/proyect/results/trimmomatic/paired/"$line"_1_paired_trim.fq.gz \
    /home/user/proyect/results/trimmomatic/unpaired/"$line"_1_unpaired_trim.fq.gz \
    /home/user/proyect/results/trimmomatic/paired/"$line"_2_paired_trim.fq.gz \
    /home/user/proyect/results/trimmomatic/unpaired/"$line"_2_unpaired_trim.fq.gz \
    ILLUMINACLIP:${adapter_file}:2:30:10:2:TRUE \
    LEADING:20 \
    SLIDINGWINDOW:4:20 \
    MINLEN:100 \
    AVGQUAL:20 \
    TRAILING:20 \
    -threads 10
    
done
```

Remove poli G artifacts
```bash
for i in {1..48}; do
    line=$(printf "eICh24_%d" "$i")

    trim_galore --cores 8 --phred33 --length 100 --stringency 3 --paired --adapter GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG --fastqc \
    -o /home/user/proyect/results/trimgalore/ \
    /home/user/proyect/results/trimmomatic/unpaired/"$line"_1_unpaired_trim.fq.gz /home/user/proyect/results/trimmomatic/unpaired/"$line"_2_unpaired_trim.fq.gz
done
```
## Step 4: Remove host and phage Phi sequence
Download the genomes and unzip them. At the time of our analyses (Jan, 2024) there is no GCF version of the Miniopterus schreibersii genome.

```bash
# Move to the folder index
cd /home/user/proyect/results/bowtie2/index/

# Download genome of *Rhinolophus ferrumequinum*
wget -P /home/user/proyect/results/bowtie2/index \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/115/265/GCF_004115265.2_mRhiFer1_v1.p/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna.gz"
gunzip /home/user/proyect/results/bowtie2/index/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna.gz
```

Repeat with:

*Myotis myotis*
```bash
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/108/235/GCF_014108235.1_mMyoMyo1.p/GCF_014108235.1_mMyoMyo1.p_genomic.fna.gz
```

*Miniopterus schreibersii*
```bash
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/146/895/GCA_964146895.1_mMinSch1.hap1.1/GCA_964146895.1_mMinSch1.hap1.1_genomic.fna.gz
```
 
*Homo sapiens*
```bash
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```
Download phage phi
```bash
if [[ ! -e "$PhiX" ]]; then
    for i in {1..10}; do
        echo "Downloading PhiX reference genome"
        esearch -db nucleotide -query "NC_001422.1" | efetch -format fasta > PhiX_NC_001422.1.fasta
        if [[ -s PhiX_NC_001422.1.fasta ]]; then
            echo "PhiX genome was downloaded"
            PhiX=$(readlink -f PhiX_NC_001422.1.fasta)
            break
        else
            echo "Try downloading again the PhiX reference genome"
        fi
    done
fi
```

Join genomes
```bash
cat \
GCA_964146895.1_mMinSch1.hap1.1_genomic.fna \
GCF_014108235.1_mMyoMyo1.p_genomic.fna \
GCF_004115265.2_mRhiFer1_v1.p_genomic.fna \
PhiX_NC_001422.1.fasta \
GCF_000001405.40_GRCh38.p14_genomic.fna \
> Mixed.fasta
```

Index the genome mix. Index construction consumes 2 to 3 times the total genome size in volatile memory.  
RAM∼30GB
```bash
bowtie2-build --large-index Mixed.fasta Mix
```

Removing unwanted sequences. This step require at least 1 TB of Storage space.
Bowtie 2 requires volatile memory proportional to the total size of the genome or concatenated genomes (∼9GB).

```bash
for i in {1..48}; do
  sample_id=$(printf "eICh24_%d" "$i")
  
  # Run Bowtie2 for the corresponding pair 
  bowtie2 \
    -x /home/user/proyect/results/bowtie2/index/Mix \
    -1 /home/user/proyect/results/trimgalore/${sample_id}_1_paired_trim_val_1.fq.gz \
    -2 /home/user/proyect/results/trimgalore/${sample_id}_2_paired_trim_val_2.fq.gz \
    -q \
    -p 10 \
    --un-conc-gz /home/user/proyect/results/bowtie2/unpaired/${sample_id}_unmapped.sam.gz \
    -S /home/user/proyect/results/bowtie2/paired/${sample_id}_mapped.sam

  echo "Processed sample: $sample_id"
done
```

Change the format of unpaired files
```bash
# Move to the folder unpaired
cd /home/user/proyect/results/bowtie2/unpaired/

# For read _1
for file in *_unmapped.1; do
  mv "$file" "${file%.1}_1.sam.gz"
done

# For read _2
for file in *_unmapped.2; do
  mv "$file" "${file%.2}_2.sam.gz"
done
```

From sam files to fastq files
```bash
for file in /home/user/proyect/results/bowtie2/unpaired/eICh24_*_unmapped_1.sam.gz; do
  sample_id=$(basename "$file" _unmapped_1.sam.gz)

  gzip -dc "/home/user/proyect/results/bowtie2/unpaired/${sample_id}_unmapped_1.sam.gz" | \
    samtools fastq - | gzip > "/home/user/proyect/results/bowtie2/unpaired/${sample_id}_unmapped_1.fastq.gz"

  gzip -dc "/home/user/proyect/results/bowtie2/unpaired/${sample_id}_unmapped_2.sam.gz" | \
    samtools fastq - | gzip > "/home/user/proyect/results/bowtie2/unpaired/${sample_id}_unmapped_2.fastq.gz"

  echo "Converted and compressed: ${sample_id}_unmapped_1 y ${sample_id}_unmapped_2"
done
```
**These files will be used in the abundance analyses**


## Step 5: Remove duplicates
These files will be used in the Humann3 analyses

```bash
for i in {1..48}; do
  sample_id=$(printf "eICh24_%d" "$i")
  
  # Temporarily unzip files
  zcat /home/user/proyect/results/bowtie2/unpaired/"$sample_id"_unmapped_1.fastq.gz > temp_1.fastq
  zcat /home/user/proyect/results/bowtie2/unpaired/"$sample_id"_unmapped_2.fastq.gz > temp_2.fastq
  
  # prinseq-lite
  prinseq-lite.pl \
    -fastq temp_1.fastq \
    -fastq2 temp_2.fastq \
    -min_len 100 \
    -min_qual_mean 20 \
    -out_format 3 \
    -derep 1 \
    -out_good /home/user/proyect/results/prinseq/filtered/"$sample_id" \
    -out_bad /home/user/proyect/results/prinseq/discarded/"$sample_id"
  
  # clean temp
  rm temp_1.fastq temp_2.fastq

  echo "Processed sample: $sample_id"
done
```

## Step6: Final quality control

Run fastqc to the samples
```bash
# Move to the folder fastqc_final
cd /home/user/proyect/results/fastqc_final

fastqc /home/user/proyect/results/prinseq/filtered/*.gz -o /home/user/proyect/results/fastqc/fastqc_final
```
