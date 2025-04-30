# Virulence analyses 

Activate the environment

```bash
conda activate bioinfo-tools
```

Continue creating the folder tree.

```bash
.
├── proyect/
│    ├── README.md
│    ├── INSTALL.md
│    ├── bioinfo-tools.yml
│    ├── checksums.md5
│    ├── adapters/
│    ├── scripts/
│    │   ├── script.py
│    │   └── script2.py
│    ├── data/
│    │   ├── eICh24_1_1.fq.gz
│    │   ├── eICh24_1_2.fq.gz
│    │   ├── eICh24_2_1.fq.gz
│    │   └── eICh24_2_2.fq.gz
│    ├── results/
│    │   └── fastqc/
│    │   │    ├── fastqc_one/
│    │   │    └── fastqc_final/
│    │   └── trimmomatic/
│    │   │    ├── paired/
│    │   │    └── unpaired/
│    │   ├── trimgalore/
│    │   ├── bowtie2/
│    │   │    ├── index/
│    │   │    ├── paired/
│    │   │    └── unpaired/
│    │   ├── prinseq/
│    │   │    ├── filtered/
│    │   │    └── discarded/
│    │   ├── centrifuge/
│    │   │    ├── index_c/
│    │   │    ├── report/
│    │   │    ├── filtered_results/
│    │   │    ├── filtered_reports/
│    │   │    └── kreport/
│    │   ├── humann/
│    │   │    ├── data_bases/
│    │   │    ├── samples_h/
│    │   │    ├── output/
│    │   │    │    ├── pathabundance/
│    │   │    │    ├── pathabundance_cpm/
│    │   │    │    ├── genefamilies/
│    │   │    │    ├── genefamilies_cpm/
│    │   │    │    └── joined/
│    │   ├── virulence/
│    │   │    ├── index_v/
│    │   │    ├── samples_v/
│    │   │    ├── results_v/
│    │   │    └── fasta_v/

```

Move to the folder indes_v and download the data base of virulence factors. This step requiere at least 1 TB of Storage space.

```bash
cd /home/user/proyect/results/virulence/index_v

# Download DB
wget https://ftp.ncbi.nlm.nih.gov/blast/db/VFDB.tar.gz

# Decompress DB
tar -xvzf VFDB.tar.gz
```

Index in BLAST format
```bash
makeblastdb -in VFDB.fasta -dbtype nucl -out VFDB
```
Sample preparation. Decompress both .gz files in memory, join them together, and recompress

```bash
cd /home/user/proyect/results/virulence/samples_v

# sample preparation
for i in {1..48}; do
  input1="/home/user/proyect/results/bowtie2/unpaired/eICh24_${i}_unmapped_1.fastq.gz"
  input2="/home/user/proyect/results/bowtie2/unpaired/eICh24_${i}_unmapped_2.fastq.gz"
  output="/home/user/proyect/results/virulence/samples_v/eICh24_${i}.fq.gz"
  zcat "$input1" "$input2" | gzip -c > "$output"
done
```
Decompress, and convert from .fq.gz to .fasta
```bash
cd /home/user/proyect/results/virulence/fasta_v

# Data conversion
for i in {1..48}; do
  line=$(printf "eICh24_%d" "$i")
  zcat "/home/user/proyect/results/virulence/samples_v/${line}.fq.gz" | sed -n '1~4s/^@/>/p;2~4p' > "/home/user/proyect/results/virulence/fasta_v/${line}.fasta"
done

```

Search for virulence factors
```bash
cd /home/user/proyect/results/virulence/results_v

# Search
for i in {1..48}; do
  line=$(printf "eICh24_%d" "$i")
  
  blastn \
    -num_threads 10 \
    -perc_identity 95 \
    -evalue 10 \
    -outfmt '6 qseqid sseqid stitle pident evalue staxids' \
    -db /home/user/proyect/results/virulence/index_v/VFDB \
    -query /home/user/proyect/results/virulence/fasta_v/"$line".fasta \
    -out "/home/user/proyect/results/virulence/results_v/$line"_VFDB.csv \
    -max_target_seqs 5
done
```
**The final table can be further analyzed using R.**
