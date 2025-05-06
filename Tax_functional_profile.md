# Taxonomic profile analyses

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

Download centrifuge database of Bacteria, Archaea, Viruses, Human
```bash
# Move to the folder index_c
cd /home/user/proyect/results/centrifuge/index_c

# Download database
wget https://genome-idx.s3.amazonaws.com/centrifuge/p%2Bh%2Bv.tar.gz
```

Decompress data base
```bash
tar -xvzf p+h+v.tar.gz
```

Run centrifige
```bash
for i in {1..48}; do
  line=$(printf "eICh24_%d" "$i")
  
  centrifuge -x /home/user/proyect/results/centrifuge/index_c/p+h+v \
    -1 /home/user/proyect/results/bowtie2/unpaired/"$line"_unmapped_1.fastq.gz \
    -2 /home/user/proyect/results/bowtie2/unpaired/"$line"_unmapped_2.fastq.gz \
    --report-file /home/user/proyect/results/centrifuge/report/report_"$line"_bacteria.txt \
    -S results_"$line"_bacteria.txt \
    --threads 10
done
```

Execute the python script to filter taxa by score scored ≤300 and had match lengths ≤ 60 bp
```bash
# Change mode to executable
chmod +x script.py

# Execute the script
/home/proyect/scripts/script.py
```

Filter the report.tsv files using the taxIDs of the previously filtered files
```bash
# Change mode to executable
chmod +x script2.py

# Execute the script
/home/proyect/scripts/script2.py
```

Transforms results.txt file to kreport.txt, for each sample
```bash
for i in {1..48}; do
  line=$(printf "eICh24_%d" "$i")

  centrifuge-kreport \
  -x /home/user/proyect/results/centrifuge/index_c/p+h+v /home/user/proyect/results/centrifuge/filtered_reports/filtered_report_"$line"_bacteria.tsv > \
  /home/user/proyect/results/centrifuge/kreport/kreport_filtered_"$line"_bacteria.txt
done
```

# Functional profiling


Download databases (ChocoPhlAn) and Uniref90

```bash
humann_databases --download chocophlan full /home/user/proyect/results/humann/data_bases
humann_databases --download uniref uniref90_diamond /home/user/proyect/results/humann/data_bases

# Available data bases
humann_databases --available
```

Concatenate paired-end FASTQ files from prinseq-filtered samples for HUMAnN 3.0 analysis

```bash
for i in {1..48}; do
    cat /home/user/proyect/results/prinseq/filtered/eICh24_${i}_1.fq.gz \
    /home/user/proyect/results/prinseq/filtered/eICh24_${i}_2.fq.gz > \
    /home/user/proyect/results/humann/samples_h/eICh24_${i}_combined.fq.gz
done
```

Run Humann. These analyses requiere around 4-6 hours per sample and ∼50-60GB of volatile memory (RAM).

```bash
for i in {1..48}; do
  line=$(printf "eICh24_%d" "$i")
  
  humann \
    --input /home/user/proyect/results/humann/samples_h/"$line"_combined.fq.gz \
    --output /home/user/proyect/results/humann/output \
    --threads 6
  echo "${line} done"
done
```
Normalization:
Humann 3.0 already adjusts for Reads Per Kilobase (RPK). It adjusts for gene length, but it needs to adjust for depth. Counts per Million (CPM) normalization adjusts abundance based on the total aligned reads in each sample.

```bash
for i in {1..48}; do
  humann_renorm_table -i /home/user/proyect/results/humann/output/pathabundance/eICh24_"$i"_combined_pathabundance.tsv \
                      -u cpm \
                      -o /home/user/proyect/results/humann/output/pathabundance_cpm/pathabundance_cpm_eICh24_"$i".tsv
  echo "${i} done"  
done


for i in {1..48}; do
  humann_renorm_table -i /home/user/proyect/results/humann/output/genefamilies/eICh24_"$i"_combined_genefamilies.tsv \
                      -u cpm \
                      -o /home/user/proyect/results/humann/output/genefamilies_cpm/genefamilies_cpm_eICh24_"$i".tsv
  echo "${i} done"  
done
```

Merge the output files (gene families and abundance) from the HUManN runs of all samples into three files
```bash
humann_join_tables --input /home/user/proyect/results/humann/output/pathabundance_cpm/ \
                   --output /home/user/proyect/results/humann/output/joined/pathabundance_cpm.tsv \
                   --file_name "pathabundance_cpm_eICh24_"


humann_join_tables --input /home/user/proyect/results/humann/output/genefamilies_cpm/ \
                   --output /home/user/proyect/results/humann/output/joined/genefamilie_cpm.tsv \
                   --file_name "genefamilies_cpm_eICh24_"
```
The final table can be further analyzed using R.


# Virulence profile

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
**With this report you can make a visualization in R**
