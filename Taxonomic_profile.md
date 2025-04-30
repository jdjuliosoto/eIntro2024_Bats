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
│    │   └── centrifuge/
│    │   │    ├── index_c/
│    │   │    ├── report/
│    │   │    ├── filtered_results/
│    │   │    ├── filtered_reports/
│    │   │    └── kreport/

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
    --threads 11
done
```

Execute the python script to filter taxa by score scored ≤300 and had match lengths ≤ 60 bp
```bash
# Change mode to executable
chmod +x script.py

# Execute the script
./script.py
```

Filter the report.tsv files using the taxIDs of the previously filtered files
```bash
# Change mode to executable
chmod +x script2.py

# Execute the script
./script2.py
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

**With this report you can make a visualization in R**
