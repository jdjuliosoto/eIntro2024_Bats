# Analysis of metabolic pathways

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

```

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
