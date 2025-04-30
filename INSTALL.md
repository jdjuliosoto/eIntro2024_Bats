# Installation Instructions

## Step 1: Prepare the setup
All the analysis were done in the latest LTS version of Ubuntu "Ubuntu 24.04.2 LTS". For the majority of bioinformatic tools we used Anaconda3 v2024-10-1.

Download and install Anaconda:
```bash
# Install Anaconda
wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
bash Anaconda3-2024.10-1-Linux-x86_64.sh
source ~/.bashrc
```
When prompted, select "Yes" to add Anaconda to the PATH automatically.


## Step 2: Install and activate bioinformatic tools
Once Conda is installed, create the environment with the bioinformatic tools, Conda chanels priority, and dependencies for the raw reads analyses. 
Use the bioinfo-tools.yml file.
```bash
conda env create -f bioinfo-tools.yml

# Check the created environment
conda env list

# List the libraries in the
conda list environment
```
This will install only the packages you mentioned.


Activate the environment
```bash
conda activate bioinfo-tools
```

## Step 3: Request data
The raw data can be requested via e-mail to the Infectious and Tropical Diseases Research Group (e-INTRO), Biomedical Research Institute of Salamanca-Research Centre for Tropical Diseases (IBSAL-CIETUS), Faculty of Pharmacy, University of Salamanca, 37008 Salamanca, Spain.

Because sequencing documents are relatively large, pre-hash comparison is recommended to ensure that downloaded files are not corrupted.

Compare hash against checksums.md5
```bash
md5sum -c checksums.md5
```
Expected result: 
```
./J1/J1_MKDL240004402-1A_22MGC7LT3_L4_1.fq.gz: OK
```
If any file is corrupt or incomplete: 
```
./J1/J1_MKDL240004402-1A_22MGC7LT3_L4_1.fq.gz: FAILED
```
