# Installation Instructions

## **Hardware and Software Setup**
All analyses were performed using a workstation equipped with:

* Processor : Intel(R) Core(TM) i9-10900F CPU @ 2.80GHz (10 cores, 20 threads) 
* Main Drive (OS + apps):  
    - Type: NVMe SSD  
    - Capacity: 238.5 GB  
    - Model: Western Digital PC SN730 NVMe WDC 256GB
* Primary Data Storage:  
    - Type: HDD  
    - Capacity: 1.8 TB  
    - Model: Seagate ST2000DM008-2FR1  
* RAM: 64 GB DDR4 (4 x 16 GB Kingston KHX2666C16D4/16GX DIMMs)
    - Speed: 2667 MT/s (configured at 2400 MT/s)
    - Type: Synchronous DDR4, Volatile memory

The operating system used was Ubuntu 24.04.2 LTS , the latest Long-Term Support (LTS) release at the time of analysis.
For most bioinformatic tools, we used Anaconda3 v2024-10-1 as the main environment manager.

## Step 1: Prepare the setup

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
