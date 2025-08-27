# Stage 0: Hackbio Internship NGS
# Project 2: Installing Bioinformatics Software on the Terminal

# 1. Activate your base conda environment
conda activate

# 2. Create a conda environment named funtools
conda create --name funtools

# 3. Activate the funtools environment
conda activate funtools

# 4. Install Figlet using conda
conda install tsnyder::figlet

# 5. Run figlet <your name> 
figlet amanda

# 6. Install bwa through the bioconda channel
conda install bioconda::bwa

# 7. Install blast through the bioconda channel
conda install bioconda::blast   

# 8. Install samtools through the bioconda channel
conda install bioconda::samtools

# 9. Install bedtools through the bioconda channel
conda install bioconda::bedtools

# 10. Install spades.py through the bioconda channel
conda install bioconda::spades

# 11. Install bcftools through the bioconda channel
conda install bioconda::bcftools

# 12. Install fastp through the bioconda channel
conda install bioconda::fastp

# 13. Install multiqc through the bioconda channel
conda install bioconda::multiqc

