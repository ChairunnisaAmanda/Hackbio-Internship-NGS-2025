# Stage 0: Hackbio Internship NGS
# Project 1: BASh Basic

# 1. Print your name
echo "Chairunnisa Amanda"

# 2. Create a folder titled your name
mkdir amanda

# 3. Create another new directory titled biocomputing and change to that directory with one line of command
mkdir biocomputing && cd biocomputing

# 4. Download these 3 files:
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# 5. Move the .fna file to the folder titled your name
mv wildtype.fna ../amanda

# 6. Delete the duplicate gbk file
rm *.gbk.1

# 7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
grep -i "tatatata" ../amanda/wildtype.fna && echo "Mutant" || echo "Wildtype"

# 8. If mutant, print all matching lines into a new file
grep -i "tatatata" ../amanda/wildtype.fna > mutant_match.txt

# 9. Count number of lines (excluding header) in the .gbk file
awk '/^ORIGIN/,0{if (!/^ORIGIN/) print}' wildtype.gbk | wc -l

# 10. Print the sequence length of the .gbk file. (Use the LOCUS tag in the first line)
awk '/^LOCUS/ {print $3}' wildtype.gbk

# 11. Print the source organism of the .gbk file. (Use the SOURCE tag in the first line)
awk '/^SOURCE/ { $1=""; print $0 }' wildtype.gbk

# 12. List all the gene names of the .gbk file. Hint {grep '/gene='}
grep '/gene=' wildtype.gbk | sed 's/.*\/gene="\([^"]*\)".*/\1/'

# 13. Clear your terminal space and print all commands used today
clear
history

# 14. List the files in the two folders
echo "amanda:"; ls amanda
echo "biocomputing: "; ls biocomputing
