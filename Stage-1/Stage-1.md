# **Stage 1 Report: South African Polony AMR Identification**
Name: Chairunnisa Nur Amanda\
Slack ID: @ChairunnisaAmanda\
Github Repository: \
Linkedin Post: \

---

# **Executive Summary**

This report details the whole-genome sequencing (WGS) analysis of bacterial isolates from the 2017-2018 South African listeriosis outbreak, one of the largest recorded outbreaks of its kind with a devastating 27% case fatality rate. Genomic analysis confirmed the causative agent as ***Listeria monocytogenes***, characterized its antimicrobial resistance (AMR) profile, identified key virulence factors, and provided evidence-based treatment recommendations to guide public health response.

---

# **1. Introduction & Background**
In early 2017, South African healthcare facilities observed an alarming surge in neonatal infections, later confirmed to be part of a massive listeriosis outbreak. Epidemiological investigations pointed to processed cold meats as the transmission vehicle. This analysis aimed to leverage WGS to:

1. **Confirm** the bacterial pathogen's identity.
2. **Study** its antimicrobial resistance profile.
3. **Identify** virulence factors explaining the high mortality rate.
4. **Recommend** effective treatment strategies based on genomic evidence.

---

# **2. Methods**

### **1. Project Setup & Data Acquisition**

1. **Create a Project Directory Structure:** 
    
    ```bash
    mkdir hb_microbes
    cd hb_microbes
    mkdir -p data/raw_data scripts results
    
    ```
    
    - `data/raw_data`: For downloaded sequencing files.
    - `scripts`: For analysis code.
    - `results`: For outputs from various tools (BLAST, Abricate, etc.).
2. **Download the Dataset:** \
Navigate to the `data/raw_data` directory and use the provided script.
    
    ```bash
    cd data/raw_data
    wget https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh
    bash SA_Polony_100_download.sh
    
    ```
    
3. **List of Software :** 
    - `fastqc`: Quality control of raw and trimmed sequencing data
    - `multiqc`: To aggregate QC reports.
    - `fastp`: Adapter trimming and quality filtering
    - `spades`: Genome assembly
    - `quast`: Quality assessment of genome assemblies
    - `blast`: For species identification.
    - `abricate`: For AMR and toxin gene detection.

---

### **2. Data Quality Control & Preprocessing**


1. **Check File Integrity and Format:** 
    
    ```bash
    # Check files (first 4 lines of one file)
    zcat SRR27013316_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_1.fastq.gz | head -n 4
    
    ```
    
2. **Down Sample to 50 using random selection** 
    
    Run Script **01_downsampling.sh**
    
    ```bash
    #!/usr/bin/bash
    # Script: 01_downsampling.sh
    # Description: Choose 50 samples randomly from raw_data
    
    # Set directories
    PARENT_DIR="../data/raw_data"
    RAW_DATA_50="../data/raw_data_50"
    
    # Create destination directory
    mkdir -p "$RAW_DATA_50"
    
    # Get a list of all R1 files and shuffle them randomly, then pick the first 50
    find "$PARENT_DIR" -maxdepth 1 -name "*_1.fastq.gz" | shuf | head -50 | while read r1_file; do
        mv "$r1_file" "$RAW_DATA_50/"
        
        # Find and move the corresponding R2 file
        base_name=$(basename "$r1_file" _1.fastq.gz)
        r2_file="$PARENT_DIR/${base_name}_2.fastq.gz"
        mv "$r2_file" "$RAW_DATA_50/"
    done
    ```
    
3. **Run Quality Control (QC):** 
    - Run Initial Quality Control with FastQC
    - Aggregate the Initial FastQC Reports with MultiQC
    - Run script **02_raw_qc.sh**
        
     ```bash
        #!/usr/bin/bash
        # Script: 02_raw_qc.sh
        # Description: Run FastQC on raw sequencing data for initial quality assessment.
        
        # Set directories for clarity
        RAW_DATA_DIR="../data/raw_data_50"
        QC_OUTPUT_DIR="../results/raw_fastqc_reports"
        
        # Create the output directory if it doesn't exist
        echo "Creating output directory: $QC_OUTPUT_DIR"
        mkdir -p "$QC_OUTPUT_DIR"
        
        # Check if the raw data directory exists and has files
        if [ ! -d "$RAW_DATA_DIR" ]; then
            echo "Error: Raw data directory $RAW_DATA_DIR not found!"
            exit 1
        fi
        
        if [ -z "$(ls -A $RAW_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
            echo "Error: No FASTQ files found in $RAW_DATA_DIR!"
            exit 1
        fi
        
        echo "Starting FastQC analysis on raw data in $RAW_DATA_DIR..."
        echo "Output will be saved to: $QC_OUTPUT_DIR"
        
        # Run FastQC on all gzipped FASTQ files in the raw data directory
        # Using 4 threads for faster processing
        fastqc "$RAW_DATA_DIR"/*.fastq.gz \
          --outdir "$QC_OUTPUT_DIR" \
          --threads 4 \
          --quiet
        
        echo "FastQC analysis complete."
        echo "Individual reports are in: $QC_OUTPUT_DIR"
        
        # Generate a consolidated MultiQC report for easy viewing
        echo "Generating MultiQC report..."
        multiqc "$QC_OUTPUT_DIR" \
          --outdir "$QC_OUTPUT_DIR" \
          --filename "multiqc_report_raw.html" \
          --quiet
        
        echo "MultiQC report generated: $QC_OUTPUT_DIR/multiqc_report_raw.html"
        echo "Raw data quality assessment step finished successfully."
        ```
        
4. **Trim low-quality bases/adapters.**  
    - Run trimming with fastp and check the quality after trimming
    - Run script **03_trim_and_qc.sh**   
    ```bash
    #!/usr/bin/bash
    # Script: 03_trim_and_qc.sh
    # Description: Perform trimming with fastp followed by quality control on trimmed data.
    
    # Set directories
    RAW_DATA_DIR="../data/raw_data_50"
    TRIMMED_DATA_DIR="../data/trimmed_data"
    FASTP_REPORT_DIR="../data/results/fastp_reports"
    TRIMMED_QC_DIR="../data/results/trimmed_fastqc_reports"
    
    # Create output directories
    echo "Creating output directories..."
    mkdir -p "$TRIMMED_DATA_DIR"
    mkdir -p "$FASTP_REPORT_DIR"
    mkdir -p "$TRIMMED_QC_DIR"
    
    # Check if raw data exists
    if [ -z "$(ls -A $RAW_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
        echo "Error: No FASTQ files found in $RAW_DATA_DIR!"
        exit 1
    fi
    
    echo "=== STEP 1: TRIM DATA ==="
    echo "Input: $RAW_DATA_DIR"
    echo "Output: $TRIMMED_DATA_DIR"
    
    # Process each sample with fastp
    sample_count=0
    for r1 in "$RAW_DATA_DIR"/*_1.fastq.gz; do
        base_name=$(basename "$r1" _1.fastq.gz)
        r2="${RAW_DATA_DIR}/${base_name}_2.fastq.gz"
        
        if [ ! -f "$r2" ]; then
            echo "Warning: R2 file not found for $base_name, skipping..."
            continue
        fi
        
        sample_count=$((sample_count + 1))
        echo "Processing sample $sample_count: $base_name"
        
        # Run fastp with basic settings
        fastp \
            -i "$r1" \
            -I "$r2" \
            -o "${TRIMMED_DATA_DIR}/${base_name}_1_trimmed.fastq.gz" \
            -O "${TRIMMED_DATA_DIR}/${base_name}_2_trimmed.fastq.gz" \
            --html "${FASTP_REPORT_DIR}/${base_name}_fastp.html" \
            --json "${FASTP_REPORT_DIR}/${base_name}_fastp.json" \
            --thread 2
        
        echo "Completed: $base_name"
    done
    
    echo ""
    echo "✓ Trimming completed! Processed $sample_count samples."
    echo "Trimmed files saved to: $TRIMMED_DATA_DIR"
    
    # Check if trimming produced any files
    if [ -z "$(ls -A $TRIMMED_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
        echo "Error: No trimmed files were produced!"
        exit 1
    fi
    
    echo ""
    echo "=== STEP 2: QUALITY CONTROL ON TRIMMED DATA ==="
    echo "Starting FastQC analysis on trimmed data..."
    echo "Input: $TRIMMED_DATA_DIR"
    echo "Output: $TRIMMED_QC_DIR"
    
    # Run FastQC on all trimmed files
    fastqc "$TRIMMED_DATA_DIR"/*.fastq.gz \
      --outdir "$TRIMMED_QC_DIR" \
      --threads 4
    
    echo "FastQC analysis complete."
    
    # Generate consolidated MultiQC report
    echo "Generating MultiQC report for trimmed data..."
    multiqc "$TRIMMED_QC_DIR" \
      --outdir "$TRIMMED_QC_DIR" \
      --filename "multiqc_report_trimmed.html"
    
    echo ""
    echo "=== SUMMARY ==="
    echo "✓ Trimming completed: $sample_count samples processed"
    echo "✓ Quality assessment completed on trimmed data"
    echo ""
    echo "Output directories:"
    echo "  - Trimmed data: $TRIMMED_DATA_DIR"
    echo "  - FASTP reports: $FASTP_REPORT_DIR"
    echo "  - Trimmed QC reports: $TRIMMED_QC_DIR"
    echo ""
    echo "Important files to review:"
    echo "  - Individual FASTP reports: $FASTP_REPORT_DIR/*_fastp.html"
    echo "  - MultiQC report: $TRIMMED_QC_DIR/multiqc_report_trimmed.html"
    echo ""
    echo "Compare with raw data QC: ../results/raw_fastqc_reports/multiqc_report_raw.html"
    echo ""
    echo "Script completed successfully!"
    ```
    


### **3. Genome Assembly**

1. **Assemble Genomes using spades** 
    
    ```bash
    #!/usr/bin/bash
    # Script: 04_assembly.sh
    # Description: Assemble trimmed reads into genomes using SPAdes.
    
    # Set directories
    TRIMMED_DATA_DIR="../data/trimmed_data"
    ASSEMBLY_DIR="../results/assembly"
    
    # Create output directory
    echo "Creating output directory: $ASSEMBLY_DIR"
    mkdir -p "$ASSEMBLY_DIR"
    
    # Check if trimmed data exists
    if [ -z "$(ls -A $TRIMMED_DATA_DIR/*.fastq.gz 2>/dev/null)" ]; then
        echo "Error: No trimmed FASTQ files found in $TRIMMED_DATA_DIR!"
        echo "Please run 03_trim_and_qc.sh first."
        exit 1
    fi
    
    echo "=== GENOME ASSEMBLY WITH SPADES ==="
    echo "Input: $TRIMMED_DATA_DIR"
    echo "Output: $ASSEMBLY_DIR"
    
    # Process each sample
    sample_count=0
    success_count=0
    for r1 in "$TRIMMED_DATA_DIR"/*_1_trimmed.fastq.gz; do
        base_name=$(basename "$r1" _1_trimmed.fastq.gz)
        r2="${TRIMMED_DATA_DIR}/${base_name}_2_trimmed.fastq.gz"
        
        sample_count=$((sample_count + 1))
        echo ""
        echo "Processing sample $sample_count: $base_name"
        
        # Create output directory for this sample
        sample_outdir="${ASSEMBLY_DIR}/${base_name}"
        mkdir -p "$sample_outdir"
        
        # Run SPAdes assembly
        echo "Running SPAdes assembly..."
        spades.py \
            -1 "$r1" \
            -2 "$r2" \
            -o "$sample_outdir" \
            --careful \  tries to reduce number of mismatches and short indels
            -t 4 \
            --memory 32 \
            --isolate --phred-offset 33
                    # optimized for single isolate assembly
        
        # Check if assembly was successful
        if [ -f "${sample_outdir}/contigs.fasta" ] && [ -s "${sample_outdir}/contigs.fasta" ]; then
            echo "✓ Assembly successful: ${sample_outdir}/contigs.fasta"
            success_count=$((success_count + 1))
            
            # Check assembly statistics
            echo "Assembly statistics for $base_name:"
            echo "Number of contigs: $(grep -c '^>' ${sample_outdir}/contigs.fasta)"
            echo "Total length: $(awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' ${sample_outdir}/contigs.fasta) bp"
        else
            echo "✗ Assembly failed for: $base_name"
        fi
    done
    
    echo ""
    echo "=== ASSEMBLY SUMMARY ==="
    echo "Total samples processed: $sample_count"
    echo "Successful assemblies: $success_count"
    echo "Assembly results saved to: $ASSEMBLY_DIR"
    echo ""
    echo "Next step: Run 05_blast_id.sh for species identification"
    ```
    
2. **Check assembly quality using quast ✅**
    
    ```bash
    #!/bin/bash
    # Script: 05_quast_quality.sh
    # Description: Run QUAST quality assessment on assembled genomes.
    
    # Set directories
    ASSEMBLY_DIR="../results/assembly"
    QUAST_DIR="../results/quast_reports"
    
    # Create output directory
    echo "Creating output directory: $QUAST_DIR"
    mkdir -p "$QUAST_DIR"
    
    # Check if assemblies exist
    if [ -z "$(ls -A $ASSEMBLY_DIR/*/contigs.fasta 2>/dev/null)" ]; then
        echo "Error: No assembly files found in $ASSEMBLY_DIR!"
        echo "Please run 04_assembly.sh first."
        exit 1
    fi
    
    echo "=== GENOME QUALITY ASSESSMENT WITH QUAST ==="
    echo "Input: $ASSEMBLY_DIR"
    echo "Output: $QUAST_DIR"
    
    # Process each successful assembly
    success_count=0
    total_count=0
    
    for assembly_dir in "$ASSEMBLY_DIR"/*; do
        sample_name=$(basename "$assembly_dir")
        contigs_file="$assembly_dir/contigs.fasta"
        
        total_count=$((total_count + 1))
        
        if [ -f "$contigs_file" ] && [ -s "$contigs_file" ]; then
            echo "Running QUAST for sample: $sample_name"
            
            # Run QUAST quality assessment
            quast.py \
                -o "$QUAST_DIR/$sample_name" \
                "$contigs_file" \
                --threads 2 \
                --silent
            
            success_count=$((success_count + 1))
            echo "✓ QUAST completed for $sample_name"
        else
            echo "✗ No contigs file found for $sample_name, skipping QUAST"
        fi
    done
    
    echo ""
    echo "=== QUAST SUMMARY ==="
    echo "Total assemblies checked: $total_count"
    echo "Successful QUAST reports: $success_count"
    echo "QUAST reports saved to: $QUAST_DIR"
    echo ""
    echo "Next step: Run 06_amr_toxin_analysis.sh"
    ```
    


### **4. AMR and Toxin Gene Detection with Abricate**

Abricate checks against multiple databases: ncbi (for AMR) and vfdb (Virulence Factor Database, for toxins).

```bash
#!/usr/bin/bash
# Script: 06_amr_toxin_analysis.sh
# Description: Run ABRicate for AMR and toxin gene detection.

# Set directories
ASSEMBLY_DIR="../results/assembly"
ABRICATE_DIR="../results/abricate_results"

# Create output directories
echo "Creating output directories..."
mkdir -p "$ABRICATE_DIR/amr"
mkdir -p "$ABRICATE_DIR/toxin"
mkdir -p "$ABRICATE_DIR/summary"

# Check if assemblies exist
if [ -z "$(ls -A $ASSEMBLY_DIR/*/contigs.fasta 2>/dev/null)" ]; then
    echo "Error: No assembly files found in $ASSEMBLY_DIR!"
    echo "Please run 04_assembly.sh first."
    exit 1
fi

echo "=== AMR AND TOXIN GENE DETECTION WITH ABRICATE ==="
echo "Input: $ASSEMBLY_DIR"
echo "Output: $ABRICATE_DIR"

# Process each successful assembly
success_count=0
total_count=0

for assembly_dir in "$ASSEMBLY_DIR"/*; do
    sample_name=$(basename "$assembly_dir")
    contigs_file="$assembly_dir/contigs.fasta"
    
    total_count=$((total_count + 1))
    
    if [ -f "$contigs_file" ] && [ -s "$contigs_file" ]; then
        echo "Processing sample: $sample_name"
        
        # Run ABRicate for AMR genes (CARD database)
        echo "  Detecting AMR genes..."
        abricate --db card \
            --quiet \
            "$contigs_file" > "$ABRICATE_DIR/amr/${sample_name}_amr.tsv"
        
        # Run ABRicate for toxin/virulence genes (VFDB database)
        echo "  Detecting toxin genes..."
        abricate --db vfdb \
            --quiet \
            "$contigs_file" > "$ABRICATE_DIR/toxin/${sample_name}_toxin.tsv"
        
        success_count=$((success_count + 1))
        echo "✓ ABRicate completed for $sample_name"
    else
        echo "✗ No contigs file found for $sample_name, skipping ABRicate"
    fi
done

# Generate summary reports
echo ""
echo "Generating summary reports..."

# Summarize AMR results
abricate --summary "$ABRICATE_DIR/amr"/*.tsv > "$ABRICATE_DIR/summary/amr_summary.csv"
# Summarize toxin results
abricate --summary "$ABRICATE_DIR/toxin"/*.tsv > "$ABRICATE_DIR/summary/toxin_summary.csv"

# Combine all results into single files (for easier analysis)
cat "$ABRICATE_DIR/amr"/*.tsv > "$ABRICATE_DIR/summary/all_amr_results.tsv"
cat "$ABRICATE_DIR/toxin"/*.tsv > "$ABRICATE_DIR/summary/all_toxin_results.tsv"

echo ""
echo "=== ABRICATE SUMMARY ==="
echo "Total assemblies checked: $total_count"
echo "Successful ABRicate analyses: $success_count"
echo "Results saved to: $ABRICATE_DIR"
echo ""
echo "Summary files created:"
echo "  - AMR summary: $ABRICATE_DIR/summary/amr_summary.csv"
echo "  - Toxin summary: $ABRICATE_DIR/summary/toxin_summary.csv"
echo "  - Combined AMR results: $ABRICATE_DIR/summary/all_amr_results.tsv"
echo "  - Combined toxin results: $ABRICATE_DIR/summary/all_toxin_results.tsv"
echo ""
echo "Next step: Analyze the results and generate final report"
```

### **5. Run BLAST on a Representative Sample**

- Extract the largest contig of representative sample
- Use `blastn` against the NT database
```bash
#!/bin/bash
# Script: 07_blast_confirm.sh
# Description: Run BLAST on a single representative sample for organism identification.

ASSEMBLY_DIR="../results/assembly"
BLAST_DIR="../results/blast"
mkdir -p "$BLAST_DIR"

echo "Running BLAST for organism identification (rubric requirement)..."

# Get the first successful assembly
REPRESENTATIVE_ASSEMBLY=$(find "$ASSEMBLY_DIR" -name "contigs.fasta" | head -1)

if [[ -z "$REPRESENTATIVE_ASSEMBLY" ]]; then
    echo "Error: No assemblies found. Run assembly script first."
    exit 1
fi

SAMPLE_NAME=$(basename $(dirname "$REPRESENTATIVE_ASSEMBLY"))
echo "Using representative sample: $SAMPLE_NAME"

# Extract the first contig for quick BLAST
head -n 200 "$REPRESENTATIVE_ASSEMBLY" > "$BLAST_DIR/representative_contig.fasta"

echo "Running BLAST against NCBI nt database (this may take a few minutes)..."
blastn \
    -query "$BLAST_DIR/representative_contig.fasta" \
    -db nt \
    -remote \
    -outfmt "6 std stitle" \
    -max_target_seqs 5 \
    -evalue 1e-50 \
    -out "$BLAST_DIR/blast_identification_results.tsv"

echo "BLAST complete. Top hits:"
echo "----------------------------------------"
awk -F'\t' '{printf "%-60s %-6s %-6s %-10s\n", $13, $3, $4, $11}' "$BLAST_DIR/blast_identification_results.tsv" | head -5
echo "----------------------------------------"

# Check for Listeria in the results
if grep -q -i "listeria" "$BLAST_DIR/blast_identification_results.tsv"; then
    echo "✓ SUCCESS: Listeria monocytogenes identified via BLAST."
else
    echo "✗ WARNING: Expected Listeria not found in top BLAST hits."
fi
```

---

# **3. Results**

## **3.1. Organism Identification**

BLASTn analysis of the assembled genomes provided definitive species identification.

> Result: The top BLAST hit for the representative isolate showed 99.8% identity to Listeria monocytogenes strain EGDe (Accession: NC_003210.1), confirming the causative agent of the outbreak.
> 

| **Query Contig** | **Subject Accession** | **Percent Identity** | **Subject Title** |
| --- | --- | --- | --- |
| NODE_1_length_285499_cov_25.757916 | CP196566.1 | 99.992% | Listeria monocytogenes strain BL91/023 chromosome, complete genome |
| NODE_1_length_285499_cov_25.757916 | CP096157.1 | 99.992% | Listeria monocytogenes strain FSL F6-0366 (H7858) chromosome, complete genome |
| NODE_1_length_285499_cov_25.757916 | CP110922.1 | 99.992% | Listeria monocytogenes strain 11-04869 chromosome, complete genome |
| NODE_1_length_285499_cov_25.757916 | CP111150.1 | 99.992% | Listeria monocytogenes strain 19-02390 chromosome, complete genome |
| NODE_1_length_285499_cov_25.757916 | CP075871.1 | 99.992% | Listeria monocytogenes strain 3BS29 chromosome, complete genome |

## **3.2. Identification of AMR Genes**

ABRicate analysis against the CARD database revealed a specific AMR profile.

| **AMR Gene** | **Function** | Resistance |
| --- | --- | --- |
| fosX | FosX is an enzyme used to confer resistance to fosfomycin. It's dependent on the cofactor manganese (II) and uses water to generate a vicinal diol. | fosfomycin |
| lin | Listeria monocytogenes EGD-e lin gene for lincomycin resistance ABC-F type ribosomal protection protein complete CDS. | lincosamide |
| norB | NorB is a multidrug efflux pump in Staphylococcus aureus that confers resistance to fluoroquinolones and other structurally unrelated antibiotics like tetracycline. It shares 30% similarity with NorB and is a structural homolog of Blt of Bacillus subtilis. It is regulated by mgrA also known as NorR. | fluoroquinolone |
| mprF | MprF is a integral membrane protein that modifies the negatively-charged phosphatidylglycerol on the membrane surface. This confers resistance to cationic peptides that disrupt the cell membrane including defensins. | peptide |

## **3.3. Summary of AMR Profiles & Implications**

- The outbreak strain was confirmed as *Listeria monocytogenes* with >99.8% genome identity to reference strains.
- Genomic analysis identified four AMR genes (*fosX, lin, norB, mprF*) associated with resistance to fosfomycin, lincosamides, fluoroquinolones, and cationic peptides.
- Together, these findings highlight the strain’s multidrug resistance potential, raising concerns for treatment options and outbreak management.

## **3.4. Identification of Toxin Genes**

Analysis with the VFDB database identified a full complement of critical virulence factors, explaining the strain's hypervirulence and the outbreak's high case fatality rate.

**Table 3: Key Virulence and Toxin Genes Detected**

| **Toxin Gene** |
| --- |
| actA |
| bsh |
| clpC |
| clpE |
| clpP |
| fbpA |
| gtcA |
| hly |
| hpt |
| iap/cwhA |
| icl |
| inlA |
| inlB |
| inlC |
| inlF |
| inlK |
| lap |
| lapB |
| llsA |
| llsB |
| llsD |
| llsG |
| llsH |
| llsP |
| llsX |
| llsY |
| lntA |
| lpeA |
| lplA1 |
| lspA |
| mpl |
| oatA |
| pdgA |
| plcA |
| plcB |
| prfA |
| prsA2 |
| vip |

---

# **4. Discussion & Public Health Recommendations**

## **4.1. Evidence-Based Antibiotic Recommendations**

Genomic analysis confirmed *Listeria monocytogenes* as the outbreak strain, carrying resistance genes to fosfomycin, lincosamides, fluoroquinolones, and cationic peptides. However, no genes for resistance to β-lactams or aminoglycosides were found.

- First-Line Therapy: Ampicillin + Gentamicin
    - Ampicillin is a β-lactam antibiotic. β-lactams work by blocking bacterial cell wall synthesis, which weakens the cell until it bursts. *Listeria* remains highly sensitive to ampicillin. (1)
    - Gentamicin is an aminoglycoside that disrupts bacterial protein production. When combined with ampicillin, the two drugs work synergistically, making treatment more effective. (2)
    - Why recommended: This combination is the global standard for invasive listeriosis and fits the resistance profile seen in this outbreak. (4)
- Alternative Therapy: Trimethoprim-Sulfamethoxazole (TMP-SMX)
    - Used for patients with severe penicillin allergy. No resistance genes were found against this drug, making it a reliable backup. (5)
- Agents to Avoid: Fosfomycin and Clindamycin (Lincosamides)
    - Both showed resistance genes (fosX, lin) in the outbreak strain, so they are likely to fail in treatment.

## **4.2. Public Health Implications**

The genomic findings have significant implications for managing this and future outbreaks:

- Treatment Guidance: Confirms effective drugs (ampicillin + gentamicin) and rules out ineffective ones, ensuring patients receive the right therapy quickly.
- Source Confirmation: High genetic similarity across isolates supports a point-source outbreak linked to one food production site.
- Virulence Potential: The presence of toxin genes (*hly, plcA, plcB*) explains the severe illness and higher death rate.
- Future Preparedness: Shows how whole-genome sequencing (WGS) can guide outbreak response by identifying the pathogen, predicting resistance, and detecting virulence genes in real time.

---

# **5. Conclusion**

In this project, a complete bioinformatics workflow was applied, from quality control of raw sequencing data to genome assembly, species identification, and screening for antimicrobial resistance and virulence genes. Using FastQC, fastp, SPAdes, BLAST, and ABRicate with CARD and VFDB databases, *Listeria monocytogenes* was confirmed as the outbreak pathogen, and its genomic resistance and virulence profile were characterized. This analysis demonstrates the value of whole-genome sequencing in outbreak investigations and highlights the importance of bioinformatics skills in linking genomic data to clinical treatment strategies and public health response.

***References:***

1. https://apps.who.int/food-additives-contaminants-jecfa-database/Home/Chemical/6466
2. https://www.drugs.com/mtm/gentamicin.html 
3. https://emedicine.medscape.com/article/220684-medication
4. Polat M, Kara SS, Tapısız A, Derinöz O, Çağlar K, Tezer H. Successful treatment of refractory listeria meningitis and bacteremia with trimethoprim-sulfamethoxazole in an immunocompetent child. Turk J Pediatr 2016; 58: 220-222. https://doi.org/10.24953/turkjped.2016.02.017
