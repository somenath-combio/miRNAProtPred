# miRNAProtPred

A Python package for miRNA protein prediction with SeqFinder and validator modules.

## Overview

miRNAProtPred is a bioinformatics tool designed to identify miRNA target sequences in DNA, RNA, or protein sequences. It uses the Boyer-Moore string matching algorithm combined with BLAST and ViennaRNA for comprehensive miRNA target prediction and validation.

## Features

- **SeqFinder**: Find miRNA target sequences in DNA, RNA, or protein sequences
  - Automatic sequence type detection (DNA/RNA/Protein)
  - BLAST integration for protein sequence analysis
  - Boyer-Moore pattern matching for efficient sequence searching
  - ViennaRNA integration for minimum free energy (MFE) calculation
  - Probability scoring (High/Medium/Low) based on MFE values
  - Results export to CSV format

- **Validator**: Validate miRNA-mRNA interactions (experimental module)
  - Example: hsa-miR-21-5p,hsa-miR-155-3p,hsa-miR-34a-5p
  - Automatically retrieves mature miRNA sequences and seed regions
  - Accepts DNA / RNA / Protein as mRNA inputg
  - Performs complementary seed matching
  - Performs complementary seed matching
  - Results export to CSV format
  

## Installation

### Prerequisites

- Python >= 3.7
- Required Python packages:
  - pandas >= 1.3.0
  - openpyxl >= 3.0.0
  - biopython
  - ViennaRNA
  - pyfiglet

### Install from source

```bash
git clone https://github.com/somenath-combio/mirnaprotpred.git
cd mirnaprotpred
pip install -e .
```

## Usage

### SeqFinder

Find miRNA target sequences in your input sequence:

```bash
SeqFinder <sequence>
```

**Examples:**

```bash
# DNA sequence
SeqFinder "ATGCATGCATGCATGC"

# RNA sequence
SeqFinder "AUGCAUGCAUGCAUGC"

# Protein sequence
SeqFinder "MKKLAVSLLLFLSSLA"
```

The tool will:
1. Automatically detect the sequence type
2. Search for miRNA seed sequences from the database
3. Calculate minimum free energy (MFE) using ViennaRNA
4. Assign probability scores (High: MFE ≤ -15, Medium: -15 < MFE ≤ -10, Low: MFE > -10)
5. Display results sorted by MFE
6. Optionally save results to CSV

### Validator

Validate miRNA-mRNA interactions:
```bash
validator <miRNA_IDs> <mRNA_sequence>
```
```bash
Examples
# Multiple miRNAs, DNA mRNA sequence
validator "hsa-miR-21-5p,hsa-miR-155-3p" "ATGCATGCATGC"

# RNA input
validator "hsa-miR-134-5p" "AUGCAUGCAUGC"

# Protein input
validator "hsa-miR-21-5p" "MKKLAVSLLLFLSSLA"
```

## Data Requirements

The SeqFinder module requires a data file located at `data/data.xlsx` containing:
- miRNA descriptions
- Human miRNA IDs
- Accession numbers
- Sequences
- Seed sequences (seed1, seed2, seed3 columns)

## Output

SeqFinder generates results with the following columns:
- **Description**: miRNA description
- **Human miRNA ID**: Identifier for the miRNA
- **Accession**: Accession number
- **Sequence**: Full miRNA sequence
- **Seed**: Matched seed sequence
- **Position**: Position of the match in the input sequence
- **CTS**: Complementary target site sequence
- **MFE**: Minimum free energy
- **Prob**: Probability score (High/Medium/Low)

Validator generates results with the following columns:

For each miRNA ID, the tool returns:
Complementary match site (if found)
- Bind = True/False
- **True** → Complementary seed-region match exists
- **False** → No valid complementary site



## How It Works

1. **Sequence Type Detection**: Automatically identifies whether the input is DNA, RNA, or protein
2. **Protein Processing**: For protein sequences, uses BLAST to retrieve the corresponding nucleotide sequence
3. **Pattern Matching**: Uses the Boyer-Moore algorithm to find miRNA seed sequences
4. **Energy Calculation**: Calculates duplex formation energy using ViennaRNA
5. **Scoring**: Assigns probability based on MFE thresholds

## License

MIT License

## Author

Somenath Dutta (somenath@pusan.ac.kr) Sudipta Sardar (sudipta@pusan.ac.kr)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use this tool in your research, please cite appropriately.
