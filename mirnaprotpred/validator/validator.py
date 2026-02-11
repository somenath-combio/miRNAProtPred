import sys
from importlib.resources import files
import pandas as pd
import re
from Bio import Blast, Entrez, SeqIO
from Bio.Blast import NCBIXML

Entrez.email = "sudipta@pusan.ac.kr"
def get_sequence_type(seq: str) -> str:
    if not seq:
        raise ValueError("Empty sequence provided")
    seq = seq.upper().strip()
    dna_chars = set('ATGCN')
    rna_chars = set('AUGCN')
    protein_chars = set('ACDEFGHIKLMNPQRSTVWY*')

    seq_chars = set(seq)
    if seq_chars.issubset(dna_chars):
        if 'T' in seq_chars:
            return 'DNA'
        elif 'U' in seq_chars:
            return 'RNA'
        else:
            return 'DNA'

    elif seq_chars.issubset(rna_chars):
        return 'RNA'

    elif seq_chars.issubset(protein_chars):
        return 'Protein'

    else:
        invalid_chars = seq_chars - (dna_chars | rna_chars | protein_chars)
        raise ValueError(f"Invalid sequence. Contains unexpected characters: {invalid_chars}")
def blast(sequence):
    print('Starting Blast')
    response = Blast.qblast(program='tblastn',
                            database='core_nt',
                            sequence=sequence,
                            format_type='XML')
    records = NCBIXML.parse(response)
    if records is None:
        print("No BLAST records found.")
        return None, None, None, None
    # if there is record, process the first one and top alignment
    record = next(records)
    if not record.alignments:
        print("No BLAST alignments found for this sequence.")
        return None, None, None, None
    top_alignment = record.alignments[0]
    ncbi_gi = top_alignment.accession
    ncbi_id = top_alignment.hit_id.split('|')[3]
    hsp = top_alignment.hsps[0]
    start_pos = hsp.sbjct_start
    end_pos = hsp.sbjct_end
    if int(start_pos) > int(end_pos):
        start_pos, end_pos = end_pos, start_pos
    return ncbi_gi, ncbi_id, start_pos, end_pos

def retrieve_seq(ncbi_id, start_pos, end_pos):
    print(f'Retrieving Sequence with NCBI ID: {ncbi_id}')
    handle = Entrez.efetch(db='nucleotide', id=ncbi_id, rettype='fasta', retmode='text')
    record = SeqIO.read(handle, 'fasta')
    seq_of_interest = record.seq[int(start_pos)-1:int(end_pos)-1]
    return seq_of_interest
def validator(mirna_ids: str, mrna: str) -> pd.DataFrame:
    data_file = files('mirnaprotpred.validator').joinpath('data').joinpath('data.xlsx')
    df = pd.read_excel(data_file, engine="openpyxl")
    mirna_id_list = mirna_ids.split(',')
    
    input_seq_type = get_sequence_type(mrna)
    if input_seq_type == "DNA":
        dna_seq = mrna.lower().strip()
    elif input_seq_type == "RNA":
        dna_seq = mrna.lower().replace('u', 't').strip()
    elif input_seq_type == "Protein":
        ncbi_gi, ncbi_id, start_pos, end_pos = blast(mrna)
        if ncbi_gi is None:
            print("Error: Could not find matching DNA sequence for the provided protein sequence.")
            print("Please try with a DNA or RNA sequence instead.")
            sys.exit(1)
        print(f'NCBI GI: {ncbi_gi}')
        dna_seq = str(retrieve_seq(ncbi_id, start_pos, end_pos)).lower().strip()
    else:
        raise ValueError("Invalid mRNA sequence type.")

    matching_rows = []
    for mirna_id in mirna_id_list:
        mirna_id = mirna_id + " "
        if mirna_id in df['Human miRNA ID'].values:
            # Get rows where 'Human miRNA ID' matches mirna_id
            rows = df[df['Human miRNA ID'] == mirna_id]
            matching_rows.append(rows)

    if matching_rows:
        result_df = pd.concat(matching_rows)
        
        for index, row in result_df.iterrows():
            seeds = row[-3:].values
            result_df.at[index, 'interaction'] = any(re.search(seed, dna_seq.lower()) for seed in seeds if pd.notna(seed))
        
        print(result_df[['Human miRNA ID', 'interaction']])

        wanna_save = input("Do you want to save this result? (y/n): ").strip().lower()
        if wanna_save == 'y':
            result_df.to_csv("validator_output.csv", index=False)
            print("Result saved to 'validator_output.csv'.")
        elif wanna_save == 'n':
            print("You chose not to save the result.")
        else:
            print("Invalid input. Result not saved.")

    else:
        print("No matching miRNA IDs found.")

def cli():
    """CLI entry point for validator"""
    if len(sys.argv) != 3:
        print("Usage: validator <miRNA_sequence> <mRNA_sequence>")
        sys.exit(1)
    validator(sys.argv[1], sys.argv[2])
