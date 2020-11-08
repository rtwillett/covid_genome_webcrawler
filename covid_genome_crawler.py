import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import requests

from Bio import Entrez
from Bio import SeqIO
Entrez.email='ryan.willett@gmail.com'


def extract_cd_features(cd):
    loc = cd['GBFeature_location']
    gene_name = [i['GBQualifier_value'] for i in cd['GBFeature_quals'] if i['GBQualifier_name'] == 'gene'][0]
    protein_seq = [i['GBQualifier_value'] for i in cd['GBFeature_quals'] if i['GBQualifier_name'] == 'translation'][0]
    return pd.DataFrame([(loc, gene_name, protein_seq)], columns = ['location', 'protein_name', 'protein_seq'])

def extract_cds(rec):
    id_ = rec['GBSeq_locus']
    genbank_url = f'https://www.ncbi.nlm.nih.gov/nuccore/{id_}?report=genbank'
    try:
        create_date = rec['GBSeq_create-date']
    except:
        create_date = ''

    try:
        rec_comment = rec['GBSeq_comment']
    except:
        rec_comment = ''

    try:
        institute = rec['GBSeq_references'][0]['GBReference_journal']
    except:
        institute = ''

    cds = [i for i in rec['GBSeq_feature-table'] if i['GBFeature_key'] == 'CDS']

    try:
        df = pd.concat([extract_cd_features(cd) for cd in cds]) #[extract_cd_features(cd) for cd in cds]
        df['seq_id'] = id_
        df['genbank_url'] = genbank_url
        df['create_date'] = create_date
        df['comment'] = rec_comment
        df['institute'] = institute
        return df
    except:
        pass

def extract_records(rec):
    seq_id = rec['GBSeq_locus']
    genbank_url = f'https://www.ncbi.nlm.nih.gov/nuccore/{seq_id}?report=genbank'
    create_date = rec['GBSeq_create-date']
    try:
        rec_comment = rec['GBSeq_comment']
    except:
        rec_comment = ''

    try:
        institute = rec['GBSeq_references'][0]['GBReference_journal']
    except:
        institute = ''

    genomic_seq = rec['GBSeq_sequence']
    return (seq_id, genbank_url, create_date, rec_comment, institute, genomic_seq)

def extract_record(accession_num):
    handle = Entrez.efetch(db="nuccore", id=accession_num, retmode="xml")
    rec = Entrez.read(handle)

    print(f"Extracting data for accession number ..... {accession_num}")

    try:
        genomic_data = [extract_records(r) for r in rec]
        genomic_df = pd.DataFrame(genomic_data, columns=['seq_id', 'genbank_url', 'created_date', 'comment', 'submission', 'genomic_seq'])
    #     genomic_df.to_csv(r'./sars_cov_2-genome.csv', mode='a')
        genomic_df.to_csv(r'./covid_genome.csv', mode='a', index=False, header=None)
    except:
        pass

    try:
        # Extract CDs from Genbank record
        cds_data = [extract_cds(r) for r in rec]
        cds_df = pd.concat(cds_data)
        cds_df.to_csv(r'./covid_cds.csv', mode='a', index=False, header=None)
    except:
        pass

    handle.close()

# Gathering list of accession numbers for SARS-CoV-2 genomes
genome_list_url = r'https://www.ncbi.nlm.nih.gov/sars-cov-2/download-nuccore-ids/'
genome_list = requests.get(genome_list_url)
accession_numbers = genome_list.text.split('\r\n')
del accession_numbers[0] # Remove the column label, which is just "ids"

try:
    df = pd.read_csv(r'./covid_genome.csv', header=None)
    df.columns = ['seq_id', 'genbank_url', 'created_date', 'comment', 'submission', 'genomic_seq']

    remaining_ids = list(set(accession_numbers) - set(df.seq_id.tolist()))
    remaining_ids = [i for i in remaining_ids if i !='']
except:
    remaining_ids = accession_numbers

# Main loop
[extract_record(a) for a in remaining_ids]
