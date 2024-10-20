from flask import Flask, render_template, request, session as flask_session, flash, redirect, url_for
from sqlalchemy import create_engine, Column, Integer, String, DateTime, func, text
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.orm import declarative_base
from contextlib import contextmanager
from datetime import datetime, timezone
from flask_mail import Mail, Message
import urllib.parse
import pandas as pd
import re
import time
import random
import string
import requests
import os
from Bio import Blast, Entrez, SeqIO

app = Flask(__name__)
#secret key for session management
app.secret_key = os.getenv('MIRNAPROTPRED_SECRET_KEY')

#Entrez EmailSetup
Entrez.email = os.getenv('MIRNAPROTPRED_EMAIL_USER')
# PostgreSQL setup
username = os.getenv('MIRNAPROTPRED_DB_USERNAME')
password = os.getenv('MIRNAPROTPRED_DB_PASSWORD')
host = os.getenv('MIRNAPROTPRED_DB_HOST')
port = os.getenv('MIRNAPROTPRED_DB_PORT')
dbname = os.getenv('MIRNAPROTPRED_DB_NAME')

# Email details
email_user = os.getenv('MIRNAPROTPRED_EMAIL_USER')
email_pass = os.getenv('MIRNAPROTPRED_EMAIL_PASS')

# SQLAlchemy setup
encoded_password = urllib.parse.quote(password)
DATABASE_URI = f'postgresql://{username}:{encoded_password}@{host}:{port}/{dbname}'
engine = create_engine(DATABASE_URI)
Session = scoped_session(sessionmaker(bind=engine))

# Configuration for Flask-Mail
app.config['MAIL_SERVER'] = 'smtp.gmail.com'
app.config['MAIL_PORT'] = 587
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USE_SSL'] = False
app.config['MAIL_USERNAME'] = email_user
app.config['MAIL_PASSWORD'] = email_pass
app.config['MAIL_DEFAULT_SENDER'] = email_user

mail = Mail(app)
Base = declarative_base()

# Define Visit model
class Visit(Base):
    __tablename__ = 'visits'
    id = Column(Integer, primary_key=True)
    ip_address = Column(String, unique=True)
    city = Column(String)
    country = Column(String)
    visit_count = Column(Integer, default=1)
    visit_time = Column(DateTime, default=datetime.now(timezone.utc))

# Genetic code dictionary
genetic_code = {
    "A": "GCU", "R": "CGU", "N": "AAU", "D": "GAU",
    "C": "UGU", "Q": "CAA", "E": "GAA", "G": "GGU",
    "H": "CAU", "I": "AUU", "L": "UUA", "K": "AAA",
    "M": "AUG", "F": "UUU", "P": "CCU", "S": "UCU",
    "T": "ACU", "W": "UGG", "Y": "UAU", "V": "GUU",
    "*": "UAA"
}

# Boyer-Moore Search Algorithm
def preprocess_bad_character_table(pattern):
    table = {}
    for i in range(len(pattern)):
        table[pattern[i]] = len(pattern) - i - 1
    return table

def boyer_moore_search(text, pattern):
    m = len(pattern)
    n = len(text)
    if m > n:
        return []

    table = preprocess_bad_character_table(pattern)
    occurrences = []

    i = m - 1
    while i < n:
        j = m - 1
        while j >= 0 and text[i] == pattern[j]:
            i -= 1
            j -= 1

        if j == -1:
            occurrences.append(i + 1)

        bad_char_skip = table.get(text[i], m)
        i += max(m - j, bad_char_skip)

    return occurrences

# Transcription and Translation
def transcription_translation(sequence):
    if all(nucleotide in "ATGC" for nucleotide in sequence):
        input_format = "DNA"
    elif all(nucleotide in "AUGC" for nucleotide in sequence):
        input_format = "RNA"
    elif all(amino_acid in "ARDNCEQGHILKMFPSTWYV-" for amino_acid in sequence):
        input_format = "Protein"
    else:
        raise ValueError("Invalid Input Sequecnce")

    if input_format == "Protein":
        ncbi_gi, ncbi_id, start_pos, end_pos = blast(sequence)
        print(f'NCBI GI: {ncbi_gi}')
        dna_sequence = retrieve_seq(ncbi_id, start_pos, end_pos)
    elif input_format == "RNA":
        dna_sequence = rna_to_dna(sequence)
    elif input_format == "DNA":
        dna_sequence = sequence
    else:
        raise ValueError("Invalid input format")

    return dna_sequence


def blast(sequence):
    print('Starting Blast')
    response = Blast.qblast(program='tblastn',
                            database='core_nt',
                            sequence=sequence,
                            format_type='XML')
    blast_record = Blast.read(response)
    ##getting the alignment of the top sequence
    alignment = blast_record[0][0].target.features[0].qualifiers.get('coded_by').split('|')
    ncbi_gi, ncbi_id, positions = alignment[1], alignment[3], alignment[4]
    start_pos, end_pos = positions.split(':')[1].split('..')
    if int(start_pos) > int(end_pos):
        start_pos, end_pos = end_pos, start_pos
    return ncbi_gi, ncbi_id, start_pos, end_pos

def retrieve_seq(ncbi_id, start_pos, end_pos):
    print(f'Retrieving Sequence with NCBI ID: {ncbi_id}')
    handle = Entrez.efetch(db='nucleotide', id=ncbi_id, rettype='fasta', retmode='text')
    record = SeqIO.read(handle, 'fasta')
    seq_of_interest = record.seq[int(start_pos)-1:int(end_pos)-1]
    return seq_of_interest

def protein_to_rna(protein_sequence):
    rna_sequence = "".join(genetic_code.get(aa, "") for aa in protein_sequence)
    return rna_sequence

def rna_to_dna(rna_sequence):
    return rna_sequence.replace("U", "T")

def fetch_rows():
    with session_scope() as session:
        rows = session.execute(text("SELECT * FROM mirdb")).fetchall()
    return rows

def search_miRNA_id_in_db(miRNA_id):
    with session_scope() as session:
        result = session.execute(text("SELECT * FROM mirdb WHERE mir_id = :miRNA_id"), {'miRNA_id': miRNA_id}).fetchall()
    return result

def get_geolocation(ip_address):
    token = os.getenv('GEOLOCATION_API_TOKEN', '39102577b08cef')
    url = f'https://ipinfo.io/{ip_address}/json?token={token}'
    response = requests.get(url)
    data = response.json()
    return {
        'city': data.get('city', 'Unknown'),
        'country': data.get('country', 'Unknown')
    }

@contextmanager
def session_scope():
    """Provide a transactional scope around a series of operations."""
    session = Session()
    try:
        yield session
        session.commit()
    except Exception as e:
        session.rollback()
        raise
    finally:
        session.close()

def get_client_ip():
    """Retrieve the client IP address, considering the possible use of proxies."""
    if request.headers.getlist("X-Forwarded-For"):
        ip = request.headers.getlist("X-Forwarded-For")[0].split(',')[0]
    else:
        ip = request.remote_addr
    return ip

@app.before_request
def track_user_location():
    if 'logged_ip' not in flask_session:
        ip_address = get_client_ip()
        app.logger.info(f"Client IP Address: {ip_address}")
        try:
            with session_scope() as db_session:
                visit = db_session.query(Visit).filter_by(ip_address=ip_address).first()
                if visit:
                    visit.visit_count += 1
                    visit.visit_time = datetime.now(timezone.utc)
                else:
                    geolocation_data = get_geolocation(ip_address)
                    visit = Visit(ip_address=ip_address, city=geolocation_data['city'], country=geolocation_data['country'])
                    db_session.add(visit)
            flask_session['logged_ip'] = True
        except Exception as e:
            app.logger.error(f"Error tracking user location: {e}")

def get_visit_counts():
    with session_scope() as session:
        visit_count = session.query(func.sum(Visit.visit_count)).scalar() or 0
        unique_visitors = session.query(Visit).count()
    return visit_count, unique_visitors

@app.route('/', methods=['GET', 'POST'])
def index():
    visit_count, unique_visitors = get_visit_counts()
    return render_template('index.html', visit_count=visit_count, unique_visitors=unique_visitors)

@app.route('/about_us', methods=['GET', 'POST'])
def about_us():
    visit_count, unique_visitors = get_visit_counts()
    return render_template('about_us.html', visit_count=visit_count, unique_visitors=unique_visitors)

@app.route('/contact_us', methods=['GET', 'POST'])
def contact_us():
    if request.method == 'POST':
        name = request.form['name']
        email = request.form['email']
        organization = request.form.get('organization', 'N/A')
        message = request.form['message']

        if not name or not email or not message:
            flash('Please fill out all required fields.')
            return redirect(url_for('contact_us'))

        try:
            msg = Message(
                subject=f"New Contact Us Message from {name}",
                sender=app.config['MAIL_DEFAULT_SENDER'],
                recipients=[app.config['MAIL_DEFAULT_SENDER']]
            )
            msg.body = f"""
            Name: {name}
            Email: {email}
            Organization: {organization}
            Message: {message}
            """
            mail.send(msg)
            flash('Your message has been sent successfully!')
            return redirect(url_for('contact_us'))

        except Exception as e:
            app.logger.error(f'Error sending email: {str(e)}')
            flash(f'An error occurred: {str(e)}')
            return redirect(url_for('contact_us'))
    visit_count, unique_visitors = get_visit_counts()
    return render_template('contact_us.html', visit_count=visit_count, unique_visitors=unique_visitors)

@app.route('/utr_prime', methods=['GET', 'POST'])
def utr_prime():
    visit_count, unique_visitors = get_visit_counts()
    if request.method == 'POST':
        try:
            sequence = request.form['sequence']
            job_id = request.form['job_id']
            sequence = re.sub(r'\s+', '', sequence).upper()
            dna_sequence = transcription_translation(sequence).lower()
            
            df = pd.DataFrame(fetch_rows(), columns=["r_id", "Description", "HumanmiRNAID", "Accession", "Sequence", "seed_1", "seed_2", "seed_3"])

            target_sequences = [seq for col in ["seed_1", "seed_2", "seed_3"] if col in df.columns for seq in df[col].dropna().tolist()]
            target_sequence_matches = [(seq, boyer_moore_search(dna_sequence, seq)) for seq in target_sequences if boyer_moore_search(dna_sequence, seq)]

            matching_rows_df = pd.DataFrame(columns=df.columns.tolist() + ["Seed", "Position"])
            added_rows = set()

            for target_sequence, occurrences in target_sequence_matches:
                matching_rows = df[df.isin([target_sequence]).any(axis=1)]
                matching_rows = matching_rows[~matching_rows.index.isin(added_rows)]
                
                if not matching_rows.empty:
                    for occurrence in occurrences:
                        occurrence_row = matching_rows.copy()
                        occurrence_row["Seed"] = target_sequence
                        occurrence_row["Position"] = occurrence + 1
                        matching_rows_df = pd.concat([matching_rows_df, occurrence_row])
                    added_rows.update(matching_rows.index)

            if matching_rows_df.empty:
                no_results_message = "No matching rows found based on your input."
                return render_template('utr_prime.html', no_results_message=no_results_message, visit_count=visit_count, unique_visitors=unique_visitors)

            display_rows = matching_rows_df.sort_values(by="Position", ascending=True)
            if len(display_rows) > 10:
                display_rows = display_rows.head(int(len(display_rows) * 0.2))

            matching_rows_df = matching_rows_df.sort_values(by="Position", ascending=True)
            matching_rows_df = matching_rows_df.iloc[:, 1:] 
            matching_rows_df = matching_rows_df.drop(["seed_1", "seed_2", "seed_3"], axis=1)

            filename = f"downloads/{job_id if job_id else 'matching_rows_' + str(int(time.time())) + '_' + ''.join(random.choices(string.ascii_letters, k=5))}.xlsx"
            matching_rows_df.to_excel(filename, index=False)

            return render_template('utr_prime.html', display_rows=display_rows.to_dict(orient='records'), excel_path=filename, job_id=job_id, visit_count=visit_count, unique_visitors=unique_visitors)
        except Exception as e:
            app.logger.error(f"Error in utr_prime: {e}")
            return render_template('utr_prime.html', error_message=str(e), visit_count=visit_count, unique_visitors=unique_visitors)
    
    return render_template('utr_prime.html', display_rows=None, visit_count=visit_count, unique_visitors=unique_visitors)

@app.route('/sequence_finder', methods=['GET', 'POST'])
def sequence_finder():
    visit_count, unique_visitors = get_visit_counts()
    if request.method == 'POST':
        try:
            sequence = request.form['sequence']
            job_id = request.form['job_id']
            sequence = re.sub(r'\s+', '', sequence).upper()
            dna_sequence = transcription_translation(sequence).lower()
            
            df = pd.DataFrame(fetch_rows(), columns=["r_id", "Description", "HumanmiRNAID", "Accession", "Sequence", "seed_1", "seed_2", "seed_3"])

            target_sequences = [seq for col in ["seed_1", "seed_2", "seed_3"] if col in df.columns for seq in df[col].dropna().tolist()]
            target_sequence_matches = [(seq, boyer_moore_search(dna_sequence, seq)) for seq in target_sequences if boyer_moore_search(dna_sequence, seq)]

            matching_rows_df = pd.DataFrame(columns=df.columns.tolist() + ["Seed", "Position"])
            added_rows = set()

            for target_sequence, occurrences in target_sequence_matches:
                matching_rows = df[df.isin([target_sequence]).any(axis=1)]
                matching_rows = matching_rows[~matching_rows.index.isin(added_rows)]
                
                if not matching_rows.empty:
                    for occurrence in occurrences:
                        occurrence_row = matching_rows.copy()
                        occurrence_row["Seed"] = target_sequence
                        occurrence_row["Position"] = occurrence + 1
                        matching_rows_df = pd.concat([matching_rows_df, occurrence_row])
                    added_rows.update(matching_rows.index)

            if matching_rows_df.empty:
                no_results_message = "No matching rows found based on your input."
                return render_template('sequence_finder.html', no_results_message=no_results_message, visit_count=visit_count, unique_visitors=unique_visitors)

            display_rows = matching_rows_df.sort_values(by="Position", ascending=True)
            if len(display_rows) > 10:
                display_rows = display_rows.head(int(len(display_rows) * 0.2))
                
            matching_rows_df = matching_rows_df.sort_values(by="Position", ascending=True)
            matching_rows_df = matching_rows_df.iloc[:, 1:] 
            matching_rows_df = matching_rows_df.drop(["seed_1", "seed_2", "seed_3"], axis=1)
            
            filename = f"downloads/{job_id if job_id else 'matching_rows_' + str(int(time.time())) + '_' + ''.join(random.choices(string.ascii_letters, k=5))}.xlsx"
            matching_rows_df.to_excel(filename, index=False)

            return render_template('sequence_finder.html', display_rows=display_rows.to_dict(orient='records'), excel_path=filename, job_id=job_id, visit_count=visit_count, unique_visitors=unique_visitors)
        except Exception as e:
            app.logger.error(f"Error in sequence_finder: {e}")
            return render_template('sequence_finder.html', error_message=str(e), visit_count=visit_count, unique_visitors=unique_visitors)
    
    return render_template('sequence_finder.html', display_rows=None, visit_count=visit_count, unique_visitors=unique_visitors)


@app.route('/validator', methods=['GET', 'POST'])
def validator():
    visit_count, unique_visitors = get_visit_counts()
    if request.method == 'POST':
        sequence = request.form['sequence']
        miRNA_ids_input = request.form['miRNA_ids']
        sequence = re.sub(r'\s+', '', sequence).upper()
        dna_sequence = transcription_translation(sequence).lower()
        
        miRNA_ids = ['hsa-' + id.strip() if not id.strip().startswith('hsa-') else id.strip() for id in miRNA_ids_input.split(',')]

        results = []
        for miRNA_id in miRNA_ids:
            search_result = search_miRNA_id_in_db(miRNA_id)
            
            if search_result:
                miRNA_info = {'id': miRNA_id, 'found': True, 'seed_match': False}
                seeds = search_result[0][-3:]
                miRNA_info['seed_match'] = any(seed in dna_sequence for seed in seeds)
                results.append(miRNA_info)
            else:
                results.append({'id': miRNA_id, 'found': False})
        
        return render_template('validator.html', results=results, visit_count=visit_count, unique_visitors=unique_visitors)
    
    return render_template('validator.html', results=[], visit_count=visit_count, unique_visitors=unique_visitors)

if __name__ == '__main__':
    app.run(debug=True)
