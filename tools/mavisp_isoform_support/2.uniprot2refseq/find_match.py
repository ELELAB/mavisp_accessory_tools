import requests
from Bio import pairwise2, SeqIO
from io import StringIO
import argparse
import csv
from time import sleep

def fetch_fasta(uniprot_id, source="uniprot"):
    if source == "uniprot":
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    else:  # use efetch for NCBI
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "protein",
            "id": uniprot_id,
            "rettype": "fasta",
            "retmode": "text"
        }
        r = requests.get(url, params=params, timeout=20)
        return r.text if r.ok and r.text.startswith(">") else None
    r = requests.get(url, timeout=20)
    return r.text if r.ok else None

def extract_sequence(fasta_text):
    try:
        fasta_io = StringIO(fasta_text)
        record = SeqIO.read(fasta_io, "fasta")
        return str(record.seq), len(record.seq)
    except Exception:
        return None, 0

def align_and_score(seq1, seq2):
    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB))
    return (matches / max(len(seq1), len(seq2))) * 100

def get_refseq_ids_from_uniprot(uac):
    url = f"https://rest.uniprot.org/uniprotkb/{uac}.json"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return [], "NA"
    data = r.json()
    entry_name = data.get("uniProtkbId", "NA")
    refseq_ids = [
        ref["id"].split(".")[0]
        for ref in data.get("uniProtKBCrossReferences", [])
        if ref["database"] == "RefSeq" and ref["id"].startswith("NP_")
    ]
    return list(set(refseq_ids)), entry_name

def main():
    parser = argparse.ArgumentParser(description="Compare UniProt canonical sequence to RefSeq NP_ sequences.")
    parser.add_argument("-u", "--uniprot", required=True, help="UniProt Accession Code (e.g. Q9UKU0)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV filename")
    args = parser.parse_args()

    uac = args.uniprot
    refseq_ids, entry_name = get_refseq_ids_from_uniprot(uac)
    if not refseq_ids:
        print(f"❌ No NP_ RefSeq IDs found for {uac}")
        return

    print(f"🔍 Retrieved {len(refseq_ids)} RefSeq ID(s) for {uac} ({entry_name})")

    uniprot_fasta = fetch_fasta(uac, source="uniprot")
    if not uniprot_fasta:
        print(f"❌ Failed to download UniProt sequence for {uac}")
        return
    uniprot_seq, uniprot_len = extract_sequence(uniprot_fasta)

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", "UniProt_AC", "UniProt_ID", "RefSeq_ID", "Sequence_Identity(%)", "UniProt_seq_length", "RefSeq_seq_length"])

        for refseq_id in refseq_ids:
            print(f"🔍 Comparing {uac} vs {refseq_id}")
            refseq_fasta = fetch_fasta(refseq_id, source="refseq")
            if not refseq_fasta:
                print(f"⚠️  Could not retrieve RefSeq sequence for {refseq_id}")
                continue
            refseq_seq, refseq_len = extract_sequence(refseq_fasta)
            if not refseq_seq:
                print(f"⚠️  RefSeq sequence for {refseq_id} was invalid or empty.")
                continue
            identity = align_and_score(uniprot_seq, refseq_seq)
            print(f"✅ Identity: {identity:.2f}% | Lengths: UniProt={uniprot_len}, RefSeq={refseq_len}")
            writer.writerow([entry_name, uac, entry_name, refseq_id, f"{identity:.2f}", uniprot_len, refseq_len])
            sleep(0.3)  # respectful pause for API

if __name__ == "__main__":
    main()

