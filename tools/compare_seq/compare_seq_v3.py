import argparse
import pandas as pd
import requests
from Bio import Entrez, SeqIO, pairwise2
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import time

Entrez.email = "your_email@example.com"  # Replace with your email

def fetch_uniprot_fasta(uniprot_id):
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            record = SeqIO.read(StringIO(response.text), "fasta")
            return str(record.seq)
    except:
        return None
    return None

def fetch_refseq_fasta(refseq_id):
    try:
        handle = Entrez.efetch(db="protein", id=refseq_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except:
        return None

def compute_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    if alignments:
        match_count = alignments[0].score
        identity = match_count / max(len(seq1), len(seq2))
        return round(identity * 100, 2)
    return 0.0

def compare_sequences(row):
    gene = row["Protein"]
    uniprot = row["Uniprot AC"]
    refseq = row["RefSeq_stripped"]

    if pd.isna(uniprot) or pd.isna(refseq):
        return ("not_found", (gene, uniprot, refseq, "missing ID"))

    print(f"🔍 Comparing {uniprot} vs {refseq} ...", end=" ")

    seq_uni = fetch_uniprot_fasta(uniprot)
    seq_ref = fetch_refseq_fasta(refseq)

    if not seq_uni or not seq_ref:
        print("❌ Sequence not found")
        return ("not_found", (gene, uniprot, refseq, "fetch error"))

    identity = compute_identity(seq_uni, seq_ref)
    print(f"✅ Identity = {identity}%")

    if identity == 100.0:
        return ("match", (gene, uniprot, refseq, identity))
    else:
        return ("mismatch", (gene, uniprot, refseq, identity))

def main():
    parser = argparse.ArgumentParser(description="Parallel comparison of UniProt vs RefSeq protein sequences.")
    parser.add_argument("-i", "--input", required=True, help="Path to MAVISp index CSV")
    parser.add_argument("-o", "--outdir", default="output", help="Directory to save results")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of parallel threads (default: 8)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input)
    df["RefSeq_stripped"] = df["RefSeq ID"].astype(str).str.extract(r"(NP_\d+)")

    matches, mismatches, not_found = [], [], []

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(compare_sequences, row) for _, row in df.iterrows()]
        for future in as_completed(futures):
            tag, result = future.result()
            if tag == "match":
                matches.append(result)
            elif tag == "mismatch":
                mismatches.append(result)
            else:
                not_found.append(result)

    pd.DataFrame(matches, columns=["Gene", "UniProt_AC", "RefSeq_ID", "Seq_Identity(%)"]).to_csv(
        f"{args.outdir}/matches.csv", index=False)
    pd.DataFrame(mismatches, columns=["Gene", "UniProt_AC", "RefSeq_ID", "Seq_Identity(%)"]).to_csv(
        f"{args.outdir}/mismatches.csv", index=False)
    pd.DataFrame(not_found, columns=["Gene", "UniProt_AC", "RefSeq_ID", "Error"]).to_csv(
        f"{args.outdir}/not_found.csv", index=False)

    print(f"\n✅ Completed. Results written to: {args.outdir}")

if __name__ == "__main__":
    main()

