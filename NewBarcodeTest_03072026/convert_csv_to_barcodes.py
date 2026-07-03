"""Build a levseq-compatible barcodes.fasta from the IDT oligo CSVs.

The IDT oligos contain a shared adapter tail after the variable barcode:
  - NB oligos end with ATCTCGATCCCGCGAAATTAATACGACTCAC (31 bp T7-promoter-side adapter)
  - RB oligos end with GCCCCAAGGGGTTATGCTAGTTATTGCTC (29 bp adapter)

LevSeq's C++ demultiplexer computes each barcode's "perfect" alignment score
from the full FASTA sequence and requires the read's front/rear window to
match >= 70% of that perfect score. If the shared adapter tail is included
in the FASTA, the perfect-score denominator is inflated by ~30 bp of
non-discriminating sequence and the 70% threshold becomes unreachable, so
every read is classified as `unclassified` (see barcoding_summary.txt).

Strip those constant tails here so only the variable 24 bp head remains,
matching the format of the bundled levseq/barcoding/minion_barcodes.fasta.
"""
import pandas as pd

# Shared adapter tails observed on every oligo. Strip if present so the
# barcodes.fasta contains only the discriminating region.
NB_TAIL = "ATCTCGATCCCGCGAAATTAATACGACTCAC"
RB_TAIL = "GCCCCAAGGGGTTATGCTAGTTATTGCTC"


def strip_tail(seq: str, tail: str) -> str:
    seq = seq.strip().upper()
    return seq[: -len(tail)] if seq.endswith(tail) else seq


df = pd.read_csv("forward_upload_IDT_100uM_DNA_oligos.csv").dropna(subset=["Name", "Sequence"])
df2 = pd.read_csv("backward_upload_IDT_100uM_DNA_oligos.csv").dropna(subset=["Name", "Sequence"])

with open("barcodes.fasta", "w") as f:
    for _, row in df.iterrows():
        f.write(f">{row['Name']}\n{strip_tail(row['Sequence'], NB_TAIL)}\n")
    for _, row in df2.iterrows():
        f.write(f">{row['Name']}\n{strip_tail(row['Sequence'], RB_TAIL)}\n")
