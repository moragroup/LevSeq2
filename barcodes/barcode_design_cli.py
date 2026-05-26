#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Norsedrunkensailor"
__license__ = "MIT"
__status__ = "beta"
__version__ = "1.0"


"""
This code use a determinantal point process approach to design distance and hybridization aware subsets of barcodes.
The DPP is an approximate optimal solution to subset selection based on kernels (which is NP-hard)

Potential cross- or self-hybridization of primers which may cause cPCR amplification problems is 
actively addressed in this subset selection. Original LevSeq primers sometimes had non-uniform amplification due to this kind of NB-RB or self-hybridization.
The method here takes the sequencing adapters/anchors and ensures that the barcodes jointly
   (a) maximise pair-wise Levinstein distance (to make demultiplexing very robust to noise) 
   (b) ensure no significant free-energy contributions between primer pairs (hybridization)
   
Example CLI command:

python barcode_design_cli.py \
    --fwd_anchor "CTCGATCCCGCGAAATTAATACG" \
    --rev_anchor "ATCCGGATATAGTTCCTCCTTTCAG" \
    --barcode_len 24 \
    --num_rb_select 12 \
    --num_nb_select 8 \
    --outdir "barcode_design_T7" \
    --name "barcode_T7_pool"

    
Where --fwd_anchor is the 5' sequencing primer to which we extent it at the 5' end with barcodes.
      --rev_anchor is the 3' sequencing primer (5'-3' direction) which will be extended in the same way with RB barcodes

      
For most applications --barcode_len 24 is perfect. Too short barcodes risk false assignmnet during demultiplexing
--num_rb_select 12 
--num_nb_select 8
This is the setup up to cover exactly one plate with 20 primers total. 8 for the rows and 12 for the columns.


All the best,
NDS
"""

# ==============================================================================

# ==============================================================================
import sys
import os
import datetime
import argparse
import csv
import os
from tqdm.auto import tqdm
import math
import re

import pandas as pd
import numpy as np

from multiprocess import Pool, cpu_count
from joblib import Parallel, delayed


import simsimd
from Bio.Align import PairwiseAligner

import RNA
from dppy.finite_dpps import FiniteDPP

# ==============================================================================
            # Trivial Functions First.
# ==============================================================================


def parse_fasta(filename):
  '''function to parse fasta file'''

  header, sequence = [],[]
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if len(line)>0:
      if line[0] == ">":
        header.append(line[1:])
        sequence.append([])
      else:
        line = line.upper()
        sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]

  return header, sequence


def mk_msa(seqs):
    '''one hot encode seq using float128
       this is for vectorised kernels:D 
    '''
    alphabet = "ATCG"
    states = len(alphabet)
    a2n = {a: n for n, a in enumerate(alphabet)}
    msa_ori = np.array([[a2n.get(aa, states-1) for aa in seq] for seq in seqs])
    return msa_ori, np.eye(states, dtype=np.float64)[msa_ori]


def reverse_complement(sequences):
    """
    Given a list of DNA sequences, returns the reverse complement for each.
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return [seq.translate(complement)[::-1] for seq in sequences]


def return_similarity_simsimd_ohe(sequencesA, sequencesB=None):
    """
    Calculates pairwise dissimilarity using SimSIMD inner product on one-hot encodings with float128.
    """
    # One-hot encode both sets
    _, msaA_ohe = mk_msa(sequencesA)
    if sequencesB is None:
        msaB_ohe = msaA_ohe
    else:
        _, msaB_ohe = mk_msa(sequencesB)
    
    nA, L, _ = msaA_ohe.shape
    nB = msaB_ohe.shape[0]

    # Flatten one-hot
    A_flat = msaA_ohe.reshape(nA, -1).astype(np.float64)
    B_flat = msaB_ohe.reshape(nB, -1).astype(np.float64)

    #inner product distance
    ip_distance_matrix = np.array(simsimd.cdist(A_flat, B_flat, metric="inner"))
    return ip_distance_matrix

# ==============================================================================
            # Alignment Kernel
# ==============================================================================


#we initialise gloablly;D
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = -1

def _pairwise_similarity_worker(seq_pair):
    seq1, seq2 = seq_pair
    num_matches = aligner.score(seq1, seq2)
    alignment_len = max(len(seq1), len(seq2))
    return num_matches / alignment_len if alignment_len > 0 else 0.0


def create_kernel_align(sequencesA, sequencesB=None, n_jobs=None):
    """
    Compute dissimilarity matrix between two sets of sequences.
    Uses symmetry speed-up if sequencesA is sequencesB.
    Args:
        sequencesA (list[str]): First list of sequences.
        sequencesB (list[str] or None): Second list of sequences. 
                                        If None, compares A vs A.
        n_jobs (int or None): Number of worker processes.
    Returns:
        np.ndarray: Dissimilarity matrix.
    """

    if sequencesB is None:
        sequencesB = sequencesA
        symmetric = True
    else:
        symmetric = sequencesA is sequencesB or sequencesA == sequencesB
    
    nA, nB = len(sequencesA), len(sequencesB)
    if n_jobs is None:
        n_jobs = cpu_count()
    
    if symmetric:  
        # only compute upper triangle (i < j)
        pairs = [(sequencesA[i], sequencesB[j]) 
                 for i in range(nA) for j in range(i+1, nB)]
    else:
        # full Cartesian product
        pairs = [(sequencesA[i], sequencesB[j]) for i in range(nA) for j in range(nB)]
    
    # Parallel execution
    with Pool(n_jobs) as pool:
        sims = list(tqdm(pool.imap(_pairwise_similarity_worker, pairs),
                         total=len(pairs), desc="Alignments"))
    
    # Fill results
    kernel = np.zeros((nA, nB))
    if symmetric:
        k = 0
        for i in range(nA):
            kernel[i, i] = 1.0
            for j in range(i+1, nB):
                kernel[i, j] = sims[k]
                kernel[j, i] = sims[k]
                k += 1
    else:
        kernel = np.array(sims).reshape(nA, nB)
    
    return kernel

# ==============================================================================
            # Free-energy Calcs
# ==============================================================================

def get_self_interaction_mfes(primer_sequences: list[str]) -> list[tuple[float, float]]:
    """Calculates MFE for hairpin and self-dimer structures."""
    mfe_values = []
    RNA.params_load_DNA_Mathews2004()
    for seq in primer_sequences:
        if not seq:
            mfe_values.append((0.0, 0.0))
            continue
        _, mfe_monomer = RNA.fold(seq)
        combined_seq = f"{seq}&{seq}"
        _, mfe_dimer = RNA.cofold(combined_seq)
        mfe_values.append((mfe_monomer, mfe_dimer))
    return mfe_values

def calculate_psd_scores(primer_sequences: list[str], k: float = 0.4) -> list[tuple[float, float]]:
    """Calculates PSD values where a higher score is better (less stable)."""
    interaction_mfes = get_self_interaction_mfes(primer_sequences)
    psd_scores = []
    for hairpin_mfe, dimer_mfe in interaction_mfes:
        psd_hairpin = math.exp(k * hairpin_mfe)
        psd_dimer = math.exp(k * dimer_mfe)
        psd_scores.append((psd_hairpin, psd_dimer))
    return psd_scores, interaction_mfes

def calculate_diagonal_terms(primer_sequences: list[str], k: float = 0.4, if_timing: bool = False) -> list[float]:
    """Calculates the final 'goodness' scores for the diagonal."""
    psd_scores, interaction_mfes = calculate_psd_scores(primer_sequences, k)
    diagonal_values = []
    
    for psd_hairpin, psd_dimer in tqdm(psd_scores, total=len(psd_scores), desc="Diagonal terms", disable=not if_timing):
        joint_score = psd_hairpin * psd_dimer
        diagonal_values.append(joint_score)
    return diagonal_values, interaction_mfes

def create_primer_kernel(primer_sequences: list[str], k: float = 0.4) -> np.ndarray:
    """
    Creates a primer interaction kernel from a list of sequences.
    - Off-diagonal elements represent the cross-dimer conflict score (higher is worse).
    - Diagonal elements represent the primer's 'goodness' score (higher is better).

    Args:
        primer_sequences: A list of DNA primer sequences.
        k: A "temperature" scaling factor -- this uses the squared exponential trick.
    Returns:
        A numpy nxn array representing the kernel matrix.
        and some other things used in dev. 
    """

    n = len(primer_sequences)
    # Start with a matrix of zeros
    kernel = np.zeros((n, n))
    kernel_raw = np.zeros((n, n))

    # --- Step 1: Pre-calculate monomer MFEs for efficiency ---
    RNA.params_load_DNA_Mathews2004()
    monomer_mfes = [RNA.fold(seq)[1] for seq in primer_sequences]

    # --- Step 2: Calculate and populate off-diagonal elements ---
    total_iters = n * (n - 1) // 2  # number of (i, j) pairs with j > i

    with tqdm(total=total_iters, desc="ViennaFold") as pbar:
        for i in range(n):
            for j in range(i + 1, n):
                seq_i, seq_j = primer_sequences[i], primer_sequences[j]

                # Calculate cross-dimer MFE
                _, cross_dimer_mfe = RNA.cofold(f"{seq_i}&{seq_j}")

                # Calculate interaction energy
                interaction_mfe = cross_dimer_mfe - (monomer_mfes[i] + monomer_mfes[j])

                kernel_raw[i,j] = interaction_mfe
                kernel_raw[j,i] = interaction_mfe

                # Convert to a score where strong interaction -> high score
                try:
                    score = 1 / (1 + math.exp(k * interaction_mfe))
                except OverflowError:
                    score = 0.0

                # Populate the symmetric elements
                kernel[i, j] = score
                kernel[j, i] = score

                pbar.update(1)

    # --- Step 3: Calculate and populate diagonal elements ---
    diagonal_scores, interaction_mfes = calculate_diagonal_terms(primer_sequences, k)
    for i in range(n):
        kernel[i, i] = diagonal_scores[i]

    return kernel, np.array(kernel_raw), np.array(diagonal_scores), np.array(interaction_mfes)


# ==============================================================================
            # sequence_levenshtein !!! The Nanopore friendly metric. 
# ==============================================================================


def levenshtein_matrix(s, t):
    """Standard Wagner–Fischer matrix for LD."""
    m, n = len(s), len(t)
    d = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        d[i][0] = i
    for j in range(n + 1):
        d[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s[i - 1] == t[j - 1] else 1
            d[i][j] = min(
                d[i - 1][j] + 1,       # deletion
                d[i][j - 1] + 1,       # insertion
                d[i - 1][j - 1] + cost # substitution
            )
    return d


def sequence_levenshtein(s, t):
    """Compute Sequence-Levenshtein Distance (SLD).
        See: PMC10614987 
             https://pmc.ncbi.nlm.nih.gov/articles/PMC10614987/
    """
    d = levenshtein_matrix(s, t)
    m, n = len(s), len(t)

    last_row_min = min(d[m])                # min over last row
    last_col_min = min(row[n] for row in d) # min over last column

    return min(last_row_min, last_col_min)


def create_kernel_levinstein(sequencesA, sequencesB=None, n_jobs=None):
    """
    Compute dissimilarity matrix between two sets of sequences.
    Uses Sequence-Levenshtein Distance (SLD).
    
    Symmetry optimization if sequencesA is sequencesB.
    Progress bar enabled.

    Args:
        sequencesA (list[str]): First list of sequences.
        sequencesB (list[str] or None): Second list of sequences.
                                        If None, compares A vs A.
        n_jobs (int or None): Number of worker processes.

    Returns:
        np.ndarray: (len(sequencesA), len(sequencesB)) dissimilarity matrix
    """

    if sequencesB is None:
        sequencesB = sequencesA
        symmetric = True
    else:
        symmetric = False

    A, B = len(sequencesA), len(sequencesB)
    kernel = np.zeros((A, B), dtype=int)

    if symmetric:
        tasks = [(i, j) for i in range(A) for j in range(i, B)]

        def compute_pair(i, j):
            return sequence_levenshtein(sequencesA[i], sequencesB[j])

        results = Parallel(n_jobs=n_jobs)(
            delayed(compute_pair)(i, j)
            for i, j in tqdm(tasks, desc="Computing SLD (symmetric)")
        )

        # fill matrix symmetrically
        k = 0
        for i, j in tasks:
            kernel[i, j] = results[k]
            kernel[j, i] = results[k]
            k += 1
    else:
        tasks = [(i, j) for i in range(A) for j in range(B)]

        def compute_pair(i, j):
            return sequence_levenshtein(sequencesA[i], sequencesB[j])

        results = Parallel(n_jobs=n_jobs)(
            delayed(compute_pair)(i, j)
            for i, j in tqdm(tasks, desc="Computing SLD (rectangular)")
        )

        kernel = np.array(results).reshape(A, B)

    return kernel


#test_seq = ["ACGT", "ACG", "GT", "GTT"]
#K = create_kernel_levinstein(test_seq, n_jobs=-1)
#print(K)


# ==============================================================================
            # Generate Psuedorandom Unique Barcodes
# ==============================================================================


def gc_content(seq: str) -> float:
    """Calculate GC content as fraction (0-1)."""
    return (seq.count("G") + seq.count("C")) / len(seq)

def has_homopolymer(seq: str, max_run: int = 3) -> bool:
    """Check if sequence has homopolymer longer than max_run. 
        This follows standard barcode design rules laid out in: Johnson et al. https://pmc.ncbi.nlm.nih.gov/articles/PMC10276077/
    """
    return bool(re.search(r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (max_run+1, max_run+1, max_run+1, max_run+1), seq))

def generate_random_barcodes_with_adapter(
    adapter: str,
    n_gen: int = 20000,
    length: int = 24,
    keep_top: int = 1000,
    k: float = 0.4,
):
    """
    Generate pseudo random DNA barcodes
        filter by GC and homopolymers, (Johnson et al. https://pmc.ncbi.nlm.nih.gov/articles/PMC10276077/)
        append adapter, score them (to find hairpins or self-primer-dimers) 
        keep highest-scoring ones. (Those without problems)
    Returns barcodes WITHOUT adapter.
    """

    alphabet = np.array(list("ATCG"))
    gen = alphabet[np.random.randint(0, 4, size=(n_gen, length))]
    random_barcodes = ["".join(seq) for seq in gen]


    #Filter by GC content and homopolymers
    filtered_barcodes = [
        b for b in random_barcodes
        if 0.2 <= gc_content(b) <= 0.6 and not has_homopolymer(b, max_run=3)
    ]
    print(f"Filtered {len(filtered_barcodes)} / {len(random_barcodes)} barcodes passed GC/homopolymer rules.")

    #Append adapter for scoring
    barcodes_with_adapter = [b + adapter for b in filtered_barcodes]

    #Score with diagonal terms MFE Vienna
    diagonal_scores, _ = calculate_diagonal_terms(barcodes_with_adapter, k=k, if_timing=True)
    diagonal_scores = np.array(diagonal_scores)

    #Keep top highest-scoring barcodes
    top_idx = np.argsort(-diagonal_scores)[:keep_top]
    top_barcodes = [filtered_barcodes[i] for i in top_idx]
    top_scores = diagonal_scores[top_idx].tolist()

    return top_barcodes, top_scores


# ==============================================================================
            # Construct All Kernel
# ==============================================================================

import numpy as np

def grab_K_levin(all_barcodes, split_idx=1000, n_jobs=-1, jitter=1e-6):
    
    # 1. Setup
    if split_idx is None:
        split_idx = len(all_barcodes) // 2

    NB_barcodes = all_barcodes[:split_idx] 
    RB_barcodes = all_barcodes[split_idx:] 
    
    n_NB = len(NB_barcodes)
    n_RB = len(RB_barcodes)
    n_total = n_NB + n_RB


    # A. Self Blocks
    D_NB = create_kernel_levinstein(NB_barcodes, n_jobs=n_jobs)
    D_RB = create_kernel_levinstein(RB_barcodes, n_jobs=n_jobs)

    # B. Cross Block (Forward vs Reverse Complement)
    # Ensure RB is reverse complemented to align orientation with NB
    D_NB_RB = create_kernel_levinstein(
        NB_barcodes, 
        sequencesB=reverse_complement(RB_barcodes), 
        n_jobs=n_jobs
    )

    # 3. Assemble the Global Distance Matrix
    D_full = np.zeros((n_total, n_total))

    # Top-left: NB-NB
    D_full[:n_NB, :n_NB] = D_NB

    # Bottom-right: RB-RB
    D_full[n_NB:, n_NB:] = D_RB

    # Off-diagonal: NB-RB / RB-NB
    D_full[:n_NB, n_NB:] = D_NB_RB        # top-right
    D_full[n_NB:, :n_NB] = D_NB_RB.T      # bottom-left

    
    off_diag_mask = ~np.eye(n_total, dtype=bool)
    sigma = np.median(D_full[off_diag_mask])
    
    # Safety check to prevent division by zero if all sequences are identical
    if sigma == 0:
        sigma = 1.0

    
    gamma = 1.0 / (2 * (sigma ** 2))
    K_levin = np.exp(-gamma * (D_full ** 2))
    
    #np.fill_diagonal(K_levin, 1.0)
    #K_levin[np.diag_indices_from(K_levin)] += jitter

    return K_levin


def grab_K_fold(sequences_with_adapter):

    K5_fold, K5_raw_fold, K5_diagonal_scores_fold, K5_interaction_mfes_fold = create_primer_kernel(sequences_with_adapter, k=0.1)

    K5_fold = K5_fold

    block_size=1000
    mask_ones = np.ones_like(K5_fold, dtype=bool)

    #Exclude top-right and bottom-left blocks from the mask
    mask_ones[:block_size, block_size:] = False
    mask_ones[block_size:, :block_size] = False

    #Apply the mask: set everything else to 0 ## ADDITIVE KERNEL -- we do not care about NB-NB self-interaction (they never see each other)
    K5_fold[mask_ones] = 1
    np.fill_diagonal(K5_fold, K5_diagonal_scores_fold)

    return K5_fold, K5_raw_fold, K5_interaction_mfes_fold


# ==============================================================================
            # DPP Prep
# ==============================================================================


def get_my_dpp_ready_kernel(anchor_fwd,anchor_rev,barcode_len=24):

    fwd_barcodes, _ = generate_random_barcodes_with_adapter(
        adapter=anchor_fwd,
        n_gen=10000,
        length=barcode_len,
        keep_top=1000,
        k=0.1
    )

    rev_barcodes, _ = generate_random_barcodes_with_adapter(
        adapter=anchor_rev,
        n_gen=10000,
        length=barcode_len,
        keep_top=1000,
        k=0.1
    )

    fwd_headers = [f"NB{str(i).zfill(2)}" for i in range(1, len(fwd_barcodes)+1)]
    rev_headers = [f"RB{str(i).zfill(2)}" for i in range(1, len(rev_barcodes)+1)]

    all_barcodes = fwd_barcodes + rev_barcodes
    all_headers = fwd_headers + rev_headers
    sequences_with_adapter = []

    for header, seq in zip(all_headers, all_barcodes):
        if "NB" in header:
            sequences_with_adapter.append(seq + anchor_fwd)
        elif "RB" in header:
            sequences_with_adapter.append(seq + anchor_rev)
        else:
            # If no NB/RB in header, leave sequence unchanged (or handle differently)
            sequences_with_adapter.append(seq)


    K_levin =  grab_K_levin(all_barcodes)
    #K5_fold_matrix, K5_raw_fold, K5_interaction_mfes = grab_K_fold(sequences_with_adapter)

    #Hadamard Product
    K =  K_levin #*K5_fold_matrix


    # Compute eigen-decomposition
    #This is needed due to our block-diagonal tricks.
    eigvals, eigvecs = np.linalg.eigh(K)


    # Clip negative eigenvalues to zero
    eigvals_clipped = np.clip(eigvals, a_min=0, a_max=None)
    K_psd = eigvecs @ np.diag(eigvals_clipped) @ eigvecs.T
    print("Min eigenvalue after clipping:", np.linalg.eigvalsh(K_psd).min())

    return K_psd, all_barcodes, all_headers, sequences_with_adapter


# ==============================================================================
            # DPP Select
# ==============================================================================


def safe_logdet(mat):
    """Return log(abs(det(mat))) or -np.inf if determinant is non-positive/invalid."""
    sign, logabsdet = np.linalg.slogdet(mat)
    return logabsdet


def sample_dpps_search(
    extended_barcodes, 
    all_headers, 
    sequences_with_adapter,
    K_psd,
    n_subset,
    num_samples_to_check=10000,
    forward_threshold=1000,
    max_determinant_for_set_init=-10000,
    show_progress=True
):
    """
    Sampling search using an exact-k DPP sampler and track the best subset
    by the log-determinant of the submatrix K_psd.

    Parameters
    ----------
    K_psd : (n, n) ndarray
        PSD kernel / L-matrix used by the DPP.
    n_subset : int
        Size of subset to sample each iteration.
    num_samples_to_check : int
        Number of samples (iterations) to draw.
    forward_threshold : int
        Indices < forward_threshold are considered NB (forward); >= -> RB (reverse).
    max_determinant_for_set_init : float
        Initial threshold for the maximum log-determinant.
    show_progress : bool
        Whether to display a tqdm progress bar.
    """


    # Initialize DPP
    dpp = FiniteDPP(kernel_type='likelihood', L=K_psd)

    best_subset = None
    max_logdet = max_determinant_for_set_init
    num_valid_subsets_found = 0

    forward_counts_in_valid = []
    reverse_counts_in_valid = []
    all_valid_logdets = []
    all_checked_logdets = []

    iterator = range(num_samples_to_check)
    if show_progress:
        iterator = tqdm(iterator, desc="Sampling DPPs")

    for _ in iterator:
        candidate_subset = dpp.sample_exact_k_dpp(size=n_subset)

        # compute logdet of the corresponding submatrix
        submat = K_psd[np.ix_(candidate_subset, candidate_subset)]
        logdet = safe_logdet(submat)
        all_checked_logdets.append(logdet)

        # Count forward/reverse barcodes
        nb_count = sum(1 for i in candidate_subset if i < forward_threshold)
        rb_count = len(candidate_subset) - nb_count

        # If logdet is finite, it's a valid subset for our purposes
        if np.isfinite(logdet):
            num_valid_subsets_found += 1
            forward_counts_in_valid.append(nb_count)
            reverse_counts_in_valid.append(rb_count)
            all_valid_logdets.append(logdet)

            # Update best subset if this logdet is better than current max
            if logdet > max_logdet:
                print(f"\nLogDet Improved: {logdet}")
                print(f"Best Subset Now: {candidate_subset}")
                print(f"NB Count: {nb_count}, RB Count: {rb_count} \n")
                max_logdet = logdet
                best_subset = list(candidate_subset)

    res = {
        "best_subset": best_subset,
        "max_logdet": max_logdet,
        "num_valid_subsets_found": num_valid_subsets_found,
        "forward_counts_in_valid": forward_counts_in_valid,
        "reverse_counts_in_valid": reverse_counts_in_valid,
        "all_valid_logdets": all_valid_logdets,
        "all_checked_logdets": all_checked_logdets,
    }


    print("Best logdet:", res["max_logdet"])
    print("Best subset:", res["best_subset"])
    print("Valid subsets found:", res["num_valid_subsets_found"])


    best_subset = res["best_subset"]
    sort_best_subset = sorted(best_subset)
    print("Selected Barcodes:", sort_best_subset)
    print("NB Counts", sum(1 for i in res["best_subset"] if i < 1000))
    print("RB Counts", len(res["best_subset"]) - sum(1 for i in res["best_subset"] if i < 1000))
    selected_barcodes = [extended_barcodes[i] for i in sort_best_subset]
    selected_headers = [all_headers[i] for i in sort_best_subset]
    selected_barcodes_adapters = [sequences_with_adapter[i] for i in sort_best_subset]
    
    
    return res, selected_headers, selected_barcodes, selected_barcodes_adapters


# ==============================================================================
            # Final Selection. 
# ==============================================================================


def process_barcodes_with_adapters(selected_headers, selected_barcodes, anchor_fwd, anchor_rev):
    
    # Check if input lists are of equal length
    if len(selected_headers) != len(selected_barcodes):
        print("Error: Header and barcode lists have different lengths.", file=sys.stderr)
        return [], [] # Return empty lists to indicate failure

    print(f"Processing Selected {len(selected_headers)} entries...")

    adapter_fwd = anchor_fwd.upper()
    adapter_rev = anchor_rev.upper()

    all_headers_selected = []
    sequences_with_adapter = []

    # Iterate through the headers and sequences together
    for header, seq in zip(selected_headers, selected_barcodes):
        if "NB" in header:
            all_headers_selected.append(header)
            sequences_with_adapter.append(seq + adapter_fwd)
        elif "RB" in header:
            all_headers_selected.append(header)
            sequences_with_adapter.append(seq + adapter_rev)
        else:
            print(f"Error: Header '{header}' contains neither 'NB' nor 'RB'. Skipping.", file=sys.stderr)

    print(f"Successfully processed {len(all_headers_selected)} entries.")
    return all_headers_selected, sequences_with_adapter



import numpy as np
import os

def analyze_primer_interactions(
    interaction_matrix,
    all_headers,
    outdir,
    plot_filename_base="selected_heatmap_ranked",
    plot_title="Interaction Heatmap: RB vs NB Primers (Ranked by Mean MFE)",
    num_rb_select=12,
    num_nb_select=8
):
    """
    Analyzes primer interactions, plots a ranked heatmap, and returns the worst-interacting primers.
    Worst interacting primers are those which are less likely to cause cPCR problems!
     
    Args:
        interaction_matrix (np.ndarray): The full, raw interaction matrix.
        all_headers (list): List of all primer headers
        outdir (str): Path to the directory where plots will be saved.
        plot_filename_base (str, optional): Base name for the saved .png and .html files.
        plot_title (str, optional): Title for the heatmap plot.
        num_rb_select (int, optional): The number of "worst" RB (reverse) primers to select.
        num_nb_select (int, optional): The number of "worst" NB (forward) primers to select.

    Returns:
        tuple: A tuple containing two lists:
            - (list): The selected NB (forward) primer names.
            - (list): The selected RB (reverse) primer names.
    """
    
    # --- Setup and Data Preparation ---
    print(f"Analyzing interactions for {len(all_headers)} primers...")
    
    # Prep
    matrix_to_plot = interaction_matrix.copy()
    np.fill_diagonal(matrix_to_plot, np.nan)
    primer_labels = [h.lstrip(">") for h in all_headers]


    # --- Separate primers ---
    rb_primers = [h.lstrip(">") for h in all_headers if "RB" in h]
    nb_primers = [h.lstrip(">") for h in all_headers if "NB" in h]

    if not rb_primers or not nb_primers:
        print("Error: Could not find both 'RB' and 'NB' primers in headers.")
        return [], []

    # Indices
    rb_indices = [primer_labels.index(p) for p in rb_primers]
    nb_indices = [primer_labels.index(p) for p in nb_primers]
    matrix_focus = matrix_to_plot[np.ix_(rb_indices, nb_indices)]

    # --- Rank RB primers (rows) ---
    mean_mfe_per_rb = np.nanmean(matrix_focus, axis=1)
    ranked_indices_rb = np.argsort(mean_mfe_per_rb)
    matrix_focus = matrix_focus[ranked_indices_rb, :]
    rb_primers_ranked = [rb_primers[i] for i in ranked_indices_rb]
    mean_mfe_per_rb_ranked = mean_mfe_per_rb[ranked_indices_rb]

    # --- Rank NB primers (columns) --- Rank NB primers by mean MFE
    mean_mfe_per_nb = np.nanmean(matrix_focus, axis=0) 
    ranked_indices_nb = np.argsort(mean_mfe_per_nb)[::-1] #to reverse
    matrix_focus = matrix_focus[:, ranked_indices_nb]

    #Ranked list of NB primers and their means
    nb_primers_ranked = [nb_primers[i] for i in ranked_indices_nb]
    mean_mfe_per_nb_ranked = mean_mfe_per_nb[ranked_indices_nb]


    # Worst is Best here ;)
    _worst_rb_primers_asc = rb_primers_ranked[-num_rb_select:]
    _worst_rb_mfes_asc = mean_mfe_per_rb_ranked[-num_rb_select:]
    worst_rb_primers_desc = _worst_rb_primers_asc[::-1]
    worst_rb_mfes_desc = _worst_rb_mfes_asc[::-1]
    worst_nb_primers_desc = nb_primers_ranked[:num_nb_select] # Already in descending order
    worst_nb_mfes_desc = mean_mfe_per_nb_ranked[:num_nb_select] # Already in descending order

    print(f"\n--- {num_rb_select} RB Primers with best Mean MFE ---")
    print("List format:", worst_rb_primers_desc)
    for i, primer in enumerate(worst_rb_primers_desc):
        print(f"{i+1:2d}: {primer} (Mean MFE: {worst_rb_mfes_desc[i]:.4f})")

    print(f"\n--- {num_nb_select} NB Primers with best Mean MFE ---")
    print("List format:", worst_nb_primers_desc)
    for i, primer in enumerate(worst_nb_primers_desc):
        print(f"{i+1:2d}: {primer} (Mean MFE: {worst_nb_mfes_desc[i]:.4f})")
    
    return worst_nb_primers_desc, worst_rb_primers_desc


# ==============================================================================
            # Final Output. 
# ==============================================================================

def resolve_and_write_barcodes_csv(
    select_list_fwd: list,
    select_list_rv: list,
    selected_headers: list,
    selected_barcodes: list,
    selected_barcodes_adapters: list,
    out_csv_file: str
) -> tuple:
    """
    Resolves requested header names against a list of selected headers,
    collects corresponding barcodes, and writes all information to a CSV file.
    """

    # --- Nested helper function ---
    def resolve_ordered(request_list, mapping):
        """Helper to resolve a request list to indices, preserving order."""
        resolved = []
        not_found = []
        # Create a case-insensitive mapping
        lower_map = {k.lower(): v for k, v in mapping.items()}

        for name in request_list:
            # Normalize requested name
            nm = name.lstrip(">").strip()
            if nm in mapping:
                resolved.append((name, mapping[nm]))
            elif nm.lower() in lower_map:
                resolved.append((name, lower_map[nm.lower()]))
            else:
                not_found.append(name)
        return resolved, not_found
    # --- End of nested helper ---

    # Build mapping: normalized header -> index (first occurrence)
    normalized_to_index = {}
    for idx, raw_h in enumerate(selected_headers):
        h = raw_h.lstrip(">").strip()
        if h not in normalized_to_index:
            normalized_to_index[h] = idx

    # Resolve forward and reverse lists
    resolved_fwd, not_found_fwd = resolve_ordered(select_list_fwd, normalized_to_index)
    resolved_rv, not_found_rv = resolve_ordered(select_list_rv, normalized_to_index)

    # Combine all resolved entries for writing
    all_resolved = [("FWD", *x) for x in resolved_fwd] + [("RV", *x) for x in resolved_rv]

    # --- Write combined CSV with all info ---
    try:
        with open(out_csv_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            # CSV Header (Corrected to match data being written)
            writer.writerow([
                "Direction", "Requested_Name", "Index", "Original_Header",
                "Barcode", "Barcode_With_Adapter"
            ])

            # Write data rows
            for direction, requested_name, idx in all_resolved:
                writer.writerow([
                    direction,
                    requested_name,
                    idx,
                    selected_headers[idx],
                    selected_barcodes[idx],
                    selected_barcodes_adapters[idx]
                ])

        print(f"Wrote combined CSV with {len(all_resolved)} entries to: {out_csv_file}")

    except IOError as e:
        print(f"Error writing to file {out_csv_file}: {e}")
    except IndexError:
        print("Error: An index was out of range. Check consistency between headers, barcodes, and adapter lists.")

    return not_found_fwd, not_found_rv


# ==============================================================================
            # Main. 
# ==============================================================================

def main():
    """
    Main function to run the primer design pipeline from the command line.
    """
    parser = argparse.ArgumentParser(description="DPP-based Barcode Design Pipeline")
    
    parser.add_argument("--fwd_anchor", type=str, required=True, 
                        help="Forward anchor sequence (e.g., CTCGATCCCGCGAAATTAATACG)")
    parser.add_argument("--rev_anchor", type=str, required=True, 
                        help="Reverse anchor sequence (e.g., ATCCGGATATAGTTCCTCCTTTCAG)")
    # --- NEW ARGUMENT ---
    parser.add_argument("--barcode_len", 
                    type=int, 
                    default=24, 
                    help="Length of the barcode sequences to generate. Default: 24 nt's")

    parser.add_argument("--num_rb_select", type=int, required=True, 
                        help="Number of row barcodes (RB) to select (e.g., 12)")
    parser.add_argument("--num_nb_select", type=int, required=True, 
                        help="Number of column/neighbor barcodes (NB) to select (e.g., 8)")
    parser.add_argument("--outdir", type=str, required=True, 
                        help="Output directory to save results and plots.")
    parser.add_argument("--name", type=str, default="my_barcodes", 
                        help="Base name for the output CSV file (default: my_barcodes).")
    
    args = parser.parse_args()

    # --- Assign args to variables ---
    anchor_fwd = args.fwd_anchor
    anchor_rev = args.rev_anchor
    barcode_len = args.barcode_len  # Get the new argument
    num_rb_select = args.num_rb_select
    num_nb_select = args.num_nb_select
    outdir = args.outdir
    name = args.name
    
    # Define file paths
    checkpoint_file = os.path.join(outdir, f"{name}_kernels_checkpoint.npz")
    out_csv_file = os.path.join(outdir, f"{name}.csv")

    # --- Start Execution ---
    os.makedirs(outdir, exist_ok=True)
    start_time = datetime.datetime.now()
    print(f"\n=== Primer Design Run Started at {start_time} ===")
    print(f"Parameters:")
    print(f"  Fwd Anchor:  {anchor_fwd}")
    print(f"  Rev Anchor:  {anchor_rev}")
    print(f"  Barcode Len: {barcode_len}") # Added for clarity
    print(f"  RB Select:   {num_rb_select}")
    print(f"  NB Select:   {num_nb_select}")
    print(f"  Out Dir:     {outdir}")
    print(f"  File Name:   {name}.csv")
    print(f"")

    # --- [1/5] Kernel Generation or Loading ---
    kernel_loaded = False
    if os.path.exists(checkpoint_file):
        print(f"Found existing kernel file: {checkpoint_file}")
        print("Attempting to load...")
        try:
            data = np.load(checkpoint_file, allow_pickle=True)
            saved_fwd = str(data['anchor_fwd'])
            saved_rev = str(data['anchor_rev'])
            saved_barcode_len = int(data['barcode_len'])

            if (saved_fwd == anchor_fwd and 
                saved_rev == anchor_rev and 
                saved_barcode_len == barcode_len):

                print("Anchors and barcode length match. Loading data from file...")
                K_psd = data['K_psd']
                all_barcodes = data['all_barcodes']
                all_headers = data['all_headers']
                sequences_with_adapter = data['sequences_with_adapter']
                kernel_loaded = True
                print("[1/5] Generating PSD kernel... SKIPPED (loaded from file)")
            else:
                print("Parameters mismatch. Regenerating kernel.")
                if saved_fwd != anchor_fwd:
                    print(f"  - Fwd Anchor mismatch: (Input: {anchor_fwd}, Saved: {saved_fwd})")
                if saved_rev != anchor_rev:
                    print(f"  - Rev Anchor mismatch: (Input: {anchor_rev}, Saved: {saved_rev})")
                if saved_barcode_len != barcode_len:
                    print(f"  - Barcode Len mismatch: (Input: {barcode_len}, Saved: {saved_barcode_len})")
        except Exception as e:
            print(f"Error loading {checkpoint_file}: {e}. Regenerating kernel.")

    if not kernel_loaded:
        print("[1/5] Generating PSD kernel...")
        K_psd, all_barcodes, all_headers, sequences_with_adapter = get_my_dpp_ready_kernel(
            anchor_fwd, anchor_rev, barcode_len
        )
        
        print(f"Saving kernel to {checkpoint_file}...")
        np.savez(
            checkpoint_file,
            K_psd=K_psd,
            sequences_with_adapter=np.array(sequences_with_adapter, dtype=object),
            all_barcodes=np.array(all_barcodes, dtype=object),
            all_headers=np.array(all_headers, dtype=object),
            anchor_fwd=np.array(anchor_fwd, dtype=object), 
            anchor_rev=np.array(anchor_rev, dtype=object),
            barcode_len=np.array(barcode_len) # Save barcode_len for checking
        )

    # --- [2/5] DPP Sampling ---
    n_rows = num_rb_select
    n_columns = num_nb_select
    n_subset = 3*n_rows + 3*n_columns

    print("[2/5] Sampling DPP subsets...")
    res, selected_headers, selected_barcodes, selected_barcodes_adapters = sample_dpps_search(
        all_barcodes,
        all_headers,
        sequences_with_adapter,
        K_psd,
        n_subset
    )

    print(f"Constructed {len(all_barcodes)} total sequences.\n")
    all_headers_selected, sequences_with_adapter_selected = process_barcodes_with_adapters(
        selected_headers, selected_barcodes, anchor_fwd, anchor_rev
    )

    # --- [3/5] Validation ---
    print("[3/5] Validating DPP Samples...")
    _, K5_raw_fold_selected, _, _ = create_primer_kernel(
        sequences_with_adapter_selected, k=0.2
    )

    # --- [4/5] QC Plotting ---
    print("[4/5] Plotting QC...")
    select_list_fwd, select_list_rv = analyze_primer_interactions(
        K5_raw_fold_selected,
        all_headers_selected,
        outdir,  # Pass the outdir here
        plot_filename_base="selected_heatmap_ranked",
        plot_title="Interaction Heatmap: RB vs NB Primers (Ranked by Mean MFE)",
        num_rb_select=num_rb_select,  # Use arg
        num_nb_select=num_nb_select   # Use arg
    )

    # --- [5/5] Output ---
    print("[5/5] Outputting...")
    not_found_f, not_found_r = resolve_and_write_barcodes_csv(
        select_list_fwd=select_list_fwd,
        select_list_rv=select_list_rv,
        selected_headers=selected_headers,
        selected_barcodes=selected_barcodes,
        selected_barcodes_adapters=selected_barcodes_adapters,
        out_csv_file=out_csv_file 
    )

    end_time = datetime.datetime.now()
    print(f"\n=== Primer Design Run Completed at {end_time} ===")
    print(f"Total runtime: {end_time - start_time}")
    print(f"Results saved to: {outdir}")


if __name__ == "__main__":
    main()