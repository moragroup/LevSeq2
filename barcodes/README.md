# Barcode design and replacement

## Install 

To run the barcode creator you'll need to install the following packages
```
pip install simsimd dppy multiprocess
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install viennarna
conda install -c plotly plotly-orca
```

## Overview
This code and readme was written by Niklas Madsen.

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

## Example

For example, to make this for the pet22b vector used in the LevSeq paper, we do the following:

python barcode_design_cli.py \
    --fwd_anchor "CTCGATCCCGCGAAATTAATACG" \
    --rev_anchor "ATCCGGATATAGTTCCTCCTTTCAG" \
    --barcode_len 24 \
    --num_rb_select 12 \
    --num_nb_select 8 \
    --outdir "barcode_design_pet22b_T7" \
    --name "barcode_T7_pet22b_pool"
