# Full Phylogenomic Pipeline with Bootstrap Analysis

Streamlit web app built with **Biopython**. Upload (or paste) a FASTA file and get an aligned,
trimmed dataset, QC statistics, a phylogenetic tree and bootstrap support values, with downloads
(Newick, PNG, CSV, FASTA).

**Live app:** `<paste your Streamlit link here>`

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
python test_pipeline.py     # quick self-test
```
MAFFT and IQ-TREE are optional locally (the app detects them). Without them, provide an
already-aligned FASTA and use Neighbor-Joining/UPGMA + bootstrap.

## Files
- `app.py` – web interface
- `pipeline.py` – analysis code (Biopython)
- `sample_data/sample_dna.fasta` – simulated demo dataset
- `requirements.txt`, `packages.txt` – Python and system dependencies for deployment
