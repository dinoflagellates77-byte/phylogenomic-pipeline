"""Core phylogenetic pipeline built on Biopython.

No Streamlit imports here, so it can be tested from the command line:
    python test_pipeline.py
Steps: parse FASTA -> (optional MAFFT alignment) -> trim gappy columns -> QC stats
       -> distance matrix -> NJ/UPGMA tree -> bootstrap support
       (or IQ-TREE ML + ultrafast bootstrap if IQ-TREE is installed).
"""
import io
import os
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Consensus import get_support
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DNA_MODELS = ["identity", "blastn"]
PROTEIN_MODELS = ["identity", "blosum62", "pam250"]
DNA_LETTERS = set("ACGT-")
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWY-")


class PipelineError(Exception):
    """Raised for user-fixable problems (bad input, missing tool)."""


# ---------------------------------------------------------------- input
def parse_fasta(text):
    records = []
    for rec in SeqIO.parse(io.StringIO(text), "fasta"):
        seq = str(rec.seq).upper().replace(".", "-").replace("*", "").replace(" ", "")
        records.append(SeqRecord(Seq(seq), id=rec.id, name=rec.id, description=""))
    if len(records) < 4:
        raise PipelineError(f"Need at least 4 sequences for a tree with support values (found {len(records)}).")
    if len(records) > 200:
        raise PipelineError("Please use 200 sequences or fewer for the web version.")
    ids = [r.id for r in records]
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        raise PipelineError(f"Duplicate sequence IDs: {', '.join(dupes[:5])}")
    return records


def detect_seq_type(records):
    chars = "".join(str(r.seq) for r in records).replace("-", "")
    if not chars:
        raise PipelineError("Sequences are empty.")
    dna_frac = sum(chars.count(c) for c in "ACGTUN") / len(chars)
    return "dna" if dna_frac > 0.9 else "protein"


def sanitize(records, seq_type):
    """Replace ambiguity codes (N, X, B, Z, ...) with gaps so distance models work."""
    allowed = DNA_LETTERS if seq_type == "dna" else PROTEIN_LETTERS
    out = []
    for r in records:
        s = str(r.seq)
        if seq_type == "dna":
            s = s.replace("U", "T")
        s = "".join(c if c in allowed else "-" for c in s)
        out.append(SeqRecord(Seq(s), id=r.id, name=r.id, description=""))
    return out


# ------------------------------------------------------------ alignment
def has_mafft():
    return shutil.which("mafft") is not None


def align_with_mafft(records):
    exe = shutil.which("mafft")
    if not exe:
        raise PipelineError("MAFFT is not installed on this server.")
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "in.fasta")
        SeqIO.write(records, inp, "fasta")
        res = subprocess.run([exe, "--auto", "--quiet", inp], capture_output=True, text=True, timeout=600)
    if res.returncode != 0:
        raise PipelineError(f"MAFFT failed: {res.stderr[-300:]}")
    aligned = list(SeqIO.parse(io.StringIO(res.stdout), "fasta"))
    by_id = {r.id: r for r in aligned}
    return [SeqRecord(Seq(str(by_id[r.id].seq).upper()), id=r.id, name=r.id, description="") for r in records]


def _to_array(aln):
    return np.array([list(str(r.seq)) for r in aln])


def _from_array(arr, ids):
    return MultipleSeqAlignment(
        [SeqRecord(Seq("".join(row)), id=i, name=i, description="") for i, row in zip(ids, arr)]
    )


def filter_columns(aln, max_gap_fraction):
    """Drop columns whose gap fraction exceeds the threshold (all-gap columns are always dropped)."""
    arr = _to_array(aln)
    frac = (arr == "-").mean(axis=0)
    keep = (frac <= max_gap_fraction) & (frac < 1.0)
    return _from_array(arr[:, keep], [r.id for r in aln])


def alignment_stats(aln):
    arr = _to_array(aln)
    n_sites = arr.shape[1]
    constant = variable = informative = 0
    for col in arr.T:
        states = col[col != "-"]
        if states.size == 0:
            continue
        vals, counts = np.unique(states, return_counts=True)
        if len(vals) == 1:
            constant += 1
        else:
            variable += 1
            if (counts >= 2).sum() >= 2:
                informative += 1
    return {
        "Sequences": arr.shape[0],
        "Alignment length (sites)": n_sites,
        "Constant sites": constant,
        "Variable sites": variable,
        "Parsimony-informative sites": informative,
        "Gaps (%)": round(100 * float((arr == "-").mean()), 2),
    }


# ----------------------------------------------------------------- trees
def distance_matrix_df(aln, model):
    dm = DistanceCalculator(model).get_distance(aln)
    names = list(dm.names)
    data = [[dm[a, b] for b in names] for a in names]
    return pd.DataFrame(data, index=names, columns=names)


def bootstrap_tree(aln, model, method, n_boot, seed, progress=None):
    """Distance tree from the original alignment + column-resampling bootstrap support (%)."""
    ctor = DistanceTreeConstructor(DistanceCalculator(model), method)
    best = ctor.build_tree(aln)
    arr = _to_array(aln)
    ids = [r.id for r in aln]
    n_cols = arr.shape[1]
    rng = np.random.default_rng(seed)
    replicates = []
    for i in range(n_boot):
        idx = rng.integers(0, n_cols, n_cols)
        replicates.append(ctor.build_tree(_from_array(arr[:, idx], ids)))
        if progress:
            progress((i + 1) / n_boot)
    tree = get_support(best, replicates)
    _normalise_support(tree, n_boot)
    return _clean_tree(tree)


def _normalise_support(tree, n_boot):
    vals = [c.confidence for c in tree.get_nonterminals() if c.confidence is not None]
    if not vals:
        return
    top = max(vals)
    scale = 1.0
    if top > 100:          # raw counts
        scale = 100.0 / n_boot
    elif top <= 1.0:       # fractions
        scale = 100.0
    if scale != 1.0:
        for c in tree.get_nonterminals():
            if c.confidence is not None:
                c.confidence = c.confidence * scale


def _clean_tree(tree):
    for c in tree.find_clades():
        if c.branch_length is not None and c.branch_length < 0:
            c.branch_length = 0.0
    for c in tree.get_nonterminals():
        c.name = None
    return tree


def iqtree_available():
    return shutil.which("iqtree2") or shutil.which("iqtree")


def run_iqtree(aln, n_boot=1000, threads=2):
    """ML tree with automatic model selection + ultrafast bootstrap (needs IQ-TREE on PATH)."""
    exe = iqtree_available()
    if not exe:
        raise PipelineError("IQ-TREE is not installed on this server.")
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout
    v1 = "version 1." in ver
    n_boot = max(1000, n_boot)  # UFBoot requires >= 1000
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, "aln.fasta")
        with open(fasta, "w") as fh:
            for r in aln:
                fh.write(f">{r.id}\n{r.seq}\n")
        cmd = [exe, "-s", fasta, "-m", "TEST" if v1 else "MFP",
               "-bb" if v1 else "-B", str(n_boot), "-bnni",
               "-nt" if v1 else "-T", str(threads), "-pre", os.path.join(d, "run"), "-redo", "-quiet"]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=1500)
        tree_path = os.path.join(d, "run.treefile")
        if res.returncode != 0 or not os.path.exists(tree_path):
            raise PipelineError(f"IQ-TREE failed: {(res.stderr or res.stdout)[-400:]}")
        tree = Phylo.read(tree_path, "newick")
        model = ""
        rep = os.path.join(d, "run.iqtree")
        if os.path.exists(rep):
            for line in open(rep):
                if line.startswith("Best-fit model"):
                    model = line.split(":", 1)[1].strip()
                    break
    return _clean_tree(tree), model


def support_table(tree):
    rows = []
    for clade in tree.get_nonterminals():
        if clade is tree.root:
            continue
        tips = sorted(t.name for t in clade.get_terminals())
        s = clade.confidence
        cat = "n/a" if s is None else ("Strong (>=95)" if s >= 95 else "Moderate (70-95)" if s >= 70 else "Weak (<70)")
        rows.append({"Clade members": ", ".join(tips), "Tips": len(tips),
                     "Support (%)": None if s is None else round(float(s), 1), "Category": cat})
    return pd.DataFrame(rows).sort_values("Tips").reset_index(drop=True)


def tree_to_newick(tree):
    handle = io.StringIO()
    Phylo.write(tree, handle, "newick")
    return handle.getvalue()


def draw_tree(tree, title=""):
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(9, max(4.0, 0.4 * n)))
    Phylo.draw(
        tree, axes=ax, do_show=False, show_confidence=False,
        label_func=lambda c: c.name if c.is_terminal() else None,
        branch_labels=lambda c: (f"{c.confidence:.0f}" if (c.confidence is not None and not c.is_terminal()) else None),
    )
    ax.set_title(title)
    ax.set_xlabel("Branch length (substitutions per site)")
    ax.set_ylabel("")
    ax.set_yticks([])
    fig.tight_layout()
    return fig


def fig_to_png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
    return buf.getvalue()


# ----------------------------------------------------------- full run
def run_analysis(text, realign=False, max_gap=0.5, method="nj", model="identity",
                 n_boot=100, seed=42, progress=None):
    records = parse_fasta(text)
    seq_type = detect_seq_type(records)
    records = sanitize(records, seq_type)
    notes = []
    unequal = len({len(r.seq) for r in records}) > 1
    if unequal or realign:
        if not has_mafft():
            raise PipelineError("Sequences are not aligned (different lengths) and MAFFT is not available. "
                                "Upload an already-aligned FASTA.")
        records = align_with_mafft(records)
        notes.append("Aligned with MAFFT (--auto).")
    aln = MultipleSeqAlignment(records)
    before = aln.get_alignment_length()
    aln = filter_columns(aln, max_gap)
    after = aln.get_alignment_length()
    notes.append(f"Column trimming (max gap fraction {max_gap}): {before} -> {after} sites.")
    if after < 20:
        raise PipelineError("Fewer than 20 usable alignment columns remain; relax the gap threshold.")

    stats = alignment_stats(aln)
    dm = distance_matrix_df(aln, model)
    if method == "iqtree":
        tree, best_model = run_iqtree(aln, n_boot)
        label = f"Maximum likelihood (IQ-TREE, {best_model or 'auto model'}), UFBoot support"
    else:
        tree = bootstrap_tree(aln, model, method, n_boot, seed, progress)
        label = f"{'Neighbor-Joining' if method == 'nj' else 'UPGMA'} ({model} distance), {n_boot} bootstrap replicates"
    fasta_out = "".join(f">{r.id}\n{r.seq}\n" for r in aln)
    return {"seq_type": seq_type, "notes": notes, "stats": stats, "distances": dm, "tree": tree,
            "label": label, "alignment_fasta": fasta_out, "support": support_table(tree)}
