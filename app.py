from pathlib import Path

import pandas as pd
import streamlit as st

import pipeline as pl

st.set_page_config(page_title="Phylogenomic Pipeline with Bootstrap", page_icon="🧬", layout="wide")
st.title("🧬 Full Phylogenomic Pipeline with Bootstrap Analysis")
st.caption("Biopython-based: FASTA → alignment → QC → tree building → bootstrap support → downloads")

_found = list(Path(__file__).parent.rglob("sample_dna.fasta"))
SAMPLE = _found[0] if _found else Path("missing")

# ------------------------------------------------------------ sidebar
with st.sidebar:
    st.header("1. Input sequences")
    source = st.radio("Source", ["Sample dataset", "Upload FASTA", "Paste FASTA"])
    text = ""
    if source == "Sample dataset":
        text = SAMPLE.read_text() if SAMPLE.exists() else ""
        st.caption("Simulated 8-taxon, 500 bp DNA alignment.")
    elif source == "Upload FASTA":
        up = st.file_uploader("FASTA file", type=["fasta", "fa", "fas", "fna", "faa", "txt"])
        if up is not None:
            text = up.getvalue().decode("utf-8", errors="ignore")
    else:
        text = st.text_area("Paste FASTA", height=200)

    seq_type = None
    if text.strip():
        try:
            seq_type = pl.detect_seq_type(pl.parse_fasta(text))
        except pl.PipelineError:
            seq_type = None
        except Exception:
            seq_type = None

    st.header("2. Settings")
    realign = st.checkbox("Re-align with MAFFT", value=False, disabled=not pl.has_mafft(),
                          help="Unaligned input is aligned automatically when MAFFT is available.")
    max_gap = st.slider("Trim columns with gap fraction above", 0.0, 1.0, 0.5, 0.05)
    methods = {"Neighbor-Joining + bootstrap": "nj", "UPGMA + bootstrap": "upgma"}
    if pl.iqtree_available():
        methods["Maximum likelihood (IQ-TREE) + UFBoot"] = "iqtree"
    method_label = st.selectbox("Tree method", list(methods))
    method = methods[method_label]
    models = pl.PROTEIN_MODELS if seq_type == "protein" else pl.DNA_MODELS
    model = st.selectbox("Distance model", models, disabled=(method == "iqtree"))
    n_boot = st.slider("Bootstrap replicates", 10, 1000, 100, 10)
    seed = st.number_input("Random seed", value=42, step=1)
    run = st.button("▶ Run pipeline", type="primary", use_container_width=True)

# --------------------------------------------------------------- run
if run:
    if not text.strip():
        st.error("Please provide sequences first.")
    else:
        bar = st.progress(0.0, text="Running bootstrap replicates...")
        try:
            with st.spinner("Running pipeline..."):
                res = pl.run_analysis(text, realign, max_gap, method, model, n_boot, int(seed),
                                      progress=lambda f: bar.progress(f, text=f"Bootstrap {f:.0%}"))
            st.session_state["res"] = res
        except pl.PipelineError as e:
            st.error(str(e))
        except Exception as e:  # keep the app alive for the evaluator
            st.error(f"Unexpected error: {e}")
        finally:
            bar.empty()

res = st.session_state.get("res")
if not res:
    st.info("Choose the sample dataset (or upload your own FASTA) in the sidebar, then click **Run pipeline**.")
    st.stop()

tab1, tab2, tab3, tab4, tab5 = st.tabs(["Tree", "Support table", "Alignment QC", "Distances", "Methods"])

with tab1:
    fig = pl.draw_tree(res["tree"], res["label"])
    st.pyplot(fig)
    c1, c2 = st.columns(2)
    c1.download_button("⬇ Tree (Newick)", pl.tree_to_newick(res["tree"]), "tree.nwk")
    c2.download_button("⬇ Tree figure (PNG)", pl.fig_to_png(fig), "tree.png", "image/png")
    st.caption("Numbers on internal branches are bootstrap support (%). ≥95 strong, 70–95 moderate, <70 weak.")

with tab2:
    df = res["support"]
    st.dataframe(df, use_container_width=True)
    st.download_button("⬇ Support table (CSV)", df.to_csv(index=False), "support.csv")
    if df["Support (%)"].notna().any():
        st.metric("Mean support", f"{df['Support (%)'].mean():.1f}%")

with tab3:
    for n in res["notes"]:
        st.write("•", n)
    st.write(f"Detected sequence type: **{res['seq_type']}**")
    st.table(pd.DataFrame(res["stats"], index=["value"]).T)
    st.download_button("⬇ Processed alignment (FASTA)", res["alignment_fasta"], "alignment.fasta")

with tab4:
    st.dataframe(res["distances"].style.format("{:.4f}"), use_container_width=True)
    st.download_button("⬇ Distance matrix (CSV)", res["distances"].to_csv(), "distances.csv")

with tab5:
    st.markdown("""
**Pipeline**
1. Parse FASTA (Biopython `SeqIO`), check IDs, detect DNA vs protein; ambiguity codes become gaps.
2. Align with MAFFT if sequences differ in length (or if requested).
3. Trim gap-rich columns; report constant / variable / parsimony-informative sites.
4. Pairwise distances (`DistanceCalculator`) and tree (`DistanceTreeConstructor`: NJ or UPGMA).
5. **Bootstrap**: resample alignment columns with replacement, rebuild the tree each time, and
   score each clade of the original tree by the % of replicates that contain it
   (`Bio.Phylo.Consensus.get_support`).
6. Optional: IQ-TREE maximum likelihood with model selection and ultrafast bootstrap (if installed).

**Notes**: trees are shown unrooted/arbitrarily rooted; support on branches next to the root
should be interpreted with care.
""")
