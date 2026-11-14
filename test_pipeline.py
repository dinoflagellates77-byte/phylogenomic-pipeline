"""Smoke test: python test_pipeline.py"""
from pathlib import Path
import pipeline as pl

text = (Path(__file__).parent / "sample_data" / "sample_dna.fasta").read_text()
res = pl.run_analysis(text, method="nj", model="identity", n_boot=20, seed=1)
tips = res["tree"].get_terminals()
assert len(tips) == 8, len(tips)
assert res["support"]["Support (%)"].notna().all()
print("OK:", res["label"])
print(res["support"].to_string())
