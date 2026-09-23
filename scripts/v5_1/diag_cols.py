"""Check pydeseq2 results_df column names."""
import pandas as pd
r = pd.read_csv(r"D:\stroke_apa\results\p1_DE_overall_astro.tsv", sep="\t", nrows=3)
print(list(r.columns))
