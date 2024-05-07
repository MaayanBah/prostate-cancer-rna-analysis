from pathlib import Path
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("transcript_to_gene_table", help="Something such as TranscriptToGene.hg38.tsv", type=Path)
parser.add_argument("salmon_main_dir",
                    help=(
                        "Where there are one (or more) subdirs named $sample, and in each, a file named `quant.sf` "
                        "(Notice that we assume that each such subdir's name is also the sample's name.)"
                    ),
                    type=Path)
parser.add_argument("out_file", help="Path to out file, containing per-sample gene-expression quantification")

args = parser.parse_args()

transcript_to_gene_table = args.transcript_to_gene_table
salmon_main_dir = args.salmon_main_dir
out_file = args.out_file

transcript_to_gene_df = pd.read_table(transcript_to_gene_table, sep="\t")
transcript_to_gene_df["Transcript"] = transcript_to_gene_df["Transcript"].str.split(".", expand=True).iloc[:, 0]

salmon_subdirs = [dir for dir in salmon_main_dir.iterdir() if dir.is_dir()]

samples_names = [subdir.name for subdir in salmon_subdirs]
quant_files = [Path(subdir, "quant.sf") for subdir in salmon_subdirs if Path(subdir, "quant.sf").exists()]
print(quant_files)

quant_dfs = [
    (
        pd.read_table(quant_file)
        .rename(columns={"Name": "Transcript"})
        .drop(["Length", "EffectiveLength"], axis=1)
        .merge(transcript_to_gene_df, how="left")
        .groupby("Gene")[["NumReads", "TPM"]].sum().reset_index()
    )
    for quant_file in quant_files
]
for quant_df, sample_name in zip(quant_dfs, samples_names):
    quant_df.insert(0, "Sample", sample_name)

quant_df = pd.concat(quant_dfs).reset_index(drop=True)

quant_df.to_csv(out_file, index=False, sep="\t")