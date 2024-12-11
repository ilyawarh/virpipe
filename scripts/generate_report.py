import pandas as pd
import sys

def parse_centrifuge(file):
    out_columns=["readID", "seqID", "taxID", "score", "hitLength"]
    df=pd.read_csv(file, sep="\t")[out_columns]
    return df[~df.apply(lambda row: row.astype(str).str.contains("unclassified").any(), axis=1)]

def parse_deepvirfinder(file):
    return pd.read_csv(file, sep="\t").drop(columns=['len'])

def parse_virsorter2(file):
    df=pd.read_csv(file, sep="\t").drop(columns=['RNA','max_score_group','length','cellular'])
    df["seqname"]=df["seqname"].str.replace(r"\|\|.*$", "", regex=True)
    return df

def merge_report(centrifuge, deepvirfinder, virsorter2):
    return centrifuge.merge(deepvirfinder, left_on="readID", right_on="name", how="left") \
                       .merge(virsorter2, left_on="readID", right_on="seqname", how="left") \
                        .rename(columns={"score_x": "score", "score_y": "dvf_score", "pvalue": "dvf_pvalue", \
                                         "max_score": "vs2_score", "viral": "vs2_viral"}).drop(columns=['name','seqname'])

if __name__ == "__main__":
    centrifuge = parse_centrifuge(snakemake.input[0])
    deepvirfinder = parse_deepvirfinder(snakemake.input[1])
    virsorter2 = parse_virsorter2(snakemake.input[2])
    report = merge_report(centrifuge, deepvirfinder, virsorter2)
    pd.DataFrame(report).to_csv(snakemake.output[0], sep="\t", index=False)

