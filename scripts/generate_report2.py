import pandas as pd

def parse_kaiju(file):
    df=pd.read_csv(file, sep="\t", usecols=[1,3], names=['C/U', 'readID', 'taxID', 'base_taxonomy'])
    return df

def parse_genomad(file):
    return pd.read_csv(file, sep="\t", usecols=['seq_name','virus_score','n_hallmarks','taxonomy'])

def parse_virsorter2(file):
    df=pd.read_csv(file, sep="\t", usecols=['seqname','max_score','hallmark','viral'])
    df["seqname"]= df["seqname"].str.replace(r"\|\|.*$", "", regex=True)
    return df

def merge_report(kaiju, virsorter2, genomad):
    return kaiju.merge(virsorter2, left_on="readID", right_on="seqname", how="left") \
                       .merge(genomad, left_on="readID", right_on="seq_name", how="left") \
                        .rename(columns={"max_score": "vs2_score", "hallmark": "vs2_hallmark", "n_hallmarks": "genomad_hallmarks", \
                                         "viral": "vs2_viral", "virus_score": "genomad_viral"}).drop(columns=['seq_name','seqname'])

if __name__ == "__main__":
    kaiju = parse_kaiju(snakemake.input[0])
    virsorter2 = parse_virsorter2(snakemake.input[1])
    genomad = parse_genomad(snakemake.input[2])
    report = merge_report(kaiju, virsorter2, genomad)
    pd.DataFrame(report).to_csv(snakemake.output[0], sep="\t", index=False)
