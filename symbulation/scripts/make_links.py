# Converts a Symbulation Interaction Snapshot file to a links.csv file for input
# into SuchTree

import pandas as pd
import sys


def make_links(df, outfile="links.csv"):
    df.columns = df.columns.str.replace(' ', '')
    hosts = df["host"].unique()

    syms = df["symbiont"].unique()
    data = [[0 for i in range(len(hosts))] for j in range(len(syms))]

    max_inter = df["count"].max()

    links_df = pd.DataFrame(data, columns=hosts)
    links_df.index = syms

    for row in df.itertuples():
        links_df.loc[row.symbiont, row.host] = row.count

    links_df.to_csv(outfile)
    return links_df


if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Usage: python make_links.py [Interaction Snapshot Filename]")
        exit(0)

    df = pd.read_csv(sys.argv[1])

    if len(sys.argv) > 1:
        make_links(df, sys.argv[2])
    else:
        make_links(df)
