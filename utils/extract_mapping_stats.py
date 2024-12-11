import pandas as pd

import os


def compile_flagstat_results(sample, reference, file, output):
    """
    read flagstat tsv output, add sample and referenc columns and write to output
    """
    output_df = pd.DataFrame()
    if os.path.exists(output):
        output_df = pd.read_csv(output, sep="\t")

    print(output_df)

    df = pd.read_csv(file, sep="\t", header=None, names=["value", "qual", "metric"])
    df["sample"] = sample
    df["reference"] = reference

    if output_df.empty:
        output_df = df
    else:
        output_df = pd.concat([output_df, df])

    output_df.to_csv(
        os.path.join(
            os.path.dirname(output),
            "_".join([sample, reference]) + "_flagstat.tsv",
        ),
        index=False,
        sep="\t",
    )


def get_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("sample", help="Sample name")
    parser.add_argument("reference", help="Reference name")
    parser.add_argument("file", help="Flagstat output file")
    parser.add_argument("output", help="Output file")
    return parser.parse_args()


def main():
    args = get_args()
    compile_flagstat_results(args.sample, args.reference, args.file, args.output)


if __name__ == "__main__":
    main()
