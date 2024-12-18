    
import pandas as pd
import os


def args_parse(parser):
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory of the pipeline')
    return parser

def main(output_dir):
    files = os.listdir(output_dir)
    files= [x for x in files if 'flag' in x]
    dfs= [pd.read_csv(os.path.join(output_dir, f),  
                      sep = "\t") for f in files]
    df= pd.concat(dfs)
    df.to_csv(
        os.path.join(output_dir, 'flagstat_summary.txt'),
    )


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser = args_parse(parser)
    args = parser.parse_args()
    main(args.output_dir)













