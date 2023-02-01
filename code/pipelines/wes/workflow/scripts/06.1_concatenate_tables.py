# -*- coding: utf-8 -*-
"""
@created: Jan 05 2022
@modified: Nov 05 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Concatenate multiple tsv tables with headers with the option of keeping the header.
"""

import argparse
import numpy as np
import os
import pandas as pd
import subprocess

# functions ============================================================================================================

def read_header(path):
    with open(path, "r") as file:
        header = [x for x in file.readlines() if x.startswith("##")]
    return header

def read_table(path):
    header = read_header(path)
    df = pd.read_table(path, skiprows=len(header), na_values=["-","."])
    return df


def save_table_with_header(df, header, output):
    output_header = os.path.join(os.path.dirname(output), "header.tsv")

    print("-writing header in file %s" % output_header)
    with open(output_header, "w") as file_header:
        for line in header:
            file_header.write(line)

    # write contents
    if output.endswith(".gz"):
        output_uncompressed = output[:-3] + ".tmp"
        output_concatenate = output[:-3]
    else:
        output_uncompressed = output + ".tmp"
        output_concatenate = output

    print("-writing contents in file %s" % output_uncompressed)
    df.to_csv(output_uncompressed, index=False, sep="\t")

    # # concat both files
    cmd = "cat %s %s >> %s" % (output_header, output_uncompressed, output_concatenate)
    print("-running the command:\n\t%s" % cmd)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    os.remove(output_header)
    os.remove(output_uncompressed)

    # compress if required
    if output_concatenate != output:
        cmd = "gzip %s" % output_concatenate
        print("-running the command:\n\t%s" % cmd)
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def main(args):
    if args.files is not None:
        files = args.files
    elif args.folder is not None:
        files = [os.path.join(args.folder, x) for x in os.listdir(args.folder)]
    else:
        raise ValueError("One of --files or --folder must be specified.")

    # load all tables 
    dfs = []
    for i, file in enumerate(files):
        df = read_table(file)
        dfs.append(df)
        if len(files)>100:
            if (i+1)%(len(files)//100)==0:
                print("-processed %d/%d files" % (i+1, len(files)), flush=True)

    # concatenate
    print("-concatenating %d tables..." % len(dfs), flush=True, end="")
    df = pd.concat(dfs, axis=0)
    print("done!", flush=True)

    # save
    print("-saving concatenated table at %s..." % args.output, flush=True, end="")
    if not args.keep_header:
        df.to_csv(args.output, index=False, sep="\t")
    else:
        header = read_header(files[0])
        save_table_with_header(df, header, args.output)
    print("done!", flush=True)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple tables.")
    parser.add_argument('--files', nargs="*", type=str, help='Paths to tables.', default=None)
    parser.add_argument('--folder', type=str, help='Path to folder with tables.', default=None)
    parser.add_argument("--keep_header", action="store_true", default=False,
                        help="If used, the header of the maf tables is preserved.")
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
