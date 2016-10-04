#!/usr/bin/env python

from argparse import ArgumentParser
from pandas import read_csv
from ZIFA import block_ZIFA
import numpy as np

def main():
    parser = ArgumentParser(description="Fit a ZIFA model on the data (with the block algorithm).")
    parser.add_argument('input_file', type=str, help="The input CSV file.")
    parser.add_argument('output_file', type=str, help="The output CSV file.")

    args = parser.parse_args()

    df = read_csv(args.input_file)
    del df['Unnamed: 0']

    lc = np.array(df)

    Z, model_params  = block_ZIFA.fitModel(lc, 2)

    np.savetxt(args.output_file, Z, delimiter=',')


if __name__ == '__main__':
        main()
