#!/usr/bin/env python

import sys
sys.path.insert(1, '/usr/local/lib/python2.7/site-packages/') ##change the path if needed

from argparse import ArgumentParser
from pandas import read_csv
from ZIFA import ZIFA,block_ZIFA
import numpy as np

def main():
    parser = ArgumentParser(description="Fit a ZIFA model on the data.")
    parser.add_argument('-b', '--block', action='store_true', default=False, help="Whether the block algorithm should be used.")
    parser.add_argument('input_file', type=str, help="The input CSV file.")
    parser.add_argument('output_file', type=str, help="The output CSV file.")

    args = parser.parse_args()

    df = read_csv(args.input_file)
    del df['Unnamed: 0']

    lc = np.array(df)
    Y = np.transpose(lc)

    if(args.block):
        Z, model_params  = block_ZIFA.fitModel(Y, 2)
    else:
        Z, model_params  = ZIFA.fitModel(Y, 2)

    np.savetxt(args.output_file, Z, delimiter=',')


if __name__ == '__main__':
        main()
