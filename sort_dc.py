#!/usr/bin/env python
import sys, argparse
import pandas as pd

def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Fit a density profile and extract valuble information')
    parser.add_argument('-f', "--filename", type=str, required=True,
                   help='file to be sorted')
    parser.add_argument("-o", "--outfilename", type=str, default='sorted.dat',
                   help='output file name, assumed to be density_profile_fit.dat')
    parser.add_argument("-c", "--collumn", type=int, default=0,
                   help='collumn to be sorted.')
    parser.add_argument('-l', "--start_line", type=int,
                   help='Set the starting line number.')

    lines = open(parser.parse_args().filename, 'r').readlines()
    output = open(parser.parse_args().outfilename, 'w')

    df = pd.read_csv(parser.parse_args().filename, sep=' ', skipinitialspace=True, header=None)
    df.sort_values(by=parser.parse_args().collumn, inplace=True)
    df.to_csv(parser.parse_args().outfilename, header=None, index=None, sep=' ')

if __name__ == '__main__':
    sys.exit(main())

