#! /usr/bin/env python
import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt

def read_peaks(peaks_file, index): 
    count = 0
    genes= [] 
    for line in open(peaks_file):
        count +=1 
        if count == 1:
            continue
        records = line.split('\t')
        if (records[index] not in genes):
                genes.append(records[index]) 
    return genes, count

def printValues(genes, outfile):
     fout = open(outfile, 'w+')
     for i in genes:  
         print(i, file = fout) 
     fout.close() 
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('peaks')
    parser.add_argument('outfile')
    parser.add_argument('index') 
    args = parser.parse_args()
    peaks_file = args.peaks  
    outfile = args.outfile 
    index = int(args.index)
    genes, count = read_peaks(peaks_file, index) 
    printValues(genes, outfile)
    print (count -1 )
if __name__ == '__main__':
    main()

   

