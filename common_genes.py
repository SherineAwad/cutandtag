#! /usr/bin/env python
import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt

def read_list(infile): 
    genes= [] 
    for line in open(infile):
       record = line.strip('\n') 
       genes.append(record)
    return(genes) 

def common(genes1, genes2): 
     count = 0 
     for i in genes1: 
          if i in genes2: 
              count +=1 
     return count 
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('infile1')
    parser.add_argument('infile2')
    args = parser.parse_args()
    infile1 = args.infile1  
    infile2 = args.infile2 
    genes1 = read_list(infile1) 
    genes2 = read_list(infile2) 
    commons = common(genes1, genes2)  
    print(commons)
if __name__ == '__main__':
    main()

   

