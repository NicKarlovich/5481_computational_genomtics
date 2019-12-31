#  Name: Nick Karlovich
#  x500: karlo015
# Class: CSCI 5481
#  Year: Fall 2019
import sys, os
import argparse
import numpy as np
from subprocess import Popen, PIPE

# arg parser, originally adapated from HW1
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='homework3.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%prog 2.0')
    parser.add_argument("-i","--input",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to sequence fasta file [required]") 
    return parser

if __name__ == '__main__':
    # -- Problem 1 Starts Here --
    parser = make_arg_parser()
    args = parser.parse_args()
    input_file = open(args.input) #open hw3.fna
    input_lines = input_file.readlines()
    query_tuples = []
    
    output_file_1 = open("hw4-var1.fna", "w+")
    begin_idx_1 = 111 # inclusive
    end_idx_1 = 219 #exclusive end of range: 1 + value in hw4-variable-regions-indicies.txt
    output_file_2 = open("hw4-var2.fna", "w+")
    begin_idx_2 = 416 # inclusive
    end_idx_2 = 463 #exclusive end of range: 1 + value in hw4-variable-regions-indicies.txt
    for x in range(0,len(input_lines),2):
        output_file_1.write(input_lines[x][:-1])
        output_file_2.write(input_lines[x][:-1])
        output_file_1.write("\n")
        output_file_2.write("\n")
        output_file_1.write(input_lines[x+1][begin_idx_1:end_idx_1])
        output_file_2.write(input_lines[x+1][begin_idx_2:end_idx_2])
        output_file_1.write("\n")
        output_file_2.write("\n")
    input_file.close()
    output_file_1.close()
    output_file_2.close()
