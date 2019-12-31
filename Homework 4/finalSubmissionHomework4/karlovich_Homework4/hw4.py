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
    parser = argparse.ArgumentParser(prog='hw4.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%prog 2.0')
    parser.add_argument("-i","--input",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to sequence fna file [required]") 
    return parser
    

if __name__ == '__main__':
    # -- Problem 1 Starts Here --
    parser = make_arg_parser()
    args = parser.parse_args()
    input_file = open(args.input) #open hw3.fna
    input_lines = input_file.readlines()
    arr_of_sequences = []
    dict_list = [dict() for x in range(1475)]
    for x in range(0,len(input_lines),2):
        arr_of_sequences.append(input_lines[x+1])
    output = open("hw4-dictionary.txt", 'w+') #problem 1 outputfile
    for line in arr_of_sequences:
        for i in range(0, len(line)):
            l_char = line[i]
            if l_char in dict_list[i]:
                dict_list[i][l_char] += 1
            else:
                dict_list[i][l_char] = 1
    idx = 1
    for entry in dict_list:
        count = 0
        if "\n" not in entry: #ignore last character which is new line
            for key in entry:
                count += entry[key]
                if(entry[key] < 1000):
                    output.write(str(key) + ": " + str(entry[key]) + " \t\t")
                else:
                    output.write(str(key) + ": " + str(entry[key]) + " \t")
            output.write("\n")
            if(count != 5088):
                print("doesn't sum!" + str(count) + " " + str(idx))
            idx += 1
    output.close()
    output = open("hw4-variability.txt", "w+")
    for entry in dict_list:
        entry.pop("-", None)
    dict_list.pop()
    #print(dict_list)
    #old method, where conserved == variability
    for dic in dict_list:
        max_corr = ("_", 0)
        for base in dic:
            #print("dic base" + str(dic[base]))
            #print("max corr[1]" + str(max_corr[1]))
            if (dic[base] > max_corr[1]):
                max_corr = (base, dic[base])
        output.write(str(max_corr[0]) + "\t" + str(round((max_corr[1] / 5088),3)) + "\n")
    



