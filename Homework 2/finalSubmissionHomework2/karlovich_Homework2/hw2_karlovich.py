#Nick Karlovich
#karlo015

import sys, os
import argparse
import numpy as np
import random as rand
from subprocess import Popen, PIPE
import multiprocessing as mp

global NUM_OF_ITERATIONS
global OUTPUT_FILE
NUM_OF_ITERATIONS = 10000
OUTPUT_FILE = "output.txt"

#Arg parsing taken from Homework 1 and edited to correct variables
def make_arg_parser():
    parser = argparse.ArgumentParser(prog='template-python3.py',
                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%prog 2.0')
    parser.add_argument("-q","--query",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to query fasta [required]")
    parser.add_argument("-r","--ref",
                      default=argparse.SUPPRESS,
                      required=True,
                      help="Path to reference fasta [required]")
    parser.add_argument("-m","--matches",
                      default=None,
                      help="whether or not to run standard Needleman-Wunsch")
    return parser


'''
In order to run multi-processing you will need to edit needleman_wunsch
to only take 1 parameter, this is done by setting the default value for seq2, to be equal to
the sequence that will not be permuted 10,000 times.

So for example the function definition of needleman_wunsch should look something like this if you're
trying to run this function:
def needleman_wunsch(seq1, seq2 = "A REFERENCE STRING THAT WILL NOT CHANGE FOR ALL 10,000 ITERATIONS", mismatch = -3, match = 1, gap = -2):

This way, we can pass needleman_wunsch along with only one paramter in the pool.map function
'''
def multi_processing(seq1):
    num_cores = mp.cpu_count()# get cores
    pool = mp.Pool(num_cores)# create pool of threads
    list_of_seq1_perm = []
    for i in range(NUM_OF_ITERATIONS):# create 10,000 permutations to be fed into map() function
        list_of_seq1_perm.append(''.join(np.random.permutation(list(seq1))))
    output = pool.map(needleman_wunsch,list_of_seq1_perm)# Do Multi-threaded map and return an array of all scores computed
    opened_file = open("output.txt", "w")# write to output.
    opened_file.write(output)

'''
A helper function to get the known matches in a sequence.
Pulls in the lines from the -m parameter and turns them into
an array of ints
'''
def get_matches(in_file):
    reader = open(in_file)
    matches_arr = []
    for line in reader:
        temp = (line.rstrip().split("\t"))# get values by delimieter "\t"
        for x in range(len(temp)):
            temp[x] = int(temp[x])
        matches_arr.append(temp)
        #an array of arrays, each item in the array is an array of 4 indexes
        #where each two numbers belong to an inclusive range of known matches
        #in Human and Fly Respectively
    #print(matches_arr)
    return matches_arr

'''
Runs needleman_wunsch on the parts of the sequence that aren't known matches
For example: if we knew matches existed from 100 to 150 in both seq1 and seq2 and
seq1 & seq2 were both 250 bases long, instead of matching all 250, it instead does
the sequences[0 to 100],then automatically puts matches in for [100 to 150] then runs
needleman_wunsch on [150 - 250].  This greatly decreases the total memory and runtime usage.
'''
def get_anchored(seq1,seq2,matches_indices):
    start_1 = 0 #first needleman wunsch starts at index 0
    start_2 = 0
    path_back =""
    score = 0
    for anchor in matches_indices: #for every set of matches
        human_end, human_start, fly_end, fly_start = anchor
        #decrease all values by one because matches is 1 indexed instead of 0th indexed
        human_end -= 1 #
        human_start -= 1
        fly_end -= 1
        fly_start -= 1

        #set up sub-sequences for needleman wunsch
        subseq_1 = seq1[start_1:human_end]
        subseq_2 = seq2[start_2:fly_end]
        (p,s) = needleman_wunsch(subseq_1,subseq_2)
        temp_path, temp_score = get_path_and_score(p,s)

        #add subsequence value to overall values
        path_back = path_back + temp_path
        score = score + temp_score

        #add values for "matched" area.
        path_back = path_back + "d"*(human_start - human_end + 1) # add one because indicies are inclusive
        score = score + (human_start - human_end) + 1 # add one because indicies are inclusive
        start_1 = human_start + 1 # add one because indicies are inclusive
        start_2 = fly_start + 1 # add one because indicies are inclusive
    #once we run out of matches run needleman_wunsch on remainder of the sequences
    subseq_1 = seq1[start_1:]
    subseq_2 = seq2[start_2:]
    (p,s) = needleman_wunsch(subseq_1,subseq_2)
    temp_path, temp_score = get_path_and_score(p,s)
    path_back = path_back + temp_path
    score = score + temp_score
    write_path_and_score(path_back, score, seq1, seq2, OUTPUT_FILE)

'''
Runs needleman_wunsch dynamic programming alignment
algorithim on two sequences.  Takes in two sequences, query then reference
Also allows for mismatch, match and gap values, defaulted to
-3, 1, and -2 respectively
'''
def needleman_wunsch(seq1, seq2, mismatch = -3, match = 1, gap = -2):
    width = len(seq1) + 1
    height = len(seq2) + 1
    #fill score and path matrix with zeroes
    score_mtx = np.zeros((height,width), dtype = int)
    path_mtx = np.zeros((height,width), dtype = str)

    #Initalize the first row and column with gap penalty values
    for j in range(width):
        score_mtx[0][j] = gap * j
        path_mtx[0][j] = "h" #previous position was horizontal
    for i in range(height):
        score_mtx[i][0] = gap * i
        path_mtx[i][0] = "v" #previous position was vertical

    #Go through entire seq1 x seq2 array
    for j in range(1,width): #go across horizontally
        for i in range(1, height): #go across vertically
        #sequences start at idx=0 so we subtract 1 to get correct value.
            if(seq1[j - 1] == seq2[i - 1]): #we have a match, and thus use match score if we go diagonally.
                diag_score = match + score_mtx[i - 1][j - 1]
            else: #we don't have a match and thus use mismatch score if we decide we want to go diagonally.
                diag_score = mismatch + score_mtx[i -1][j - 1]
            #calculate scores if we were to pull from horizontal and vertical positions
            horz_score = gap + score_mtx[i][j - 1]
            vert_score = gap + score_mtx[i - 1][j]

            #the order of the if statements means that this particaular algorithim will decide to go horizontal rather than vertically
            #if horizontal and vertical are the same values.  Doesn't really make a difference but a point of note.
            if(diag_score >= horz_score and diag_score >= vert_score): #if diagonal score is the greatest
                #we do greater OR EQUALS because if they're equal it doesn't matter
                #and we'd rather have matches than mismatches and gaps for
                #arbitrary reasons that I decided on now.
                score_mtx[i][j] = diag_score
                path_mtx[i][j] = "d"
            elif (horz_score >= diag_score and horz_score >= vert_score): #if horizontal score is the greatest
                score_mtx[i][j] = horz_score
                path_mtx[i][j] = "h"
            else: #only other case could be vertical score is the greatest.
                score_mtx[i][j] = vert_score
                path_mtx[i][j] = "v"
    return (path_mtx,score_mtx)



'''
Retuns the path back in terms of movements ("d" diagonal, "h" horizontal, "v" vertical)
Also returns the int score value of the path
'''
def get_path_and_score(path_matrix, score_matrix):
    path_back = ""
    width = len(path_matrix[0]) - 1
    height = len(path_matrix) - 1
    score = score_matrix[height][width] #get value in bottom rigth corner, which is the score
    while(height != 0 or width != 0): #until we reach the top left corner
        curr_move = path_matrix[height][width]
        path_back = curr_move + path_back
        #if it's diagonal we decrement x and y indicies by 1
        if(curr_move == "d"):
            height-=1
            width-=1
        elif(curr_move == "h"): #if its horizontal decrement x by 1
            width-=1
        else: # if its vertical decrement y by 1
            height-=1
    return (path_back,score) #return tuple of the string of the order back and the integer score of the alignment.


'''
takes both input sequences and the
path_back which consists of "directions" ("d","h","v")
and converts it into the actual bases and gaps of the original sequences.
'''
def convert_directions_to_sequences(seq1, seq2, path_back):
    seq_1_output = ""
    seq_2_output = ""
    w_counter = 0
    h_counter = 0
    for s in path_back:
        if(s == "d"):
            seq_1_output += seq1[w_counter] #in the diag case, seq1 & seq2  and w_counter & h_counter are identical
            seq_2_output += seq2[h_counter] #in the diag case, seq1 & seq2  and w_counter & h_counter are identical
            w_counter+=1
            h_counter+=1
        elif(s == "h"): #if we go horizontal, include space, assuming we're aliging seq2 to seq1
            seq_1_output += seq1[w_counter]
            seq_2_output += "_"
            w_counter+=1
        else:
            seq_1_output += "_"
            seq_2_output += seq2[h_counter]
            h_counter+=1
    return (seq_1_output, seq_2_output)


'''
Seq1 is the query sequence running along the top
Seq2 is the reference sequence running along the side
Used for debugging in cmd line.
'''
def print_path_and_score(path_back, score, seq1, seq2):
    print("Score: " + str(score))
    print("Path Back: " + path_back)
    (seq_1_output, seq_2_output) = convert_directions_to_sequences(seq1, seq2, path_back)
    print("-------------------------")
    print("Sequence 1: " + str(seq1))
    print("Sequence 2: " + str(seq2))
    print("--- Aligned Sequences ---")
    print("Seq1: " + str(seq_1_output))
    print("Seq2: " + str(seq_2_output))

def write_path_and_score(path_back, score, seq1, seq2, out_file = "output.txt"):
    output = open(out_file, 'w')
    output.write("Score: " + str(score) + "\n")
    output.write("Path Back: " + path_back + "\n")
    (seq_1_output, seq_2_output) = convert_directions_to_sequences(seq1, seq2, path_back)
    output.write("\n-------------------------")
    output.write("\nSequence 1: " + str(seq1))
    output.write("\nSequence 2: " + str(seq2))
    output.write("\n--- Aligned Sequences ---")
    output.write("\nSeq1: " + str(seq_1_output))
    output.write("\nSeq2: " + str(seq_2_output))


'''
used for multi-processing when we don't care about the path, only the alignment score
'''
def get_score(score_mtx):
    width = len(score_mtx[0]) - 1
    height = len(score_mtx) - 1
    score = score_mtx[height][width]
    return score


'''main function'''
if __name__ == '__main__':
    #setup arg parsers
    parser = make_arg_parser()
    args = parser.parse_args()

    #setup input fiels
    query_file = open(args.query)
    ref_file = open(args.ref)
    seq1 = ""
    seq2 = ""
    #skip the first line in the file, couldn't get file.next() to work
    first = True
    for line in query_file:
        if(not first):
            seq1 = seq1 + line.rstrip() # rstrip() removes all "\n,\t" and such
        else:
            first = False
    first = True
    for line in ref_file:
        if(not first):
            seq2 = seq2 + line.rstrip() # rstrip() removes all "\n,\t" and such
        else:
            first = False
    query_file.close()
    ref_file.close()
    if(args.matches == None): #no -m argument means that we just run default needleman_wunsch
        #TO ENABLE MULTI-PROCESSING, UNCOMMENT THE LINE BELOW and comment lines below it in if statement
        #multi_processing(seq1)

        # Comment out three lines below to enable multi-processing
        (p,s) = needleman_wunsch(seq1,seq2)
        (pb,sc) = get_path_and_score(p,s)
        write_path_and_score(pb,sc,seq1,seq2)

    else: #with a -m matches argument means we run the anchored version.
        matches = get_matches(args.matches)
        get_anchored(seq1,seq2,matches)
