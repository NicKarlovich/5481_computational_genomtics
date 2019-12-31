#  Name: Nick Karlovich
#  x500: karlo015
# Class: CSCI 5481
#  Year: Fall 2019
import sys, os
import argparse
import numpy as np
from subprocess import Popen, PIPE

# the master tree where we store all the nodes
tree = []

# Node class
class Node:
    def __init__(self, key, val, node = [], dist = 0.0, p = None):
        self.k = key # we assign (based on ntips)
        self.v = val # id from file
        self.nodes = node
        self.dist_to_parent = dist
        self.parent = p
    
    def add_child(self, a_node):
        self.nodes.append(a_node)
    
    def set_dist(self, fdist): #set distance to parent
        self.dist_to_parent = fdist

    def set_parent(self, parent):
        self.parent = parent
    
    # a test function for printing out status of a node
    def out(self):
        if(self.parent == None):
            print("key: " + str(self.k) + "  val: " + str(self.v) + "  parent: None  dist_to_parent: 0")
        else:
            print("key: " + str(self.k) + "  val: " + str(self.v) + "  parent: " + str(self.parent.k) + "  dist_to_parent: " + str(self.dist_to_parent))

# helper functions
def preorder_traversal(root):
    root.out()
    for node in root.nodes:
        if(node != None):
            preorder_traversal(node)

def postorder_traversal(root):
    for node in root.nodes:
        if(node != None):
            postorder_traversal(node)
    root.out()

# a function to format output into newick style
def newick_output(root, file):
    if(len(root.nodes) == 0):
        file.write(str(root.v) + ":" + str(root.dist_to_parent))
    else:
        file.write("(")
        for x in range(0,len(root.nodes)):
            newick_output(root.nodes[x], file)
            if(x != len(root.nodes) - 1):
                file.write(",")
        file.write(")")
        if(root.dist_to_parent != 0.0):
            file.write(":" + str(root.dist_to_parent))
    
# formatting edges for Problem 1 file
def edges_print(root, child):
    return str(root.k) + "\t" + str(child.k) + "\t" + str(child.dist_to_parent) + "\n"

# outputting for Problem 1 file
def edges_output(root, file):
    for node in root.nodes:
        file.write(edges_print(root, node))
        if(node != None):
            edges_output(node, file)
    return 0

# helper functions
def print_all_nodes():
    for x in tree:
        x.out()

def find_root():
    for x in tree:
        if(x.parent == None):
            return x
    return None

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
    
# Dissimilarity calculation is sometimes a tiny bit imprecise but good enough for our implementation according to TA.
def find_dissim(seq1, seq2):
    assert(len(seq1) == len(seq2))
    if(seq1 == seq2): # if they are the same don't do the computation
        return 0
    numerator = 0
    denominator = len(seq1)
    for idx in range(len(seq1)):
        if(seq1[idx] == seq2[idx]):
            numerator += 1

    # return (1 - (numerator / denominator)),15)
    return round(1 - (numerator / float(denominator)),15)

# when calculating a,b,c we know they are the last three
# values in the dist_mtx so we use absolute references
# rather than variables like i or j, also more constants
# Calculates distance between all three remaining nodes
def calc_a_b_c_u(dist_mtx):
    distance_i = 0.0
    distance_j = 0.0
    for a in range(1, 4):
        distance_i += float(dist_mtx[1][a])
        distance_j += float(dist_mtx[a][2])
    
    dist_a_b = float(dist_mtx[1][2])
    dist_a_c = float(dist_mtx[1][3])
    a = (0.5 * dist_a_b) + (0.5) * (distance_i - distance_j)
    b = dist_a_b - a
    c = dist_a_c - a
    return (a,b,c)

# calculates distance between two nodes (a,b) and their combined output (u)
def calc_a_b_u(i, j, dist_mtx):
    n = (len(dist_mtx[0])) - 1 # how many taxa
    distance_i = 0.0
    distance_j = 0.0
    for a in range(1, n + 1):
        distance_i += float(dist_mtx[i + 1][a])
        distance_j += float(dist_mtx[a][j + 1])
    
    dist_a_b = float(dist_mtx[i + 1][j + 1])
    a = (0.5 * dist_a_b) + (1 / (2 * (n - 2))) * (distance_i - distance_j)
    b = dist_a_b - a
    return (a,b)

# creating the Q-matrix
def calc_q(i,j,distance_mtx):
    n = (len(distance_mtx[0])) - 1 # how many taxa
    distance_i = 0.0
    distance_j = 0.0
    for a in range(1,n + 1):
        distance_i += float(distance_mtx[i + 1][a])
        distance_j += float(distance_mtx[a][j + 1])
    return round((n - 2.0) * float(distance_mtx[i + 1][j + 1]) - distance_i - distance_j,15)

#returns the matrix q and the max value in the matrix
def build_q(distance_mtx):
    min = (sys.maxsize, 0, 0)
    n = len(distance_mtx[0]) - 1
    q = np.zeros((n,n), dtype = float)
    for i in range(0, n):
        for j in range(0, n):
            if(i != j):
                q[i][j] = calc_q(i,j,distance_mtx)
                if(q[i][j] < min[0]):
                    min = (q[i][j], i, j)
    return (q,min)

# We decided we want to remove entries at indexes
# i and j from our distance matrix and replace it 
# with a new node u, here we are removing i and j
# from the distance matrix.
def reconstruct_dist(dist_mtx,i_r, j_r):
    new_dist = []
    for i in range(0,len(dist_mtx)):
        new_row = []
        for j in range(0,len(dist_mtx[0])):
            if((i != i_r and i != j_r) and (j != j_r and j != i_r)):
                new_row.append(str(dist_mtx[i][j]))
        if(len(new_row) > 0):
            new_dist.append(new_row)
    return new_dist

# Here is where we calculate the 1 x n matrix that 
# we will be adding to the updated dist_mtx
# plan is to add new entry to the end of the array which is why
# the last value in the array is 0 because you are 0 dist away from yourself.
# the other values in 1 x n array is the distance from each current node in 
# dist_mtx to the new one were adding here.
def create_new_entry(i, j, old_dist, new_dist, new_val):
    output = []
    output.append(new_val[0])
    keys = old_dist[0]
    for entry in new_dist[0]:
        if entry != '0':
            k = keys.index(entry)
            output.append(str(dist_update(i,j,k,old_dist)))
    output.append('0')
    return output

# helper function
def dist_update(i, j, k, dist_mtx):
    return 0.5 * (float(dist_mtx[i][k]) + float(dist_mtx[k][j]) - float(dist_mtx[i][j]))

# combine the outputs of reconstruct_dist and create_new_entry
# into the final dist_mtx that we will use in the next iteration of nei_saitou
def combine_new_entry_new_dist(ne, nd):
    for r in range(0,len(nd)):
        nd[r].append(ne[r]) 
    nd.append(ne)
    return nd

# helper
def get_node_by_key(key):
    for node in tree:
        if (node.k == key):
            return node
    return None

#returns list of child nodes
def get_child_nodes(keys):
    out = []
    for x in keys:
        out.append(get_node_by_key(x))
    return out

def nei_saitou(dist_mtx):
    n = 0
    while(len(dist_mtx) > 4): # we do nei_saitou until we have 3 distinct nodes left then special case we handle outside of loop.
        (q,(min_val, i, j)) = build_q(dist_mtx)
        n = len(tree)
        c1 = str(dist_mtx[0][i + 1])
        c2 = str(dist_mtx[0][j + 1])
        c = [c1,c2]
        [c1, c2] = get_child_nodes(c)
        new_node = Node(str(n + 1), None, [c1,c2])

        #setting parents and distances for new nodes
        (a_d, b_d) = calc_a_b_u(i, j, dist_mtx)
        c1.set_parent(new_node)
        c1.set_dist(a_d)
        c2.set_parent(new_node)
        c2.set_dist(b_d)
        tree.append(new_node)

        #reconstructing the distance matrix to account for 2 removed entries and 1 new entry.
        new_val = (str(n + 1), (str(dist_mtx[0][i+1]) + "-" + str(dist_mtx[j + 1][0])))
        new_dist_mtx = reconstruct_dist(dist_mtx, i + 1, j + 1)
        new_entry = create_new_entry(i + 1, j + 1, dist_mtx, new_dist_mtx, new_val)
        dist_mtx = combine_new_entry_new_dist(new_entry, new_dist_mtx)
    
    # Special case for top of tree where tree has 3 child nodes not 2.
    n += 1
    [c1, c2, c3] = get_child_nodes([str(dist_mtx[0][1]), str(dist_mtx[0][2]), str(dist_mtx[0][3])])
    new_node = Node(str(n + 1), None, [c1, c2, c3])
    (d1, d2, d3) = calc_a_b_c_u(dist_mtx)

    # Updating parents and distances
    c1.set_parent(new_node)
    c1.set_dist(d1)
    c2.set_parent(new_node)
    c2.set_dist(d2)
    c3.set_parent(new_node)
    c3.set_dist(d3)
    tree.append(new_node)

if __name__ == '__main__':
    # -- Problem 1 Starts Here --
    parser = make_arg_parser()
    args = parser.parse_args()
    input_file = open(args.input) #open hw3.fna
    input_lines = input_file.readlines()
    query_tuples = []   
    for x in range(0,len(input_lines),2):
        query_tuples.append((input_lines[x][1:-1], input_lines[x+1][:-1]))
    #output = open("hw4-var4-genetic-distances.txt", 'w+') #problem 1 outputfile
    col_headers = ""
    file_col_headers = ""
    dist_mtx = []
    count = 1
    for x in query_tuples:
        (y,z) = x
        col_headers = col_headers + "\t" + str(count)
        file_col_headers = file_col_headers + "\t" + str(y) #we need a different header for the file output that what we store in our internal nei_saitou matrix.
        tree.append(Node(str(count),str(y)))
        count += 1
    #output.write(file_col_headers)
    # -- Problem 1 Done Here --
    # -- Problem 2 Starts Here --
    dist_mtx.append(col_headers.split("\t"))
    dist_mtx[0][0] = '0'
    for row_iter in query_tuples:
        (a,b) = row_iter
        col_info = str(a) + "\t"
        for col_iter in query_tuples:
            (c,d) = col_iter
            col_info = col_info + str(round(find_dissim(b,d),11)) + "\t"
        #output.write("\n" + col_info)
        dist_mtx.append(col_info.split("\t")[:-1])
    #output.close()
    nei_saitou(dist_mtx)
    # -- Problem 2 Done Here --
    # -- Problem 3 Starts Here --
    root = find_root()
    num_of_var_from_filename = args.input[-5:-4]
    if(args.input == "hw4-100-sequences.fna"):
        edges_filename = "edges4.txt"
    else:
        edges_filename = "edges4-var" + str(num_of_var_from_filename) + ".txt"
    edges = open(edges_filename, "w+")
    edges_output(root, edges)
    edges.close()
    if(args.input == "hw4-100-sequences.fna"):
        tree_filename = "tree4.tre"
    else:
        tree_filename = "tree4-var" + str(num_of_var_from_filename) + ".tre"
    tree_file = open(tree_filename, "w+")
    newick_output(root, tree_file)
    tree_file.write(";")
    tree_file.close()
    # -- Problem 3 Done Here --


