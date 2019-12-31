from Bio import SeqIO

from graph import Graph
from shotgun_to_kmer import ShotgunToKmer
from de_bruin import DeBruin
from simulated_shotgun import SimulatedShotgun
from eulerian_circuit import EulerianCircuit
from fuzzywuzzy import fuzz

import random, copy, time

NUM_ITERATIONS = 100

if __name__ == '__main__':
    """Script to test if scores go down as larger segments are used"""
    chromosome_22 = SeqIO.read("chr22.fa", "fasta")
    chromosome_22 = str(chromosome_22.seq).lower()

    # Write headers to file.
    with open("scriptb_result.txt", "w") as file:
        file.write("DeBruin and Eulerian Run on the different KMERS 15 times.\n")
        file.write("|Sequence Length|\t|Score|\t|Total Time|\n")

    # Test 10, 100, 1000, 10000, 100000
    steps = [25, 50, 100, 1000, 10000, 25000, 50000, 75000, 100000]
    for num_seq in steps:
        min_score = 100
        max_score = 0
        total = 0
        total_time = 0

        for count in xrange(NUM_ITERATIONS):
            start = time.time()
            rand_index = random.randint(0, (len(chromosome_22)-num_seq)-1)
            sub_chromosome = copy.deepcopy(chromosome_22[rand_index:rand_index+num_seq])

            len_shotgun = min(len(sub_chromosome)/2, 150)

            shotgun_inst = SimulatedShotgun(sub_chromosome, shotgun_length = len_shotgun, coverage = 10)
            shotguns = shotgun_inst.shotguns

            len_kmer = min(len_shotgun/2, 31)
            kmer_inst = ShotgunToKmer(shotguns, kmer_length = len_kmer)
            kmers = kmer_inst.kmers

            de_bruin_inst = DeBruin(kmers, True)
            graph = de_bruin_inst.graph.get_graph()

            eulerian_circuit_inst = EulerianCircuit(graph)

            score = fuzz.ratio(eulerian_circuit_inst.assembled_seq, sub_chromosome)
            total += score
            total_time += (time.time() - start)

            min_score = min(min_score, score)
            max_score = max(max_score, score)

            print num_seq, count

        with open("scriptb_result.txt", "a") as file:
            file.write("{}\t{}\t{}\n".format(num_seq, total / float(NUM_ITERATIONS), total_time))
