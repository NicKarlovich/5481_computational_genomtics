from Bio import SeqIO

from graph import Graph
from shotgun_to_kmer import ShotgunToKmer
from de_bruin import DeBruin
from simulated_shotgun import SimulatedShotgun
from eulerian_circuit import EulerianCircuit

import random, copy

from matplotlib import pyplot as plt
from numpy.polynomial.polynomial import polyfit
import numpy as np

def plot_scores(score, rotated_score, errors):
    """Plot scores vs error"""
    score = np.array(score)
    rotated_score = np.array(rotated_score)
    errors = np.array(errors)

    plt.title("Error vs Score")
    plt.ylabel('Score')
    plt.xlabel('Percent Error')

    plt.scatter(errors, score, c = 'b', marker = '+')
    plt.scatter(errors, rotated_score, c = 'g', marker =  'x')

    plt.savefig('Error_Vs_Score.png', bbox_inches='tight')

    plt.show()

if __name__ == '__main__':
    """Script to test simmulation error"""
    chromosome_22 = SeqIO.read("chr22.fa", "fasta")
    chromosome_22 = str(chromosome_22.seq).lower()

    # Write headers to file.
    with open("scriptc_result.txt", "w") as file:
        file.write("Everything run on same sub_chromosome 5 times. Different error introduced.\n")
        file.write("|Sequence Length|\t|Error|\t|Score|\t|Best Rotated Score|\n")

    test_errors = []
    for i in range(101):
        test_errors += [i]

    sequence_length = 1000

    rand_index = random.randint(0, (len(chromosome_22)-sequence_length)-1)
    sub_chromosome = copy.deepcopy(chromosome_22[rand_index:rand_index+sequence_length])

    score_lst = []
    rotated_score_lst = []

    for error in test_errors:

        total_score = 0
        total_best_score = 0
        num_iterations = 5

        for i in range(num_iterations):
            shotgun_inst = SimulatedShotgun(sub_chromosome, shotgun_length = 150, coverage = 10, error_percent = error)
            shotguns = shotgun_inst.shotguns

            kmer_inst = ShotgunToKmer(shotguns, kmer_length = 31)
            kmers = kmer_inst.kmers

            de_bruin_inst = DeBruin(kmers, True)
            graph = de_bruin_inst.graph.get_graph()

            eulerian_circuit_inst = EulerianCircuit(graph)

            score, best_score, idx = eulerian_circuit_inst.compare_string(sub_chromosome)

            total_score += score
            total_best_score += best_score

        score = total_score/num_iterations
        rotated_score = total_best_score/num_iterations
        score_lst += [score]
        rotated_score_lst += [rotated_score]

        with open("scriptc_result.txt", "a") as file:
            file.write("{}\t{}\t{}\t{}\n".format(sequence_length, error, score, rotated_score))

    plot_scores(score_lst, rotated_score_lst, test_errors) 
