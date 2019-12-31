from Bio import SeqIO

from graph import Graph
from shotgun_to_kmer import ShotgunToKmer
from de_bruin import DeBruin
from simulated_shotgun import SimulatedShotgun
from eulerian_circuit import EulerianCircuit


if __name__ == '__main__':
    """Script to test different length strings"""
    chromosome_22 = SeqIO.read("chr22.fa", "fasta")
    chromosome_22 = str(chromosome_22.seq).lower()

    # Write headers to file.
    with open("scripta_result.txt", "w") as file:
        file.write("DeBruin and Eulerian Run on the same KMERS 5 times. Tests different start nodes.\n")
        file.write("|Sequence Length|\t|Score|\t|Best Rotated Score|\t|Number of Rotations|\n")

    # Test 10, 100, 1000, 10000
    for i in xrange(4):
        num_seq = 10 ** (i+1)

        sub_chromosome = chromosome_22[:num_seq]

        len_shotgun = min(len(sub_chromosome)/2, 150)

        shotgun_inst = SimulatedShotgun(sub_chromosome, shotgun_length = len_shotgun, coverage = 10)
        shotguns = shotgun_inst.shotguns

        len_kmer = min(len_shotgun/2, 31)
        kmer_inst = ShotgunToKmer(shotguns, kmer_length = len_kmer)
        kmers = kmer_inst.kmers

        # Run 5 times.
        for i in range(5):
            de_bruin_inst = DeBruin(kmers, True)
            graph = de_bruin_inst.graph.get_graph()

            eulerian_circuit_inst = EulerianCircuit(graph)

            score, best_score, idx = eulerian_circuit_inst.compare_string(sub_chromosome)

            with open("scripta_result.txt", "a") as file:
                file.write("{}\t{}\t{}\t{}\n".format(num_seq, score, best_score, idx))

            print num_seq, score, best_score, idx
