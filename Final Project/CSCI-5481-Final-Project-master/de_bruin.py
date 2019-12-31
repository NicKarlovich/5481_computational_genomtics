from threading import Thread, Lock
import string
from graph import Graph
from shotgun_to_kmer import ShotgunToKmer

THREADS = 3
OVERLAP_THRESH = 1

# create a de bruin graph from a a list of kmers


class DeBruin:
    def __init__(self, kmer_dict, eularian=False):
        self.graph = Graph()
        self.result_lock = Lock()
        self.distinct_kmers = kmer_dict.keys()
        self.kmer_dict = kmer_dict

        threads = []
        num_kmers_per_thread = (len(self.distinct_kmers) / THREADS)
        for i in xrange(THREADS):
            start_index = i * num_kmers_per_thread

            if i == (THREADS - 1):  # add the remaining kmers to last thread
                end_index = len(self.distinct_kmers)
            else:
                end_index = (i+1) * num_kmers_per_thread

            if (eularian):
                thread = Thread(target=self.create_eularian_graph,
                                args=(start_index, end_index))
            else:
                thread = Thread(target=self.create_hamiltonaian_graph,
                                args=(start_index, end_index))

            threads.append(thread)
            thread.start()

        for thread in threads:
            thread.join()

    # create a graph were edges repsent kmers, nodes are k-1 mers
    def create_eularian_graph(self, kmer_start_index, kmer_end_index):

        list_of_edges = []
        # for every distinct kmer in given range
        for i in xrange(kmer_start_index, kmer_end_index):
            current = self.distinct_kmers[i]
            number_of_current = self.kmer_dict[current]

            # for repeat kmers
            for j in xrange(number_of_current):
                # add edge between left and right k - 1 mers
                l, r = self.left_and_right_kmers(current)
                list_of_edges += [(l, r)]

        # add edges to graph
        self.result_lock.acquire()
        for (l, r) in list_of_edges:
            self.graph.create_edge(l, r)

        self.result_lock.release()

    # create graph with kmers as nodes and edges with weights of overlap with other nodes
    def create_hamiltonaian_graph(self, kmer_start_index, kmer_end_index):

        list_of_edges = []
        length = len(self.distinct_kmers)
        # for every distinct kmer in given range
        for i in xrange(kmer_start_index, kmer_end_index):
            current = self.distinct_kmers[i]
            # for every other kmer
            for j in xrange(length):
                if (i != j):
                    # create an edge if overlap between kmers is > than threshold
                    other_kmer = self.distinct_kmers[j]
                    overlap = self.overlap(current, other_kmer)
                    if (overlap >= OVERLAP_THRESH):
                        list_of_edges += [(current, other_kmer, overlap)]

        # add edges to graph
        self.result_lock.acquire()
        for (k1, k2, w) in list_of_edges:
            self.graph.create_weighted_edge(k1, k2, w)

        self.result_lock.release()

    # split kmer into left and right k - 1 mers

    def left_and_right_kmers(self, kmer):
        return kmer[:-1], kmer[1:]

    # determines the amount of overlap between suffix of kmer1 and prefix of kmer2
    def overlap(self, kmer1, kmer2):
        length = len(kmer1)
        overlap = length

        suffix = kmer1[length - overlap:length]
        prefix = kmer2[0:overlap]

        while(prefix != suffix and overlap > 0):
            overlap -= 1
            suffix = kmer1[length - overlap:length]
            prefix = kmer2[0:overlap]

        return overlap


if __name__ == "__main__":
    # eularian debruin graph with k-1 mers
    # example from lecture on oct 9 at 1 hour in
    shotgun_inst = ShotgunToKmer(["AAABBBA"], kmer_length=3)
    print(shotgun_inst.kmers)
    de_bruin_inst = DeBruin(shotgun_inst.kmers, True)
    de_bruin_inst.graph.print_dict()

    # hamiltonian debruin graph with k mers
    # example from lecture on oct 7 at 1 hour in
    shotgun_inst = ShotgunToKmer(["ATCCAGT"], kmer_length=3)
    print(shotgun_inst.kmers)
    de_bruin_inst = DeBruin(shotgun_inst.kmers, False)
    de_bruin_inst.graph.print_dict()
