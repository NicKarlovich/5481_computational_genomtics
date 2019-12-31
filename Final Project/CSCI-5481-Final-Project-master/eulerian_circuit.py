import random
from collections import Counter

from graph import Graph
from shotgun_to_kmer import ShotgunToKmer
from de_bruin import DeBruin
from simulated_shotgun import SimulatedShotgun

from Bio import SeqIO

from fuzzywuzzy import fuzz

class EulerianCircuit:
    def __init__(self, graph):
        self.graph = graph
        self.vertex_order = []
        self.assembled_seq = ""

        # Check if there are any imbalanced nodes.
        self.imbalanced_nodes = self.get_imbalanced_nodes(self.graph)

        if len(self.imbalanced_nodes) != 0:
            # Iterate through all imbalanced nodes to visit all edges.
            # TODO: The order of nodes to visit should be made differently.
            for start_node in self.imbalanced_nodes:
                self.visit_all_edges(start_node)

        else:
            # Make a random selection for the first node.
            # NOTE: Because of the randomness, each run will produce a different assembly.
            # TODO: This choice is arbitrary right now. We should either improve
            #       the method for selecting the first node, or iterate through all nodes.
            start_node = random.choice(list(self.graph.keys()))

            # Visit all Edges and Assemble the sequences.
            self.visit_all_edges(start_node)

        # Assemble the Sequence.
        self.assemble_sequence()

    def assemble_sequence(self):
        """Assembled the sequence from the Eulerian Circuit"""
        num_vertex = len(self.vertex_order)
        for i, v in enumerate(self.vertex_order):
            if i+1 != num_vertex:
                self.assembled_seq += v[0]
            else:
                self.assembled_seq += v

    def print_assembled_sequence(self):
        """Print the assembled sequence"""
        print "\nAssemebled Sequence:\n", self.assembled_seq

    def visit_all_edges(self, start_node):
        """Visit all the edges to form the Eulerian Circuit."""
        # Iterate through each node in the graph.
        stack = [start_node]
        while len(stack):
            node = stack.pop()
            for n in self.graph[node]:
                # Explore the path of each node, and record the order.
                if not n.get_visited():
                    n.set_visited(True)
                    node2 = n.get_vertex()
                    self.vertex_order.append(self.get_edge(node, node2))
                    stack.append(node2)
                    break

    def print_vertex_order(self):
        """Print the order or the vertexes"""
        for i, n in enumerate(self.vertex_order):
            spaces = " "*i
            print(spaces + n)

    def get_edge(self, node1, node2):
        """Get the edge name between two given nodes"""
        return node1 + node2[-1]

    def get_imbalanced_nodes(self, graph):
        """Get all nodes that are imbalanced"""
        # Get all destination nodes and flatten out into a single list.
        dest_nodes = graph.values()
        dest_nodes = [item.get_vertex() for sublist in dest_nodes for item in sublist]
        dest_nodes_count = Counter(dest_nodes)

        imbalanced_nodes = []
        # Iterate through each node and determine if it is imbalanced.
        for node in graph:
            if len(graph[node]) != dest_nodes_count[node]:
                 imbalanced_nodes.append(node)

        # Randomize before returning.
        random.shuffle(imbalanced_nodes)
        return imbalanced_nodes

    def compare_string(self, original):
        """Compare the original String to the assembled sequence and test its rotations"""
        # Align the two sequences.

        def right_rotate(string, num):
            """Right Rotate a string"""
            front = string[:num]
            end = string[num:]
            return end + front

        def binary_search(string, original):
            """Use binary search to find the maximum Fuzzy Score - O(LOG(N))"""
            # Initialize indexes
            idx_start = 0
            idx_middle = len(string) / 2
            idx_end = len(string) - 1

            # Initialize scores
            score_start, score_middle, score_end = [0]*3

            best_score = 0
            best_idx = None
            max_found = False

            # Search for the best rotation score.
            while not max_found:
                # Rotate all strings.
                str_start = right_rotate(string, idx_start)
                str_middle = right_rotate(string, idx_middle)
                str_end = right_rotate(string, idx_end)

                # Calculate all scores.
                # Optimized to only calculate neccessary scores. (fuzz.ratio takes a long time for long strings)
                if score_end > score_start:
                    score_middle = fuzz.ratio(str_middle, original)
                    score_end = fuzz.ratio(str_end, original)
                elif score_end < score_start:
                    score_start = fuzz.ratio(str_start, original)
                    score_middle = fuzz.ratio(str_middle, original)
                else:
                    score_start = fuzz.ratio(str_start, original)
                    score_middle = fuzz.ratio(str_middle, original)
                    score_end = fuzz.ratio(str_end, original)

                if score_start == 100 or score_middle == 100 or score_end == 100:
                    # Best score found. Stop searching.
                    best_score = 100
                    if score_start == 100:
                        best_idx = idx_start
                    elif score_middle == 100:
                        best_idx = idx_middle
                    else:
                        best_idx = idx_end
                    max_found = True

                if idx_start == idx_middle:
                    # Converging, pick the best score.
                    best_score = max(score_start, score_middle)
                    if score_start > score_middle:
                        best_idx = idx_start
                    else:
                        best_idx = idx_middle
                    max_found = True

                else:
                    # Adjust indexes and use already calculated scores.
                    if score_end > score_start:
                        idx_start = idx_middle
                        score_start = score_middle
                    else:
                        idx_end = idx_middle
                        score_end = score_middle

                    idx_middle = (idx_end+idx_start)/2

            return (best_score, best_idx)

        # Compute the score and the best score.
        score = fuzz.ratio(self.assembled_seq, original)
        best_score, best_idx = binary_search(self.assembled_seq, original)

        return (score, best_score, best_idx)

if __name__ == '__main__':
    """Run `python eulerian_circuit.py` to test"""
    tst_str = "AABBCCDD"

    shotgun_inst = SimulatedShotgun(tst_str, shotgun_length=6, coverage = 10)
    shotguns = shotgun_inst.shotguns

    kmer_inst = ShotgunToKmer(shotguns, kmer_length=3)

    de_bruin_inst = DeBruin(kmer_inst.kmers, True)

    de_bruin_inst.graph.order_dict()
    de_bruin_inst.graph.print_dict()

    eulerian_circuit_inst = EulerianCircuit(de_bruin_inst.graph.get_graph())

    eulerian_circuit_inst.print_vertex_order()
    eulerian_circuit_inst.print_assembled_sequence()
    print "\nOriginal Sequence:\n", tst_str

    score, best_score, idx = eulerian_circuit_inst.compare_string(tst_str)

    print score, best_score, idx
