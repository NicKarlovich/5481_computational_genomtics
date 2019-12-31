# Graphs as dictionary https://www.python.org/doc/essays/graphs/

# A representation of a directed graph as a dictionary with vertices and list of edges
from vertex import vertex

class Graph:

    def __init__(self):
        self.graph = {}

    # add a vertex to the graph, if not already present in graph
    def add_vertex(self, vert, edges=[]):
        if (self.graph.has_key(vert)):
            return False
        else:
            self.graph[vert] = edges
            return True

    # adds an edge between two existing vertcies in the graph
    def add_edge(self, vert1, vert2):
        # create edge if both verticies are present
        if (self.graph.has_key(vert1) and self.graph.has_key(vert2)):
            self.graph[vert1] = self.graph[vert1] + [vert2]
            return True
        else:
            return False

    # adds an edge between 2 vertcies, creates verticies if not in graph
    def create_edge(self, vert1, vert2):
        # create vert2 if its not in the graph
        if (not self.graph.has_key(vert2)):
            self.graph[vert2] = []

        # create edge between vert1 and vert2
        if (self.graph.has_key(vert1)):
            self.graph[vert1] = self.graph[vert1] + [vertex(vert2)]
        else:
            self.graph[vert1] = [vertex(vert2)]

    # adds an edge between 2 vertcies, creates verticies if not in graph
    def create_weighted_edge(self, vert1, vert2, w):
        # create vert2 if its not in the graph
        if (not self.graph.has_key(vert2)):
            self.graph[vert2] = []

        # create edge between vert1 and vert2
        if (self.graph.has_key(vert1)):
            self.graph[vert1] = self.graph[vert1] + [vertex(vert2, w)]
        else:
            self.graph[vert1] = [vertex(vert2, w)]

    # gets all outbound edges from a vertex
    def get_edges(self, vert):
        if(self.graph.has_key(vert)):
            return self.graph[vert], True
        return [], False

    def print_dict(self):
        """Print the graph in an easy to view format"""
        for vertex in self.graph:
            tmp = []
            for sub_vertex in self.graph[vertex]:
                weight = sub_vertex.get_weight()
                if weight == None:
                    tmp.append(sub_vertex.get_vertex())
                else:
                    vert_weight = [sub_vertex.get_vertex(), weight]
                    tmp.append(vert_weight)
            print("{} -> {}".format(vertex, tmp))

    def order_dict(self):
        """If a vertex is in pointing to itself, put that vertex at the beginning"""
        for vertex in self.graph:
            tmp = []
            for sub_vertex in self.graph[vertex]:
                if sub_vertex.get_vertex() == vertex:
                    tmp = [sub_vertex] + tmp
                else:
                    tmp += [sub_vertex]
            self.graph[vertex] = tmp

    def get_graph(self):
        """Return the graph"""
        return self.graph


if __name__ == "__main__":
    test = Graph()

    # add some verticies
    test.add_vertex("A")
    test.add_vertex("B")

    # add an edge
    test.add_edge("A", "B")

    # create some edges and verticies
    test.create_edge("B", "C")
    test.create_edge("B", "D")
    test.create_edge("C", "D")

    # edge not added because D does not exist
    test.add_edge("E", "F")
    # edge not added because F does not exist
    test.add_edge("D", "F")

    print("Edges from A:" + str(test.get_edges("A")))
    print("Edges from B:" + str(test.get_edges("B")))
    print("Edges from C:" + str(test.get_edges("C")))
    print("Edges from D:" + str(test.get_edges("D")))
    print("Edges from E:" + str(test.get_edges("E")))
    print("Edges from F:" + str(test.get_edges("F")))

    test.print_dict()

#    A --------------> B
#                     | |
#           +---------+ +--------+
#           |                    |
#           |                    |
#       C <-+                    +-> D
#       |                            ^
#       |                            |
#       |                            |
#       +----------------------------+
