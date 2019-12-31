class vertex:
    def __init__(self, vertex, weight = None):
        self.vertex = vertex
        self.visited = False
        self.weight = weight

    def get_vertex(self):
        """Return the vertex"""
        return self.vertex

    def get_visited(self):
        """Return whether the vertex was visited"""
        return self.visited

    def get_weight(self):
        """Return the weight"""
        return self.weight

    def set_visited(self, visit):
        """Set visited to visit"""
        self.visited = visit
