import math


class DeBruijnGraph:

    def __init__(self, nodes, k=21):
        """DeBruijnGraph Class Constructor."""
        self.map_size = 7
        self.nb_of_nodes = 0
        self.kmers_graph = [None] * self.map_size
        self.k = k

        for i, kmers in enumerate(nodes):
            self.add(kmers)

    def __contains__(self, node: str) -> bool:
        """Determines whether the DeBruijnGraph contains a given node."""
        index = self.index(node)
        if index < 0 or index >= self.map_size or self[index] != node:
            return False
        return True

    def __iter__(self):
        for node in self.kmers_graph:
            if node is not None and node != "AVAIL":
                yield node

    def __getitem__(self, item):
        if item < 0 or item >= self.map_size:
            raise ValueError("Specified index out of bound.")
        return self.kmers_graph[item]

    def __setitem__(self, key, value):
        if key < 0 or key >= self.map_size:
            raise ValueError("Specified index out of bound.")
        self.kmers_graph[key] = value

    def _hash_func(self, node):
        """DNA Hashing function that returns the index for a given node."""
        symbols = {'A': "1", 'T': "5", 'C': "7", 'G': "9"}
        coded_node = ""

        for strand in node:
            coded_node += symbols[strand]

        return int(coded_node) % self.map_size

    def _load_factor(self) -> float:
        """Returns the load factor of the DeBruijnGraph's underlying hash map."""
        return self.nb_of_nodes / self.map_size

    def _resize_map(self):
        """Resize hash map"""
        self.nb_of_nodes = 0
        graph_nodes = [i for i in self if i is str]
        self.map_size = self._next_size()
        self.kmers_graph = [None] * self.map_size

        for i in graph_nodes:
            self.add(i)

    def index(self, node: str):
        """Returns the index for a node"""
        index = self._get_new_index(node)
        value = self[index]
        if value is None or value is "AVAIL":
            return -1
        return index

    def _get_new_index(self, node: str):
        """Returns the index for a new node or the current index if the graph already contains the node."""
        current_index = self._hash_func(node)
        current_node = self[current_index]

        while current_node != node:
            if current_node is None or current_node == "AVAIL":
                break

            # Linear hashing
            current_index += 1
            if current_index == self.map_size:
                current_index = 0

            current_node = self[current_index]

        return current_index

    def _next_size(self):
        """Returns a prime number more than twice the size of the current hash map."""
        candidate_size = 2 * self.map_size + 1
        is_prime = False

        while not is_prime:
            is_prime = True
            for i in range(2, int(math.sqrt(candidate_size)) + 1):
                if candidate_size % i == 0:
                    is_prime = False
                    candidate_size += 2
                    break

        return candidate_size

    def add(self, node: str):
        """Adds the given node to the DeBruijnGraph."""
        if self._load_factor() >= 0.75:
            self._resize_map()

        potential_index = self._get_new_index(node)
        value_at_index = self[potential_index]

        if value_at_index != node:
            self[potential_index] = node
            self.nb_of_nodes += 1

    def remove(self, node: str):
        """Removes the given node from the DeBruijnGraph."""
        node_index = self.index(node)
        if node_index >= 0:
            self[node_index] = "AVAIL"
            self.nb_of_nodes -= 1

            # if self.load_factor < 0.1:
            #    resizeDown()

    def get_nodes(self):
        """Gets the nodes."""
        list_nodes = []
        for i in self:
            list_nodes.append(i)
        return list_nodes

    def predecessors(self, node: str):
        """Returns the predecessors of the given node."""
        predecessors = []
        for symbol in ['A', 'T', 'C', 'G']:
            candidate = symbol + node[:self.k - 1]
            if candidate in self:
                predecessors.append(candidate)
        return predecessors

    def successors(self, node: str):
        """Returns the successors of the given node."""
        successors = []
        for symbol in ['A', 'T', 'C', 'G']:
            candidate = node[1:] + symbol
            if candidate in self:
                successors.append(candidate)
        return successors

    def contigs_depth_first_search(self, start):
        """Depth first search of the hash map, returns an array of contigs."""
        contigs = []
        stepped_nodes = set()
        current_node = start
        current_contig = start
        is_backtracking = False

        while True:
            neighboring_nodes = self.successors(current_node)
            new_neighboring_nodes = [node for node in neighboring_nodes if node not in stepped_nodes]

            if not new_neighboring_nodes:
                if current_node == start:
                    break
                if not is_backtracking:
                    contigs.append(current_contig)
                    is_backtracking = True
                current_contig = current_contig[:-1]
                current_node = current_contig[len(current_contig) - self.k:]
                continue

            is_backtracking = False
            chosen_node = new_neighboring_nodes[0]

            current_contig += chosen_node[-1]
            stepped_nodes.add(chosen_node)
            current_node = chosen_node

        return contigs
