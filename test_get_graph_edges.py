from stoich import get_substance_adjacency, get_graph_stoich
from graph import get_graph_edges
from parse_mechanism  import get_network_from_mechanism
import numpy as np

alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/test_mechanism.txt', 4)
G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
edges = get_graph_edges(G)
