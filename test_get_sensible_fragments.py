from stoich import get_substance_adjacency, get_graph_stoich
from fragments import get_sensible_fragments
from parse_mechanism  import get_network_from_mechanism
import numpy as np

#alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/test_mechanism.txt', 4)
#alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/double_layer_mapk_mechanism.txt', 12)
alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/single_layer_mapk_mechanism.txt', 9)

G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

fragments = get_sensible_fragments(G, stoich_rank)
