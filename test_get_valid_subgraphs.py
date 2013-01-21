from stoich import get_substance_adjacency, get_graph_stoich
from fragments import get_sensible_fragments
from subgraphs import *
from parse_mechanism  import get_network_from_mechanism
import numpy as np

#alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/test_mechanism.txt', 4)
#alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/single_layer_mapk_mechanism.txt', 9)
alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/double_layer_mapk_mechanism.txt', 12)
G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

fragments = get_sensible_fragments(G, stoich_rank)

f = fragments[0]
sc = get_subgraph_components(G, f)
all_subgraphs = get_all_subgraphs(sc, stoich_rank)
sensible_subgraphs = get_sensible_subgraphs(sc)

#valid_subgraphs = validate_subgraphs(all_subgraphs, sc, f)
#f, sc, valid_subgraphs = get_all_valid_subgraphs(G, stoich_rank, f)
