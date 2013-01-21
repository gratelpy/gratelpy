from stoich import get_substance_adjacency
from parse_mechanism  import get_network_from_mechanism
import numpy as np

alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('../mechanisms/test_mechanism.txt', 4)

subs_adj = get_substance_adjacency(alpha, beta)
