import unittest

import os
path = os.path.split(os.path.realpath(__file__))[0]

from multiprocessing import freeze_support
from gratelpy.gensg import gen_valid_subgraphs_mps
from gratelpy.parse_mechanism import get_network_from_mechanism
from gratelpy.stoich import get_graph_stoich
from gratelpy.fragments import get_sensible_fragments
from gratelpy.pyin import resource_path

def main():
    alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, \
        dict_constants_reverse = get_network_from_mechanism(
                resource_path(os.path.join(path, '..', 'mechanisms'),
                              'reversible_substrate_inhibition.txt'), 
                4)
			
    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)
    valid_fragments = get_sensible_fragments(G, stoich_rank)

    gen_valid_subgraphs_mps(G, valid_fragments, 4)

if __name__ == '__main__':
    # Required on Windows 
    # http://docs.python.org/2/library/multiprocessing.html#multiprocessing.freeze_support
    freeze_support()
    main()