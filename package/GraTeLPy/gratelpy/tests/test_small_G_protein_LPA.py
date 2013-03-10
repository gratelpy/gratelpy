import numpy as np
import networkx as nx

from fragments import get_valid_fragments, score_fragment, get_sensible_fragments, get_edges_of_sensible_fragments
from subgraphs import get_all_valid_subgraphs
from stoich import get_graph_stoich
from parse_mechanism import get_network_from_mechanism
from gensg import gen_valid_subgraphs, gen_valid_subgraphs_mp, gen_valid_subgraphs_mps
from graph import get_lpa_alpha_beta

import cPickle as pickle

mechanism_file = '../mechanisms/small_G_protein_mechanism.txt'
num_complexes = 11

alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism(mechanism_file, num_complexes)

slow_names = ['[A]','[A.GEF]','[A.GAP]','[B]','[B.GAP]','[B.GEF]','[B.GDI]']
slow_indices = [dict_complexes[name] for name in slow_names]

alpha_lpa, beta_lpa = get_lpa_alpha_beta(alpha, beta, slow_indices)

G, stoich, stoich_rank = get_graph_stoich(alpha_lpa, beta_lpa)

edges_in_fragments = get_edges_of_sensible_fragments(G, stoich_rank)
