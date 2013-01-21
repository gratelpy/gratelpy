import numpy as np
import networkx as nx

from fragments import get_valid_fragments, score_fragment, get_sensible_fragments, get_edges_of_sensible_fragments, get_sensible_sc, get_sscc
from subgraphs import get_all_valid_subgraphs, get_sensible_subgraphs, get_subgraph_components
from stoich import get_graph_stoich
from parse_mechanism import get_network_from_mechanism
from gensg import gen_valid_subgraphs, gen_valid_subgraphs_mp, gen_valid_subgraphs_mps
from graph import get_lpa_alpha_beta, get_path_graph
from matplotlib import pylab as pl
import itertools as it

import copy
import cPickle as pickle

from matplotlib import pylab as pl

#mechanism_file = '../mechanisms/small_G_protein_mechanism.txt'
#num_complexes = 11
#mechanism_file = '../mechanisms/testing_sscc.txt'
#num_complexes = 4
mechanism_file = '../mechanisms/small_G_protein_protein_X_mechanism.txt'
num_complexes = 15

alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism(mechanism_file, num_complexes)

#G, stoich, stoich_rank = get_graph_stoich(alpha, beta, complex_names = dict_complexes_reverse, constant_names=dict_constants_reverse)
G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

#stoich_rank = 4

fragments = get_sensible_fragments(G, stoich_rank)

f = fragments[0]
sc = get_subgraph_components(G, f)
#sensible_sc = get_sensible_sc(G, f)

path_graph = get_path_graph(sc)

all_edges_subgraphs, mixed_edges_cycles_subgraphs, all_cycles_subgraphs = get_sensible_subgraphs(sc)

#sg_sensible_sc = get_sensible_subgraphs(sensible_sc)

# sscc, substance_graph = get_sscc(sc)

# cycles_in_sscc = []
# for an_sscc in sscc:
#         cycles_in_sscc.append(nx.simple_cycles(substance_graph.subgraph(an_sscc)))

# path_graph = get_path_graph(sc)
# cycles_in_path_graph = nx.simple_cycles(path_graph)
# for cycle in cycles_in_path_graph:
# 	del cycle[-1]

# cycles_annotated_in_path_graph = []
# for cycle in cycles_in_path_graph:
# 	cycle_annotated = []
	
# 	substance_beginnings = [path[0] for path in cycle]
# 	if len(substance_beginnings)>len(set(substance_beginnings)):
# 		continue # ignore current cycle since at least one substance node beginning of more than one path
	
# 	cycle_annotated.append(substance_beginnings)
# 	cycle_annotated.append(cycle)

# 	cycles_annotated_in_path_graph.append(cycle_annotated)

# all_edges = [sc[key]['edges'] for key in sc.keys()]
# all_edges_subgraphs = it.product(*all_edges)

# cycle_combs_iters = []
# for repeat in range(1,len(cycles_annotated_in_path_graph)+1):
#     cycle_combs_iters.append(it.combinations(cycles_annotated_in_path_graph, repeat))

# cycle_combs = []
# for iterator in cycle_combs_iters:
#     for cycle_comb in iterator:
#         path_beginnings = [set(cycle[0]) for cycle in cycle_comb]
# #	print path_beginnings
# 	pairwise_disjoint = all(len(set.intersection(*set_pair))==0 for set_pair in it.product(path_beginnings, path_beginnings) if set_pair[0]!=set_pair[1])
#         if pairwise_disjoint:
#             cycle_combs.append({'substances': set.union(*path_beginnings), 'cycles': [cycle[1] for cycle in cycle_comb]})

# mixed_edges_cycles_subgraphs = []
# for cycle_comb in cycle_combs:
# 	for edges_tuple_of_subgraph in it.product(*[sc[key]['edges'] for key in sc.keys() if key not in cycle_comb['substances']]):
# 		curr_subgraph = [edge for edge in edges_tuple_of_subgraph]
# 		curr_subgraph = curr_subgraph + [path for cycle in cycle_comb['cycles'] for path in cycle]

# 		mixed_edges_cycles_subgraphs.append(curr_subgraph)

# all_cycles_subgraphs = []
# for cycle_comb in cycle_combs:
# 	if len(cycle_comb['substances'])==stoich_rank:
# 		all_cycles_subgraphs.append([path for cycle in cycle_comb['cycles'] for path in cycle])

# # any all-cycle subgraph?
# all_cycles_subgraphs_present = any(len(cycle)==stoich_rank for comb in cycle_combs for cycle in comb)

# record for each cycle the substance nodes used in it
# combine all cycles with non-overlapping substance node sets
# combine all combinations of cycles with non-overlapping substance node sets with all remaining edges

# paths_in_scc = get_paths_in_scc(f, sc)

# sc_amended = copy.deepcopy(sc)
# for substance in sc_amended:
#     for path in sc_amended[substance]['n_paths']:
#         if path not in paths_in_scc['n_paths']:
#             # path is not part of an scc of length > 1
#             # remove path from subgraph components only if it does not form a cycle by itself
#             if path[0] != path[2]:
#                 sc_amended[substance]['n_paths'].remove(path)
# #                print 'removed '+str(path)
#     for path in sc_amended[substance]['p_paths']:
#         if path not in paths_in_scc['p_paths']:
#             # path is not part of an scc of length > 1
#             # remove path from subgraph components only if it does not form a cycle by itself
#             if path[0] != path[2]:
#                 sc_amended[substance]['p_paths'].remove(path)
# #            print 'removed '+str(path)




# slow_names = ['[A]','[A.GEF]','[A.GAP]','[B]','[B.GAP]','[B.GEF]','[B.GDI]']
# slow_indices = [dict_complexes[name] for name in slow_names]

# alpha_lpa, beta_lpa = get_lpa_alpha_beta(alpha, beta, slow_indices)

# G, stoich, stoich_rank = get_graph_stoich(alpha_lpa, beta_lpa)

# edges_in_fragments = get_edges_of_sensible_fragments(G, stoich_rank)
