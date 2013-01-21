import networkx as nx
from graph import get_all_cliques
from subgraphs import get_sensible_subgraphs_paths_in_cycles_using_cycle_graph
frag = (('s9', 's8', 's2', 's1', 's5', 's4'), ('w12', 'w9', 'w4', 'w1', 'w6', 'w10'))
sc = {'s9': {'p_paths': [('s9', 'w12', 's1')], 'edges': [('s9', 'w12')], 'n_paths': []}, 's8': {'p_paths': [('s8', 'w9', 's4')], 'edges': [('s8', 'w9')], 'n_paths': []}, 's2': {'p_paths': [('s2', 'w4', 's5')], 'edges': [('s2', 'w4'), ('s2', 'w1')], 'n_paths': [('s2', 'w4', 's4'), ('s2', 'w1', 's1')]}, 's1': {'p_paths': [], 'edges': [('s1', 'w1')], 'n_paths': [('s1', 'w1', 's2')]}, 's5': {'p_paths': [('s5', 'w6', 's2')], 'edges': [('s5', 'w6')], 'n_paths': []}, 's4': {'p_paths': [('s4', 'w4', 's5'), ('s4', 'w10', 's9')], 'edges': [('s4', 'w4'), ('s4', 'w10')], 'n_paths': [('s4', 'w4', 's2')]}}

sg = get_sensible_subgraphs_paths_in_cycles_using_cycle_graph(sc)
