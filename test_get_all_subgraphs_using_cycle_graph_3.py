import numpy as np
import networkx as nx

from fragments import get_valid_fragments, score_fragment, get_sensible_fragments, get_edges_of_sensible_fragments, get_sensible_sc, get_sscc, score_fragment
from subgraphs import get_all_valid_subgraphs, get_sensible_subgraphs, validate_subgraphs
from stoich import get_graph_stoich
from parse_mechanism import get_network_from_mechanism
from gensg import gen_valid_subgraphs, gen_valid_subgraphs_mp, gen_valid_subgraphs_mps
from graph import get_lpa_alpha_beta, get_path_graph, get_valid_path_graph_cycles
from matplotlib import pylab as pl
import itertools as it

import copy
import cPickle as pickle

from matplotlib import pylab as pl

def get_sensible_subgraphs_paths_in_cycles_using_cycle_graph(sc):
    # code from test_get_all_subgraphs_using_cycle_graph_2.py

    # defintion of cliques / complete subgraphs
    # http://en.wikipedia.org/wiki/Clique_(graph_theory)
    # NOTE: in this method, the word 'subgraph' or phrase 'complete subgraph' most likely refers to the above graph-theoretic concept of cliques
    # and not the concept used by Mincheva et al.

    # get path graph
    path_graph = get_path_graph(sc)

    #cycles_in_path_graph = nx.simple_cycles(path_graph)
    cycles_in_path_graph = list(get_valid_path_graph_cycles(path_graph))

    # create cycle graph
    # nodes are cycles
    # edge(c1,c2) exists iff c1 and c2 have at least one substance as beginning of path in common
    cycle_graph = nx.Graph()
    for cycle_1 in cycles_in_path_graph:
        for cycle_2 in cycles_in_path_graph:
            if cycle_1 is not cycle_2 and any(path_in_cycle_1[0] in [path_in_cycle_2[0] for path_in_cycle_2 in cycle_2] for path_in_cycle_1 in cycle_1):
                cycle_graph.add_edge(tuple(cycle_1),tuple(cycle_2)) # nodes must be immutable objects, hence create tuple from list

    # cycles_in_path_graph not used hereafter
    del(cycles_in_path_graph)

    # look for all complete subgraphs in complement of cycle_graph
    # complete subgraphs in complement are those cycles that don't share any complexes as beginning of path
    cycle_graph = nx.complement(cycle_graph)
    
    # collect all complete subgraphs in complement of cycle graph
    complete_subgraphs = []
    # every node by itself is a complete subgraph
    for n in cycle_graph.nodes():
        if set(n) not in complete_subgraphs:
            #print [n]
            complete_subgraphs.append([n])
    # every edge in cycle graph is a complete subgraph (of cardinality 2)
    print cycle_graph.edges()
    for e in cycle_graph.edges():
        # TODO: below code a bit messy with switching between different types
        if set([list([cycle for edge_endpoint in e for cycle in edge_endpoint])]) not in complete_subgraphs:
            complete_subgraphs.append([list([cycle for edge_endpoint in e for cycle in edge_endpoint])])
    k = 3
    k_cliques = list(nx.k_clique_communities(cycle_graph,k))
    while len(k_cliques) > 0:
        for k_clique in k_cliques:
            if set(k_clique) not in complete_subgraphs:
                complete_subgraphs.append(tuple(k_clique))
        k = k+1
        k_cliques = list(nx.k_clique_communities(cycle_graph,k))

    # all-edges subgraphs
    all_edges = [sc[key]['edges'] for key in sc.keys()]
    all_edges_subgraphs = it.product(*all_edges)

    # mixed edges-cycles and all-cycles subgraphs
    mixed_edges_cycles_subgraphs = []
    for cycle_comb in complete_subgraphs:
        # IMPORTANT: cycle_comb is expected to be a list of tuples (each tuple being one cycle)
        if type(cycle_comb) != type(list()):
            #print "get_sensible_subgraphs: before raise: cycle_comb = "+str(cycle_comb)
            cycle_comb = list(cycle_comb)
            

        unused_edges_per_substance = [sc[key]['edges'] for key in sc.keys() if key not in [path[0] for cycle in cycle_comb for path in cycle]]
#        print 'cycle combination: '+str(cycle_comb)
#        print 'substance beginnings: '+str([path[0] for cycle in cycle_comb for path in cycle])
#        print [path for cycle in cycle_comb for path in cycle]
#        print 'unused edges: '+str(unused_edges_per_substance)

        # now generate all combinations of unused edges (combined in a way that each substance node is beginning of exactly one edge per combination) and concatenate each edges combination with current cycle_comb
        # note that if cycle_comb already uses all substances as beginning of exactly one path then unused_edges_per_substances is empty and the generated subgraph is an all-cycles subgraph
        for edges_combination in it.product(*unused_edges_per_substance):
            mixed_edges_cycles_subgraphs.append([edge for edge in edges_combination]+[path for cycle in cycle_comb for path in cycle])

    # return list of subgraphs
    return [tuple(sg) for sg in all_edges_subgraphs]+[tuple(sg) for sg in mixed_edges_cycles_subgraphs]


# first critical fragment of single-layer MAPK network
#frag = (('s9', 's8', 's2', 's1', 's5', 's4'),('w12', 'w8', 'w4', 'w1', 'w6', 'w10'))
#sc = {'s9': {'p_paths': [('s9', 'w12', 's1')], 'edges': [('s9', 'w12')], 'n_paths': []}, 's8': {'p_paths': [], 'edges': [('s8', 'w8')], 'n_paths': []}, 's2': {'p_paths': [('s2', 'w4', 's5')], 'edges': [('s2', 'w4'), ('s2', 'w1')], 'n_paths': [('s2', 'w4', 's4'), ('s2', 'w1', 's1')]}, 's1': {'p_paths': [], 'edges': [('s1', 'w1')], 'n_paths': [('s1', 'w1', 's2')]}, 's5': {'p_paths': [('s5', 'w6', 's2')], 'edges': [('s5', 'w6')], 'n_paths': []}, 's4': {'p_paths': [('s4', 'w4', 's5'), ('s4', 'w10', 's9')], 'edges': [('s4', 'w4'), ('s4', 'w10')], 'n_paths': [('s4', 'w4', 's2')]}}

# frag = (('s2', 's1', 's7', 's6', 's5', 's4'),('w1', 'w1', 'w10', 'w7', 'w5', 'w4'))
# sc = {'s2': {'p_paths': [('s2', 'w4', 's5')], 'edges': [('s2', 'w4'), ('s2', 'w1')], 'n_paths': [('s2', 'w4', 's4'), ('s2', 'w1', 's1')]}, 's1': {'p_paths': [], 'edges': [('s1', 'w1')], 'n_paths': [('s1', 'w1', 's2')]}, 's7': {'p_paths': [], 'edges': [('s7', 'w7'), ('s7', 'w10')], 'n_paths': [('s7', 'w7', 's6'), ('s7', 'w10', 's4')]}, 's6': {'p_paths': [], 'edges': [('s6', 'w7')], 'n_paths': [('s6', 'w7', 's7')]}, 's5': {'p_paths': [('s5', 'w5', 's2'), ('s5', 'w5', 's4')], 'edges': [('s5', 'w5')], 'n_paths': []}, 's4': {'p_paths': [('s4', 'w4', 's5')], 'edges': [('s4', 'w4'), ('s4', 'w10')], 'n_paths': [('s4', 'w4', 's2'), ('s4', 'w10', 's7')]}}

frag = (('s2', 's1', 's7', 's6', 's5', 's4'), ('w1', 'w1', 'w10', 'w7', 'w5', 'w4'))
sc = {'s2': {'p_paths': [('s2', 'w4', 's5')], 'edges': [('s2', 'w4'), ('s2', 'w1')], 'n_paths': [('s2', 'w4', 's4'), ('s2', 'w1', 's1')]}, 's1': {'p_paths': [], 'edges': [('s1', 'w1')], 'n_paths': [('s1', 'w1', 's2')]}, 's7': {'p_paths': [], 'edges': [('s7', 'w7'), ('s7', 'w10')], 'n_paths': [('s7', 'w7', 's6'), ('s7', 'w10', 's4')]}, 's6': {'p_paths': [], 'edges': [('s6', 'w7')], 'n_paths': [('s6', 'w7', 's7')]}, 's5': {'p_paths': [('s5', 'w5', 's2'), ('s5', 'w5', 's4')], 'edges': [('s5', 'w5')], 'n_paths': []}, 's4': {'p_paths': [('s4', 'w4', 's5')], 'edges': [('s4', 'w4'), ('s4', 'w10')], 'n_paths': [('s4', 'w4', 's2'), ('s4', 'w10', 's7')]}}

sg = get_sensible_subgraphs_paths_in_cycles_using_cycle_graph(sc)
valid_subgraphs = validate_subgraphs(sg, sc, frag)
fs = score_fragment(valid_subgraphs, sc, frag)
