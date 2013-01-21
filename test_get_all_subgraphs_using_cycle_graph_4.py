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
    # pl.figure()
    # nx.draw(path_graph)

    #cycles_in_path_graph = nx.simple_cycles(path_graph)
    cycles_in_path_graph = list(get_valid_path_graph_cycles(path_graph))

    # create cycle graph
    # nodes are cycles
    # edge(c1,c2) exists iff c1 and c2 have at least one substance as beginning of path in common
    cycle_graph = nx.Graph()
    for cycle in cycles_in_path_graph:
        cycle_graph.add_node(tuple(cycle))

    for cycle_1 in cycles_in_path_graph:
        for cycle_2 in cycles_in_path_graph:
            if cycle_1 is not cycle_2 and any(path_in_cycle_1[0] in [path_in_cycle_2[0] for path_in_cycle_2 in cycle_2] for path_in_cycle_1 in cycle_1):
                cycle_graph.add_edge(tuple(cycle_1),tuple(cycle_2)) # nodes must be immutable objects, hence create tuple from list

    # pl.figure()
    # nx.draw(cycle_graph)

    # cycles_in_path_graph not used hereafter
    del(cycles_in_path_graph)

    # look for all complete subgraphs in complement of cycle_graph
    # complete subgraphs in complement are those cycles that don't share any complexes as beginning of path
    cycle_graph = nx.complement(cycle_graph)
    # pl.figure()
    # nx.draw(cycle_graph)
    # pl.show()
    
    # collect all complete subgraphs in complement of cycle graph
    complete_subgraphs = []
    # every node by itself is a complete subgraph
    for n in cycle_graph.nodes():
        if set(n) not in complete_subgraphs:
            #print [n]
            complete_subgraphs.append([n])
    # every edge in cycle graph is a complete subgraph (of cardinality 2)
#    print cycle_graph.edges()
    for e in cycle_graph.edges():
        # TODO: below code a bit messy with switching between different types
        if set(tuple([cycle for edge_endpoint in e for cycle in edge_endpoint])) not in complete_subgraphs:
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


frag = (('s9', 's8', 's3', 's12', 's1', 's10', 's6', 's5', 's4'),('w16', 'w10', 'w2', 'w18', 'w1', 'w12', 'w13', 'w6', 'w4'))
sc = {'s1': {'edges': [('s1', 'w1')],
        'n_paths': [],
        'p_paths': [('s1', 'w1', 's3')]},
 's10': {'edges': [('s10', 'w12')],
         'n_paths': [],
         'p_paths': [('s10', 'w12', 's9'), ('s10', 'w12', 's6')]},
 's12': {'edges': [('s12', 'w18')],
         'n_paths': [],
         'p_paths': [('s12', 'w18', 's1'), ('s12', 'w18', 's9')]},
 's3': {'edges': [('s3', 'w2')],
        'n_paths': [],
        'p_paths': [('s3', 'w2', 's1')]},
 's4': {'edges': [('s4', 'w16'), ('s4', 'w4')],
        'n_paths': [('s4', 'w16', 's9')],
        'p_paths': [('s4', 'w16', 's12'), ('s4', 'w4', 's5')]},
 's5': {'edges': [('s5', 'w6')],
        'n_paths': [],
        'p_paths': [('s5', 'w6', 's6')]},
 's6': {'edges': [('s6', 'w13')],
        'n_paths': [('s6', 'w13', 's9')],
        'p_paths': []},
 's8': {'edges': [('s8', 'w10')],
        'n_paths': [('s8', 'w10', 's9')],
        'p_paths': [('s8', 'w10', 's10')]},
 's9': {'edges': [('s9', 'w16'), ('s9', 'w13'), ('s9', 'w10')],
        'n_paths': [('s9', 'w16', 's4'),
                    ('s9', 'w13', 's6'),
                    ('s9', 'w10', 's8')],
        'p_paths': [('s9', 'w16', 's12'), ('s9', 'w10', 's10')]}}

missing_sgs = [(('s9', 'w16', 's12'), ('s8', 'w10'), ('s3', 'w2', 's1'), ('s12', 'w18', 's9'), ('s1', 'w1', 's3'), ('s10', 'w12'), ('s6', 'w13'), ('s5', 'w6'), ('s4', 'w4')),
(('s9', 'w16'), ('s8', 'w10'), ('s3', 'w2', 's1'), ('s12', 'w18'), ('s1', 'w1', 's3'), ('s10', 'w12'), ('s6', 'w13'), ('s5', 'w6'), ('s4', 'w4')),
(('s9', 'w16', 's4'), ('s8', 'w10'), ('s3', 'w2', 's1'), ('s12', 'w18'), ('s1', 'w1', 's3'), ('s10', 'w12'), ('s6', 'w13', 's9'), ('s5', 'w6', 's6'), ('s4', 'w4', 's5'))]
missing_sgs_set = [set(s) for s in missing_sgs]

sg = get_sensible_subgraphs_paths_in_cycles_using_cycle_graph(sc)
sg_set = [set(s) for s in sg]
valid_subgraphs = validate_subgraphs(sg, sc, frag)
fs = score_fragment(valid_subgraphs, sc, frag)
