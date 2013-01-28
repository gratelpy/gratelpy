#!/usr/bin/env python
"""check_data.py

Usage: check_data.py data_file"""

import cPickle as pickle
import sys
from collections import Counter
import errno
from graph import get_path_graph, get_valid_path_graph_cycles
from subgraphs import get_subgraph_motifs, score_subgraph
from fragments import score_fragment

import networkx as nx
import itertools as it

import math

# data format
frag_i = 0
sc_i = 1
sg_i = 2
ks_i = 3

# global variables available to interactive python sessions
filename = None
data = None

fragments = None
fragment_tags = None
fragment_counter = None
subgraph_components = None
subgraphs = None 
ks = None

try:
    filename = sys.argv[1]
except IndexError:
    print __doc__
    sys.exit(2)

try:
    print 'loading your data from',filename
    data = pickle.load(open(filename))
except IOError, e:
    if e.errno == errno.ENOENT: # no such file or directory
        raise Exception('unable to open your file \'',filename,'\'')
    else:
        raise Exception('unknown error')

fragment_tags = {s[0]: (tuple(sorted(s[frag_i][0])), tuple(sorted(s[frag_i][1]))) for s in data}
fragments = [fragment_tags[s[0]] for s in data] 

fragment_counter = Counter(fragments)
duplicate_fragments = [f for f in fragment_counter.keys() if fragment_counter[f]>1]
if len(duplicate_fragments) > 0:
    print str(len(duplicate_fragments)),'duplicate fragments detected'
    print 'on average, duplicate fragments occured',str(float(sum([fragment_counter[f] for f in duplicate_fragments]))/len(duplicate_fragments)),'times'

subgraph_components = {fragment_tags[s[frag_i]]: s[sc_i] for s in data}
subgraphs = {fragment_tags[s[frag_i]]: s[sg_i] for s in data}
edges_only_subgraphs = {}
for f in fragments:
    edges_only_subgraphs[f] = []
    for sg in subgraphs[f]:
        if all(len(el)==2 for el in sg):
            edges_only_subgraphs[f].append(sg)

ks = {fragment_tags[s[frag_i]]: s[ks_i] for s in data}

if len(duplicate_fragments) > 0:
    print str(len(fragments)),'fragments reported'
    print str(len(fragments)-len(duplicate_fragments)),'unique fragments'
    print str(len(subgraph_components)),'subgraph component dictionary entries'
    print str(len(subgraphs)),'subgraph dictionary entries'
    print str(len(ks)),'ks dictionary entries'

critical_fragments_reported = [f for f in fragments if ks[f]<0]
critical_fragments_unique = [f for f in ks.keys() if ks[f]<0]

print str(len(critical_fragments_reported)),'critical fragments reported'
print str(len(critical_fragments_unique)),'unique critical fragments'

# for every critical fragment, check:
# (i) that the multiplicity of each reaction and complex node is correct in every subgraph
# (ii) that every path in every subgraph is part of a cycle
# (iii) that every subgraph is scored correctly, hence the overall fragment score is correct too

# check multiplicity
for f in critical_fragments_unique:
    complex_multiplicity = Counter(f[0])
    reaction_multiplicty = Counter(f[1])
    for sg in subgraphs[f]:
        if complex_multiplicity != Counter([el[0] for el in sg]) or reaction_multiplicty != Counter([el[1] for el in sg]):
            print 'fragment',str(f)
            print 'subgraph with mismatched multiplicity:'
            print sg
            raise

# check that every path in every subgraph is in a cycle
subgraphs_tested = 0
subgraphs_with_all_paths_in_cycles = 0
no_subgraphs = sum([1 for frag in critical_fragments_unique for sg in subgraphs[frag]])
print 'checking if all paths in every subgraph are witihin a cycle ...'
for frag in critical_fragments_unique:
#frag = critical_fragments_unique[0]
    path_graph = get_path_graph(subgraph_components[frag])
    cycles = get_valid_path_graph_cycles(path_graph)
    sc = subgraph_components[frag]
    
    # create cycle graph
    cycles_in_path_graph_1 = get_valid_path_graph_cycles(path_graph)
    cycles_in_path_graph_2 = None
    cycle_graph = nx.Graph()
    cycles_in_path_graph = get_valid_path_graph_cycles(path_graph)
    for cycle in cycles_in_path_graph:
        cycle_graph.add_node(tuple(cycle))
    for cycle_1 in cycles_in_path_graph_1:
        cycles_in_path_graph_2 = get_valid_path_graph_cycles(path_graph)
        for cycle_2 in cycles_in_path_graph_2:
            if cycle_1 is not cycle_2 and any(path_in_cycle_1[0] in [path_in_cycle_2[0] for path_in_cycle_2 in cycle_2] for path_in_cycle_1 in cycle_1):
                cycle_graph.add_edge(tuple(cycle_1),tuple(cycle_2)) # nodes must be immutable objects, hence create tuple from list
    cycle_graph = nx.complement(cycle_graph)

    sg_motifs = get_subgraph_motifs(subgraph_components[frag])
    for s in subgraphs[frag]:
        subgraphs_tested += 1
        s_sg_motifs = sg_motifs[frozenset(s)]
        all_paths_in_cycles = all([s_el in [path for cycle in s_sg_motifs['cycles'] for path in cycle] for s_el in s if len(s_el)==3])
        if all_paths_in_cycles:
            subgraphs_with_all_paths_in_cycles += 1
        else:
            print 'fragment',str(frag)
            print 'subgraph',str(s)
            print 'subgraph motifs',str(s_sg_motifs)
            raise Exception('not all paths in cycles')

        if subgraphs_tested%100==0:
            print str(subgraphs_tested),'of',str(no_subgraphs),'tested -',100*float(subgraphs_tested)/float(no_subgraphs),'%'

print str(subgraphs_tested),'subgraphs tested and',str(subgraphs_with_all_paths_in_cycles),'of these have all paths in cycles'

# check that fragment score K_S is correct for every fragment
# we read out the K_S from data, generate K_S anew using the current library method, and use our new method using sg_motifs
fragments_checked = 0
k_s_ok = 0
for frag in critical_fragments_unique:
#frag = critical_fragments_unique[0]
    fragments_checked += 1
    sc = subgraph_components[frag]
    sg_motifs = get_subgraph_motifs(subgraph_components[frag])
    k_s_frag_this_method = 0
    k_s_frag_library_method = score_fragment(subgraphs[frag], sc, frag)
    k_s_frag_library_method_redone = 0

    for sg in subgraphs[frag]:
        k_g = float(1)
        t_g = len(sg_motifs[frozenset(sg)]['cycles'])
        k_g = k_g * math.pow(-1, t_g)
        for edge in sg_motifs[frozenset(sg)]['edges']:
            # NOTE: assuming that all stoichiometric coefficients == 1
            k_g = k_g * 1
        for cycle in sg_motifs[frozenset(sg)]['cycles']:
            for path in cycle:
                if path in sc[path[0]]['n_paths'] and path not in sc[path[0]]['p_paths']:
                    # NOTE: assuming that all stoichiometric coefficients == 1
                    k_g = k_g * (-1)
                elif path not in sc[path[0]]['n_paths'] and path in sc[path[0]]['p_paths']:
                    # NOTE: assuming that all stoichiometric coefficients == 1
                    k_g = k_g * 1
                else:
                    raise

        k_g_library_method = score_subgraph([tuple(sg), sc])
        if k_g != k_g_library_method:
            print 'k_g this method =',str(k_g)
            print 'k_g library method =',str(k_g_library_method)
            print 'subgraph =',str(sg)
            raise

        k_s_frag_this_method += k_g
        k_s_frag_library_method_redone += k_g_library_method

    different_K_S = [k_s_frag_this_method, k_s_frag_library_method[-1], k_s_frag_library_method_redone, ks[frag]]
    if not all(ks_1 == ks_2 for ks_1, ks_2 in zip(different_K_S, different_K_S)):
        print 'fragment',str(frag) 
        print 'K_S this method =',str(k_s_frag_this_method)
        print 'K_S library method =',str(k_s_frag_library_method[-1])
        print 'K_S library method redone =',str(k_s_frag_library_method_redone)
        print 'K_S in saved data =',str(ks[frag])
    else:
        k_s_ok += 1

    if fragments_checked % 10 == 0:
        print str(fragments_checked),'of',str(len(critical_fragments_unique)),'fragments checked -',str(100*float(fragments_checked)/float(len(critical_fragments_unique))),'%'

print str(fragments_checked),'fragments checked of total',str(len(critical_fragments_unique)),'fragments. k_s scores ok:',str(k_s_ok)
