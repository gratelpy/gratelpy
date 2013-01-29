import functools
import itertools as it
import networkx as nx
from multiprocessing import Pool
from networkx.algorithms import bipartite
from subgraphs import score_subgraph
from stoich import get_substance_adjacency
from graph import get_graph_edges, get_bipartite_sets
from decorators import deprecated

def pretty_print(valid_subgraphs):
    # pretty printing for fragments:
    # rearrange substance and reaction nodes in fragment tuples such that
    # substance and reaction nodes are paired up as in the fragment's 'all-edges' subgraph
    frag_print = None
    pretty_printed = False
    for sg in valid_subgraphs:
        if all(len(sg_el)==2 for sg_el in sg):
            if pretty_printed:
                # there should only be one all-edges subgraph per fragment
                # waltherg Jan 29: this is not true
                pass
            else:
                pretty_printed = True

            frag_print_sub = []
            frag_print_rxn = []
            for sg_el in sg:
                frag_print_sub.append(sg_el[0])
                frag_print_rxn.append(sg_el[1])
            frag_print = (tuple(frag_print_sub), tuple(frag_print_rxn))

    return frag_print

def validate_fragments(G, stoich_rank, f):
    indgr = G.subgraph([item for sublist in f for item in sublist])
    edges = []
    non_edges = [] # those arcs that go from reaction node to complex node
    for arc in indgr.edges():
        if indgr.node[arc[0]]['bipartite']==0:
            edges.append(arc)
        else:
            non_edges.append(arc)
                
    v1_nodes = []
    for e in edges:
        if not e[0] in v1_nodes:
            v1_nodes.append(e[0])

    f_valid = False
    if len(v1_nodes) == stoich_rank:
        f_valid = True # assume f is valid since v1 count is fine
        for a in non_edges: 
            # all we need to do here is catch those non-edge arcs that aren't part of a positive path - no need for fancy cycle detection at this point
            if a[1] not in v1_nodes:
                f_valid = False
            
    if f_valid:
        return f
    else:
        return None

def get_unique_fragments(some_frags):
    # some_frags: list of fragments
    # ideally, pass a list of validated fragments
    
    # first, sort current fragments internally
    sorted_frags = []
    for f in some_frags:
        sorted_frags.append([sorted(f[0]), sorted(f[1])])

    # go through list of sorted fragments and append to list of unique fragments if not present yet
    unique_frags = []
    for f in sorted_frags:
        if f not in unique_frags:
            unique_frags.append(f)

    # code following current call to get_unique_fragments may expect frags to be a list of tuples of tuples
    for f_i in range(len(unique_frags)):
        unique_frags[f_i] = (tuple(unique_frags[f_i][0]), tuple(unique_frags[f_i][1]))

    return unique_frags

@deprecated('get_sensible_fragments')
def get_valid_fragments(G, stoich_rank):
    #reactions, complexes = bipartite.sets(G)
    complexes, reactions = bipartite.sets(G)

    complexes = list(complexes)
    reactions = list(reactions)

    if 'w1' not in complexes and 'w1' not in reactions:
        raise Exception('my hack to resolve this unexpected behavior shown by bipartite.sets assumes that reaction nodes are named \'w1\', \'w2\', ...')
    
    if 'w1' in complexes:
        complexes, reactions = reactions, complexes

    if not ('w1' in reactions and 's1' in complexes):
        raise Exception('Something went wrong generating the lists of complexes of reactions.')

    complex_perms = list(it.combinations(complexes,stoich_rank))
    reaction_perms = list(it.combinations_with_replacement(reactions,stoich_rank))
    fragments = list(it.product(complex_perms, reaction_perms))

    valid_fragments = []
    
    pool = Pool()
    chunksize = 100

    myval = functools.partial(validate_fragments, G, stoich_rank)
    
    fragment_list = pool.imap(myval, fragments, chunksize)
    valid_fragments = [f for f in fragment_list if f is not None]

    return get_unique_fragments(valid_fragments)

def score_fragment(vs, sc, frag):
    # vs: a list of validated subgraphs, one list per fragment
    # sc: dictionary of corresponding subgraph components
    # frags: validated fragment
    # vs = args[0]
    # sc = args[1]
    # frag = args[2]

    K_S = 0
    for sub_graph in vs:
        # get K_g of 'sub_graph'
        K_g = score_subgraph([sub_graph, sc])
        # update K_S
        K_S = K_S + K_g

    return frag, sc, K_S

def get_all_substance_combinations_with_cycles(alpha, beta):
    # alpha, beta are stoichiometry matrices as used throughout code

    # number of reactions = number of columns of alpha
    no_rxn = alpha.shape[1]
    # number of substance = number of rows of alpha
    no_sub = alpha.shape[0]

    # check
    if no_rxn != beta.shape[1] or no_sub != beta.shape[0]:
        raise

    # get substance adjacency matrix
    subs_adj = get_substance_adjacency(alpha, beta)

    # get directed substance graph
    subs_G = nx.from_numpy_matrix(subs_adj, create_using=nx.DiGraph())

    # get cycles in substance graph
    subs_cycles  = nx.simple_cycles(subs_G)
    # remove substance index repetitions
    for c_i in range(len(subs_cycles)):
        subs_cycles[c_i] = list(set(subs_cycles[c_i]))
    
def get_sensible_fragments(G, stoich_rank):
    # G: the graph
    # stoich_rank

    # edges: dictionary of substance nodes with list of edges each induces
    edges = get_graph_edges(G)
#    print 'len(list(edges)) = '+str(len(list(edges)))
#    print edges
    # get list of complex nodes
    complexes, reactions = get_bipartite_sets(G)
    complexes = list(complexes)

    # get all complex combinations
    complex_perms = it.combinations(complexes,stoich_rank)
#    print 'len(complex_perms) = '+str(len(complex_perms))

    # for each complex combination, create all sensible fragments
    # realize that every fragment must contain a unique all-edges subgraph (unique b/c each substance node can only be beginning of exactly one edge in an all-edges subgraph)
    fragments = []
    for perm in complex_perms:
        edges_in_perm = []
        for c in perm:
            edges_in_perm.append(edges[c])

        edge_combinations = it.product(*edges_in_perm)
        #print 'len(list(edge_combinations)) = '+str(len(list(edge_combinations)))
        for edge_comb in edge_combinations:
            reaction_nodes = [e[1] for e in edge_comb]
            fragments.append(tuple([tuple(perm),tuple(reaction_nodes)]))
                             
    # return
    return fragments

def has_large_scc_in_substance_graph(f, sc):
    # f: fragment
    # sc: subgraph components in fragment 'f'

    # tells you if fragment 'f' has strongly-connected component of length greater than one in corresponding directed graph induced by paths in sc
    # this should tell you if fragment f contains subgraphs that have at least one cycle of length greater than one
    
    substance_edges_in_paths = set()
    for key in sc.keys():
        for path in sc[key]['n_paths']+sc[key]['p_paths']:
            # 
            curr_edge = (path[0],path[2])
#            if curr_edge not in substance_edges_in_paths:
            substance_edges_in_paths.add(curr_edge)
    substance_edges_in_paths = list(substance_edges_in_paths)

    substance_graph = nx.DiGraph()
    substance_graph.add_edges_from(substance_edges_in_paths)

    strong_comp = nx.strongly_connected_components(substance_graph)

    if max([len(el) for el in strong_comp]) > 1:
        # strongly-connected component of size > 1 detected in directed substance graph
        return True
    else:
        # all strongly-connected components are individual substance nodes
        return False
def get_substance_strongly_connected_components(sc):
    substance_edges_in_paths = {}
    for key in sc.keys():
        for path in sc[key]['n_paths']:
            curr_edge = (path[0],path[2])
            if curr_edge not in substance_edges_in_paths.keys():
                substance_edges_in_paths[curr_edge] = {'n_paths': [], 'p_paths': []}
                substance_edges_in_paths[curr_edge]['n_paths'].append(path)
            else:
                substance_edges_in_paths[curr_edge]['n_paths'].append(path)
        for path in sc[key]['p_paths']:
            curr_edge = (path[0],path[2])
            if curr_edge not in substance_edges_in_paths.keys():
                substance_edges_in_paths[curr_edge] = {'n_paths': [], 'p_paths': []}
                substance_edges_in_paths[curr_edge]['p_paths'].append(path)
            else:
                substance_edges_in_paths[curr_edge]['p_paths'].append(path)
        # for path in sc[key]['n_paths']+sc[key]['p_paths']:
        #     curr_edge = (path[0],path[2])
        #     if curr_edge not in substance_edges_in_paths.keys():
        #         substance_edges_in_paths[curr_edge] = [path]
        #     else:
        #         substance_edges_in_paths[curr_edge].append(path)
#    substance_edges_in_paths = list(substance_edges_in_paths)

    substance_graph = nx.DiGraph()
    substance_graph.add_edges_from([(edge[0],edge[1],{'n_paths': substance_edges_in_paths[edge]['n_paths'],'p_paths': substance_edges_in_paths[edge]['p_paths']}) for edge in substance_edges_in_paths.keys()])

    strong_comp = nx.strongly_connected_components(substance_graph)

    # if a substance node in sc is not beginning of a path, it is not considered in the above analysis of strongly-connected components
    # however, for how 'strong_comp' will be used in remainder of code, those substance nodes that are beginning of edges only in 'sc' will also be considered scc of length 1
    substance_nodes_used = [node for scc in strong_comp for node in scc]
    for key in sc.keys():
        if key not in substance_nodes_used:
            strong_comp.append([key])
            substance_graph.add_node(key)

    return strong_comp, substance_graph

def get_sscc(sc):
    return get_substance_strongly_connected_components(sc)

def get_paths_in_scc(f, sc):
    # f: fragment
    # sc: subgraph components in fragment 'f'

    # returns list of those paths that are part of at least one 'large' strongly-connected component
    # i.e. only those paths are returned that are part of an scc of length > 1

    substance_edges_in_paths = set()
    for key in sc.keys():
        for path in sc[key]['n_paths']+sc[key]['p_paths']:
            curr_edge = (path[0],path[2])
#            if curr_edge not in substance_edges_in_paths:
            substance_edges_in_paths.add(curr_edge)
    substance_edges_in_paths = list(substance_edges_in_paths)

    substance_graph = nx.DiGraph()
    substance_graph.add_edges_from(substance_edges_in_paths)

    strong_comp = nx.strongly_connected_components(substance_graph)

    # collect paths for return
    n_paths = set()
    p_paths = set()
    count = 0
    for scc in strong_comp:
        if len(scc) > 1:
            for substance in scc:
                for path in sc[substance]['n_paths']:
                    n_paths.add(path)
                    count = count + 1
                for path in sc[substance]['p_paths']:
                    p_paths.add(path)
                    count = count + 1
    # return
    return {'n_paths': list(n_paths), 'p_paths': list(p_paths), 'count': count}


def get_sensible_subgraph_components(G, f):
    return get_sensible_sc(G, f)

def get_sensible_sc(G, f):
    sc = get_subgraph_components(G, f)
    sc = sc_remove_unreasonable_paths(f, sc)

    return sc

def sc_remove_unreasonable_paths(f, sc):
    paths_in_scc = get_paths_in_scc(f, sc)

    # amend sc to remove those paths that neither form a cycle by themselves nor are part of a strongly-connected component of length > 1
    for substance in sc:
        for path in sc[substance]['n_paths']:
            if path not in paths_in_scc['n_paths']:
                # path is not part of an scc of length > 1
                # remove path from subgraph components only if it does not form a cycle by itself
                if path[0] != path[2]:
                    sc[substance]['n_paths'].remove(path)
    #                print 'removed '+str(path)
        for path in sc[substance]['p_paths']:
            if path not in paths_in_scc['p_paths']:
                # path is not part of an scc of length > 1
                # remove path from subgraph components only if it does not form a cycle by itself
                if path[0] != path[2]:
                    sc[substance]['p_paths'].remove(path)
        #            print 'removed '+str(path)

    return sc

def get_edges_of_sensible_fragments(G, stoich_rank):
    # G: the graph
    # stoich_rank

    # edges: dictionary of substance nodes with list of edges each induces
    edges = get_graph_edges(G)
#    print 'len(list(edges)) = '+str(len(list(edges)))
#    print edges
    # get list of complex nodes
    complexes, reactions = get_bipartite_sets(G)
    complexes = list(complexes)

    # get all complex combinations
#    complex_perms = list(it.combinations(complexes,stoich_rank))
    complex_perms = it.combinations(complexes,stoich_rank)
#    print 'len(complex_perms) = '+str(len(complex_perms))

    # for each complex combination, create all sensible fragments
    edges_in_fragments = []
    for perm in complex_perms:
        edges_in_perm = []
        for c in perm:
            edges_in_perm.append(edges[c])

        edge_combinations = it.product(*edges_in_perm)
        #print 'len(list(edge_combinations)) = '+str(len(list(edge_combinations)))
        for edge_comb in edge_combinations:
#            reaction_nodes = [e[1] for e in edge_comb]
#            fragments.append(tuple([tuple(perm),tuple(reaction_nodes)]))
            edges_in_fragments.append(edge_comb)
                             
    # return
    return edges_in_fragments
