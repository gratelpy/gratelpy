import networkx as nx
import numpy as np
from networkx.algorithms import bipartite
from collections import defaultdict

def get_valid_path_graph_cycles(path_graph):
    # path_graph: digraph with nodes = paths & edge n1 -> n2 => endpoint(n1) = beginning(n2)
    # returns all valid cycles in path_graph,
    # valid meaning that throughout each cycle every substance is beginning of exactly one path

    cycles = simple_cycles_unique_complexes(path_graph) # this returns a list!
    #valid_cycles = []
    for cycle in cycles:
        del cycle[-1]
        v1_beginnings = []
        cycle_valid = True
        for path in cycle:
            if path[0] not in v1_beginnings:
                v1_beginnings.append(path[0])
            else:
                cycle_valid = False
                break

        if cycle_valid:
            #valid_cycles.append(cycle)
            yield cycle
        else:
            continue

        # tidy up after each 'cycle'
        del(v1_beginnings)

    # tidy up
    #del(cycles)
        
    #return valid_cycles

    
def get_graph_edges(G):
    # G: networkx directed bipartite graph
    # returns dictionary of substance nodes and corresponding edges that they induce

    # get complex and reaction nodes
    complexes, reactions = get_bipartite_sets(G)

    # get complexes list
    complexes = list(complexes)
    
    # go over all complexes and collect all edges that they induce
    edges = {}
    for c in complexes:
        edges[c] = []
        for arc in G.edges():
            if arc[0] == c and G.node[arc[0]]['bipartite']==0:
                edges[c].append(arc)

    # return dictionary of edges
    return edges

def get_bipartite_sets(G):
    complexes, reactions = bipartite.sets(G)

    # some unexpected behaviour that networkx shows:
    # sometimes 'complexes' and 'reactions' get swapped around by bipartite.sets
    # this seems to happen for larger reaction networks

    if 'w1' not in complexes and 'w1' not in reactions:
        raise Exception('my hack to resolve this unexpected behavior shown by bipartite.sets assumes that reaction nodes are named \'w1\', \'w2\', ...')

    if 'w1' in complexes:
        return reactions, complexes
    else:
        return complexes, reactions

def get_lpa_alpha_beta(alpha, beta, slow_indices, complex_dict=None, constant_dict=None):
    # alpha: reactant stoichiometric coefficients of well-mixed system
    # beta: product stoichiometric coefficients of well-mixed system
    # slow_indices: indices of slow species in system; indices are row indices of alpha, beta

    # construct reverse dictionaries
    if complex_dict is not None:
        complex_dict_reverse = {v:k for k,v in complex_dict.iteritems()}
    else:
        complex_dict_reverse = None
    if constant_dict is not None:
        constant_dict_reverse = {v:k for k,v in constant_dict.iteritems()}
    else:
        constant_dict_reverse = None

    # construct LPA dictionaries
    complex_dict_lpa_reverse = {}
    constant_dict_lpa_reverse = {}

    if not type(alpha) is type(np.array([])) or not type(beta) is type(np.array([])):
        raise Exception('alpha and beta are expected to be NumPy arrays')

    # number of reactions = number of columns of stoich
    no_rxn = alpha.shape[1]
    if no_rxn != beta.shape[1]:
        raise Exception('alpha and beta need to have the same number of columns')
    # number of substance = number of rows of stoich
    no_sub = alpha.shape[0]
    if no_sub != beta.shape[0]:
        raise Exception('alpha and beta need to have the same number of rows')
    
    # fast species
    fast_indices = list(set(range(no_sub))-set(slow_indices))
    #fast_indices = [index if index not in slow_indices for index in range(len(no_sub))] 
    if len(fast_indices)==0:
        raise Exception('no fast species present in system: fast species are required to connect the global and local slow subgraphs')

    # how many slow and fast species do we have
    no_slow = len(slow_indices)
    no_fast = len(fast_indices)

    # row indices of slow (global and local) and fast variables in resulting LPA alpha and LPA beta matrices
    #slow_global_indices = range(0,no_slow-1,1)
    #fast_indices = range(no_slow, no_slow+no_fast-1,1)
    #slow_local_indices = range(no_slow+no_fast, 2*no_slow+no_fast-1,1)

    # fill in LPA dictionary of complexes with global variables
    if complex_dict_reverse is not None:
        complex_dict_lpa_reverse = {}
        complex_dict_lpa_key_i = 0
    else:
        complex_dict_lpa_reverse = None

    # fill in LPA dictionary of constants
    # TODO: check if following assumption correct: 
    # reactions 0 & 0+no_reactions have same rate constant
    constant_dict_lpa_reverse = {}
    for rxn_i in range(no_rxn):
        constant_dict_lpa_reverse[rxn_i] = constant_dict_reverse[rxn_i]
        constant_dict_lpa_reverse[rxn_i+no_rxn] = constant_dict_reverse[rxn_i]

    # lpa alpha matrix
    alpha_lpa = []
    # slow global variables
    for slow_i in slow_indices:
        alpha_lpa_row = list(alpha[slow_i])
        alpha_lpa_row = alpha_lpa_row + [0 for r in range(no_rxn)]
        alpha_lpa.append(alpha_lpa_row)

        if complex_dict_reverse is not None:
            complex_dict_lpa_reverse[complex_dict_lpa_key_i] = complex_dict_reverse[slow_i]+'_G'
            complex_dict_lpa_key_i += 1

        #print "global slow_i = "+str(slow_i)
        #print alpha_lpa_row
    # fast variables
    for fast_i in fast_indices:
        alpha_lpa_row = list(alpha[fast_i])
        alpha_lpa_row = alpha_lpa_row + list(alpha[fast_i])
        alpha_lpa.append(alpha_lpa_row)

        if complex_dict_reverse is not None:
            complex_dict_lpa_reverse[complex_dict_lpa_key_i] = complex_dict_reverse[fast_i]+'_G'
            complex_dict_lpa_key_i += 1

        #print "fast_i = "+str(fast_i)
        #print alpha_lpa_row
    # slow local variables
#    for slow_i in reversed(slow_indices): # not sure why I'm using reversed here ... ?
    for slow_i in slow_indices:
        alpha_lpa_row = [0 for r in range(no_rxn)]
        alpha_lpa_row = alpha_lpa_row + list(alpha[slow_i])
        alpha_lpa.append(alpha_lpa_row)

        if complex_dict_reverse is not None:
            complex_dict_lpa_reverse[complex_dict_lpa_key_i] = complex_dict_reverse[slow_i]+'_L'
            complex_dict_lpa_key_i += 1

        #print "local slow_i = "+str(slow_i)
        #print alpha_lpa_row

    # lpa beta matrix
    beta_lpa = []
    # global slow variables
    for slow_i in slow_indices:
        beta_lpa_row = list(beta[slow_i])
        beta_lpa_row = beta_lpa_row + [0 for r in range(no_rxn)]
        beta_lpa.append(beta_lpa_row)
    # fast variables
    for fast_i in fast_indices:
        # the part of beta_lpa pertaining to global slow variables is unaltered
        beta_lpa_row = list(beta[fast_i])
        # beta_lpa associated with local slow variables needs to be modified
        # nothing ever flows from slow local variables to fast variables, and
        # fast variables are never consumed in producing slow local variables
        beta_lpa_row_local = list(beta[fast_i])
        for rxn_i in range(no_rxn):
            if beta[fast_i, rxn_i] > 0 and any(coeff > 0 for coeff in list(alpha[slow_indices,rxn_i])):
                # if coeff > 0 and beta[fast_i, rxn_i] > 0 then mass flows from global slow variable to fast variable
                # this may not happen for the corresponding local slow variable
                beta_lpa_row_local[rxn_i] = 0
            if alpha[fast_i, rxn_i] > 0 and any(coeff > 0 for coeff in list(beta[slow_indices, rxn_i])):
                # if coeff > 0 and alpha[fast_i, rxn_i] > 0 then mass flows from fast variable to global slow variable
                # there is no net flow of mass from fast variable to local slow variables but the rate of formation
                # of local slow mass still is still proportional to the amount of the fast variable present, hence fast variable acts like a catalyst in these reactions
                beta_lpa_row_local[rxn_i] = alpha[fast_i, rxn_i] # mass of fast going into reaction, equals mass of fast coming out of it (=> like a catalyst)
        beta_lpa_row = beta_lpa_row + beta_lpa_row_local
        beta_lpa.append(beta_lpa_row)
    # local slow variables
    for slow_i in slow_indices:
        beta_lpa_row = [0 for r in range(no_rxn)]
        beta_lpa_row = beta_lpa_row + list(beta[slow_i])
        beta_lpa.append(beta_lpa_row)

    if complex_dict_reverse is not None:
        complex_dict_lpa = {v:k for k,v in complex_dict_lpa_reverse.iteritems()}
    else:
        complex_dict_lpa = None
    if constant_dict_reverse is not None:
        constant_dict_lpa = {v:k for k,v in constant_dict_lpa_reverse.iteritems()}
    else:
        constant_dict_lpa = None

    # NOTE: keys in constant_dict_lpa are not unique!!

    # convert to expected numpy.array type
    if constant_dict_lpa is not None and complex_dict_lpa is not None:
        return np.array(alpha_lpa), np.array(beta_lpa), complex_dict_lpa_reverse, constant_dict_lpa_reverse
    else:
        return np.array(alpha_lpa), np.array(beta_lpa)

def get_path_graph(sc):
    # sc: subgraph components

    # see #413 (in lab notebook) for description
    # generates path graph for fragment encoded in 'sc'
    
    # each negative path is converted to two "positive" paths
    # since negative path (A,k,B) != negative path (B,k,A) (i.e. these are two distinct paths)
    paths = set()
    for substance in sc.keys():
        for path in sc[substance]['p_paths']:
            paths.add((path[0],path[1],path[2],path[3]))
        for path in sc[substance]['n_paths']:
            paths.add((path[0],path[1],path[2],path[3]))
            paths.add((path[2],path[1],path[0],path[3]))
    paths = list(paths)

    # collect all starting and end points of paths
    path_starts = {key:[] for key in sc.keys()}#dict.fromkeys(sc.keys(),[])
 #   path_ends = dict.fromkeys(sc.keys(),[])

    for path in paths:
        path_starts[path[0]].append(path)
#        path_ends[path[2]].append(path)

    # paths in object 'paths' are nodes of path graph
    path_graph = nx.DiGraph()
    #path_graph.add_nodes_from(paths)

    # go through list of paths and add appropriate edges
    for path in paths:
        end_of_path = path[2]
        for target_path in path_starts[end_of_path]:
            path_graph.add_edge(path, target_path)
#            print 'added edge '+str(path)+' -> '+str(target_path)

    # tidy up
    #del(paths)
    #del(path_starts)

    # return
    return path_graph
        

#@not_implemented_for('directed')
def get_all_cliques(G):
    """Returns all cliques in an undirected graph.

    This method returns cliques of size (cardinality) k = 1, 2, 3, ..., maxDegree - 1.
    Where maxDegree is the maximal degree of any node in the graph.

    Keyword arguments
    -----------------
    G: undirected graph

    Returns
    -------
    generator of lists: generator of list for each clique.

    Notes
    -----
    To obtain a list of all cliques, use list(get_all_cliques(G)).

    Based on the algorithm published by Zhang et al. (2005) [1]_ and adapted to output all cliques discovered.
    
    This algorithm is not suitable for directed graphs.

    This algorithm ignores self-loops and parallel edges as
    clique is not conventionally defined with such edges.

    There are often many cliques in graphs. This algorithm however, hopefully, does not run out of memory
    since it only keeps candidate sublists in memory and continuously removes exhausted sublists.

    References
    ----------
    .. [1] Yun Zhang, Abu-Khzam, F.N., Baldwin, N.E., Chesler, E.J., Langston, M.A., Samatova, N.F., 
       Genome-Scale Computational Approaches to Memory-Intensive Applications in Systems Biology 
       Supercomputing, 2005. Proceedings of the ACM/IEEE SC 2005 Conference , vol., no., pp. 12, 12-18 Nov. 2005
       doi: 10.1109/SC.2005.29
       http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1559964&isnumber=33129
    """

    def greater_neighbors(G, a_node):
        """Helper method used in get_all_cliques"""
        nodes_sorted = sorted(G.nodes())
        a_node_index = nodes_sorted.index(a_node)
        
        neighbors_of_a_node = []

        for another_node_index, another_node in enumerate(nodes_sorted):
            if another_node_index > a_node_index and another_node in G.neighbors(a_node):
                neighbors_of_a_node.append(another_node)
    
        return tuple(neighbors_of_a_node)

    # sorted list of nodes in graph
    nodes_sorted = sorted(G.nodes())

    # starting point: build all 2-clique sublists
    clique_sublists = []
    for a_node_index, a_node in enumerate(nodes_sorted):
        clique_sublist = {}
        # sublist base, sb
        clique_sublist['sb'] = [tuple(a_node)]
        # common neighbors, cn
        clique_sublist['cn'] = greater_neighbors(G, a_node)
        clique_sublists.append(clique_sublist)

    # output cliques of size k = 1
    for node in nodes_sorted:
        yield [node]

    # output cliques of size k >= 2
    while clique_sublists:
        a_sublist = clique_sublists.pop(0)
        for node_added in a_sublist['cn']:
            neighbors_of_node_added = greater_neighbors(G, node_added)

            current_sublist_base = [] + a_sublist['sb'] + [tuple(node_added)]
            current_sublist_cn = tuple(sorted(set(neighbors_of_node_added).intersection(a_sublist['cn'])))

            #print 'clique: '+str(current_sublist_base)
            yield [node for node in current_sublist_base]

            for node in current_sublist_cn:
                new_sublist_base = [] + current_sublist_base 
                new_sublist_base.append(tuple(node))
                #print 'current_sublist_based =',str(current_sublist_base)
                #print 'new_sublist_base =',str(new_sublist_base)
                new_sublist_cn = tuple(sorted(set(current_sublist_cn).intersection(greater_neighbors(G, node))))
    
                if len(new_sublist_cn) == 0:
                    #print 'clique: '+str(new_sublist_base)
                    yield [n for n in new_sublist_base]
                elif len(new_sublist_cn) == 1:
                    #print 'clique: '+str(new_sublist_base)
                    #print 'new_sublist_base + list(new_sublist_cn):',new_sublist_base+list(new_sublist_cn)
                    yield [n for n in new_sublist_base]
                    #print 'clique: '+str(new_sublist_base+new_sublist_cn)
                    
                    yield [n for n in new_sublist_base + list(new_sublist_cn)]
                else:
                    #print 'candidate sublist: '+str([new_sublist_base, new_sublist_cn])
                    clique_sublists.append({'sb': new_sublist_base, 'cn': new_sublist_cn})

def simple_cycles_unique_complexes(G):
    """Find simple cycles (elementary circuits) of a directed graph and do not revisit complexes.

    This method is a slightly modified clone of the method networkx.algorithms.cycles.simple_cycles provided by NetworkX.
    The original method can be found at https://github.com/networkx/networkx/blob/master/networkx/algorithms/cycles.py#L98.
    All credit goes to the authors of the original method [2].

    An simple cycle, or elementary circuit, is a closed path where no
    node appears twice, except that the first and last node are the same.
    Two elementary circuits are distinct if they are not cyclic permutations
    of each other.

    Parameters
    ----------
    G : NetworkX DiGraph
       A directed graph

    Returns
    -------
    A list of circuits, where each circuit is a list of nodes, with the first
    and last node being the same.

    Example:
    >>> G = nx.DiGraph([(0, 0), (0, 1), (0, 2), (1, 2), (2, 0), (2, 1), (2, 2)])
    >>> nx.simple_cycles(G)
    [[0, 0], [0, 1, 2, 0], [0, 2, 0], [1, 2, 1], [2, 2]]

    See Also
    --------
    cycle_basis (for undirected graphs)

    Notes
    -----
    The implementation follows pp. 79-80 in [1]_.

    The time complexity is O((n+e)(c+1)) for n nodes, e edges and c
    elementary circuits.

    References
    ----------
    .. [1] Finding all the elementary circuits of a directed graph.
       D. B. Johnson, SIAM Journal on Computing 4, no. 1, 77-84, 1975.
       http://dx.doi.org/10.1137/0204007
    .. [2] Exploring network structure, dynamics, and function using {NetworkX}.
       Aric A. Hagberg and Daniel A. Schult and Pieter J. Swart. Proceedings of the 7th Python in Science Conference (SciPy2008).
       http://math.lanl.gov/~hagberg/Papers/hagberg-2008-exploring.pdf

    See Also
    --------
    cycle_basis
    """
    # Jon Olav Vik, 2010-08-09
    def _unblock(thisnode):
        """Recursively unblock and remove nodes from B[thisnode]."""
        if blocked[thisnode]:
            blocked[thisnode] = False
            while B[thisnode]:
                _unblock(B[thisnode].pop())

    def circuit(thisnode, startnode, component):
        closed = False # set to True if elementary path is closed
        path.append(thisnode)
        blocked[thisnode] = True
        for nextnode in component[thisnode]: # direct successors of thisnode
            if nextnode == startnode:
                path_to_be_appended = path + [startnode]

                if len(set([p[0] for p in path_to_be_appended])) == len(path_to_be_appended[0:-1]): 
                    # length of path minus last element since that one is the first one again
                    # only add those paths to result list that have equal numbers of nodes (paths in Mincheva sense) and complexes visited
                    # number of complexes == number of paths, hence each complex visited exactly once
                    result.append(path_to_be_appended)
                closed = True
            elif not blocked[nextnode]:
                if circuit(nextnode, startnode, component):
                    closed = True
        if closed:
            _unblock(thisnode)
        else:
            for nextnode in component[thisnode]:
                if thisnode not in B[nextnode]: # TODO: use set for speedup?
                    B[nextnode].append(thisnode)
        path.pop() # remove thisnode from path
        return closed

    path = [] # stack of nodes in current path
    blocked = defaultdict(bool) # vertex: blocked from search?
    B = defaultdict(list) # graph portions that yield no elementary circuit
    result = [] # list to accumulate the circuits found
    # Johnson's algorithm requires some ordering of the nodes.
    # They might not be sortable so we assign an arbitrary ordering.
    ordering=dict(zip(G,range(len(G))))
    for s in ordering:
        # Build the subgraph induced by s and following nodes in the ordering
        subgraph = G.subgraph(node for node in G
                              if ordering[node] >= ordering[s])
        # Find the strongly connected component in the subgraph
        # that contains the least node according to the ordering
        strongcomp = nx.strongly_connected_components(subgraph)
        mincomp=min(strongcomp,
                    key=lambda nodes: min(ordering[n] for n in nodes))
        component = G.subgraph(mincomp)
        if component:
            # smallest node in the component according to the ordering
            startnode = min(component,key=ordering.__getitem__)
            for node in component:
                blocked[node] = False
                B[node][:] = []
            dummy=circuit(startnode, startnode, component)

    return result
