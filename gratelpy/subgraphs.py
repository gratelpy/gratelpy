import itertools as it
import networkx as nx
import math
from graph import get_path_graph, get_valid_path_graph_cycles, get_all_cliques

def get_subgraph_components(G, f):
    species = [spec for spec in f[0]]
    edges = [(e0, e1) for e0, e1 in zip(f[0],f[1])]
    if len(edges) != len(f[0]):
        print 'get_subgraph_components: fragment',fragment
        print 'get_subgraph_components: edges found',edges
        raise('More edges discovered than number of edges indicated by fragment order')

    subgraph_components = {}
    for n in f[0]:
        subgraph_components[n]={'edges' : [], 'p_paths':[], 'n_paths':[]}
    
    # add edges starting in each complex node
    for e in edges:
        if e not in subgraph_components[e[0]]['edges']:
            subgraph_components[e[0]]['edges'].append(e)

    # add all positive paths
    # NOTE: it doesn't matter if a positive path is part of a cycle. here you want to add ALL positive paths
    # whether or not every positive path contributes to a sensible subgraph is decided later in validate_subgraph
    # for e in cycle_edges:
    #     for n in cycle_non_edges:
    for e in edges:
        for sp in species:
            if (e[1], sp) in G.edges():
                p_edge = (e[0],e[1],sp,'p')
                if p_edge not in subgraph_components[e[0]]['p_paths']:
                    subgraph_components[e[0]]['p_paths'].append(p_edge)

    # add all negative paths
    # NOTE: it doesn't matter if a negative path is part of a cycle. here you want to add ALL negative paths
    # whether or not every negative path contributes to a sensible subgraph is decided later in validate_subgraph
    # for e1 in cycle_edges:
    #     for e2 in cycle_edges:
    for e1 in edges:
        for sp in species:
            if e1[0] != sp and (sp, e1[1]) in G.edges():
                #print (e1[0],e1[1],e2[0])
                n_edge = (e1[0],e1[1],sp,'n')
                if n_edge not in subgraph_components[e1[0]]['n_paths']:
                    subgraph_components[e1[0]]['n_paths'].append(n_edge)

    #print 'get_subgraph_components: len(subgraph_components) = '+str(len(subgraph_components))

    return subgraph_components

def is_cycle(paths):
    # paths: tuple of paths, i.e. ordered set of paths

    paths = list(paths) # resulting paths list is different from paths tuple, i.e. popping elements won't alter global 'paths' object

    try:
        path = paths.pop()
    except IndexError:
        return False

    first_v1 = path[0]
    last_v1 = path[2]

    v2_seen = [path[1]]

    while paths:
        path = paths.pop()
        if last_v1 == path[0]:
            # new path connects up with old path
            last_v1 = path[2]
        else:
            # paths don't connect
            return False

        # positive path consisting of two negative edges (s1, w1, s2) + (s2, w1, s1) is valid
        # hence the following check was incorrect:

        # if not path[1] in v2_seen:
        #     v2_seen.append(path[1])
        # else:
        #     # we don't want to traverse the same edge multiple times
        #     # if there is a reaction node we already traversed we are about to traverse that reaction node again with the current path
        #     return False

    if last_v1 == first_v1:
        return True
    else:
        return False

def get_all_subgraphs(sc, stoich_rank):
    # sc: subgraph components dictionary
    gr_el = []
    for key in sc:
        for c_type in sc[key]:
            for el in sc[key][c_type]:
                if el not in gr_el:
                    gr_el.append(el)
#    return sc, list(it.combinations(gr_el, stoich_rank))

    all_subgraphs = it.combinations(gr_el, stoich_rank)
    
    # DEBUG
#    print 'get_all_subgraphs: number of subgraphs generated = '+str(len(all_subgraphs))

    return all_subgraphs

def get_sensible_subgraphs(sc):
    #path_graph = get_path_graph(sc)

    # if len(path_graph.edges()) <= 50:
    #     return get_sensible_subgraphs_paths_in_cycles(sc)
    # else:
    #     return get_sensible_subgraphs_unique_beginnings(sc)
    
    #return get_sensible_subgraphs_paths_in_cycles(sc)

    return get_sensible_subgraphs_paths_in_cycles_using_cycle_graph(sc)


def get_subgraph_motifs(sc):
    # defintion of cliques / complete subgraphs
    # http://en.wikipedia.org/wiki/Clique_(graph_theory)
    # NOTE: in this method, the word 'subgraph' or phrase 'complete subgraph' most likely refers to the above graph-theoretic concept of cliques
    # and not the concept used by Mincheva et al.

    # get path graph
    path_graph = get_path_graph(sc)

    #cycles_in_path_graph = nx.simple_cycles(path_graph)
    cycles_in_path_graph_1 = get_valid_path_graph_cycles(path_graph)
    
    cycles_in_path_graph_2 = None

    # create cycle graph
    # nodes are cycles
    # edge(c1,c2) exists iff c1 and c2 have at least one substance as beginning of path in common
    cycle_graph = nx.Graph()

    # add all cycles (nodes) to cycle graph
    cycles_in_path_graph = get_valid_path_graph_cycles(path_graph)
    for cycle in cycles_in_path_graph:
        cycle_graph.add_node(tuple(cycle))

    # NOTE: code below only makes certain that nodes (cycles) connected to at least one other node (cycle) are added to the graph
    # nodes (cycles) that are entirely disconnected from remainder of graph would be left out erroneously! Hence, we now add every cycle to the graph above before connecting all cycles.
    for cycle_1 in cycles_in_path_graph_1:
        # get_valid_path_graph_cycles returns a generator so receiving object of this method is 'empty' after each full loop across the list of valid path graph cycles
        # we could also turn the returned generator into a list but here we choose a loss in performance (keeping the generator nature of returned object) over increased memory usage (conversion to list type)
        cycles_in_path_graph_2 = get_valid_path_graph_cycles(path_graph)
        for cycle_2 in cycles_in_path_graph_2:
            if cycle_1 is not cycle_2 and any(path_in_cycle_1[0] in [path_in_cycle_2[0] for path_in_cycle_2 in cycle_2] for path_in_cycle_1 in cycle_1):
                cycle_graph.add_edge(tuple(cycle_1),tuple(cycle_2)) # nodes must be immutable objects, hence create tuple from list


    # look for all complete subgraphs in complement of cycle_graph
    # complete subgraphs in complement are those cycles that don't share any complexes as beginning of path
    cycle_graph = nx.complement(cycle_graph)

    # subgraph motifs dictionary
    sg_motifs = {}

    # all-edges subgraphs
    all_edges = [sc[key]['edges'] for key in sc.keys()]
    all_edges_subgraphs = it.product(*all_edges)

    #for cycle_comb in complete_subgraphs:
    for cycle_comb in get_all_cliques(cycle_graph):
        # IMPORTANT: cycle_comb is expected to be a list of tuples (each tuple being one cycle)
        if type(cycle_comb) != type(list()):
            #print "get_sensible_subgraphs: before raise: cycle_comb = "+str(cycle_comb)
            cycle_comb = list(cycle_comb)
            

        unused_edges_per_substance = [sc[key]['edges'] for key in sc.keys() if key not in [path[0] for cycle in cycle_comb for path in cycle]]
        for edges_combination in it.product(*unused_edges_per_substance):
            curr_key = frozenset([edge for edge in edges_combination]+[path for cycle in cycle_comb for path in cycle])
            if curr_key not in sg_motifs.keys():
                sg_motifs[curr_key] = {}
                sg_motifs[curr_key]['edges'] = [edge for edge in edges_combination]
                sg_motifs[curr_key]['cycles'] = [cycle for cycle in cycle_comb]
            else:
                raise

    # add all-edges subgraphs
    for sg in all_edges_subgraphs:
        curr_key = frozenset(sg)
        if curr_key not in sg_motifs.keys():
            sg_motifs[curr_key] = {}
            sg_motifs[curr_key]['edges'] = [edge for edge in sg]
            sg_motifs[curr_key]['cycles'] = []
        else:
            raise

    # return subgraph motifs dictionary
    return sg_motifs

def get_sensible_subgraphs_paths_in_cycles_using_cycle_graph(sc):
    # code from test_get_all_subgraphs_using_cycle_graph_2.py

    # defintion of cliques / complete subgraphs
    # http://en.wikipedia.org/wiki/Clique_(graph_theory)
    # NOTE: in this method, the word 'subgraph' or phrase 'complete subgraph' most likely refers to the above graph-theoretic concept of cliques
    # and not the concept used by Mincheva et al.

    # get path graph
    path_graph = get_path_graph(sc)

    #cycles_in_path_graph = nx.simple_cycles(path_graph)
    cycles_in_path_graph_1 = get_valid_path_graph_cycles(path_graph)
    
    cycles_in_path_graph_2 = None

    # create cycle graph
    # nodes are cycles
    # edge(c1,c2) exists iff c1 and c2 have at least one substance as beginning of path in common
    cycle_graph = nx.Graph()

    # add all cycles (nodes) to cycle graph
    cycles_in_path_graph = get_valid_path_graph_cycles(path_graph)
    for cycle in cycles_in_path_graph:
        cycle_graph.add_node(tuple(cycle))

    # NOTE: code below only makes certain that nodes (cycles) connected to at least one other node (cycle) are added to the graph
    # nodes (cycles) that are entirely disconnected from remainder of graph would be left out erroneously! Hence, we now add every cycle to the graph above before connecting all cycles.
    for cycle_1 in cycles_in_path_graph_1:
        # get_valid_path_graph_cycles returns a generator so receiving object of this method is 'empty' after each full loop across the list of valid path graph cycles
        # we could also turn the returned generator into a list but here we choose a loss in performance (keeping the generator nature of returned object) over increased memory usage (conversion to list type)
        cycles_in_path_graph_2 = get_valid_path_graph_cycles(path_graph)
        for cycle_2 in cycles_in_path_graph_2:
            if cycle_1 is not cycle_2 and any(path_in_cycle_1[0] in [path_in_cycle_2[0] for path_in_cycle_2 in cycle_2] for path_in_cycle_1 in cycle_1):
                cycle_graph.add_edge(tuple(cycle_1),tuple(cycle_2)) # nodes must be immutable objects, hence create tuple from list


    # look for all complete subgraphs in complement of cycle_graph
    # complete subgraphs in complement are those cycles that don't share any complexes as beginning of path
    cycle_graph = nx.complement(cycle_graph)

    # all-edges subgraphs
    all_edges = [sc[key]['edges'] for key in sc.keys()]
    all_edges_subgraphs = it.product(*all_edges)

    # mixed edges-cycles and all-cycles subgraphs
    mixed_edges_cycles_subgraphs = []
    #for cycle_comb in complete_subgraphs:
    for cycle_comb in get_all_cliques(cycle_graph):
        # IMPORTANT: cycle_comb is expected to be a list of tuples (each tuple being one cycle)
        if type(cycle_comb) != type(list()):
            #print "get_sensible_subgraphs: before raise: cycle_comb = "+str(cycle_comb)
            cycle_comb = list(cycle_comb)
            

        unused_edges_per_substance = [sc[key]['edges'] for key in sc.keys() if key not in [path[0] for cycle in cycle_comb for path in cycle]]
        for edges_combination in it.product(*unused_edges_per_substance):
            mixed_edges_cycles_subgraphs.append([edge for edge in edges_combination]+[path for cycle in cycle_comb for path in cycle])

    # return list of subgraphs
    return [tuple(sg) for sg in all_edges_subgraphs]+[tuple(sg) for sg in mixed_edges_cycles_subgraphs]

def get_sensible_subgraphs_paths_in_cycles(sc):
    # get path graph
    path_graph = get_path_graph(sc)

    #cycles_in_path_graph = nx.simple_cycles(path_graph)
    cycles_in_path_graph = get_valid_path_graph_cycles(path_graph)

    # for cycle in cycles_in_path_graph:
    #     del cycle[-1]

    cycles_annotated_in_path_graph = []
    for cycle in cycles_in_path_graph:
        cycle_annotated = []
	
        substance_beginnings = [path[0] for path in cycle]
	if len(substance_beginnings)>len(set(substance_beginnings)):
            continue # ignore current cycle since at least one substance node beginning of more than one path
	
	cycle_annotated.append(substance_beginnings)
	cycle_annotated.append(cycle)

	cycles_annotated_in_path_graph.append(cycle_annotated)

    all_edges = [sc[key]['edges'] for key in sc.keys()]
    all_edges_subgraphs = it.product(*all_edges)

    cycle_combs_iters = []
    for repeat in range(1,len(cycles_annotated_in_path_graph)+1):
        cycle_combs_iters.append(it.combinations(cycles_annotated_in_path_graph, repeat))

    cycle_combs = []
    for iterator in cycle_combs_iters:
        for cycle_comb in iterator:
            path_beginnings = [set(cycle[0]) for cycle in cycle_comb]
        #	print path_beginnings
            pairwise_disjoint = all(len(set.intersection(*set_pair))==0 for set_pair in it.product(path_beginnings, path_beginnings) if set_pair[0]!=set_pair[1])
            if pairwise_disjoint:
                cycle_combs.append({'substances': set.union(*path_beginnings), 'cycles': [cycle[1] for cycle in cycle_comb]})

    mixed_edges_cycles_subgraphs = []
    for cycle_comb in cycle_combs:
	for edges_tuple_of_subgraph in it.product(*[sc[key]['edges'] for key in sc.keys() if key not in cycle_comb['substances']]):
            curr_subgraph = [edge for edge in edges_tuple_of_subgraph]
            curr_subgraph = curr_subgraph + [path for cycle in cycle_comb['cycles'] for path in cycle]

            mixed_edges_cycles_subgraphs.append(tuple(curr_subgraph))

    all_cycles_subgraphs = []
    # NOTE: all-cycles subgraphs probably covered by mixed edges-cycles subgraphs!!
    for cycle_comb in cycle_combs:
	if len(cycle_comb['substances'])==len(sc.keys()):
            all_cycles_subgraphs.append(tuple([path for cycle in cycle_comb['cycles'] for path in cycle]))

    # making sure that we don't return the same subgraph twice ... there must be a more efficient way of doing this altogether
#    return list(set(list(all_edges_subgraphs)+mixed_edges_cycles_subgraphs+all_cycles_subgraphs))
#    return list(set([tuple(sorted(sg)) for sg in all_edges_subgraphs]+[tuple(sorted(sg)) for sg in mixed_edges_cycles_subgraphs]+[tuple(sorted(sg)) for sg in all_cycles_subgraphs]))
    #return list(set([tuple(sorted(sg)) for sg in all_edges_subgraphs]+[tuple(sorted(sg)) for sg in mixed_edges_cycles_subgraphs]))
    return [tuple(sg) for sg in all_edges_subgraphs]+[tuple(sg) for sg in mixed_edges_cycles_subgraphs]

def get_sensible_subgraphs_unique_beginnings(sc):
    # sc: dictionary of dictionaries of subgraph component types per complex node in fragment
    # 'sensible' means that in a valid subgraph, every complex node is the beginning of exactly one edge or path
    # this method constructs subgraphs that meets this requirement (get_all_subgraphs returns many subgraphs where complex nodes are the beginning of multiple edges and/or paths)

    # this method for getting sensible subgraphs merely ensures that in every subgraph generated, each substance node is the beginning of exactly one edge or path
    list_of_sc = [] # list of lists of subgraph components
    for comp in sc.keys():
        sc_of_curr_comp = []
        for component_type in sc[comp]:
            for component in sc[comp][component_type]:
                sc_of_curr_comp.append(component)

        list_of_sc.append(sc_of_curr_comp)

    sensible_subgraphs = it.product(*list_of_sc)

    # # this method for getting sensible subgraphs ensures that in every subgraph generated, each substance node is the beginning of exactly one edge or path and that only those paths which form a cycle of length > 1 or single-path cycles are used
    # sscc = get_sscc(sc)
    # scc_list = []
    # for scc in sscc:
    #     all_edges = [] # all-edges and single-path cycles
    #     all_paths = [] # all-paths
    #     for substance in scc:
    #         substance_edges = []
    #         substance_paths = []
    #         for path in sc[substance]['n_paths']+sc[substance]['p_paths']:
    #             if path[0] != path[2] and len(scc)>1:
    #                 substance_paths.append(path)
    #             elif path[0] == path[2]:
    #                 substance_edges.append(path)
    #         for edge in sc[substance]['edges']:
    #             substance_edges.append(edge)
    #         all_edges.append(substance_edges)
    #         all_paths.append(substance_paths)
    #     scc_list.append([all_edges,all_paths])
        
    # # scc_list = [for each strongly-connected component:[[all-edges and single-path cycles], [all paths]] ]
    # # get all combinations of scc: 
    # # scc1 all-edges + scc2 all-edges + ...
    # # scc1 all-edges + scc2 all-paths + ...
    # scc_combinations = it.product(*scc_list)

    # # for each scc combination, create sensible subgraphs (i.e. those where each substance node is beginning of exactly one edge or path)
    # # ignore those scc combinations that have an empty entry for one of the substance nodes (this happens when e.g. one substance node is not beginning of any paths)
    # sensible_subgraphs = []
    # for comb in scc_combinations:
    #     pass # TODO: continue writing here
    
    return sensible_subgraphs

#def validate_subgraph(element_of_all_subgraphs):
def validate_subgraphs(all_subs, sc, f):
    # element_of_all_subgraphs: generator returned by get_all_subgraphs call, holds subgraph components dictionary as first element
#    sc = element_of_all_subgraphs[0]
#    all_subs = element_of_all_subgraphs[1:][0]

    # CHECK: make sure that sc-all_subs pairing is correct, i.e. every preliminary subgraph created contains elements found in dictionary of all subgraph components, sc
    # gr_el = []
    # for key in sc:
    #     for c_type in sc[key]:
    #         for el in sc[key][c_type]:
    #             if el not in gr_el:
    #                 gr_el.append(el)

    # NOTE: don't do this if all_subs is an iterator!
    # for item in all_subs:
    #     for el in item:
    #         if el not in gr_el:
    #             raise Exception('validate_subgraph: sc - all_subs pairing incorrect')
    # END OF CHECK

    # v2 nodes used in present fragment
    v2_nodes_in_fragment = list(f[1])
#    print 'fragment '+str(f)
#    print 'v2 in fragment '+str(v2_nodes_in_fragment)

    valid_subs = []
    for sub in all_subs:
        sub_valid = True
        v1_beginning_in_sub = [] # list of v1 nodes that start a 'sub' element
        for el in sub:
            if el[0] in v1_beginning_in_sub:
                sub_valid = False # if v1 node starts two 'sub' elements, then subgraph 'sub' is invalid
                break
            else:
                v1_beginning_in_sub.append(el[0])

        if not sub_valid:
            continue

        # make sure every v1 node in current fragment is beginning of exactly subgraph element
        # I think this check is redundant
        if len(v1_beginning_in_sub)!=len(sc):
            sub_valid = False

        if not sub_valid:
            continue
            
        # collect v2 nodes
        v2_used = []
        for el in sub:
            #if el[1] not in v2_used:
            # must collect all v2 nodes at their right number!
            v2_used.append(el[1])

        # make sure the right number of v2 nodes is used
        if len(v2_used) != len(v2_nodes_in_fragment):
            sub_valid = False

        if not sub_valid:
            continue

#        if len(v2_used) == len(v2_nodes_in_fragment):
#        print 'sub graph: '+str(sub)
#        print 'v2_used: '+str(v2_used)
#        print 'v2_in_fragment: '+str(v2_nodes_in_fragment)

        # make sure that each v2 node (of present fragment!) is used exactly as often as it appears in fragment listing
        if len(v2_used) == len(v2_nodes_in_fragment):
            for a_v2_n in v2_nodes_in_fragment:
                #print 'trying to remove '+str(a_v2_n)
                try:
                    v2_used.remove(a_v2_n)
                 #   print 'removed '+str(a_v2_n)
                except ValueError:
                    # if a_v2_n is not present (anymore) in v2_nodes_in_fragment then .remove throws this exception
                    # this means that either a_v2_n is not used in present subgraph 'sub' or
                    # a_v2_n is used more often than it should be
                  #  print 'ValueError raised'
                    sub_valid = False
                  #  print 'failed to remove '+str(a_v2_n)
                    break

        if not sub_valid:
            continue
        
        # to validate a subgraph we only want to make sure that every path of the subgraph is part of a cycle
        # all those paths with identical starting and end points form a shortest possible cycle by themselves and can't be part of a bigger cycle
        paths_to_check = []
        for el in sub:
            if len(el)==4 and el[0]!=el[2]:
                # collect all paths in subgraph 'sub'
                # omit those (presumably!) positive paths that form a positive cycle (sk, wi, sk)
                paths_to_check.append(el)

        paths_in_cycle = set()
        for l in range(len(paths_to_check)+1): # REMEMBER: range(3)=[0,1,2]!
#            paths_perms = list(it.combinations(paths_to_check,l))
            #paths_perms = list(it.permutations(paths_to_check,l))  
            
            # itertools.permutations is the correct choice (over combinations) since it returns all possible orderings of permutations of length l, combinations returns just one sorted combination
            paths_perms = it.permutations(paths_to_check,l)

            for pp in paths_perms:
                if is_cycle(pp):
                    paths_in_cycle = paths_in_cycle.union(set(pp))
        if len(paths_in_cycle) != len(paths_to_check):
            sub_valid = False

        if not sub_valid:
            continue

        # while paths_in_sub:
        #     if len(paths_in_sub) == 1:
        #         # if I am about to pop the last remaining (non-self-circular (sk,wi,sk)) path, then I know that this one path won't be able to form a complete cycle, hence subgraph is invalid
        #         sub_valid = False
        #     path = paths_in_sub.pop()
        #     first_v1_visited = path[0]
        #     last_v1_visited=path[2]
        #     for p in paths_in_sub:
        #         if path[2]==p[0]:
        #             last_v1_visited=p[2]
        #             path=p
        #             paths_in_sub.remove(p)

        #         if last_v1_visited==first_v1_visited:
        #             # cycle found, continue to next 'while' iteration
        #             break # this breaks smallest enclosing 'for' loop, i.e. jumps back to 'while' loop

        #     # if all paths in paths_in_sub were removed by last for loop, then we either used up all paths to build a cycle, or we used up all paths chasing down a non-cyclic chain of paths
        #     if len(paths_in_sub)==0 and first_v1_visited!=last_v1_visited:
        #         sub_valid = False
             
        if sub_valid:
            valid_subs.append(sub)

    #print 'validate_subgraphs: len(valid_subs) = '+str(len(valid_subs))

    #return sc, valid_subs
    return valid_subs

def get_all_valid_subgraphs(G, stoich_rank, f):
    # f: validated fragment

    new_variant = True

    if new_variant:
        # new variant: scoring done for all subgraphs on the fly and sum appended to output
        return get_all_valid_subgraphs_with_scoring(G, stoich_rank, f)
    else:
        # old variant: scoring done seperately
        return get_all_valid_subgraphs_without_scoring(G, stoich_rank, f)

def get_all_valid_subgraphs_with_scoring(G, stoich_rank, f):
    # get all subgraph components, i.e. for every substance node all edges, p_paths, n_paths of fragment 'f'
    sc = get_subgraph_components(G, f)

    # generate all possible subgraphs (combinatorically)
    #all_subgraphs = get_all_subgraphs(sc, stoich_rank)
    #all_subgraphs = get_sensible_subgraphs(sc)
    all_subgraphs = get_sensible_subgraphs(sc)
    
    # prune combinatoric list of subgraphs, i.e. get those subgraphs that make sense for present reaction network
    valid_subgraphs = validate_subgraphs(all_subgraphs, sc, f)
    
    # score of fragment f
    K_S = 0
    for sub_graph in valid_subgraphs:
        # get K_g of 'sub_graph'
        K_g = score_subgraph([sub_graph, sc])
        # update K_S
        K_S = K_S + K_g

    if len(valid_subgraphs)==0:
        return None
    else:
        return f, sc, valid_subgraphs, K_S

def get_all_valid_subgraphs_without_scoring(G, stoich_rank, f):
    # get all subgraph components, i.e. for every substance node all edges, p_paths, n_paths of fragment 'f'
    sc = get_subgraph_components(G, f)

    # generate all possible subgraphs (combinatorically)
    #all_subgraphs = get_all_subgraphs(sc, stoich_rank)
    #all_subgraphs = get_sensible_subgraphs(sc)
    all_subgraphs = get_sensible_subgraphs(sc)
    
    # prune combinatoric list of subgraphs, i.e. get those subgraphs that make sense for present reaction network
    valid_subgraphs = validate_subgraphs(all_subgraphs, sc, f)
    
    if len(valid_subgraphs)==0:
        return None
    else:
        return f, sc, valid_subgraphs

def score_subgraph(args):
    sub_graph = args[0] # preferably a previously validated subgraph!
    sc = args[1] # subgraph components found in present fragment

    K_g = 1
    paths = []
    edges = []
    for el in sub_graph:
        if len(el)==2:
            edges.append(el)
        else:
            paths.append(el)
        
    for edge in edges:
        K_g = K_g * 1 # assuming all stoich coefficients = 1!

    # here we need to work out the exact number of unique cycles
    # i.e. (s1, w1, s2)+(s2,w2,s1) and (s2,w1,s1)+(s1,w1,s2) are the same cycle!
    # this is important since t_g = no. of unique cycles in subgraph
    cycles_in_sub_graph = []
    # using it.permutations, we'll get the same unique cycle multiple times
    # keep track of which subgraph components have already been used to avoid counting the same cycle more than once
    sc_used = []
    for l in range(len(paths)+1): # REMEMBER: range(3)=[0,1,2]!
#        el_perms = list(it.combinations(paths,l))
#        el_perms = list(it.permutations(paths,l))
        el_perms = it.permutations(paths,l)
        for elp in el_perms:
            if is_cycle(elp):
                cycle_is_unique=True
                for elp_sc in elp:
                    if not elp_sc in sc_used:
                        sc_used.append(elp_sc)
                    else:
                        cycle_is_unique = False

                if cycle_is_unique:
                    cycles_in_sub_graph.append(elp)

    t_g = len(cycles_in_sub_graph)
    K_C = 1
    for c in cycles_in_sub_graph:
        for p in c:
            if p[3] == 'p':
                K_C = K_C*1
            elif p[3] == 'n':
                K_C = K_C * -1
            else:
                raise

    K_g = K_g * math.pow(-1, t_g) * K_C

    return K_g
