# Created by Georg Walther (georg.walther@jic.ac.uk)
# Bugfixes and improvements by Matthew Hartley

import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import itertools as it
import numpy as np
import sys
from multiprocessing import Pool
import math
import re
import datetime
import cPickle as pickle

# First, identify all fragments in graph. Fragments as defined on p. 74 of [1]
# Then look at each fragment and identify all subgraphs (p. 73, [1]) induced by nodes of fragment
# For each fragment, throw all edges and cycles (of the corresponding undirected graph - induced by nodes of fragment) in a bag and create all permutations (of all cardinalities) with mutually exclusive V1 (complex / substance nodes) from this bag - V2 nodes (reaction nodes) can be visited multiple times

# [1] M. Mincheva, M.R. Roussel (2007): Graph-theoretic methods for the analysis of chemical and biochemical networks. I. ...

def on_draw(event):
    """Auto-wraps all text objects in a figure at draw-time"""
    import matplotlib as mpl
    fig = event.canvas.figure

    # Cycle through all artists in all the axes in the figure
    for ax in fig.axes:
        for artist in ax.get_children():
            # If it's a text artist, wrap it...
            if isinstance(artist, mpl.text.Text):
                autowrap_text(artist, event.renderer)

    # Temporarily disconnect any callbacks to the draw event...
    # (To avoid recursion)
    func_handles = fig.canvas.callbacks.callbacks[event.name]
    fig.canvas.callbacks.callbacks[event.name] = {}
    # Re-draw the figure..
    fig.canvas.draw()
    # Reset the draw event callbacks
    fig.canvas.callbacks.callbacks[event.name] = func_handles

def autowrap_text(textobj, renderer):
    """Wraps the given matplotlib text object so that it exceed the boundaries
    of the axis it is plotted in."""
    import textwrap
    # Get the starting position of the text in pixels...
    x0, y0 = textobj.get_transform().transform(textobj.get_position())
    # Get the extents of the current axis in pixels...
    clip = textobj.get_axes().get_window_extent()
    # Set the text to rotate about the left edge (doesn't make sense otherwise)
    textobj.set_rotation_mode('anchor')

    # Get the amount of space in the direction of rotation to the left and 
    # right of x0, y0 (left and right are relative to the rotation, as well)
    rotation = textobj.get_rotation()
    right_space = min_dist_inside((x0, y0), rotation, clip)
    left_space = min_dist_inside((x0, y0), rotation - 180, clip)

    # Use either the left or right distance depending on the horiz alignment.
    alignment = textobj.get_horizontalalignment()
    if alignment is 'left':
        new_width = right_space 
    elif alignment is 'right':
        new_width = left_space
    else:
        new_width = 2 * min(left_space, right_space)

    # Estimate the width of the new size in characters...
    aspect_ratio = 0.5 # This varies with the font!! 
    fontsize = textobj.get_size()
    pixels_per_char = aspect_ratio * renderer.points_to_pixels(fontsize)

    # If wrap_width is < 1, just make it 1 character
    wrap_width = max(1, new_width // pixels_per_char)
    try:
        wrapped_text = textwrap.fill(textobj.get_text(), wrap_width)
    except TypeError:
        # This appears to be a single word
        wrapped_text = textobj.get_text()
    textobj.set_text(wrapped_text)

def min_dist_inside(point, rotation, box):
    """Gets the space in a given direction from "point" to the boundaries of
    "box" (where box is an object with x0, y0, x1, & y1 attributes, point is a
    tuple of x,y, and rotation is the angle in degrees)"""
    from math import sin, cos, radians
    x0, y0 = point
    rotation = radians(rotation)
    distances = []
    threshold = 0.0001 
    if cos(rotation) > threshold: 
        # Intersects the right axis
        distances.append((box.x1 - x0) / cos(rotation))
    if cos(rotation) < -threshold: 
        # Intersects the left axis
        distances.append((box.x0 - x0) / cos(rotation))
    if sin(rotation) > threshold: 
        # Intersects the top axis
        distances.append((box.y1 - y0) / sin(rotation))
    if sin(rotation) < -threshold: 
        # Intersects the bottom axis
        distances.append((box.y0 - y0) / sin(rotation))
    return min(distances)

def is_cycle(paths):
    # paths: tuple of paths, i.e. ordered set of paths
#    assert(type(paths)==type(tuple()))

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

# def validate_fragments(f):
#     indgr = G.subgraph([item for sublist in f for item in sublist])
#     edges = []
#     non_edges = [] # those arcs that go from reaction node to complex node
#     for arc in indgr.edges():
#         if indgr.node[arc[0]]['bipartite']==0:
#             edges.append(arc)
#         else:
#             non_edges.append(arc)
                
#     cycles_undir = nx.cycle_basis(nx.Graph(indgr))
#     v1_nodes = []
#     for e in edges:
#         if not e[0] in v1_nodes:
#             v1_nodes.append(e[0])
#     for c in cycles_undir:
#         for n in c:
#             if 's' in n and n not in v1_nodes:
#                 v1_nodes.append(n)
#     if len(v1_nodes) == stoich_rank:
#         f_valid = True # assume f is valid since v1 count is fine
#         if len(non_edges) > 0:
#             # now check if f is really valid or if it has any ('w','s') arcs that are not part of a positive path within a cycle!
#             cycles_undir_indgr=G.subgraph([item for sublist in cycles_undir for item in sublist])
#             for arc in non_edges:
#                 if arc not in cycles_undir_indgr.edges():
#                     # the non_edge arc being part of an undirected cycle is necessary
#                     f_valid = False
#                 else:
#                     # for each non-edge arc ('wx','si') of a cycle, there must be an edge ('sk','wx'), where k == i is allowed
#                     # hence, extract all edges arcs from cycles_undir_indgr and compare current non-edge arc 'arc' to each edge
#                     edges_in_undir_cycle = []
#                     for cycle_arc in cycles_undir_indgr.edges():
#                         if cycles_undir_indgr.node[cycle_arc[0]]['bipartite'] == 0:
#                             edges_in_undir_cycle.append(cycle_arc)
                    
#                     arc_in_positive_path = False
#                     for cycle_edge in edges_in_undir_cycle:
#                         if arc[0] == cycle_edge[1]:
#                             arc_in_positive_path = True
                    
#                     # as soon as one 'arc' is not in a positive path, then the entire fragment 'f' is invalid
#                     f_valid = f_valid and arc_in_positive_path

#         if f_valid:
#             return f
#             #valid_fragments.append(f)

#     return None
#     #return valid_fragments

def validate_fragments(f):
    indgr = G.subgraph([item for sublist in f for item in sublist])
    edges = []
    non_edges = [] # those arcs that go from reaction node to complex node
    for arc in indgr.edges():
        if indgr.node[arc[0]]['bipartite']==0:
            edges.append(arc)
        else:
            non_edges.append(arc)
                
 #   print 'fragment: '+str(f)
  #  print 'edges: '+str(edges)
   # print 'non-edges: '+str(non_edges)

    v1_nodes = []
    for e in edges:
        if not e[0] in v1_nodes:
            v1_nodes.append(e[0])

    f_valid = False
    if len(v1_nodes) == stoich_rank:
        f_valid = True # assume f is valid since v1 count is fine
#        print 'initially: '+str(f_valid)
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

def get_subgraph_components(f):
    indgr = G.subgraph([item for sublist in f for item in sublist])
    cycles_undir = nx.cycle_basis(nx.Graph(indgr))
    cycles_undir_indgr=G.subgraph([item for sublist in cycles_undir for item in sublist])
    
    edges = []
    non_edges = [] # those arcs that go from reaction node to complex node
    for arc in indgr.edges():
        if indgr.node[arc[0]]['bipartite']==0:
            edges.append(arc)
        else:
            non_edges.append(arc)
    
    # cycle_edges = []
    # cycle_non_edges = []
    # for arc in cycles_undir_indgr.edges():
    #     if cycles_undir_indgr.node[arc[0]]['bipartite']==0:
    #         cycle_edges.append(arc)
    #     else:
    #         cycle_non_edges.append(arc)
    
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
        for n in non_edges:
            if e[1]==n[0]:
                #print (e[0],e[1],n[1])
                p_edge = (e[0],e[1],n[1])
                if p_edge not in subgraph_components[e[0]]['p_paths']:
                    subgraph_components[e[0]]['p_paths'].append(p_edge)

    # add all negative paths
    # NOTE: it doesn't matter if a negative path is part of a cycle. here you want to add ALL negative paths
    # whether or not every negative path contributes to a sensible subgraph is decided later in validate_subgraph
    # for e1 in cycle_edges:
    #     for e2 in cycle_edges:
    for e1 in edges:
        for e2 in edges:
            if e1[0] != e2[0] and e1[1]==e2[1]:
                #print (e1[0],e1[1],e2[0])
                n_edge = (e1[0],e1[1],e2[0])
                if n_edge not in subgraph_components[e1[0]]['n_paths']:
                    subgraph_components[e1[0]]['n_paths'].append(n_edge)

    return subgraph_components

def get_all_subgraphs(sc):
    # sc: subgraph components dictionary
    gr_el = []
    for key in sc:
        for c_type in sc[key]:
            for el in sc[key][c_type]:
                if el not in gr_el:
                    gr_el.append(el)
#    return sc, list(it.combinations(gr_el, stoich_rank))
    return list(it.combinations(gr_el, stoich_rank))

#def validate_subgraph(element_of_all_subgraphs):
def validate_subgraphs(all_subs, sc, f):
    # element_of_all_subgraphs: generator returned by get_all_subgraphs call, holds subgraph components dictionary as first element
#    sc = element_of_all_subgraphs[0]
#    all_subs = element_of_all_subgraphs[1:][0]

    # CHECK: make sure that sc-all_subs pairing is correct, i.e. every preliminary subgraph created contains elements found in dictionary of all subgraph components, sc
    gr_el = []
    for key in sc:
        for c_type in sc[key]:
            for el in sc[key][c_type]:
                if el not in gr_el:
                    gr_el.append(el)

    for item in all_subs:
        for el in item:
            if el not in gr_el:
                raise Exception('validate_subgraph: sc - all_subs pairing incorrect')
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
            else:
                v1_beginning_in_sub.append(el[0])

        # make sure every v1 node in current fragment is beginning of exactly subgraph element
        # I think this check is redundant
        if len(v1_beginning_in_sub)!=len(sc):
            sub_valid = False
            
        # collect v2 nodes
        v2_used = []
        for el in sub:
            #if el[1] not in v2_used:
            # must collect all v2 nodes at their right number!
            v2_used.append(el[1])

        # make sure the right number of v2 nodes is used
        if len(v2_used) != len(v2_nodes_in_fragment):
            sub_valid = False

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
        
        # to validate a subgraph we only want to make sure that every path of the subgraph is part of a cycle
        # all those paths with identical starting and end points form a shortest possible cycle by themselves and can't be part of a bigger cycle
        paths_to_check = []
        for el in sub:
            if len(el)==3 and el[0]!=el[2]:
                # collect all paths in subgraph 'sub'
                # omit those (presumably!) positive paths that form a positive cycle (sk, wi, sk)
                paths_to_check.append(el)

        paths_in_cycle = set()
        for l in range(len(paths_to_check)+1): # REMEMBER: range(3)=[0,1,2]!
#            paths_perms = list(it.combinations(paths_to_check,l))
            paths_perms = list(it.permutations(paths_to_check,l))  # itertools.permutations is the correct choice (over combinations) since it returns all possible orderings of permutations of length l, combinations returns just one sorted combination
            for pp in paths_perms:
                if is_cycle(pp):
                    paths_in_cycle = paths_in_cycle.union(set(pp))
        if len(paths_in_cycle) != len(paths_to_check):
            sub_valid = False

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

    #return sc, valid_subs
    return valid_subs

def get_all_valid_subgraphs(f):
    # f: validated fragment

    # get all subgraph components, i.e. for every substance node all edges, p_paths, n_paths of fragment 'f'
    sc = get_subgraph_components(f)

    # generate all possible subgraphs (combinatorically)
    all_subgraphs = get_all_subgraphs(sc)

    # prune combinatoric list of subgraphs, i.e. get those subgraphs that make sense for present reaction network
    valid_subgraphs = validate_subgraphs(all_subgraphs, sc, f)

    if len(valid_subgraphs)==0:
        return None
    else:
        return f, sc, valid_subgraphs
        
def score_fragment(args):
    # vs: a list of validated subgraphs, one list per fragment
    # sc: dictionary of corresponding subgraph components
    # frags: validated fragment
    vs = args[0]
    sc = args[1]
    frag = args[2]

    K_S = 0
    for sub_graph in vs:
        # get K_g of 'sub_graph'
        K_g = score_subgraph([sub_graph, sc])
        # update K_S
        K_S = K_S + K_g

        # K_g = 1
#         paths = []
#         edges = []
#         for el in sub_graph:
#             if len(el)==2:
#                 edges.append(el)
#             else:
#                 paths.append(el)
        
#         for edge in edges:
#             K_g = K_g * 1 # assuming all stoich coefficients == 1!

#         cycles_in_sub_graph = []
#         for l in range(len(paths)+1): # REMEMBER: range(3)=[0,1,2]!
# #            el_perms = list(it.permutations(paths,l))
#             el_perms = list(it.combinations(paths,l))
#             for elp in el_perms:
#                 if is_cycle(elp):
#                     #if len(elp)>1:
#                     cycles_in_sub_graph.append(elp)

#         t_g = len(cycles_in_sub_graph)
#         K_C = 1
#         for c in cycles_in_sub_graph:
#             for p in c:
#                 #if p in sc[p[0]]['p_paths'] and p in sc[p[0]]['n_paths']:
#                 if False:
#                     raise Exception('score_fragments: path is both positive and negative')
#                 elif p in sc[p[0]]['p_paths']:
#                     K_C = K_C * 1
#                 elif p in sc[p[0]]['n_paths']:
#                     K_C = K_C * -1
#                 else:
#                     raise Exception('score_fragments: path is neither positive nor negative')

        # while paths:
        #     path = paths.pop()
        #     first_v1 = path[0]
        #     last_v1 = path[2]
            
        #     # Update K_C depending on whether 'path' is positive or negative
        #     # We know we're looking at validated subgraphs, i.e. every path is part of a cycle so
        #     # we're good updating K_C at this point
        #     if path in sc[path[0]]['p_paths']:
        #         K_C = K_C * 1
        #     elif path in sc[path[0]]['n_paths']:
        #         K_C = K_C * (-1)
        #     else:
        #         raise Exception('score_fragments: path is neither positive nor negative')

        #     if first_v1 == last_v1:
        #         # single-path cycle detected
        #         t_g = t_g + 1

        #         # make sure single-path cycle is a positive path!
        #         if not path in sc[path[0]]['p_paths']:
        #             raise Exception('score_fragments: single-path cycle is not a positive path')

        #         # skip remainder of this 'while' iteration
        #         continue

        #     # 'path' doesn't induce a single-path cycle, hence let's construct and score the sequence of paths making up the current cycle until we revisit 'first_v1'
        #     for p in paths:
        #         if path[2]==p[0]:
        #             # update last_v1
        #             last_v1 = p[2]
        #             # update 'path'
        #             path = p
        #             # update K_C
        #             if path in sc[path[0]]['p_paths']:
        #                 K_C = K_C * 1
        #             elif path in sc[path[0]]['n_paths']:
        #                 K_C = K_C * (-1)
        #             else:
        #                 raise Exception('score_fragments: path is neither positive nor negative')
        #             # remove 'p' from 'paths' so this isn't revisited
        #             paths.remove(p)
        #             # check if cycle completed
        #             if first_v1 == last_v1:
        #                 t_g = t_g + 1

        #                 # if cycle completed, skip to next 'while' iteration
        #                 continue

    # return K_S
    return frag, sc, K_S

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
        el_perms = list(it.permutations(paths,l))
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
            if False:
                raise Exception('score_fragments: path is both positive and negative')
            elif p in sc[p[0]]['p_paths']:
                K_C = K_C * 1
            elif p in sc[p[0]]['n_paths']:
                K_C = K_C * -1
            else:
                raise Exception('score_fragments: path is neither positive nor negative')

    K_g = K_g * math.pow(-1, t_g) * K_C

    return K_g

def get_graph_stoich(alpha, beta):

    # stoichiometry matrix
    stoich = beta - alpha
    # rank of stoichiometry matrix
    stoich_rank = np.linalg.matrix_rank(stoich)
    # number of reactions = number of columns of stoich
    no_rxn = stoich.shape[1]
    # number of substance = number of rows of stoich
    no_sub = stoich.shape[0]

    # make directed graph
    G = nx.DiGraph()

    # add nodes to graph: add all reactions 'w' and substances 's'
    G.add_nodes_from(['w'+str(n) for n in range(1,no_rxn+1)], bipartite=1)
    G.add_nodes_from(['s'+str(n) for n in range(1,no_sub+1)], bipartite=0)

    # add ('s','w') edges, i.e. positive element in alpha
    for row_i in range(no_sub):
        curr_row = alpha[row_i]
        for col_i in range(no_rxn):
            curr_coeff = curr_row[col_i]
            if curr_coeff > 0:
                G.add_edge('s'+str(row_i+1), 'w'+str(col_i+1), coeff=curr_coeff)

    # add ('w','s') edges, i.e. positive element in beta
    for row_i in range(no_sub):
        curr_row = beta[row_i]
        for col_i in range(no_rxn):
            curr_coeff = curr_row[col_i]
            if curr_coeff > 0:
                G.add_edge('w'+str(col_i+1), 's'+str(row_i+1), coeff=curr_coeff)

    return G, stoich, stoich_rank
def get_lpa_alpha_beta(alpha, beta, slow_indices):
    # alpha: reactant stoichiometric coefficients of well-mixed system
    # beta: product stoichiometric coefficients of well-mixed system
    # slow_indices: indices of slow species in system; indices are row indices of alpha, beta

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

    # lpa alpha matrix
    alpha_lpa = []
    # slow global variables
    for slow_i in slow_indices:
        alpha_lpa_row = list(alpha[slow_i])
        alpha_lpa_row = alpha_lpa_row + [0 for r in range(no_rxn)]
        alpha_lpa.append(alpha_lpa_row)
        #print "global slow_i = "+str(slow_i)
        #print alpha_lpa_row
    # fast variables
    for fast_i in fast_indices:
        alpha_lpa_row = list(alpha[fast_i])
        alpha_lpa_row = alpha_lpa_row + list(alpha[fast_i])
        alpha_lpa.append(alpha_lpa_row)
        #print "fast_i = "+str(fast_i)
        #print alpha_lpa_row
    # slow local variables
#    for slow_i in reversed(slow_indices): # not sure why I'm using reversed here ... ?
    for slow_i in slow_indices: # not sure why I'm using reversed here ... ?
        alpha_lpa_row = [0 for r in range(no_rxn)]
        alpha_lpa_row = alpha_lpa_row + list(alpha[slow_i])
        alpha_lpa.append(alpha_lpa_row)
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

    # convert to expected numpy.array type
    return np.array(alpha_lpa), np.array(beta_lpa)

# def network_from_gratelpy(filename):
#     # read reactions from gratelpy file and create network
#     # returns alpha and beta matrices that define reaction network

#     # open file
#     f = open(filename)

#     # rate constants
#     rates = []
#     # species
#     species = []

#     # go through file and collect all rate constants and species
#     for l in f:
#         rates_in_line = re.findall('k_\{{1}.*?\}{1}',l)
#         species_in_line = re.findall('\[{1}.*?\]{1}',l.split("=")[1])
        
#         for rate_in_line in rates_in_line:
#             for not rate_in_line in rates:
#                 rates.append(rate_in_line)

def get_network_from_mechanism(filename, no_complexes):
    # reads mechanism file: each line is of the form alpha_1 [A] -> beta_1 [B]; k_rat
    # NOTE: mechanism file is read case-sensitively

    # returns dictionaries of all detected complexes and rate constants as well as alpha and beta matrices

    # open file
    mechanism_file = open(filename)

    # go through file line-by-line and create columns of alpha and beta
    complexes = {}
    complexes_reverse = {}
    constants = {}
    constants_reverse = {}
    alpha_transposed = []
    beta_transposed = []

    # count number of lines = number of rate constants
    no_constants = len([1 for el in enumerate(mechanism_file)])
    mechanism_file.seek(0) # go back to beginning of file

    for line in mechanism_file:
        # tidy up current line
        curr_line = line.strip() # removes trailing and leading whitespaces
        curr_line = curr_line.replace(' ', '') # removes all remaining whitespaces on line

        # get left, reactant side
        curr_left = curr_line.split('->')[0]
        # get right, product side
        curr_right = curr_line.split('->')[1]
        curr_right_split = curr_right.split(';')
        curr_right = curr_right_split[0]
        # get current constant name and add to dictionary
        curr_const = curr_right_split[1]
        if not curr_const in constants.keys():
            constants[curr_const] = len(constants)
            constants_reverse[len(constants_reverse)] = curr_const

        # create current alpha column and add to alpha_transposed
        curr_left = re.split('(\+)', curr_left)
        curr_left = [el for el in curr_left if len(el)>0] # remove empty elements from list
        alpha_column = [0 for el in range(no_complexes)]
        for el in curr_left:
            if el != '+':
                # detect alpha coefficient in 'el': 'curr_coeff'
                curr_coeff_match = re.findall('^\d+', el)
                curr_coeff = -1
                if len(curr_coeff_match) == 0:
                    curr_coeff = 1
                elif len(curr_coeff_match) == 1:
                    curr_coeff = int(curr_coeff_match[0])
                else:
                    raise Exception('too many coefficients detected: '+ el)
                    
                # detect complex in 'el': 'curr_complex'
                curr_complex_match = re.findall('\[{1}.*?\]{1}', el)
                curr_complex = 'not_a_complex'
                if len(curr_complex_match)==0:
                    raise Exception('no complex detected: '+ el)
                elif len(curr_complex_match)==1:
                    curr_complex = curr_complex_match[0]
                else:
                    raise Exception('more than one complex detected: '+ el)

                # add current complex to dictionary if not yet recorded
                if not curr_complex in complexes.keys():
                    complexes[curr_complex] = len(complexes)
                    complexes_reverse[len(complexes_reverse)] = curr_complex

                # modify alpha_column accordingly
                alpha_column[complexes[curr_complex]] = curr_coeff

        # add alpha_column to alpha_transposed
        alpha_transposed.append(alpha_column)

        # create current beta column and add to beta_transposed
        curr_right = re.split('(\+)', curr_right)
        curr_right = [el for el in curr_right if len(el)>0] # remove empty elements from list
        beta_column = [0 for el in range(no_complexes)]
        for el in curr_right:
            if el != '+':
                # detect beta coefficient in 'el': 'curr_coeff'
                curr_coeff_match = re.findall('^\d+', el)
                curr_coeff = -1
                if len(curr_coeff_match) == 0:
                    curr_coeff = 1
                elif len(curr_coeff_match) == 1:
                    curr_coeff = int(curr_coeff_match[0])
                else:
                    raise Exception('too many coefficients detected: '+ el)
                    
                # detect complex in 'el': 'curr_complex'
                curr_complex_match = re.findall('\[{1}.*?\]{1}', el)
                curr_complex = 'not_a_complex'
                if len(curr_complex_match)==0:
                    raise Exception('no complex detected: '+ el)
                elif len(curr_complex_match)==1:
                    curr_complex = curr_complex_match[0]
                else:
                    raise Exception('more than one complex detected: '+ el)

                # add current complex to dictionary if not yet recorded
                if not curr_complex in complexes.keys():
                    complexes[curr_complex] = len(complexes)
                    complexes_reverse[len(complexes_reverse)] = curr_complex

                # modify beta_column accordingly
                beta_column[complexes[curr_complex]] = curr_coeff

        # add alpha_column to alpha_transposed
        beta_transposed.append(beta_column)

    # close file
    mechanism_file.close()

    # check if we detected all complexes
    if len(complexes) != no_complexes:
        raise Exception('number of complexes indicated: '+str(no_complexes) + ' - complexes detected: '+str(complexes.keys()))
    # check if we detected all rate constants
    if len(constants) != no_constants:
        raise Exception('number of rows in mechanism file counted: '+str(no_constants) + ' - constants detected: '+str(constants.keys()))

    # create alpha and beta
    alpha_transposed = np.array(alpha_transposed)
    beta_transposed = np.array(beta_transposed)
    alpha = np.transpose(alpha_transposed)
    beta = np.transpose(beta_transposed)
    
    return alpha, beta, complexes, constants, complexes_reverse, constants_reverse

# no_substances = 2
# no_reactions = 3

# row_base = list(it.product([0,1], repeat=no_reactions))
# matrices = list(it.product(row_base, repeat=no_substances))
# alpha_beta_combinations=list(it.product(matrices, repeat=2))

alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism('test_mechanism.txt', 4)
#alpha, beta, complexes, constants, complexes_reverse, constants_reverse = get_network_from_mechanism('test_mechanism.txt', 4)

ct = datetime.datetime.now()
timestamp = str(ct.year)+'_'+str(ct.month)+'_'+str(ct.day)+'_'+str(ct.hour)+'_'+str(ct.minute)+'_'+str(ct.second)
alpha_beta_combinations=[[alpha,beta]]

print str(len(alpha_beta_combinations))+' networks tested'

# graph index, for interesting graphs to .pdf
graph_i = 0

for comb in alpha_beta_combinations:
    alpha = np.array(comb[0])
    beta = np.array(comb[1])

    G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

    print 'rank of stoichiometry matrix = ' + str(stoich_rank)

    complexes, reactions = bipartite.sets(G)

    # bipartite = 1: reaction, w, nodes ('bottom nodes')
    # bipartite = 0: complex, s, nodes ('top nodes')
    # bottom_nodes, top_nodes = bipartite.sets(G)
#    reactions, complexes = bipartite.sets(G)

    complexes = list(complexes)
    reactions = list(reactions)

    complex_perms = list(it.combinations(complexes,stoich_rank))
    reaction_perms = list(it.combinations_with_replacement(reactions,stoich_rank))
    fragments = list(it.product(complex_perms, reaction_perms))

    #print 'number of fragments: ' + str(len(fragments))

    # super efficient: let's go through all fragments and only keep those that are valid (i.e. have as many complex nodes as the rank of the stoichiometric matrix)
    valid_fragments = []

    
    pool = Pool()
    chunksize = 100

    
    # fragment_list = pool.imap(validate_fragments, fragments, chunksize)
    # valid_fragments = [f for f in fragment_list if f is not None]
    # print 'valid fragments received from \'validate_fragments\': '+str(len(valid_fragments))
    # #for f in valid_fragments:
    # #    print f
    # valid_fragments = get_unique_fragments(valid_fragments)
    # print 'valid fragments received from \'get_unique_fragments\': '+str(len(valid_fragments))

    # pickle.dump(valid_fragments, open('valid_fragments.p', 'wb'))
    
    valid_fragments = pickle.load(open('valid_fragments.p'))

    res_get_all_valid_subgraphs = pool.imap(get_all_valid_subgraphs, valid_fragments, chunksize)
    valid_fragments = []
    sc = []
    valid_subs = []

    pickle.dump(sc, open('sc.p', 'wb'))
    pickle.dump(valid_subs, open('valid_subs.p', 'wb'))

    for item in res_get_all_valid_subgraphs:
        print item
        if item is not None:
            valid_fragments.append(item[0])
            sc.append(item[1])
            valid_subs.append(item[2:][0])

    pickle.dump(sc, open('sc.p', 'wb'))
    pickle.dump(valid_subs, open('valid_subs.p', 'wb'))

        #print item[0]
        #print item[2:][0]

    #print valid_fragments
    #print sc
    #print valid_subs
    
    # for f in valid_fragments:
    #     print f

    # valid_fragments = get_unique_fragments(valid_fragments)
    # print len(valid_fragments)
    # print len(get_unique_fragments(valid_fragments))

    # for f in valid_fragments:
    #     print f

    # #print 'number of valid fragments: '+str(len(valid_fragments))

    # subgraph_components = pool.imap(get_subgraph_components, valid_fragments, chunksize)

    # valid_fragments = []
    # subgraph_comps = []
    # for item in subgraph_components:
    #     valid_fragments.append(item[0])
    #     subgraph_comps.append(item[1:][0])
    
    #print 'number of fragments we got subgraph components of: ' + str(len(subgraph_comps))


    # for each valid fragment, construct all possible constituent subgraphs
    # brute force, combinatorical way of doing this:
    # for each complex node, enumerate the edges, positive paths, and negative paths (the latter two embedded in cycles) that you can start with them and create a bag for each complex node
    # construct all subgraphs combinatorically

    # all_subgraphs = pool.imap(get_all_subgraphs, subgraph_comps, chunksize)
    # validated_subgraphs = pool.imap(validate_subgraph, [item for item in all_subgraphs], chunksize)

    # sc = []
    # valid_subs = []
    # for item in validated_subgraphs:
    #     if len(item[1:][0])>0: # validate_subgraphs returns empty list if no valid subgraphs supplied
    #         sc.append(item[0])
    #         valid_subs.append(item[1:][0])

#print type(sc) # sc is a list of dictionaries

    fragment_scores = pool.imap(score_fragment, [[valid_subs[i], sc[i], valid_fragments[i]] for i in range(len(sc))], 10)
    sc = []
    fs = []
    frags = []
    for item in fragment_scores:
        frags.append(item[0])
        sc.append(item[1])
        fs.append(item[2:][0])

    # frags = []
    # for d in sc:
    #     curr_v1 = []
    #     curr_v2 = []
    #     for key in d:
    #         for g_type in d[key]:
    #             for g_el in d[key][g_type]:
    #                 for n in g_el:
    #                     if G.node[n]['bipartite']==0 and n not in curr_v1:
    #                         curr_v1.append(n)
    #                     if G.node[n]['bipartite']==1 and n not in curr_v2:
    #                         curr_v2.append(n)

    #     frags.append((tuple(curr_v1),tuple(curr_v2)))

    # if len(fs)==0:
    #     print G.nodes()
    #     print G.edges()
        
        
    if len(fs)>0:
        substance_nodes = []
        reaction_nodes = []
        for n in G.nodes():
            if G.node[n]['bipartite']==0 and n not in substance_nodes:
                substance_nodes.append(n)
            if G.node[n]['bipartite']==1 and n not in reaction_nodes:
                reaction_nodes.append(n)
        substance_graph = nx.DiGraph()
        substance_graph.add_nodes_from(substance_nodes)
        reaction_graph = nx.DiGraph()
        reaction_graph.add_nodes_from(reaction_nodes)

        # print G.nodes()
        # print G.edges()
        # print sc
        # print frags
        # print fs
        
        # number of reactions = number of columns of stoich
        no_rxn = stoich.shape[1]
        # number of substance = number of rows of stoich
        no_sub = stoich.shape[0]

        # pos = {
        #     'w6': [0,0], 's4': [1,1], 's3': [3,0],
        #     'w5': [2,2], 'w4': [3,2],
        #     'w2': [0,4], 's1': [1,4], 'w3': [2,4], 's2': [3,4],
        #     'w1': [1,6]
        #     }
        
        pos = nx.spring_layout(G)

        labels={}
        for n in G.nodes():
                labels[n]=str(n)


        # page number
        page_i = 0

        # number of graph drawings per page
        no_draw_rows = 2
        no_draw_columns = 2

        # draw full graph with red nodes
        ctr = 1
        plt.figure()
        plt.subplot(no_draw_rows,no_draw_columns,ctr)
#        plt.suptitle('graph '+str(graph_i)+', page '+str(page_i))
        plt.suptitle('page '+str(page_i+1))

        plt.axis('off')
        nx.draw_networkx_nodes(substance_graph, pos, node_shape='o', node_size=400, node_color='r')
        nx.draw_networkx_nodes(reaction_graph, pos, node_shape='s', node_size=600, node_color='r')
        nx.draw_networkx_edges(G, pos, width=2)
        nx.draw_networkx_labels(G, pos, labels, font_size=18)

        # cycle over all fragments in graph G
        for frag_i in range(len(frags)):
            
            print 'fragment '+str(frag_i)
            print frags[frag_i]
            # if frag_i == 3:
            #     print 'subgraphs in fragment 3:'
            #     for sg in valid_subs[frag_i]:
            #         print sg
            #     print 'subgraph components:'
            #     print sc[frag_i]

            # draw fragment nodes in green
            ctr = ctr + 1
            # get all nodes in current fragment
            nodes_in_frag = [item for sublist in frags[frag_i] for item in sublist]
            # determine which are substance nodes and which are reaction nodes and create corresponding graphs
            substance_nodes = []
            reaction_nodes = []

            for n in nodes_in_frag:
                if G.node[n]['bipartite']==0 and n not in substance_nodes:
                    substance_nodes.append(n)
                if G.node[n]['bipartite']==1 and n not in reaction_nodes:
                    reaction_nodes.append(n)

            substance_graph = nx.DiGraph()
            substance_graph.add_nodes_from(substance_nodes)
            reaction_graph = nx.DiGraph()
            reaction_graph.add_nodes_from(reaction_nodes)

            # get all arcs in fragment
            arcs_in_frag = []
            for n in nodes_in_frag:
                for key in G.edge[n]:
                    an_arc = (n,key)
                    if an_arc not in arcs_in_frag:
                        arcs_in_frag.append(an_arc)

            frag_graph = nx.DiGraph()
            frag_graph.add_edges_from(arcs_in_frag)

            #substance_graph.add_edges_from(arcs_in_frag)
            #reaction_graph.add_edges_from(arcs_in_frag)

            plt.subplot(no_draw_rows,no_draw_columns,ctr)
            plt.title('frag. '+str(frag_i+1)+', K_S = '+str(fs[frag_i]))
            plt.axis('off')
            plt.suptitle('page '+str(page_i+1))
            plt.annotate('S_'+str(stoich_rank)+' '+str(frags[frag_i]), (0,0), (0, -3), xycoords='axes fraction', textcoords='offset points', va='top', fontsize=8)
            plt.subplots_adjust(hspace=0.2, top=.9, bottom=.1)

            nx.draw_networkx_nodes(substance_graph, pos, node_shape='o', node_size=400, node_color='g', alpha=.5)
            nx.draw_networkx_nodes(reaction_graph, pos, node_shape='s', node_size=600, node_color='g', alpha=.5)
#            nx.draw_networkx_edges(frag_graph, pos)
            #nx.draw_networkx_edges(reaction_graph, pos)
            labels = {}
            for n in substance_nodes:
                labels[n]=str(n)
            for n in reaction_nodes:
#                labels[n]=str(n)+':'+str(frags[frag_i][1].count(n))
                labels[n]=str(n)
            nx.draw_networkx_labels(G, pos, labels, font_size=18)

            # check if ctr == 6, if it does save current plot to .pdf with graph number and page number and set ctr = 1 again
            if ctr == no_draw_rows*no_draw_columns:
                plt.savefig('graph_'+str(graph_i).zfill(4)+'_page_'+str(page_i).zfill(4)+'.pdf')
                plt.clf()
                ctr = 0 # for both frag and subgraph drawing we increase ctr first and then draw, so want ctr == 1 for first drawing
                page_i = page_i + 1

            # draw every valid subgraph found in current fragment
            sg_ctr = -1
            for sg in valid_subs[frag_i]:
                sg_ctr = sg_ctr + 1 # subgraph counter, for subgraph labelling in output

                arcs_in_sg = []
                nodes_in_sg = [item for sublist in sg for item in sublist]
                # find all edges in subgraph
                for el in sg:
                    if len(el)==2 and el not in arcs_in_sg:
                        arcs_in_sg.append(el)
                    elif len(el)==3:
                        if el in sc[frag_i][el[0]]['p_paths']:
                            arc_1_in_path = [el[0],el[1]]
                            arc_2_in_path = [el[1],el[2]]
                        elif el in sc[frag_i][el[0]]['n_paths']:
                            arc_1_in_path = [el[0], el[1]]
                            arc_2_in_path = [el[2], el[1]]
                        if arc_1_in_path not in arcs_in_sg:
                            arcs_in_sg.append(arc_1_in_path)
                        if arc_2_in_path not in arcs_in_sg:
                            arcs_in_sg.append(arc_2_in_path)
                sg_graph = nx.DiGraph()
                sg_graph.add_edges_from(arcs_in_sg)


                # collect all nodes in subgraph
                substance_nodes = []
                reaction_nodes = []
                for n in nodes_in_sg:
                    if G.node[n]['bipartite']==0 and n not in substance_nodes:
                        substance_nodes.append(n)
                    if G.node[n]['bipartite']==1 and n not in reaction_nodes:
                        reaction_nodes.append(n)
                substance_graph = nx.DiGraph()
                substance_graph.add_nodes_from(substance_nodes)
                reaction_graph = nx.DiGraph()
                reaction_graph.add_nodes_from(reaction_nodes)

                ctr = ctr + 1
                plt.subplot(no_draw_rows,no_draw_columns,ctr)
                plt.suptitle('page '+str(page_i+1))
#                plt.suptitle('graph '+str(graph_i)+', page '+str(page_i))
                plt.title('frag. '+str(frag_i+1)+', s.gr. '+str(sg_ctr+1)+', K_g = '+str(score_subgraph([sg,sc[frag_i]])))
                plt.axis('off')
                subgraph_text = str(sg)
                plt.annotate(subgraph_text, (0,0), (0, -3), xycoords='axes fraction', textcoords='offset points', va='top', fontsize=8)
                plt.subplots_adjust(hspace=0.2, top=.9, bottom=.1)

                nx.draw_networkx_nodes(substance_graph, pos, node_shape='o', node_size=400, node_color='y')
                nx.draw_networkx_nodes(reaction_graph, pos, node_shape='s', node_size=600, node_color='y')
                nx.draw_networkx_edges(sg_graph, pos)
                labels = {}
                for n in substance_nodes:
                    labels[n]=str(n)
                for n in reaction_nodes:
                    labels[n]=str(n)
                nx.draw_networkx_labels(G, pos, labels, font_size=18)

#                ax1.axis([0.,1.,.1,1.])
                
            
                # check if ctr == 6, if it does save current plot to .pdf with graph number and page number and set ctr = 1 again
                if ctr == no_draw_rows*no_draw_columns:
                    plt.savefig('graph_'+str(graph_i).zfill(4)+'_page_'+str(page_i).zfill(4)+'.pdf')
                    plt.clf()
                    ctr = 0 # for both frag and subgraph drawing we increase ctr first and then draw, so want ctr == 1 for first drawing
                    page_i = page_i + 1

        plt.savefig('graph_'+str(graph_i).zfill(4)+'_page_'+str(page_i).zfill(4)+'.pdf')
            
        # increase graph_i since we just finished one interesting graph
        graph_i = graph_i+1

    pool.close()
