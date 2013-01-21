import networkx as nx
from matplotlib import pylab as pl

def neighbors_bits(G, node):
    nodes_sorted = sorted(G.nodes())
    
    neighbors = [0 for n in nodes_sorted]
    for index, n in enumerate(nodes_sorted):
        if n in G.neighbors(node):
            neighbors[index] = 1
    return neighbors

def greater_neighbors_bits(G, a_node):
    nodes_sorted = sorted(G.nodes())
    a_node_index = nodes_sorted.index(a_node)

    neighbors_of_a_node = [0 for n in nodes_sorted]

    for another_node_index, another_node in enumerate(nodes_sorted):
        # only set bits to one if index of another_node is greater than index of a_node
        # thus a_node isn't indicated as its own neighbour
        if another_node_index > a_node_index and another_node in G.neighbors(a_node):
            neighbors_of_a_node[another_node_index] = 1

    return neighbors_of_a_node

def greater_neighbors_strings(G, a_node):
    nodes_sorted = sorted(G.nodes())
    a_node_index = nodes_sorted.index(a_node)

    neighbors_of_a_node = []

    for another_node_index, another_node in enumerate(nodes_sorted):
        if another_node_index > a_node_index and another_node in G.neighbors(a_node):
            neighbors_of_a_node.append(another_node)
    
    return tuple(neighbors_of_a_node)

def neighbors_bits_to_string(G, neighbors_bits):
    nodes_sorted = sorted(G.nodes())
    
    return [nodes_sorted[n_index] for n_index, n_bit in enumerate(neighbors_bits) if n_bit == 1]

G = nx.Graph()
edges_fig_4 = [('a','b'),('a','c'),('a','d'),('a','e'),
 ('b','c'),('b','d'),('b','e'),
 ('c','d'),('c','e'),
 ('d','e'),
 ('f','b'),('f','c'),('f','g'),
 ('g','f'),('g','c'),('g','d'),('g','e')]
pos = {
'b': (0,4), 'f': (0,0),
'c': (2,2),
'a': (4,6),
'd': (6,2),
'e': (8,4), 'g': (8,0)
}

G.add_edges_from(edges_fig_4)

#pl.figure()
#nx.draw(G, pos=pos)
#pl.show()

# sorted list of nodes in graph
nodes_sorted = sorted(G.nodes())

# starting point: build all 2-clique sublists
clique_sublists = []
for a_node_index, a_node in enumerate(nodes_sorted):
    clique_sublist = {}
    clique_sublist['sublist_base'] = clique_sublist['sb'] = tuple(a_node)
    clique_sublist['common_neighbors'] = clique_sublist['cn'] = greater_neighbors_strings(G, a_node)
    clique_sublists.append(clique_sublist)

# now, for each 2-clique sublist, build all constituent 2-cliques in turn
# for each 2-clique you build, check if you can combine that 2-clique with another 2-clique from the same 2-clique sublist to form a 3-clique

# if you can combine successfuly then you found a 3-clique which is either:
# (i) non-maximal: combining the former 2-clique with the latter 2-clique results in a 3-clique sublist with non-empty common neighbors list
# (ii) maximal: the resultant 3-clique sublist has empty common neighbours list

#for a_sublist in clique_sublists:
while clique_sublists:
    a_sublist = clique_sublists.pop(0)
#a_sublist = clique_sublists[0]
    for node_added in a_sublist['cn_strings']:
#node_added = a_sublist['cn_strings'][0]
        neighbors_of_node_added = greater_neighbors_strings(G, node_added)

        current_sublist_base = a_sublist['sb']+tuple(node_added)
        current_sublist_cn = tuple(sorted(set(neighbors_of_node_added).intersection(a_sublist['cn_strings'])))

        print 'clique: '+str(current_sublist_base)

        for node in current_sublist_cn:
            new_sublist_base = current_sublist_base+tuple(node)
            new_sublist_cn = tuple(sorted(set(current_sublist_cn).intersection(greater_neighbors_strings(G, node))))
    
            if len(new_sublist_cn) == 0:
                print 'clique: '+str(new_sublist_base)
            elif len(new_sublist_cn) == 1:
                print 'clique: '+str(new_sublist_base)
                print 'clique: '+str(new_sublist_base+new_sublist_cn)
            else:
                print 'candidate sublist: '+str([new_sublist_base, new_sublist_cn])
                clique_sublists.append({'sb': new_sublist_base, 'cn_strings': new_sublist_cn})
                
