import networkx as nx
from matplotlib import pylab as pl

def neighbors_bits(G, node):
    nodes_sorted = sorted(G.nodes())
    
    neighbors = [0 for n in nodes_sorted]
    for index, n in enumerate(nodes_sorted):
        if n in G.neighbors(node):
            neighbors[index] = 1
    return neighbors

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
clique_sublists = {}
for a_node_index, a_node in enumerate(nodes_sorted):
    # note: here each neighbour of a_node is a common neighbour of the 1-clique 'a_node'
    #clique_sublists[(a_node)] = [neighbour for neighbour in sorted(G.neighbors(a_node)) if neighbour > a_node]
    common_neighbours = [0 for node in G.nodes()]
    neighbours_of_a_node = G.neighbors(a_node)

    for another_node_index, another_node in enumerate(nodes_sorted):
        if (a_node is not another_node) and (another_node in neighbours_of_a_node) and (another_node_index > a_node_index):
            common_neighbours[another_node_index] = 1

    clique_sublists[(a_node)] = common_neighbours
            

# now, for each 2-clique sublist, build all constituent 2-cliques in turn
# for each 2-clique you build, check if you can combine that 2-clique with another 2-clique from the same 2-clique sublist to form a 3-clique

# if you can combine successfuly then you found a 3-clique which is either:
# (i) non-maximal: combining the former 2-clique with the latter 2-clique results in a 3-clique sublist with non-empty common neighbors list
# (ii) maximal: the resultant 3-clique sublist has empty common neighbours list

# let's build 3-clique sublists
new_clique_sublists = {}
#for a_sublist_index, a_sublist in enumerate(clique_sublist.keys()):
#a_sublist = clique_sublists[0]
#common_neighbours_in_a_sublist = tuple([cn for cn ])

# look at first sublist only
a_sublist = clique_sublists.keys()[0]
nb = neighbors_bits(G, a_sublist)
ns = neighbors_bits_to_string(G,nb)

for k_clique_kth_node in ns[0:-1]:
    print 'adding ',k_clique_kth_node
    k_clique = [a_sublist, k_clique_kth_node]

    cn_k_clique = [0 for index in range(len(nodes_sorted))]
    for index in range(len(nodes_sorted)):
        if index > nodes_sorted.index(k_clique[-1]):
            cn_k_clique[index] = clique_sublists[a_sublist][index]&neighbors_bits(G,k_clique_kth_node)[index]

    print k_clique, cn_k_clique

    new_sublist = []
    
    # construct tentative (k+1)-cliques based on current k_clique and check if that tentative clique
    # has more common neighbors (i.e. nodes you can add in future) or if that (k+1)-clique is end of road
    for k_p_1_clique_k_p_1_node in ns[1:]:
        if k_p_1_clique_k_p_1_node in G.neighbors(k_clique[-1]):
            cn_k_p_1_clique = [0 for index in range(len(nodes_sorted))]
            for index in range(len(nodes_sorted)):
                if index > nodes_sorted.index(k_p_1_clique_k_p_1_node):
                    cn_k_p_1_clique[index] = cn_k_clique[index]&neighbors_bits(G, k_p_1_clique_k_p_1_node)[index]
                    print cn_k_p_1_clique
                    if any([bit==1 for bit in cn_k_p_1_clique]):
                        print 'adding node '+k_p_1_clique_k_p_1_node+' to sublist extending from '+str(k_clique)
                        
