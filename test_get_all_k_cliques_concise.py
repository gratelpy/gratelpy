import networkx as nx
from matplotlib import pylab as pl

def greater_neighbors(G, a_node):
    nodes_sorted = sorted(G.nodes())
    a_node_index = nodes_sorted.index(a_node)

    neighbors_of_a_node = []

    for another_node_index, another_node in enumerate(nodes_sorted):
        if another_node_index > a_node_index and another_node in G.neighbors(a_node):
            neighbors_of_a_node.append(another_node)
    
    return tuple(neighbors_of_a_node)

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
    # sublist base, sb
    clique_sublist['sb'] = tuple(a_node)
    # common neighbors, cn
    clique_sublist['cn'] = greater_neighbors(G, a_node)
    clique_sublists.append(clique_sublist)


while clique_sublists:
    a_sublist = clique_sublists.pop(0)
    for node_added in a_sublist['cn']:
        neighbors_of_node_added = greater_neighbors(G, node_added)

        current_sublist_base = a_sublist['sb']+tuple(node_added)
        current_sublist_cn = tuple(sorted(set(neighbors_of_node_added).intersection(a_sublist['cn'])))

        print 'clique: '+str(current_sublist_base)

        for node in current_sublist_cn:
            new_sublist_base = current_sublist_base+tuple(node)
            new_sublist_cn = tuple(sorted(set(current_sublist_cn).intersection(greater_neighbors(G, node))))
    
            if len(new_sublist_cn) == 0:
                print 'clique: '+str(new_sublist_base)
            elif len(new_sublist_cn) == 1:
                print 'clique: '+str(new_sublist_base)
                print 'clique: '+str(new_sublist_base+new_sublist_cn)
            else:
                print 'candidate sublist: '+str([new_sublist_base, new_sublist_cn])
                clique_sublists.append({'sb': new_sublist_base, 'cn': new_sublist_cn})
