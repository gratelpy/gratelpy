import networkx as nx
from matplotlib import pylab as pl
from graph import greater_neighbors, get_all_cliques

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

for clique in get_all_cliques(G):
    print clique
