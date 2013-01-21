import networkx as nx
from networkx_modified import get_all_cliques
from matplotlib import pylab as pl

G = nx.Graph()

G.add_edge('A','B')
G.add_edge('A','C')
G.add_edge('A','D')
G.add_edge('D','B')
G.add_edge('D','C')
G.add_edge('C','B')

for cl in get_all_cliques(G):
    print cl

pl.figure()
nx.draw(G)
pl.show()
