import numpy as np
import networkx as nx
from matplotlib import pylab as pl

G = nx.Graph()
G_edges = [
    ('A','B'),('A','D'),
    ('B','A'),('B','E'),('B','C'),
    ('C','B'),('C','G'),
    ('D','A'),('D','F'),('D','E'),
    ('E','F'),('E','D'),('E','B'),('E','G'),('E','F'),('E','H'),
    ('F','D'),('F','E'),
    ('G','C'),('G','E'),('G','I'),('G','H'),
    ('H','E'),('H','I'),('H','G'),
    ('I','E'),('I','G'),('I','H')
    ]

G.add_edges_from(G_edges)
#nx.draw(G)
#pl.show()

complete_subgraphs = []
for n in G.nodes():
    if set(n) not in complete_subgraphs:
        complete_subgraphs.append(set(n))
for e in G.edges():
    if set(e) not in complete_subgraphs:
        complete_subgraphs.append(set(e))
k = 3
k_cliques = list(nx.k_clique_communities(G,k))
while len(k_cliques) > 0:
    for k_clique in k_cliques:
        if set(k_clique) not in complete_subgraphs:
            complete_subgraphs.append(set(k_clique))
    k = k+1
    k_cliques = list(nx.k_clique_communities(G,k))
