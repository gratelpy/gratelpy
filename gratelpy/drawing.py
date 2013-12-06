import networkx as nx
from networkx.algorithms import bipartite

try:
    from matplotlib import pyplot as plt
    import matplotlib
except ImportError:
    print 'The drawing module requires matplotlib. Please install matplotlib.'
    raise

def gratelpy_draw(G, positions=None, dictionary_complexes=None, dictionary_reactions=None, filename=None, subgraph=None, rnsize = 1600, cnsize = 1400):
    # draws entire graph or subgraph (subgraph being a fragment)
    # squares for reaction nodes, circles for complex nodes
    # inscribes s, w (complex and reaction nodes) labels into nodes

    # positions dictionary expected of the form {'si': [x_pos, y_pos], 'sj': ...}
    # supplied dictionary_complexes expected of the form {i: 'descriptive name of s(i+1)'}
    # supplied dictionary_reactions expected of the form {i: 'descriptive name of w(i+1)'}
    # filename expected as string including prefix

    font = {'family' : 'sans-serif','sans-serif':['Helvetica'],
            'weight' : 'normal',
            'size'   : 26}

    matplotlib.rc('font', **font)

    # generate positions of nodes if not supplied
    if positions is None:
        positions = nx.spring_layout(G)

    # generate and modify figure
    fig = plt.figure(num=None, figsize=(20,10), dpi=80, facecolor='w', edgecolor='k')
    fig_axis = fig.add_subplot(111)
    fig_axis.axis('off')

    # generate separate graphs for both complexes and reactions so we can draw them differently
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

    # if drawing subgraph, then generate specifically edges that are to be displayed
    if subgraph is not None:
        edges_graph = nx.DiGraph()
        edges_graph.add_nodes_from(substance_nodes+reaction_nodes)
        for el in subgraph:
            if len(el) == 2:
                # edge
                edges_graph.add_edge(el[0], el[1])
            else:
                if el[-1] == 'p':
                    edges_graph.add_edge(el[0], el[1])
                    edges_graph.add_edge(el[1], el[2])
                elif el[-1] == 'n':
                    edges_graph.add_edge(el[0], el[1])
                    edges_graph.add_edge(el[2], el[1])
                else:
                    raise
    else:
        edges_graph = None

    # generate complex labels
    if dictionary_complexes is None:
        complex_labels = {}
        for n in substance_graph.nodes():
            complex_labels[n] = str(n)
    else:
        complex_labels = {}
        for n in substance_graph.nodes():
            complex_labels[n] = dictionary_complexes[int(n[1:])-1].translate(None, '[]')

    # generate reaction labels
    if dictionary_reactions is None:
        reaction_labels = {}
        for n in reaction_graph.nodes():
            reaction_labels[n] = str(n)
    else:
        reaction_labels = {}
        for n in reaction_graph.nodes():
            reaction_labels[n] = dictionary_reactions[int(n[1:])-1].translate(None, '[]')
    
    # draw substance and reaction nodes
    nx.draw_networkx_nodes(substance_graph, positions, node_shape='o', node_size=cnsize, node_color='white')
    nx.draw_networkx_nodes(reaction_graph, positions, node_shape='s', node_size=rnsize, node_color='white')

    if subgraph is None:
        nx.draw_networkx_edges(G, positions, width=2)
    else:
        nx.draw_networkx_edges(edges_graph, positions, width=2)

    nx.draw_networkx_labels(substance_graph, positions, complex_labels, font_size=26)
    nx.draw_networkx_labels(reaction_graph, positions, reaction_labels, font_size=26)

    return fig
