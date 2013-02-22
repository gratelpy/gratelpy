from stoich import get_graph_stoich
from parse_mechanism import get_network_from_mechanism
from drawing import gratelpy_draw
from matplotlib import pyplot as plt
from fragments import get_sensible_fragments
from subgraphs import get_all_valid_subgraphs

#mechanism_file = '/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/mechanisms/reversible_substrate_inhibition.txt'
mechanism_file = '/home/waltherg/Dropbox/GratelPy/mechanisms/reversible_substrate_inhibition.txt'

alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism(mechanism_file, 4)

G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

pos = {
    'w6': [0,0], 's4': [1,1], 's3': [3,0],
    'w5': [2,2], 'w4': [3,2],
    'w2': [0,4], 's1': [1,4], 'w3': [2,4], 's2': [3,4],
    'w1': [1,6]
    }

fig = gratelpy_draw(G, positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
plt.draw()
plt.savefig('example_4_graph.pdf', bbox_inches='tight')

valid_fragments = get_sensible_fragments(G, stoich_rank)
critical_fragments = []
for frag in valid_fragments:
    res = get_all_valid_subgraphs(G, stoich_rank, frag)
    if res[-1] < 0:
        critical_fragments.append(res)
print 'number of critical fragments:',str(len(critical_fragments))

frag_unfurled = [n for li in critical_fragments[0][0] for n in li]
subgraphs = critical_fragments[0][2]
subgraphs_unfurled = []
for sg in subgraphs:
    all_nodes = set(n for el in sg for n in el)
    subgraphs_unfurled.append(list(all_nodes))

gratelpy_draw(G.subgraph(frag_unfurled), positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
plt.draw()
plt.savefig('example_4_critical_fragment.pdf', bbox_inches='tight')

for sg_i, (sg_unfurled, sg) in enumerate(zip(subgraphs_unfurled, subgraphs)):
    gratelpy_draw(G.subgraph(sg_unfurled), subgraph=sg, positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
    plt.draw()
    plt.savefig('example_4_critical_fragment_sg_'+str(sg_i)+'.pdf', bbox_inches='tight')
plt.show()
