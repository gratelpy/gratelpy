from stoich import get_graph_stoich
from parse_mechanism import get_network_from_mechanism
from drawing import gratelpy_draw
from matplotlib import pyplot as plt
from fragments import get_sensible_fragments
from subgraphs import get_all_valid_subgraphs

#mechanism_file = '/usr/users/cbu/waltherg/JIC/Dropbox/GratelPy/mechanisms/reversible_substrate_inhibition.txt'
mechanism_file = '/home/waltherg/Dropbox/GratelPy/mechanisms/glycolysis_mechanism.txt'

base_name = 'glycolysis'

alpha, beta, dict_complexes, dict_constants, dict_complexes_reverse, dict_constants_reverse = get_network_from_mechanism(mechanism_file, 7)
print dict_complexes
print dict_constants
G, stoich, stoich_rank = get_graph_stoich(alpha, beta)

dx = 1
pos_mech = {
    'k10': [1*dx,7], '[A4]': [4*dx,7], '[A2]': [8*dx,7],
    'k4': [4*dx,6], 'k3': [6*dx,6],
    'k7': [0*dx,5], '[A3]': [1*dx,5], 'k6': [2*dx,5], '[A6]': [3*dx,5], '[A1]': [5*dx,5], 'k2': [9*dx,5],
    'k1': [6*dx,4],
    '[A5]': [3*dx,3], 'k5': [4*dx,4], '[A7]': [8*dx,3],
    'k8': [6*dx,2], 
    'k9': [9*dx,0]
    }

pos = {}
for s in dict_complexes:
    pos['s'+str(dict_complexes[s]+1)] = pos_mech[s]
for w in dict_constants:
    pos['w'+str(dict_constants[w]+1)] = pos_mech[w]

fig = gratelpy_draw(G, positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse, rnsize=3600, cnsize=2000)
plt.draw()
plt.savefig(base_name+'_graph.pdf', bbox_inches='tight')

# generate rank 2 critical fragment
stoich_rank = 2

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
plt.savefig(base_name+'_critical_fragment_'+str(stoich_rank)+'.pdf', bbox_inches='tight')

for sg_i, (sg_unfurled, sg) in enumerate(zip(subgraphs_unfurled, subgraphs)):
    gratelpy_draw(G.subgraph(sg_unfurled), subgraph=sg, positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
    plt.draw()
    plt.savefig(base_name+'_critical_fragment_sg_'+str(sg_i)+'_'+str(stoich_rank)+'.pdf', bbox_inches='tight')

# rank 3
stoich_rank = 3

valid_fragments = get_sensible_fragments(G, stoich_rank)
critical_fragments = []
for frag in valid_fragments:
    res = get_all_valid_subgraphs(G, stoich_rank, frag)
    if res[-1] < 0:
        critical_fragments.append(res)
print 'number of critical fragments:',str(len(critical_fragments))
for f_i, f in enumerate(critical_fragments):
    print f_i, [dict_complexes_reverse[int(s[1:])-1] for s in f[0][0]], [dict_constants_reverse[int(w[1:])-1] for w in f[0][1]]

frag_i = 4
frag_unfurled = [n for li in critical_fragments[frag_i][0] for n in li]
subgraphs = critical_fragments[frag_i][2]
subgraphs_unfurled = []
for sg in subgraphs:
    all_nodes = set(n for el in sg for n in el)
    subgraphs_unfurled.append(list(all_nodes))

gratelpy_draw(G.subgraph(frag_unfurled), positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
plt.draw()
plt.savefig(base_name+'_critical_fragment_'+str(frag_i)+'_'+str(stoich_rank)+'.pdf', bbox_inches='tight')

for sg_i, (sg_unfurled, sg) in enumerate(zip(subgraphs_unfurled, subgraphs)):
    gratelpy_draw(G.subgraph(sg_unfurled), subgraph=sg, positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
    plt.draw()
    plt.savefig(base_name+'_critical_fragment_'+str(frag_i)+'_sg_'+str(sg_i)+'_'+str(stoich_rank)+'.pdf', bbox_inches='tight')

frag_i = 7
frag_unfurled = [n for li in critical_fragments[frag_i][0] for n in li]
subgraphs = critical_fragments[frag_i][2]
subgraphs_unfurled = []
for sg in subgraphs:
    all_nodes = set(n for el in sg for n in el)
    subgraphs_unfurled.append(list(all_nodes))

gratelpy_draw(G.subgraph(frag_unfurled), positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
plt.draw()
plt.savefig(base_name+'_critical_fragment_'+str(frag_i)+'_'+str(stoich_rank)+'.pdf', bbox_inches='tight')

for sg_i, (sg_unfurled, sg) in enumerate(zip(subgraphs_unfurled, subgraphs)):
    gratelpy_draw(G.subgraph(sg_unfurled), subgraph=sg, positions=pos, dictionary_complexes=dict_complexes_reverse, dictionary_reactions=dict_constants_reverse)
    plt.draw()
    plt.savefig(base_name+'_critical_fragment_'+str(frag_i)+'_sg_'+str(sg_i)+'_'+str(stoich_rank)+'.pdf', bbox_inches='tight')

plt.show()
