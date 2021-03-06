def invert_dict(a_dict):
    rev_dict = {}
    for key in a_dict.keys():
        rev_dict[a_dict[key]] = key
    return rev_dict

ks_index = -1
frag_i = 0
sc_i = 1
sg_i = 2

spec_i = 0
rxn_i = 1


def result_get_fragment(result):
    return result[0]


def result_get_ks(result):
    return result[-1]


def result_get_sc(result):
    return result[1]


def result_get_sg(result):
    return result[2]


def fragment_get_species(fragment):
    return fragment[0]


def fragment_get_reactions(fragment):
    return fragment[1]


def subgraph_get_species(sg):
    return [s[spec_i] for s in sg]


def subgraph_get_reactions(sg):
    return [s[rxn_i] for s in sg]


# index in alpha, beta, and stoichiometry matrices
def species_get_index(species):
    return int(species[1:])-1


def reaction_get_index(reaction):
    return int(reaction[1:])-1


def edge_get_species(edge):
    return edge[0]
