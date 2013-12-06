import networkx as nx
import linalg
import sys

def get_substance_adjacency(alpha, beta):
    # alpa:
    # beta:
    
    # number of reactions = number of columns of alpha
    no_rxn = len(alpha[0])
    # number of substance = number of rows of alpha
    no_sub = len(alpha)

    # check
    if no_rxn != len(beta) or no_sub != len(beta[0]):
        raise

    # substance adjacency matrix
    subs_adj = []

    # build subs_adj matrix row-by-row
    for sub_i in range(no_sub):
        subs_adj_row = [0 for i in range(no_sub)]
        for rxn_i in range(no_rxn):
            if alpha[sub_i][rxn_i] > 0:
                for adj_sub_i in range(no_sub):
                    if beta[adj_sub_i][rxn_i] > 0:
                        subs_adj_row[adj_sub_i] = 1
        subs_adj.append(subs_adj_row)

    # return adjacency matrix
    return subs_adj

def get_random_alpha_beta(no_complexes, no_reactions, no_times_complexes_tested, remove_empty_reactions = True):
    try:
        import numpy as np
    except ImportError:
        print 'This method requires NumPy. Please install NumPy first.'
        raise

    # build alpha
    alpha_as_list = []
    for s_i in range(no_complexes):
        s_alpha = [0 for w_i in range(no_reactions)]
        test_indices = np.random.randint(0, no_reactions-1, size=no_times_complexes_tested)
        for index in test_indices:
            if np.random.rand() <= 0.5:
                s_alpha[index] = 1
        alpha_as_list.append(s_alpha)
    alpha = np.array(alpha_as_list)

    # build beta
    beta_as_list = []
    for s_i in range(no_complexes):
        s_beta = [0 for w_i in range(no_reactions)]
        test_indices = np.random.randint(0, no_reactions-1, size=no_times_complexes_tested)
        for index in test_indices:
            if np.random.rand() <= 0.5:
                s_beta[index] = 1
        beta_as_list.append(s_beta)
    beta = np.array(beta_as_list)

    if remove_empty_reactions:
        reaction_indexes = []
        for rxn_i in range(no_reactions):
            if all(entries == 0 for entries in alpha[:,rxn_i]) and all(entries == 0 for entries in beta[:,rxn_i]):
                reaction_indexes.append(rxn_i)
        alpha = np.delete(alpha, np.s_[reaction_indexes], 1)
        beta = np.delete(beta, np.s_[reaction_indexes], 1)

    # make sure nothing funny happened
    assert alpha.shape == beta.shape

    # remove s1 -> s1 type of reactions
    reaction_indexes = []
    for rxn_i in range(alpha.shape[1]):
        if sum(alpha[:,rxn_i]) == 1 and sum(beta[:,rxn_i]) == 1:
            if all(a == b for a,b in zip(alpha[:,rxn_i], beta[:,rxn_i])):
                reaction_indexes.append(rxn_i)
    alpha = np.delete(alpha, np.s_[reaction_indexes], 1)
    beta = np.delete(beta, np.s_[reaction_indexes], 1)

    # make sure nothing funny happened
    assert alpha.shape == beta.shape

    # remove empty rows (unused complexes)
    complex_indexes = []
    for cmpl_i in range(alpha.shape[0]):
        if all(entry == 0 for entry in alpha[cmpl_i,:]) and all(entry == 0 for entry in beta[cmpl_i,:]):
            complex_indexes.append(cmpl_i)
    alpha = np.delete(alpha, np.s_[complex_indexes], 0)
    beta = np.delete(beta, np.s_[complex_indexes], 0)

    # make sure nothing funny happened
    assert alpha.shape == beta.shape

    # remove trimolcular and higher order reactions
    reaction_indexes = []
    for rxn_i in range(alpha.shape[1]):
        if sum(alpha[:,rxn_i]) >= 3:
            reaction_indexes.append(rxn_i)

    if min(alpha.shape) != 0:
        alpha = np.delete(alpha, np.s_[reaction_indexes], 1)
        beta = np.delete(beta, np.s_[reaction_indexes], 1)

    # make sure nothing funny happened
    assert alpha.shape == beta.shape

    return alpha, beta

def get_graph_stoich(alpha, beta, complex_names = None, constant_names = None):

    # number of reactions = number of columns of alpha
    no_rxn = len(alpha[0])
    # number of substance = number of rows of alpha
    no_sub = len(alpha)

    # stoichiometry matrix
    stoich = []
    for row_i in range(no_sub):
        stoich_row = []
        for col_i in range(no_rxn):
            stoich_row.append(beta[row_i][col_i] - alpha[row_i][col_i])
        stoich.append(stoich_row)

    # rank of stoichiometry matrix
    if no_sub < no_rxn:
        # need to transpose stoich for svd call
        stoich_tr = []
        for i in range(no_rxn):
            stoich_tr_row = []
            for j in range(no_sub):
                stoich_tr_row.append(stoich[j][i])
            stoich_tr.append(stoich_tr_row)
        matrix_U, sing_vals, matrix_V_t = linalg.svd(stoich_tr)
    else:
        matrix_U, sing_vals, matrix_V_t = linalg.svd(stoich)
    # tol computation lent from NumPy
    # https://github.com/numpy/numpy/blob/85b83e6938fa6f5176eaab8e8fd1652b27d53aa0/numpy/linalg/linalg.py#L1513
    tol = max(sing_vals) * max(no_rxn, no_sub) * sys.float_info.epsilon
    stoich_rank = sum([1 for s in sing_vals if s > tol])

    # make directed graph
    G = nx.DiGraph()

    # add nodes to graph: add all reactions 'w' and substances 's'
    if constant_names == None:
        G.add_nodes_from(['w'+str(n) for n in range(1,no_rxn+1)], bipartite=1)
    else:
        G.add_nodes_from([constant_names[n-1] for n in range(1,no_rxn+1)], bipartite=1)

    if complex_names == None:
        G.add_nodes_from(['s'+str(n) for n in range(1,no_sub+1)], bipartite=0)
    else:
        G.add_nodes_from([complex_names[n-1] for n in range(1,no_sub+1)], bipartite=0)

    # add ('s','w') edges, i.e. positive element in alpha
    for row_i in range(no_sub):
        curr_row = alpha[row_i]
        for col_i in range(no_rxn):
            curr_coeff = curr_row[col_i]
            if curr_coeff > 0:
                if constant_names == None and complex_names == None:
                    G.add_edge('s'+str(row_i+1), 'w'+str(col_i+1), coeff=curr_coeff)
                elif constant_names == None and complex_names != None:
                    G.add_edge(complex_names[row_i], 'w'+str(col_i+1), coeff=curr_coeff)
                elif constant_names != None and complex_names == None:
                    G.add_edge('s'+str(row_i+1), constant_names[col_i], coeff=curr_coeff)
                elif constant_names != None and complex_names != None:
                    G.add_edge(complex_names[row_i], constant_names[col_i], coeff=curr_coeff)
                else:
                    raise

    # add ('w','s') edges, i.e. positive element in beta
    for row_i in range(no_sub):
        curr_row = beta[row_i]
        for col_i in range(no_rxn):
            curr_coeff = curr_row[col_i]
            if curr_coeff > 0:
                if constant_names == None and complex_names == None:
                    G.add_edge('w'+str(col_i+1), 's'+str(row_i+1), coeff=curr_coeff)
                elif constant_names == None and complex_names != None:
                    G.add_edge('w'+str(col_i+1), complex_names[row_i], coeff=curr_coeff)
                elif constant_names != None and complex_names == None:
                    G.add_edge(constant_names[col_i], 's'+str(row_i+1), coeff=curr_coeff)
                elif constant_names != None and complex_names != None:
                    G.add_edge(constant_names[col_i], complex_names[row_i], coeff=curr_coeff)
                else:
                    raise

    return G, stoich, stoich_rank
