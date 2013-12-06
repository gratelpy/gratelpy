import re

def fgsl(v):
    return format_gsl(v)

def format_gsl(v):
    v = v.replace('[', '__')
    v = v.replace(']', '__')
    v = v.replace('.', '_')
    v = v.replace('{', '')
    v = v.replace('}', '')

    return v

def print_jac_gsl_footer(mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None):
    print_gsl_footer('jac', mechanism_file, alpha, beta, complex_dict, constant_dict, gsl_conservation_rules)

def print_ode_gsl_footer(mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None):
    print_gsl_footer('ode', mechanism_file, alpha, beta, complex_dict, constant_dict, gsl_conservation_rules)

def print_gsl_footer(footer_type, mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None):

    no_complexes, no_reactions = alpha.shape
    
    # check for conservation rules
    gsl_conserved = []
    if gsl_conservation_rules is not None:
        gsl_conserved = [a_rule[0] for a_rule in gsl_conservation_rules]
    no_conserved = len(gsl_conserved)

    if footer_type == 'jac':
        for complex_i in range(no_complexes-no_conserved):
            mechanism_file.write('\tdfdt['+str(complex_i)+'] = 0.0;\n')
    
    mechanism_file.write('\n\n')
    mechanism_file.write('\treturn GSL_SUCCESS;\n')
    mechanism_file.write('}')
    
def print_jac_gsl_header(mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None):
    print_gsl_header('jac', mechanism_file, alpha, beta, complex_dict, constant_dict, gsl_conservation_rules)

def print_ode_gsl_header(mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None):
    print_gsl_header('ode', mechanism_file, alpha, beta, complex_dict, constant_dict, gsl_conservation_rules)
 
def print_gsl_header(header_type, mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None):

    no_complexes, no_reactions = alpha.shape

    # check for conservation rules
    gsl_conserved = []
    if gsl_conservation_rules is not None:
        gsl_conserved = [a_rule[0] for a_rule in gsl_conservation_rules]
    no_conserved = len(gsl_conserved)

    gsl_header = ''
    if header_type == 'jac':
        gsl_header += 'int jac (double t, const double y[], double *dfdy, double dfdt[], void *void_params){\n'
        gsl_header += '\tdouble * params = (double*)void_params;\n'
        gsl_header += '\tgsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, '+str(no_complexes-no_conserved)+','+str(no_complexes-no_conserved)+');\n'
        gsl_header += '\tgsl_matrix * m = &dfdy_mat.matrix;\n'
    elif header_type == 'ode':
        gsl_header += 'int func (double t, const double y[], double f[], void *void_params){\n'
        gsl_header += '\tdouble * params = (double*)void_params;\n'
    else:
        raise

    complex_print_i = 0
    for complex_i in range(no_complexes):
        if complex_dict is not None:
            if complex_dict[complex_i] not in gsl_conserved:
                gsl_header += '\tdouble '+fgsl(complex_dict[complex_i])+' = y['+str(complex_print_i)+'];\n'
                complex_print_i += 1
        else:
            if 's'+str(complex_i+1) not in gsl_conserved:
                gsl_header += '\tdouble s'+str(complex_print_i+1)+' = y['+str(complex_print_i)+'];\n'
                complex_print_i += 1

    constants_printed = []
    for reaction_i in range(no_reactions):
        if constant_dict is not None:
            if fgsl(constant_dict[reaction_i]) not in constants_printed:
                gsl_header += '\tdouble '+fgsl(constant_dict[reaction_i])+' = params['+str(reaction_i)+'];\n'
                constants_printed.append(fgsl(constant_dict[reaction_i]))
        else:
            if 'k'+str(reaction_i+1) not in constants_printed:
                gsl_header += '\tdouble k'+str(reaction_i+1)+' = params['+str(reaction_i)+'];\n'
                constants_printed.append('k'+str(reaction_i+1))
    no_constants_printed = len(constants_printed)

    # conservation rules if any
    if gsl_conservation_rules is not None:
        gsl_header += '\n'
        for cons_rule_i in range(len(gsl_conservation_rules)):
            gsl_header += '\tdouble '+fgsl(gsl_conservation_rules[cons_rule_i][1])+' = params['+str(no_constants_printed+cons_rule_i)+'];\n'
        gsl_header += '\n'
        for cons_rule_i in range(len(gsl_conservation_rules)):
            gsl_header += '\tdouble '+fgsl(gsl_conservation_rules[cons_rule_i][0])+' = '+fgsl(gsl_conservation_rules[cons_rule_i][1])
            for el in gsl_conservation_rules[cons_rule_i][2:]:
                gsl_header += '-'+fgsl(el)

            gsl_header += ';\n'


    gsl_header += '\n\n'

    mechanism_file.write(gsl_header)


def print_jac_from_alpha_beta(basename, alpha, beta, complex_dict=None, constant_dict=None, gsl=False, gsl_conservation_rules = None):
    # conservation rules expected as
    # [variable_to_be_replaced, new_constant, independent_variable_1, independent_variable_2, ...]
    # variable_to_be_replaced = new_constant - independent_variable_1 - independent_variable_2 - ...
    
    # check if we were passed a 'reverse' dict and reverse it if needed
    if complex_dict is not None:
        if type(complex_dict.keys()[0]) is not type(int()):
            complex_dict = {v:k for k,v in complex_dict.iteritems()}

    if constant_dict is not None:
        if type(constant_dict.keys()[0]) is not type(int()):
            constant_dict = {v:k for k,v in constant_dict.iteritems()}

    # open file
    mechanism_file_name = basename+'.jac'
    mechanism_file = open(mechanism_file_name, 'wb')

    # no complexes and reactions
    no_complexes, no_reactions = alpha.shape
    if (no_complexes, no_reactions) != tuple(beta.shape):
        raise

    # check for conservation rules
    gsl_conserved = []
    if gsl_conservation_rules is not None:
        gsl_conserved = [a_rule[0] for a_rule in gsl_conservation_rules]
    
    # net stoichiometric matrix
    gamma = beta - alpha

    # if printing GSL matrix, print header now
    if gsl:
        print_jac_gsl_header(mechanism_file, alpha, beta, complex_dict=complex_dict, constant_dict=constant_dict, gsl_conservation_rules=gsl_conservation_rules)
    
    jac_row_print_i = 0
    for jac_row_i in range(no_complexes):
        jac_column_print_i = 0
        for jac_col_i in range(no_complexes):
            # print right-hand side for each complex
            jac_entry = ''
            for rxn_i in range(no_reactions):
                # rxn_i relevant for rhs of complex jac_row_i iff jac_row_i either produced or consumed in that reaction
                if alpha[jac_col_i,rxn_i] != 1:
                    jac_entry_complexes = '0.0'
                    
                else:
                    jac_entry_complexes = ''

                    if alpha[jac_row_i, rxn_i] == 1 or beta[jac_row_i, rxn_i] == 1:
                        for cmp_i in range(no_complexes):
                            if complex_dict is None:
                                if alpha[cmp_i,rxn_i] == 1 and alpha[jac_col_i,rxn_i] == 1 and cmp_i != jac_col_i:
                                    jac_entry_complexes += '*[s'+str(cmp_i+1)+']'#+' + '
                                elif alpha[cmp_i,rxn_i] > 1:
                                    raise Exception('stoichiometric coefficients > 1 not implemented yet')
                            else:
                                if alpha[cmp_i,rxn_i] == 1 and alpha[jac_col_i,rxn_i] == 1 and cmp_i != jac_col_i:
                                    jac_entry_complexes += '*' +complex_dict[cmp_i]# + ' + '
                                elif alpha[cmp_i,rxn_i] > 1:
                                    raise Exception('stoichiometric coefficients > 1 not implemented yet')

                if constant_dict is None:
                    if jac_entry_complexes == '0.0':
                        pass
                    elif gamma[jac_row_i, rxn_i] == 0:
                        pass
                    elif gamma[jac_row_i, rxn_i] == 1:
                        jac_entry += '+' + 'k'+str(rxn_i+1) + jac_entry_complexes
                    elif gamma[jac_row_i, rxn_i] == -1:
                        jac_entry += '-' + 'k'+str(rxn_i+1) + jac_entry_complexes
                else:
                    if jac_entry_complexes == '0.0':
                        pass
                    elif gamma[jac_row_i, rxn_i] == 0:
                        pass
                    elif gamma[jac_row_i, rxn_i] == 1:
                        jac_entry += '+' + constant_dict[rxn_i] + jac_entry_complexes
                    elif gamma[jac_row_i, rxn_i] == -1:
                        jac_entry += '-' + constant_dict[rxn_i] + jac_entry_complexes
         

            if gsl:
                if (complex_dict[jac_row_i] not in gsl_conserved if complex_dict is not None else 's'+str(jac_row_i+1) not in gsl_conserved) and (complex_dict[jac_col_i] not in gsl_conserved if complex_dict is not None else 's'+str(jac_col_i+1) not in gsl_conserved):
                    if len(jac_entry)==0:
                        mechanism_file.write('\tgsl_matrix_set(m, '+str(jac_row_print_i)+', '+str(jac_column_print_i)+', 0.0);')
                        mechanism_file.write('\n')
                    else:
                        mechanism_file.write('\tgsl_matrix_set(m, '+str(jac_row_print_i)+', '+str(jac_column_print_i)+', '+ fgsl(jac_entry)+');')
                        mechanism_file.write('\n')
                        
                    jac_column_print_i += 1
                    if jac_column_print_i == no_complexes - len(gsl_conserved):
                        jac_row_print_i += 1
            else:
                if len(jac_entry)==0:
                    mechanism_file.write('jac['+str(jac_row_i)+']['+str(jac_col_i)+'] = 0.0')
                    mechanism_file.write('\n')
                else:
                    mechanism_file.write('jac['+str(jac_row_i)+']['+str(jac_col_i)+'] = '+ jac_entry)
                    mechanism_file.write('\n')

            

    # if printing GSL matrix, print footer now
    if gsl:
        print_jac_gsl_footer(mechanism_file, alpha, beta, complex_dict=complex_dict, constant_dict=constant_dict, gsl_conservation_rules=gsl_conservation_rules)

    mechanism_file.close()

    
def print_ode_from_alpha_beta(basename, alpha, beta, complex_dict=None, constant_dict=None, gsl=False, gsl_conservation_rules=None):
    # conservation rules expected as
    # [variable_to_be_replaced, new_constant, independent_variable_1, independent_variable_2, ...]
    # variable_to_be_replaced = new_constant - independent_variable_1 - independent_variable_2 - ...
    
    # check if we were passed a 'reverse' dict and reverse it if needed
    if complex_dict is not None:
        if type(complex_dict.keys()[0]) is not type(int()):
            complex_dict = {v:k for k,v in complex_dict.iteritems()}

    if constant_dict is not None:
        if type(constant_dict.keys()[0]) is not type(int()):
            constant_dict = {v:k for k,v in constant_dict.iteritems()}
    
    # open file
    mechanism_file_name = basename+'.ode'
    mechanism_file = open(mechanism_file_name, 'wb')

    # number of complexes and number of reactions
    no_complexes, no_reactions = alpha.shape
    if (no_complexes, no_reactions) != tuple(beta.shape):
        raise

    # check for conservation rules
    gsl_conserved = []
    if gsl_conservation_rules is not None:
        gsl_conserved = [a_rule[0] for a_rule in gsl_conservation_rules]
    
    # net stoichiometric matrix
    gamma = beta - alpha
    
    # if printing GSL ode, print header now
    if gsl:
        print_ode_gsl_header(mechanism_file, alpha, beta, complex_dict=complex_dict, constant_dict=constant_dict, gsl_conservation_rules=gsl_conservation_rules)
    
    ode_row_print_i = 0
    for ode_cmp_i in range(no_complexes):
        # print right-hand side for each complex
        rhs = ''
        for rxn_i in range(no_reactions):
            # rxn_i relevant for rhs of complex ode_cmp_i iff ode_cmp_i either produced or consumed in that reaction
            reactant_side = ''
            if alpha[ode_cmp_i, rxn_i] == 1 or beta[ode_cmp_i, rxn_i] == 1:
                for cmp_i in range(no_complexes):
                    if complex_dict is None:
                        if alpha[cmp_i,rxn_i] == 1:
                            reactant_side = reactant_side + '*[s'+str(cmp_i+1)+']'#+' + '
                        elif alpha[cmp_i,rxn_i] > 1:
                            raise Exception('stoichiometric coefficients > 1 not implemented yet')
                    else:
                        if alpha[cmp_i,rxn_i] == 1:
                            reactant_side = reactant_side + '*' +complex_dict[cmp_i]# + ' + '
                        elif alpha[cmp_i,rxn_i] > 1:
                            raise Exception('stoichiometric coefficients > 1 not implemented yet')
                
                # prune last '+'
                # reactant_side = reactant_side.rstrip()
                # reactant_side = reactant_side.rstrip('+')
                # reactant_side = reactant_side.rstrip()

            if constant_dict is None:
                if gamma[ode_cmp_i, rxn_i] == 0:
                    pass
                elif gamma[ode_cmp_i, rxn_i] == 1:
                    rhs += '+' + 'k'+str(rxn_i+1) + reactant_side
                elif gamma[ode_cmp_i, rxn_i] == -1:
                    rhs += '-' + 'k'+str(rxn_i+1) + reactant_side
            else:
                if gamma[ode_cmp_i, rxn_i] == 0:
                    pass
                elif gamma[ode_cmp_i, rxn_i] == 1:
                    rhs += '+' + constant_dict[rxn_i] + reactant_side
                elif gamma[ode_cmp_i, rxn_i] == -1:
                    rhs += '-' + constant_dict[rxn_i] + reactant_side
         
        # write kinetics of cmp_i to file
        if gsl:
            if (complex_dict[ode_cmp_i] not in gsl_conserved if complex_dict is not None else 's'+str(ode_cmp_i+1) not in gsl_conserved):
                mechanism_file.write('\tf['+str(ode_row_print_i)+'] = '+fgsl(rhs)+';\n')
                ode_row_print_i += 1
        else:
            if complex_dict is None:
                mechanism_file.write('[s'+str(ode_cmp_i+1)+']\' = '+rhs)
                mechanism_file.write('\n')
            else:
                mechanism_file.write(complex_dict[ode_cmp_i] + '\' = ' +rhs)
                mechanism_file.write('\n')

    if gsl:
        print_ode_gsl_footer(mechanism_file, alpha, beta, complex_dict=None, constant_dict=None, gsl_conservation_rules=None)

    mechanism_file.close()

def print_mechanism_from_alpha_beta(basename, alpha, beta):
    # open file
    mechanism_file_name = basename+'.mechanism'
    mechanism_file = open(mechanism_file_name, 'wb')

    no_complexes, no_reactions = alpha.shape
    if (no_complexes, no_reactions) != tuple(beta.shape):
        raise
    
    for rxn_i in range(no_reactions):
        # print reactants
        reactant_side = ''
        for cmp_i in range(no_complexes):
            if alpha[cmp_i,rxn_i] == 1:
                reactant_side = reactant_side +'[s'+str(cmp_i+1)+']'+' + '
            elif alpha[cmp_i,rxn_i] > 1:
                raise Exception('stoichiometric coefficients > 1 not implemented yet')
        # prune last '+'
        reactant_side = reactant_side.rstrip()
        reactant_side = reactant_side.rstrip('+')
        reactant_side = reactant_side.rstrip()

        # write to file
        mechanism_file.write(reactant_side)

        # add '->'
        mechanism_file.write(' -> ')

        # print products
        product_side = ''
        for cmp_i in range(no_complexes):
            if beta[cmp_i,rxn_i] == 1:
                product_side = product_side +'[s'+str(cmp_i+1)+']'+' + '
            elif beta[cmp_i,rxn_i] > 1:
                raise Exception('stoichiometric coefficients > 1 not implemented yet')
        # prune last '+'
        product_side = product_side.rstrip()
        product_side = product_side.rstrip('+')
        product_side = product_side.rstrip()

        # write to file
        mechanism_file.write(product_side)

        # add constant
        mechanism_file.write('; k'+str(rxn_i+1))

        # add newline
        mechanism_file.write('\n')
    
    # close file
    mechanism_file.close()
            
def get_network_from_mechanism(filename, no_complexes):
    # reads mechanism file: each line is of the form alpha_1 [A] -> beta_1 [B]; k_rat
    # NOTE: mechanism file is read case-sensitively

    # returns dictionaries of all detected complexes and rate constants as well as alpha and beta matrices

    # open file
    mechanism_file = open(filename)

    # go through file line-by-line and create columns of alpha and beta
    complexes = {}
    complexes_reverse = {}
    constants = {}
    constants_reverse = {}
    alpha_transposed = []
    beta_transposed = []

    # count number of lines = number of rate constants
    no_constants = len([1 for el in enumerate(mechanism_file)])
    mechanism_file.seek(0) # go back to beginning of file

    for line in mechanism_file:
        # tidy up current line
        curr_line = line.strip() # removes trailing and leading whitespaces
        curr_line = curr_line.replace(' ', '') # removes all remaining whitespaces on line

        # check if this line is a comment
        if curr_line[0] == '#':
            no_constants = no_constants - 1 # no_constants is taken as number of non-comment lines in mechanism file
            continue

        # get left, reactant side
        curr_left = curr_line.split('->')[0]
        # get right, product side
        curr_right = curr_line.split('->')[1]
        curr_right_split = curr_right.split(';')
        curr_right = curr_right_split[0]
        # get current constant name and add to dictionary
        curr_const = curr_right_split[1]
        if not curr_const in constants.keys():
            constants[curr_const] = len(constants)
            constants_reverse[len(constants_reverse)] = curr_const

        # create current alpha column and add to alpha_transposed
        curr_left = re.split('(\+)', curr_left)
        curr_left = [el for el in curr_left if len(el)>0] # remove empty elements from list
        alpha_column = [0 for el in range(no_complexes)]
        for el in curr_left:
            if el != '+':
                # detect alpha coefficient in 'el': 'curr_coeff'
                curr_coeff_match = re.findall('^\d+', el)
                curr_coeff = -1
                if len(curr_coeff_match) == 0:
                    curr_coeff = 1
                elif len(curr_coeff_match) == 1:
                    curr_coeff = int(curr_coeff_match[0])
                else:
                    raise Exception('too many coefficients detected: '+ el)
                    
                # detect complex in 'el': 'curr_complex'
                curr_complex_match = re.findall('\[{1}.*?\]{1}', el)
                curr_complex = 'not_a_complex'
                if len(curr_complex_match)==0:
                    raise Exception('no complex detected: '+ el)
                elif len(curr_complex_match)==1:
                    curr_complex = curr_complex_match[0]
                else:
                    raise Exception('more than one complex detected: '+ el)

                # add current complex to dictionary if not yet recorded
                if not curr_complex in complexes.keys():
                    complexes[curr_complex] = len(complexes)
                    complexes_reverse[len(complexes_reverse)] = curr_complex

                # modify alpha_column accordingly
                alpha_column[complexes[curr_complex]] = curr_coeff

        # add alpha_column to alpha_transposed
        alpha_transposed.append(alpha_column)

        # create current beta column and add to beta_transposed
        curr_right = re.split('(\+)', curr_right)
        curr_right = [el for el in curr_right if len(el)>0] # remove empty elements from list
        beta_column = [0 for el in range(no_complexes)]
        for el in curr_right:
            if el != '+':
                # detect beta coefficient in 'el': 'curr_coeff'
                curr_coeff_match = re.findall('^\d+', el)
                curr_coeff = -1
                if len(curr_coeff_match) == 0:
                    curr_coeff = 1
                elif len(curr_coeff_match) == 1:
                    curr_coeff = int(curr_coeff_match[0])
                else:
                    raise Exception('too many coefficients detected: '+ el)
                    
                # detect complex in 'el': 'curr_complex'
                curr_complex_match = re.findall('\[{1}.*?\]{1}', el)
                curr_complex = 'not_a_complex'
                if len(curr_complex_match)==0:
                    raise Exception('no complex detected: '+ el)
                elif len(curr_complex_match)==1:
                    curr_complex = curr_complex_match[0]
                else:
                    raise Exception('more than one complex detected: '+ el)

                # add current complex to dictionary if not yet recorded
                if not curr_complex in complexes.keys():
                    complexes[curr_complex] = len(complexes)
                    complexes_reverse[len(complexes_reverse)] = curr_complex

                # modify beta_column accordingly
                beta_column[complexes[curr_complex]] = curr_coeff

        # add alpha_column to alpha_transposed
        beta_transposed.append(beta_column)

    # close file
    mechanism_file.close()

    # check if we detected all complexes
    if len(complexes) != no_complexes:
        raise Exception('number of complexes indicated: '+str(no_complexes) + ' - complexes detected: '+str(complexes.keys()))
    # check if we detected all rate constants
    if len(constants) != no_constants:
        raise Exception('number of rows in mechanism file counted: '+str(no_constants) + ' - constants detected: '+str(constants.keys()))
    
    # create alpha and beta
    # we do not need NumPy to transpose these
    no_rows = len(alpha_transposed)
    no_cols = len(alpha_transposed[0])
    alpha = []
    beta = []
    for i in range(no_cols):
        alpha_row = []
        beta_row = []
        for j in range(no_rows):
            alpha_row.append(alpha_transposed[j][i])
            beta_row.append(beta_transposed[j][i])
        alpha.append(alpha_row)
        beta.append(beta_row)
    
    return alpha, beta, complexes, constants, complexes_reverse, constants_reverse


