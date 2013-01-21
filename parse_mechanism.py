import re

import numpy as np

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
    alpha_transposed = np.array(alpha_transposed)
    beta_transposed = np.array(beta_transposed)
    alpha = np.transpose(alpha_transposed)
    beta = np.transpose(beta_transposed)
    
    return alpha, beta, complexes, constants, complexes_reverse, constants_reverse


