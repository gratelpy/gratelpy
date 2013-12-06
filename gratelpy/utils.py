def invert_dict(a_dict):
    rev_dict = {}
    for key in a_dict.keys():
        rev_dict[a_dict[key]] = key
    return rev_dict
