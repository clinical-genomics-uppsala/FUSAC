def base_check(read, record_pos):
    try:
        # Gets the positions the sequence maps to in the reference
        # Full length with soft clips is required for the index selection to be correct
        read_pos = read.get_reference_positions(full_length=True)
        # Gets the index of the position the sequence maps to
        ind_pos = read_pos.index(record_pos)
        # Obtains the query sequence (the sequence as it were read)
        read_seq = read.query_sequence
        read_seq = list(read_seq)
        # Gets the base present at the index position of the variant
        read_base = read_seq[ind_pos]
        return read_base
    except ValueError:
        pass


def base_count(read_base):
    base_lst = ["A", "T", "G", "C", "N", "-"]
    base_cl = [0]*6
    # Checks the read_base against the base_lst and increases the corresponding counter if a match is found
    for base in base_lst:
        if base == read_base:
            base_cl[base_lst.index(base)] += 1
    # Returns a dict with the counts for each entry in the base_lst
    pos_dict = {"A": base_cl[0], "T": base_cl[1], "G": base_cl[2], "C": base_cl[3], "N": base_cl[4], "-": base_cl[5]}
    return pos_dict


def ffpe_finder(base_res, ref_var, ref_base):
    n_sym = "N"
    del_sym = "-"
    ffpe_dict = {}
    ffpe_ind = 0
    mut_dict = {}
    mut_ind = 0
    ref_dict = {}
    ref_ind = 0
    n_ind = 0
    n_dict = {}
    del_ind = 0
    del_dict = {}

    try:
        # Retrieves the forward and reverse base hits for the key value
        f_hits = base_res["String_1_Hits"]
        r_hits = base_res["String_2_Hits"]
        # Iterates through each variant called by the program
        # First checks if the variant is present in both forward and reverse molecule, if so it is a Mutation
        # If both positions are deemed to be the reference, it is added to the ref_dict
        if f_hits == r_hits and f_hits == ref_base:
            ref_dict = {"String_1": f_hits, "String_2": r_hits}
            ref_ind += 1
        # If both positions contains the variant, it is added to the mut_dict
        elif f_hits == r_hits and f_hits == ref_var:
            mut_dict = {"String_1": f_hits, "String_2": r_hits}
            mut_ind += 1
        # If either of the positions are equal to an N or a deletion, they are deemed as a N or Del respectively
        elif f_hits == n_sym or r_hits == n_sym:
            n_dict = {"String_1": f_hits, "String_2": r_hits}
            n_ind += 1
        elif f_hits == del_sym or r_hits == del_sym:
            del_dict = {"String_1": f_hits, "String_2": r_hits}
            del_ind += 1
        elif f_hits != n_sym and f_hits != del_sym and r_hits != n_sym and r_hits != del_sym and f_hits != r_hits:
            ffpe_dict = {"String_1": f_hits, "String_2": r_hits}
            ffpe_ind += 1
    except KeyError as e:
        print("No match for: " + str(e) + " found, comparison not possible")
    # Returns all dictionaries in a list along with their total counts
    var_dict = {"Reference_Hits": ref_dict, "Mutation_Hits": mut_dict, "FFPE_Hits": ffpe_dict,
                "N_Hits": n_dict, "Del_Hits": del_dict, "Reference_Support": ref_ind, "Mutation_Support": mut_ind,
                "FFPE_Support": ffpe_ind, "N_Support": n_ind, "Del_Support": del_ind}
    return var_dict
