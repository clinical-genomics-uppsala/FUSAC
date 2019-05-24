def base_check(read, rec_pos):
    """ Checks the variant-record position against the supplemented read and extracts the base belonging to this
    postion in the read.
    Args:
        :param read: Input read
        :param rec_pos: Position of the variant in the reference genome
    Returns:
        :return: Returns the base found at read matching the position of the variant in the reference genome
    Raises:
        :raises ValueError: If a ValueError is found, the function returns nothing
    """
    try:
        # Gets the positions the sequence maps to in the reference
        # Full length with soft clips is required for the index selection to be correct, returns the reverse complement
        # for any reverse strand
        read_pos = read.get_reference_positions(full_length=True)
        # Gets the index of the position the sequence maps to
        ind_pos = read_pos.index(rec_pos)
        # Obtains the query sequence (the sequence as it were read), reverse complemented if in reverse strand
        read_seq = read.query_sequence
        read_seq = list(read_seq)
        # Gets the base present at the index position of the variant
        read_base = read_seq[ind_pos]
        return read_base
    except ValueError:
        pass


def ffpe_finder(base_dict, ref_var, ref_base, ffpe_b):
    """ Function for classying the UMI variant-type. Alltogether the UMI and its variantrecord position can be
    classifies as: No mutation, Mutation, FFPE-artefact, Unknown (N-symbol) or Deletion (-). The function takes the
    str1 and str2 base and compares these to one another to classify the position, and add to the appropriate index.
    Args:
        :param base_dict: Dict containing Str1 and Str2 with their representative bases
        :param ref_var: The variant-record variant base, this is always aligned to the positive strand
        :param ref_base: The variant-record reference genome base, always aligned to the  positive strand
        :param ffpe_b: Parameter determining if all bases should be included for FFPE-classification, or just
        C>T:G>A
    Returns:
        :return: Returns a dict with separate dicts for every possible variant type, which in turns contains the
        String_1 and String_2 base for this position. Also returns an index for the variant type.
    Raises:
        :raises KeyError: Raises a key-error if somehow, there are not two read-bases to compare
    """
    var_dict = {}
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
    n_sym = "N"
    del_sym = "-"
    try:
        # Retrieves the forward and reverse base hits for the key value
        str1_bas = base_dict["String_1_Hits"]
        str2_bas = base_dict["String_2_Hits"]

        if str1_bas == str2_bas and str1_bas == ref_base:
            ref_dict = {"String_1": str1_bas, "String_2": str2_bas}
            ref_ind += 1
        # If both positions contains the variant, it is added to the mut_dict
        elif str1_bas == str2_bas and str1_bas == ref_var:
            mut_dict = {"String_1": str1_bas, "String_2": str2_bas}
            mut_ind += 1
        # If either of the positions are equal to an N or a deletion, they are deemed as a N or Del respectively
        elif str1_bas == n_sym or str2_bas == n_sym:
            n_dict = {"String_1": str1_bas, "String_2": str2_bas}
            n_ind += 1
        elif str1_bas == del_sym or str2_bas == del_sym:
            del_dict = {"String_1": str1_bas, "String_2": str2_bas}
            del_ind += 1
        elif str1_bas != str2_bas:
            if ffpe_b == 'standard':
                # Checks the b_trans parameter, if standard only "FFPE" instances are classified as FFPE
                # In order:
                # Hydrolytic deamination causes a C>T or G>A change due to the loss of the amino group
                # PySAM's get_forward_sequence returns the reverse complement for any reverse reads
                # Deamination of C-G [C-C] > T-G [T-C] , returned as T-C, as G is returned as C by PySAM
                # Deamination of G-C [G-G] > G-T [G-A] , returned as G-A, as T is returned as A by PySAM
                # Deamination of C-G [C-C] > C-A [C-T], returned as C-T, as A is returned as T by PySAM
                # Deamination of G-C [G-G] > A-C [A-G], returned as A-G, as C is returned as G by PySAM
                if str1_bas == 'T' and str2_bas == 'C' or \
                        str1_bas == 'G' and str2_bas == 'A' or \
                        str1_bas == 'C' and str2_bas == 'T' or \
                        str1_bas == 'A' and str2_bas == 'G':
                    print("New UMI")
                    print("Ref: " + ref_base)
                    print("Var: " + ref_var)
                    print("str1: " + str1_bas)
                    print("str2: " + str2_bas)
                    # Checks to see if any of the bases are equal to the variant call
                    if str1_bas == ref_var or str2_bas == ref_var:
                        print("New FFPE")
                        print("Ref: " + ref_base)
                        print("Var: " + ref_var)
                        print("str1: " + str1_bas)
                        print("str2: " + str2_bas)
                        ffpe_dict = {"String_1": str1_bas, "String_2": str2_bas}
                        ffpe_ind += 1
                    else:
                        mut_dict = {"String_1": str1_bas, "String_2": str2_bas}
                        mut_ind += 1
                else:
                    mut_dict = {"String_1": str1_bas, "String_2": str2_bas}
                    mut_ind += 1
            elif ffpe_b == 'all':
                    ffpe_dict = {"String_1": str1_bas, "String_2": str2_bas}
                    ffpe_ind += 1
        var_dict = {"Reference_Hits": ref_dict, "Mutation_Hits": mut_dict, "FFPE_Hits": ffpe_dict,
                    "N_Hits": n_dict, "Del_Hits": del_dict, "Reference_Support": ref_ind, "Mutation_Support": mut_ind,
                    "FFPE_Support": ffpe_ind, "N_Support": n_ind, "Del_Support": del_ind}
    except KeyError as e:
        print("No match for: " + str(e) + " found, comparison not possible")
    # Returns all dictionaries in a list along with their total counts
    return var_dict