def nuc_check(read, rec_pos):
    """ The nuc\_check function checks the variant-record position against the supplemented read, and then extracts
    the nucleotide belonging to this position in the read. Returns the nucleotide in the read mapping against
    the variant-record position

    Args:
        :param read: Input read
        :param rec_pos: The position of the called variant in the reference genome

    Returns:
        :return: Returns the nucleotide in the read mapping against the variant-record position

    Raises:
        :raises ValueError: If a ValueError is found, the function returns nothing
    """
    try:
        # Gets the positions the sequence maps to in the reference
        # Full length with soft clips is required for the index selection to be correct, nucleotides are always returned
        # as they would be on the plus strand
        # for any reverse strand
        read_pos = read.get_reference_positions(full_length=True)
        # Gets the index of the position the sequence maps to
        ind_pos = read_pos.index(rec_pos)
        # Obtains the query sequence (the sequence as it were read), reverse complemented if in reverse strand
        read_seq = read.query_sequence
        read_seq = list(read_seq)
        # Gets the nucleotide present at the index position of the variant
        read_nuc = read_seq[ind_pos]
        return read_nuc
    except ValueError:
        pass


def ffpe_finder(cons_dict, var_nuc, ref_nuc, ffpe_n):
    """ The ffpe\_finder function is made to classify the variant type for paired UMI-reads. All-together the UMI and
    its variant-record position can be classified as: No mutation, Mutation, FFPE-artefact, Unknown (N) or Deletion (-).
    The function uses a dict of paired reads containing their consensus nucleotides categorized through their UMI-tag,
    which is then iterated through for every UMI. It then uses the consensus nucleotide originating from string 1 and
    string 2 to classify the UMI through comparing these to one another. If the two consensus nucleotides are equal
    to one another and furthermore equal to the nucleotide in the reference genome, the UMI is determined to be
    "Reference". If the two consensus nucleotides are equal to one another and furthermore equal to the variant in the
    variant-record, they are instead deemed to be a "True Variant". In default mode, a FFPE classification only occurs
    if there is a mismatch between the two consensus nucleotides, if one of the consensus nucleotides is equal to the
    variant in the variant record, and finally if the mismatch is of a C:T, T:C or G:A, A:G type. Alternatively, if the
    flag "ffpe_n" has been called with the input "all", the function instead classifies any mismatch between the
    consensus nucleotides for the positive and negative strand as a FFPE-artefact. If any of the consensus nucleotides
    are equal to N or -. the UMi is instead deemed to be "Unknown" or "Deleteion" respectively.
    After a UMI is classified, a counter is added to, and the UMi is stored within a dict named after the variant type.
    Once the algorithm has iterated through every UMI within the cons_dict, it creates a new dict containing all
    variant-type dicts as well as their molecular support, which is then returned

    Args:
        :param cons_dict: Dict containing the consensus nucleotides for String 1 and String 2,
        classified through their UMI-tag
        Example dict for a FFPE artefact:
        cons_dict = {Pos_Str_Hits: C, Neg_Str_Hits: T}
        :param var_nuc: The nucleotide called in the variant-record
        :param ref_nuc: The nucleotide found in the reference genome at the variant-call position
        :param ffpe_n: Parameter determining if all nucleotides should be included for FFPE-classification, or just
        C>T:G>A

    Returns:
        :return: Returns a dict with separate dicts for every possible variant type, which in turns contains the
        positive and negative strand nucleotide for this position. Also returns an index for the variant type.
        Example dict for a FFPE-artefact:
        var_dict = "Reference_Hits": {}, True_Variant_Hits": {}, "FFPE_Hits": {"Pos_Str": C, "Neg_Str": T},
        "N_Hits": {}, "Del_Hits": {}, "Reference_Support": 0, "True_Variant_Support": 0, "FFPE_Support": 1,
        "N_Support": 0, "Del_Support": 0}

    Raises:
        :raises KeyError: Raises a key-error if there are not two nucleotides to compare
    """
    var_dict = {}
    ffpe_dict = {}
    ffpe_ind = 0
    true_var_dict = {}
    true_var_ind = 0
    ref_dict = {}
    ref_ind = 0
    n_ind = 0
    n_dict = {}
    del_ind = 0
    del_dict = {}
    n_sym = "N"
    del_sym = "-"
    try:
        # Retrieves the positive and negative strand nucleotide hits for the key value
        pos_str_nuc = cons_dict["Pos_Str_Hits"]
        neg_str_nuc = cons_dict["Neg_Str_Hits"]

        if pos_str_nuc == ref_nuc and neg_str_nuc == pos_str_nuc:
            ref_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
            ref_ind += 1
        # If both positions contains the variant, it is added to the true_var_dict
        elif pos_str_nuc == var_nuc and neg_str_nuc == pos_str_nuc:
            true_var_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
            true_var_ind += 1
        # If either of the positions are equal to an N or a deletion, they are deemed as a N or Del respectively
        elif pos_str_nuc == n_sym or neg_str_nuc == n_sym:
            n_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
            n_ind += 1
        elif pos_str_nuc == del_sym or neg_str_nuc == del_sym:
            del_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
            del_ind += 1
        elif pos_str_nuc != neg_str_nuc:
            if ffpe_n == 'standard':
                # Checks the b_trans parameter, if standard only "FFPE" instances are classified as FFPE
                # Hydrolytic deamination causes C:G>T:A changes due to hydrolytic deamination of cytosine
                # PySAM's get_forward_sequence returns the reverse complement for any reverse reads

                # Deamination of C to U [U-G] > [U-A], [C-G] > [T-A], [C,G],
                # returned as [T-T],[C,C] as A is returned as T by PySAM, ie: T vs. C

                # Deamination of C to U [G-U] > [G-C], [A-U] > [G-C], [A-T],
                # returned as [G,G], [A-A] as T is returned as A by PySAM, ie G vs. A

                # Deamination of methylated C to T [C-G] > [T-G] > [T-A],[C-G],
                # returned as [T-T],[C-C], as A is returned as T by PySAM, ie: T vs.C

                # Deamination of methylated C to T [G-C] > [G-T] > [A-T],[G-C],
                # returned as [A-A],[G-G] as T is returned as A by PySAM, ie: G vs. A

                if pos_str_nuc == 'T' and neg_str_nuc == 'C' or \
                        pos_str_nuc == 'G' and neg_str_nuc == 'A':
                    # Checks to see if any of the nucleotides are equal to the variant call
                    if pos_str_nuc == var_nuc or neg_str_nuc == var_nuc:
                        ffpe_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
                        ffpe_ind += 1
                    else:
                        true_var_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
                        true_var_ind += 1
                else:
                    true_var_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
                    true_var_ind += 1
            elif ffpe_n == 'all':
                    ffpe_dict = {"Pos_Str": pos_str_nuc, "Neg_Str": neg_str_nuc}
                    ffpe_ind += 1
        var_dict = {"Reference_Hits": ref_dict, "True_Variant_Hits": true_var_dict, "FFPE_Hits": ffpe_dict,
                    "N_Hits": n_dict, "Del_Hits": del_dict, "Reference_Support": ref_ind,
                    "True_Variant_Support": true_var_ind, "FFPE_Support": ffpe_ind, "N_Support": n_ind,
                    "Del_Support": del_ind}
    except KeyError as e:
        print("No match for: " + str(e) + " found, comparison not possible")
    # Returns all dictionaries in a list along with their total counts
    return var_dict
