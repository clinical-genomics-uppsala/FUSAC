import pos_function
import base_function
import count_function


def vcf_extract(record, bam_file, ffpe_b, ext_fun, spl_fun, spl_cha):
    """ Uses the supplemented variant-record to extract all reads in the BAM-file overlapping with its position. This
    newly generated list is used for the pos_checker function to return molecular data. The output from pos_checker is
    then subsequently used in the inf-builder function. Finally, the output from inf_builder is added to the
    input_record.

    Args:
        :param record: Variant-record of interest
        :param bam_file: Path to the BAM-file of interest
        :param ffpe_b: Parameter to determine if all mismatches should be classified as ffpe, or solely C:G>T:A
        :param ext_fun: Function for extracting the UMI-tag from a read
        :param spl_fun: Function used for splitting the UMI-tag in a read
        :param spl_cha: The character used for splitting the UMI-tag

    Returns:
        :return: Returns a copy of the variant-record modified by the inf_builder output. More specifically, adds a
        a list containing the support for each variant-type, as well as the support for the reference
        and variant call for str1 and str2 to the "samples" field. Furthermore, if any record has support for containing
        an FFPE-artefact, the "filter" tag will be modified to say "FFPE"
    """
    bam_lst = []
    n_ref = record.ref
    n_ref = ''.join(n_ref)
    list(n_ref)
    n_alt = record.alts
    n_alt = ''.join(n_alt)
    list(n_alt)
    # Checks so that the length of the list is not greater then 1 (temporary solution for handling SNVs only)
    if len(n_ref) > 1 or len(n_alt) > 1:
        return
    # Copies the record information
    n_cop = record.copy()
    n_pos = record.pos
    rec_chr = str(record.chrom)
    # The position that is returned to Python is 0 - based, NOT 1 - based as in the VCF file.
    n_pos = (n_pos - 1)

    # Use the record position to fetch all reads matching it, then append to the list
    for read in bam_file.fetch(rec_chr, n_pos, n_pos+1):
        bam_lst.append(read)

    # Calls the pos_checker function to obtain ffpe_data
    mpd_data, unmpd_data = var_extract(bam_lst, n_pos, n_alt, n_ref, ffpe_b, ext_fun, spl_fun, spl_cha)

    mpd_inf = inf_builder(mpd_data, n_alt, n_ref)
    unmpd_inf = inf_builder(unmpd_data, n_alt, n_ref)

    for sample in n_cop.samples:
        # n_cop.samples[sample]['UMI'] = str(mpd_inf[0][0]) + ";" + str(mpd_inf[0][1]) + ";" + str(mpd_inf[0][2]) + \
        #                                    ";" + str(mpd_inf[0][3]) + ";" + str(mpd_inf[0][4]) + ";" + mpd_inf[1] + \
        #                                    ";" + mpd_inf[2] + ";" + mpd_inf[3] + ";" + mpd_inf[4]
        n_cop.samples[sample]['UUMI'] = str(unmpd_inf[0][0]) + ";" + str(unmpd_inf[0][1]) + \
                                        ";" + str(unmpd_inf[0][2]) + ";" + str(unmpd_inf[0][3]) + \
                                        ";" + str(unmpd_inf[0][4]) + ";" + unmpd_inf[1] + ";" + unmpd_inf[2] \
                                        + ";" + unmpd_inf[3] + ";" + unmpd_inf[4]

        n_cop.samples[sample]['UMI'] = "{No_Mutation};{True_Mutation};{FFPE_Artefact};{Unknown};{Deletion};" \
                                       "{Ref_Paired};{Var_Paired};{Ref_Single};{Var_Single}"\
            .format(No_Mutation=mpd_inf[0][0], True_Mutation=mpd_inf[0][1], FFPE_Artefact=mpd_inf[0][2],
                    Unknown=mpd_inf[0][3], Deletion=mpd_inf[0][4], Ref_Paired=mpd_inf[1], Var_Paired=mpd_inf[2],
                    Ref_Single=mpd_inf[3], Var_Single=mpd_inf[4])

        # Checks if any record in the returned dict indicates an FFPE, if so updates the n_fil parameter
    for umi_key in mpd_data:
        if mpd_data[umi_key]["Mate_Hits"]:
            if mpd_data[umi_key]["Mate_Hits"]["FFPE_Hits"]:
                n_cop.filter.add("FFPE")
                break
    return n_cop


def var_extract(bam_lst, rec_pos, var_nuc, ref_nuc, ffpe_b, ext_fun, spl_fun, spl_cha):
    """ Function with the purpose of creating a dict based on the directionality and umi-tags of the supplemented
    reads in the bam_lst. Then using said dict to call the pos_hits and ffpe_finder functions to return a dict with
    data regarding positional data and variant types for the variant-record position and the reads aligning to it.
    Args:
        :param bam_lst: Input list of BAM-reads aligning to the variant call
        :param rec_pos: The position of the variant in the reference genome
        :param var_nuc: The nucleotide called in the variant-record
        :param ref_nuc: The nucleotide found in the reference genome at the variant-call position
        :param ffpe_b: Optional input argument controlling which mismatches to consider for FFPE-classification
        :param ext_fun: Function for extracting the UMI-tag from a read
        :param spl_fun: Function used for splitting the UMI-tag in a read
        :param spl_cha: Character used for splitting the UMI-tag

    Returns:
        :return: Returns a dict for mapped and unmapped reads. Each of these dicts containing a single-hits and a
        mate-hits dict. The mate-hits dict in turn contains data regarding if the variant is a mutation, no mutation,
        FFPE-artfefact deletion or N-call. Whereas the single-hits dict contains positional data for reads with no mate.

    Raises:
        :raises KeyError: If a umi is not found within the umi_dict, adds said umi_id to the umi_dict as well as two
        empty dicts for String_1 and String_2 for the umi
        :raises KeyError: Raises a key-error if the requested read/dict_key does not exist
    """
    mpd_dict = {}
    unmpd_dict = {}
    mpd_b_dict = {}
    unmpd_b_dict = {}
    umi_dict = {}
    mpd_res = {}
    unmpd_res = {}
    try:
        for read in bam_lst:
            umi = ext_fun(read)
            splt_umi = spl_fun(umi, spl_cha)
            umi_res = pos_function.umi_maker(read, splt_umi)
            qr_nm = umi_res[0]
            strand = umi_res[1]
            umi_id = umi_res[2]

            try:
                umi_dict[umi_id][strand][qr_nm].append(read)
            except KeyError:
                if umi_id not in umi_dict:
                    umi_dict[umi_id] = {"String_1": dict(), "String_2": dict()}
                umi_dict[umi_id][strand][qr_nm] = [read]
        umi_dict = {k: v for k, v in umi_dict.items() if v}

        str1_ms_hits = {}
        str2_ms_hits = {}
        str1_us_hits = {}
        str2_us_hits = {}
        # Iterates through every UMI-key in the dict
        for umi_key in umi_dict.keys():
            # Retrieves the forward and reverse molecule hits from said UMI-key
            str1_lst = umi_dict[umi_key]["String_1"]
            str2_lst = umi_dict[umi_key]["String_2"]
            # If any of the keys have an empty list entry, separately calculates the pos_hits and stores it
            if str1_lst:
                if str2_lst:
                    mpd_dict["String_1_Hits"] = pos_function.pos_hits(str1_lst, rec_pos)[0]
                    mpd_dict["String_2_Hits"] = pos_function.pos_hits(str2_lst, rec_pos)[0]
                    unmpd_dict["String_1_Hits"] = pos_function.pos_hits(str1_lst, rec_pos)[1]
                    unmpd_dict["String_2_Hits"] = pos_function.pos_hits(str2_lst, rec_pos)[1]
                    mpd_b_dict = base_function.ffpe_finder(mpd_dict, var_nuc, ref_nuc, ffpe_b)
                    unmpd_b_dict = base_function.ffpe_finder(unmpd_dict, var_nuc, ref_nuc, ffpe_b)
                else:
                    str1_ms_hits = pos_function.pos_hits(str1_lst, rec_pos)[0]
                    str1_us_hits = pos_function.pos_hits(str1_lst, rec_pos)[1]
            elif str2_lst:
                if not str1_lst:
                    str2_ms_hits = pos_function.pos_hits(str2_lst, rec_pos)[0]
                    str2_us_hits = pos_function.pos_hits(str2_lst, rec_pos)[1]
            sing_m_dict = {"String_1_Single": str1_ms_hits, "String_2_Single": str2_ms_hits}
            sing_u_dict = {"String_1_Single": str1_us_hits, "String_2_Single": str2_us_hits}
            mpd_res[umi_key] = {"Single_Hits": sing_m_dict, "Mate_Hits": mpd_b_dict}
            unmpd_res[umi_key] = {"Single_Hits": sing_u_dict, "Mate_Hits": unmpd_b_dict}
    except KeyError as e:
        print("ERROR: The requested key " + str(e) + " does not exist")
    return [mpd_res, unmpd_res]


def inf_builder(inp_dict, n_ref, n_alt):
    """ Uses the output from pos_checker to generate a list containing strings representing the data found for each
    record, more specifically support for each variant-type, as well as the support for the reference and variant
    call for str1 and str2.
    The inp_dict is meant to be output from the var_extract function,and is divided into two dicts named
    Single Hits and Mate Hits. The Single-dict contains the molecular support for the reference genome nucleotide
    and the variant nucleotide based on all reads without a mate.
    The Mate-Hits dicts instead contaisn data regarding the variant-classification, the support for each variant type,
    and the molecular support for the reference genome nucleotide and the variant nucleotide based on reads with a mate.

    Args:
        :param inp_dict: Input dict dict for mapped and unmapped reads. Each of these dicts containing a single-hits
        and a mate-hits dict. The mate-hits dict in turn contains data regarding if the variant is a mutation,
        no mutation, FFPE-artfefact deletion or N-call. Whereas the single-hits dict contains positional data for reads
        with no mate.
        :param n_ref: Nucletoide in the reference genome for the variant-record variant position
        :param n_alt: Variant for the variant-record variant position

    Returns:
        :return: Returns a list containing the support for each variant-type, as well as the support for the reference
        and variant call for str1 and str2
    """
    # Calls the sup_count function to obtain the support for each nucleotide called
    ref_sup = count_function.nuc_count(inp_dict, n_ref)
    alt_sup = count_function.nuc_count(inp_dict, n_alt)
    type_sup = count_function.mol_count(inp_dict)
    ref_p = None
    ref_s = None
    var_p = None
    var_s = None

    for var in n_alt:
        var_p = str(alt_sup["Paired"]["String_1"][var]) + ";" + str(alt_sup["Paired"]["String_2"][var])
        var_s = str(alt_sup["String_1_Single"][var]) + ";" + str(alt_sup["String_2_Single"][var])

    for ref in n_ref:
        ref_p = str(ref_sup["Paired"]["String_1"][ref]) + ";" + str(ref_sup["Paired"]["String_2"][ref])
        ref_s = str(ref_sup["String_1_Single"][ref]) + ";" + str(ref_sup["String_2_Single"][ref])

    return [type_sup, ref_p, var_p, ref_s, var_s]