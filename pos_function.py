import base_function
from collections import Counter


def pos_checker(bam_lst, record_pos, ref_var, ref_base, b_trans):
    mpd_dict = {}
    unmpd_dict = {}
    mpd_b_dict = {}
    unmpd_b_dict = {}
    umi_dict = {}
    mpd_res = {}
    unmpd_res = {}
    try:
        for read in bam_lst:
            # Gets the name of the sequence
            qn = read.query_name
            # Splits the name based on "_" character
            qn_spl = qn.split("_")
            # Selects the last element of the vector
            umi = qn_spl[-1]
            # Obtains the barcode's based on splitting the vector  using the "+" symbol
            brc = umi.split("+")
            umi_l = brc[0]
            umi_r = brc[1]

            # Checks read polarity to see whether or not it is reverse paired
            strand = "String_1"
            qr_nm = read.query_name
            if read.is_read1:
                if read.is_reverse:
                    # Forward direction on the reverse molecule
                    strand = "String_2"
                    umi_id = umi_r + "_" + umi_l
                else:
                    # Forward direction on the forward molecule
                    umi_id = umi_l + "_" + umi_r
            else:
                if read.is_reverse:
                    # Reverse direction on the reverse molecule
                    umi_id = umi_l + "_" + umi_r
                else:
                    # Forward direction on the reverse molecule
                    strand = "String_2"
                    umi_id = umi_r + "_" + umi_l

            # Adds every read to a list corresponding to its polarity, the lists in turn are part of a dict represented
            # By its unique id-tag created from its leftmost position and UMIs
            try:
                umi_dict[umi_id][strand][qr_nm].append(read)
            except KeyError:
                if umi_id not in umi_dict:
                    umi_dict[umi_id] = {"String_1": dict(), "String_2": dict()}
                umi_dict[umi_id][strand][qr_nm] = [read]

        str1_ms_hits = {}
        str2_ms_hits = {}
        str1_us_hits = {}
        str2_us_hits = {}
        umi_dict = {k: v for k, v in umi_dict.items() if v}
        # Iterates through every UMI-key in the dict
        for umi_key in umi_dict.keys():
            # Retrieves the forward and reverse molecule hits from said UMI-key
            str1_lst = umi_dict[umi_key]["String_1"]
            str2_lst = umi_dict[umi_key]["String_2"]
            # If any of the keys have an empty list entry, separately calculates the pos_hits and stores it
            if str1_lst:
                if str2_lst:
                    mpd_dict["String_1_Hits"] = pos_hits(str1_lst, record_pos)[0]
                    mpd_dict["String_2_Hits"] = pos_hits(str2_lst, record_pos)[0]
                    unmpd_dict["String_1_Hits"] = pos_hits(str1_lst, record_pos)[1]
                    unmpd_dict["String_2_Hits"] = pos_hits(str2_lst, record_pos)[1]
                    mpd_b_dict = base_function.ffpe_finder(mpd_dict, ref_var, ref_base, b_trans)
                    unmpd_b_dict = base_function.ffpe_finder(unmpd_dict, ref_var, ref_base, b_trans)
                else:
                    str1_ms_hits = pos_hits(str1_lst, record_pos)[0]
                    str1_us_hits = pos_hits(str1_lst, record_pos)[1]
            elif str2_lst:
                if not str1_lst:
                    str2_ms_hits = pos_hits(str2_lst, record_pos)[0]
                    str2_us_hits = pos_hits(str2_lst, record_pos)[1]
            sing_m_dict = {"String_1_Single": str1_ms_hits, "String_2_Single": str2_ms_hits}
            sing_u_dict = {"String_1_Single": str1_us_hits, "String_2_Single": str2_us_hits}
            mpd_res[umi_key] = {"Single_Hits": sing_m_dict, "Mate_Hits": mpd_b_dict}
            unmpd_res[umi_key] = {"Single_Hits": sing_u_dict, "Mate_Hits": unmpd_b_dict}
    except KeyError as e:
        print("ERROR: The requested key " + str(e) + " does not exist")
    return [mpd_res, unmpd_res]


def pos_hits(inp_lst, record_pos):
    # Loops through each key and its respective reads to extract their variant position data and then counts
    # The no. hits for each respective letter for this position
    uread_base = None
    read_base = None
    mpd_dict = Counter({"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "-": 0})
    unmpd_dict = Counter({"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "-": 0})

    # Iterates through every query_name entry within the given UMI-key for the direction
    for query_name, read in inp_lst.items():
        # If two entries exist for the same query name (ie: forward and reverse strand), calls the base_check function
        # For both reads, then checks if these are identical or if they are different. If they are different, skips the
        # Query-name pair as this  indicates some form of error, as the same molecule should have a identical base
        if len(read) == 2:
            read_base = base_function.base_check(read[0], record_pos)
            r2b = base_function.base_check(read[1], record_pos)
            # If the read mates do not agree on which base the position represents, the data cannot be used
            if read_base != r2b:
                continue
        elif len(read) == 1:
            # Checks if the read simply is alone as its mate did not align,
            # If it lacks a mate it is categorized as unmapped.
            if read[0].mate_is_unmapped:
                uread_base = base_function.base_check(read[0], record_pos)
            else:
                read_base = base_function.base_check(read[0], record_pos)
        else:
            raise Exception("No. entries exceeds 2, cannot handle the no. reads: {}".format(len(read)))
        if read_base:
            tmp_m_dict = Counter(base_function.base_count(read_base))
            mpd_dict = mpd_dict + tmp_m_dict
            mpd_dict.update(base_function.base_count(read_base))
        if uread_base:
            tmp_u_dict = Counter(base_function.base_count(uread_base))
            unmpd_dict = unmpd_dict + tmp_u_dict
    # Selects the most prominent base in the unmapped/mapped dict
    mpd_bas = max(mpd_dict, key=mpd_dict.get)
    unmpd_bas = max(unmpd_dict, key=unmpd_dict.get)
    # Returns a list of the mapped and unmapped most prominent base
    bas_lst = [mpd_bas, unmpd_bas]
    return bas_lst
