import base_function
from collections import Counter


def pos_checker(bam_lst, rec_pos, ref_var, ref_base, ffpe_b, ext_fun, splt_fun, spl_cha):
    """ Function with the purpose of creating a dict based on the directionality and umi-tags of the supplemented
    reads in the bam_lst. Then using said dict to call the pos_hits and ffpe_finder functions to return a dict with
    data regarding positional data and variant types for the variant-record position and the reads aligning to it.

    Args:
        :param bam_lst: Input list of BAM-reads aligning to the variant call
        :param rec_pos: The position of the variant in the reference genome
        :param ref_var: The variant (nucleotide) in the variant call
        :param ref_base: The reference nucleotide in the reference genome for the variant call
        :param ffpe_b: Optional input argument for the FFPE base-types to consider, standard = C>T:G>A
        :param ext_fun: Extraction-function based on data type, used to extract the UMI from a read of interest
        :param splt_fun: Function used to split the UMI-tag for read of interest. Default: Query-name based,
        Optional: RX-tag based
        :param spl_cha: Character used to split the query-name tag or RX-tag to obtain the UMI-sequence halfs. Default
        for query name based :"+", for RX-tag based: ""

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
            splt_umi = splt_fun(umi, spl_cha)
            umi_res = umi_maker(read, splt_umi)
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
                    mpd_dict["String_1_Hits"] = pos_hits(str1_lst, rec_pos)[0]
                    mpd_dict["String_2_Hits"] = pos_hits(str2_lst, rec_pos)[0]
                    unmpd_dict["String_1_Hits"] = pos_hits(str1_lst, rec_pos)[1]
                    unmpd_dict["String_2_Hits"] = pos_hits(str2_lst, rec_pos)[1]
                    mpd_b_dict = base_function.ffpe_finder(mpd_dict, ref_var, ref_base, ffpe_b)
                    unmpd_b_dict = base_function.ffpe_finder(unmpd_dict, ref_var, ref_base, ffpe_b)
                else:
                    str1_ms_hits = pos_hits(str1_lst, rec_pos)[0]
                    str1_us_hits = pos_hits(str1_lst, rec_pos)[1]
            elif str2_lst:
                if not str1_lst:
                    str2_ms_hits = pos_hits(str2_lst, rec_pos)[0]
                    str2_us_hits = pos_hits(str2_lst, rec_pos)[1]
            sing_m_dict = {"String_1_Single": str1_ms_hits, "String_2_Single": str2_ms_hits}
            sing_u_dict = {"String_1_Single": str1_us_hits, "String_2_Single": str2_us_hits}
            mpd_res[umi_key] = {"Single_Hits": sing_m_dict, "Mate_Hits": mpd_b_dict}
            unmpd_res[umi_key] = {"Single_Hits": sing_u_dict, "Mate_Hits": unmpd_b_dict}
    except KeyError as e:
        print("ERROR: The requested key " + str(e) + " does not exist")
    return [mpd_res, unmpd_res]


def qrn_ext(read):
    """ Function for extracting the umi-key, based on the key being present as the last item in the query-name.

    Args:
        :param read: Input bam-read from the umi_dict

    Returns:
        :return: A dict containing the umi-tag
    """
    return read.query_name.split("_")[-1]


def rx_ext(read):
    """ Function for extracting the umi-key, based on the key being present in the RX-tag

    Args:
        :param read: Input bam-read from the umi_dict

    Returns:
        :return: A dict containing the umi-tag
    """
    return str(read.get_tag("RX"))


def cha_splt(umi_str, char):
    """ Split function for the query-name, splits the umi_string based on the split-character argument

        Args:
            :param umi_str: A string representing the umi-tag to be split
            :param char: Character to split the umi-string by

        Returns:
            :return: A dict containing the umi-tag split in two
        """
    return umi_str.split(char)


def hlf_splt(umi_str, char):
    """ Split function for splitting the umi-tag in half based on half its length

    Args:
        :param umi_str: A string representing the umi-tag to be split
        :param char: Not used in this function, however needs to be supplemented as the function call for the query-name
        based split function requires a split character

    Returns:
        :return: A dict containing the umi-tag split in two
    """
    return umi_str[:len(umi_str) // 2], umi_str[len(umi_str) // 2:]


def umi_maker(read, splt_umi):
    """ Changes the directionality of the UMI-tag based on if the read is read_1 or read_2 in combination with its
    directionality. Returns the query-name of the read, the strand it belongs to, and the corrected UMI-sequence

    Args:
        :param read: Input bam-read
        :param splt_umi: UMI-tag for the read, split into a list in half
    Returns:

        :return: Returns a dict with the query-name, the correct string, and the corrected umi-sequence belonging to the
        read
    """
    umi_l = splt_umi[0]
    umi_r = splt_umi[1]

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

    return [qr_nm, strand, umi_id]


def pos_hits(inp_lst, record_pos):
    """ Function for determining the most prominent base for a UMI of interest. Determines if the query-name has a mate
    or not, then uses the base matching the record-pos and counts this for each read, adding to a dict which is then
    used to determine the most prominent base for the UMI belonging to the String in the supplemented input list .

    Args:
        :param inp_lst: Input list of reads aligning to one directionality of the string, categorized by their read-name
        :param record_pos: Position in the reference  genome of the variant

    Returns:
        :return: Returns a dict with mapped and unmapped reads with the most prominent base for the
        reads belonging to a UMI

    Raises:
        :raises Warning: Raises a warning if a query-name has more then 2 reads belonging to it. If this happens it
        is assumed to be a software error and these reads are ignored.
    """
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
            raise Warning("Warning! No. reads belonging to: " + str(query_name) + " exceeds 2, skipping these")

        # Adds to the count of bases belonging to the query-name to the mpd/unmpd dict for the inp_lst belonging to the
        # umi_key
        if read_base:
            for base in mpd_dict.keys():
                if base == read_base:
                    mpd_dict[base] += 1
        elif uread_base:
            for base in unmpd_dict.keys():
                if base == uread_base:
                    unmpd_dict[base] += 1
    # Selects the most prominent base in the unmapped/mapped dict
    mpd_bas = max(mpd_dict, key=mpd_dict.get)
    unmpd_bas = max(unmpd_dict, key=unmpd_dict.get)
    # Returns a list of the mapped and unmapped most prominent base
    bas_lst = [mpd_bas, unmpd_bas]
    return bas_lst
