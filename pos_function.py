import nuc_function
import warnings
from collections import Counter


def umi_maker(read, splt_umi):
    """ The umi_maker function rearranges the UMI-tag belonging to a read, based on if the read is read 1 or read 2
    in combination with its directionality. To extract the UMI from the read the ext_fun function is used call either
    qrn_ext or rx_ext based on user input. The UMI is then transformed into a string and used as input for the spl_fun
    function. Returns the query-name of the read, the strand it belongs to, and the adjusted UMI-sequence

    Args:
        :param read: Read of interest
        :param splt_umi: UMI-tag for the read split into two strings

    Returns:
        :return: Returns a dict with the query-name, the correct string, and the corrected umi-sequence belonging
        to the read
    """
    umi_l = splt_umi[0]
    umi_r = splt_umi[1]

    # The umi-maker will assign the read as from Strand 1 or Strand 2 of the original molecule dependent on its
    # directionality and if its read 1 or read 2.

    # Checks read polarity to see whether or not it is reverse paired
    strand = "Pos_Str"
    qr_nm = read.query_name
    if read.is_read1:
        if read.is_reverse:
            # Forward direction on the reverse molecule
            strand = "Neg_Str"
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
            strand = "Neg_Str"
            umi_id = umi_r + "_" + umi_l

    # Adds every read to a list corresponding to its polarity, the lists in turn are part of a dict represented
    # By its unique id-tag created from its leftmost position and UMIs

    return [qr_nm, strand, umi_id]


def qrn_ext(read, q_spl_cha):
    """ The qrn_ext function extracts a umi-tag from a read, based on the key being present as the last item in the
    query-name. Returns the umi-tag.

    Args:
        :param read: Read from which the umi-tag is to be extracted from
        :param q_spl_cha: Character from which to split the UMI-tag from the query-name

    Returns:
        :return: A dict containing the umi-tag
    """
    return read.query_name.split(q_spl_cha)[-1]


def rx_ext(read, q_spl_cha):
    """ The rx_ext function extracts a umi-key from a read, based on the key being present in the RX-tag.
    Returns the umi-tag.

    Args:
        :param read: Read from which the umi-tag is to be extracted from
        :param q_spl_cha: Character from which to split the UMI-tag from the query-name (required by the function call)

    Returns:
        :return: A dict containing the umi-tag
    """
    return str(read.get_tag("RX"))


def cha_splt(umi_str, u_spl_cha):
    """ The cha\_splt function splits the umi\_string based on the split-character argument. Returns a list containing
     the umi-tag split into two components.

    >>> cha_splt("ACTACTA+GCTGCTG", "+")
    '[ACTACTA, GCTGCTG]'
    >>> cha_splt("ACTACTA_GCTGCTG", "_")
    '[ACTACTA, GCTGCTG]'

        Args:
            :param umi_str: A string representing the umi-tag to be split
            :param u_spl_cha: Character to split the umi-string by

        Returns:
            :return: A dict containing the umi-tag split in two
        """
    return umi_str.split(u_spl_cha)


def hlf_splt(umi_str, u_spl_cha):
    """ The hlf_splt function splits the umi_string in half based on its length. Returns a list containing the umi-tag
    split into two components.

     >>> hlf_splt("ACTACTAGCTGCTG", "")
    '[ACTACTA, GCTGCTG]'

    Args:
        :param umi_str: A string representing the umi-tag to be split
        :param u_spl_cha: Character to split the umi-string by (required by the function call)

    Returns:
        :return: A dict containing the umi-tag split in two
    """
    tgg = umi_str[:len(umi_str) // 2], umi_str[len(umi_str) // 2:]
    return list(tgg)


def pos_hits(inp_dict, rec_pos):
    """ The pos\_hits function selects the most prominent nuc for a UMI of interest. The function works through
    iterating through all query-names in the input list and determines if the query-name has a mate or not.
    The function then calls the nuc\_check function to retrieve the nuc matching the variant-record position for each
     read. In the next step the retrieved nuc is matched against a dict, and depending on the outcome adds to a
    counter representative of the nuc. This process is repeated for each query name and its subsequent read, and the
    resulting dict is then used to determine the most prominent nucleotide for the UMI, effectively collapsing all
    reads belonging to a UMI. Returns the consensus nucleotide

    Args:
        :param inp_dict: Input dict of reads categorized by their query-name
        Example dict:
        input_dict = {example_name_UMI_ACTGCA+ACTGCA: {read1, read2}, example_2_name_UMI_TGACGT+TGACGT: {read2}}
        :param rec_pos: The position of the called variant in the reference genome

    Returns:
        :return: Returns a dict with mapped and unmapped reads with the most prominent nuc for the
        reads belonging to a UMI
        Example dict:
        cons_dict = {UMI_tag_1: C, UMI_tag_2: T}

    Raises:
        :raises Warning: Raises a warning if a query-name has more then 2 reads belonging to it. If this happens it
        is assumed to be a software error and these reads are ignored.
    """
    # Loops through each key and its respective reads to extract their variant position data and then counts
    # The no. hits for each respective letter for this position
    read_nuc = None
    cons_nuc = None
    singleton_nuc = None
    mate_dict = Counter({"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "-": 0})
    singleton_dict = Counter({"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "-": 0})

    # Iterates through every query_name entry within the given UMI-key for the direction
    for query_name, read in inp_dict.items():
        # If two entries exist for the same query name (ie: forward and reverse strand), calls the nuc_check function
        # For both reads, then checks if these are identical or if they are different. If they are different, skips the
        # Query-name pair as this  indicates some form of error, as the same molecule should have a identical nuc
        if len(read) == 2:
            read_nuc = nuc_function.nuc_check(read[0], rec_pos)
            r2b = nuc_function.nuc_check(read[1], rec_pos)
            # If the read mates do not agree on which nuc the position represents, the data cannot be used
            if read_nuc != r2b:
                continue
        elif len(read) == 1:
            # Checks if the read simply is alone as its mate did not align,
            # If it lacks a mate it is categorized as unmapped.
            if read[0].mate_is_unmapped:
                singleton_nuc = nuc_function.nuc_check(read[0], rec_pos)
            else:
                read_nuc = nuc_function.nuc_check(read[0], rec_pos)
        else:
            warnings.warn("Warning! No. reads belonging to: " + str(query_name) + " exceeds 2, skipping these")

        # Adds to the count of nucleotides belonging to the query-name to the mpd/unmpd dict for the
        # inp_lst belonging to the umi_key
        if read_nuc:
            for nuc in mate_dict.keys():
                if nuc == read_nuc:
                    mate_dict[nuc] += 1
        elif singleton_nuc:
            for nuc in singleton_dict.keys():
                if nuc == singleton_nuc:
                    singleton_dict[nuc] += 1
    # Selects the most prominent nuc (the consensus nucleotide) in the unmapped/mapped dict if the dict have any values
    if max(mate_dict.values()) > 0:
        cons_nuc = max(mate_dict, key=mate_dict.get)
    if max(singleton_dict.values()) > 0:
        singleton_nuc = max(singleton_dict, key=singleton_dict.get)
    # Returns a list of the mapped and unmapped most prominent nuc
    cons_lst = [cons_nuc, singleton_nuc]
    return cons_lst
