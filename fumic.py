# First attempt to read in BAM as SAM
# !/usr/bin/env python3

# FUMIC - FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data
# By Hugo Swenson, with assistance from Patrik Smeds and Claes Edenvall
# Made for Klinisk Genetik, Uppsala Akademiska Sjukhus 2019

# Imports modules
import pysam


def vcf_extract(vcf_file, bam_file):
    # parser = argparse.ArgumentParser(description='FUMIC - FFPE-artefact UMI-based Mapper for Imputation in '
    #                                              'Cancer-sample tissue data')
    #
    # parser.add_argument('-i', '--inputfile', help='Input txt file (Required)', required=True)
    # parser.add_argument('-o', '--outputfile', help='Output txt file (Required)', required=True)
    #
    # args = vars(parser.parse_args())
    # infile = open(args['inputfile'], "rb")
    # outfile = open(args['outputfile'], "rb")

    vcf_head = vcf_file.header
    vcf_head.filters.add('FFPE', None, None, 'FFPE Artefact')
    vcf_head.formats.add("UMI", ".", "String", "UMI information for reference,variant Paired:SForward:SReverse")
    n_vcf = pysam.VariantFile("output_vcf.vcf", mode='w', header=vcf_head)

    # Initialize null lists
    # Retrieve the position from every record in the VCF file
    for record in vcf_file.fetch():
        bam_lst = []
        var_sd = []
        ref_sd = []
        # Copies the record information
        n_cop = record.copy()
        n_pos = record.pos
        # The position that is returned to Python is 0 - based, NOT 1 - based as in the VCF file.
        record_pos = (n_pos - 1)
        # Converts the tuple to a string, splits it into individual characters in a list, then removes duplicates
        n_ref = record.ref
        n_ref = ''.join(n_ref)
        list(n_ref)
        n_ref = set(n_ref)
        n_alt = record.alts
        n_alt = ''.join(n_alt)
        list(n_alt)
        n_alt = set(n_alt)

        # Checks so that the length of the list is not greater then 1 (temporary solution for handling SNVs only)
        if len(n_ref) > 1 or len(n_alt) > 1:
            continue

        # # Removes reference bases also present in the variant call
        # for alt in n_alt:
        #     n_ref = [r_base for r_base in n_ref if r_base != alt]

        # # Use the record position to fetch all reads matching it, then append to the list
        for read in bam_file.fetch('chr8', record_pos, record_pos+1):
            bam_lst.append(read)

        # Calls the pos_checker function to obtain ffpe_data
        ffpe_data = (pos_checker(bam_lst, record_pos, n_alt, n_ref))

        # Calls the sup_count function to obtain the support for each nucleotide called
        ref_sup = sup_count(ffpe_data, n_ref)
        alt_sup = sup_count(ffpe_data, n_alt)
        pos_sup = {"Variant": alt_sup, "Reference": ref_sup}

        for var in n_alt:
            var_sd = str(pos_sup["Paired"]["String_1"][var]) + ":" + str(pos_sup["Paired"]["String_2"][var]) + \
                     ":" + str(pos_sup["String_1_Single"][var]) + ":" + str(pos_sup["String_2_Single"][var])
        for ref in n_ref:
            ref_sd = str(pos_sup["Paired"]["String_1"][ref]) + ":" + str(pos_sup["Paired"]["String_2"][ref]) + \
                     ":" + str(pos_sup["String_1_Single"][ref]) + ":" + str(pos_sup["String_2_Single"][ref])
        for sample in n_cop.samples:
            n_cop.samples[sample]['UMI'] = var_sd + ref_sd
        # umi_lst = []
        # for ref in pos_sup["Reference"]:
        #     ref_c = pos_sup["Variant"][ref]
        #     umi_lst.append(ref_c + ":")
        # for var in pos_sup["Variant"]:
        #     umi_lst.append(pos_sup["Variant"][var])

        # Checks if any record in the returned dict indicates an FFPE, if so updates the n_fil parameter
        for umi_key in ffpe_data:
            if ffpe_data[umi_key]["Mate_Hits"]:
                if ffpe_data[umi_key]["Mate_Hits"]["FFPE_Hits"]:
                    n_cop.filter.add("FFPE")
                    break
        # Adds a record to the new VCF-file
        n_vcf.write(n_cop)


def pos_checker(bam_lst, record_pos, ref_var, ref_base):
    umi_dict = {}
    base_res = {}
    var_dict = {}
    pos_res = {}
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

        str1_s_hits = {}
        str2_s_hits = {}
        umi_dict = {k: v for k, v in umi_dict.items() if v}
        counter = 0
        # Iterates through every UMI-key in the dict
        for umi_key in umi_dict.keys():
            # Retrieves the forward and reverse molecule hits from said UMI-key
            str1_lst = umi_dict[umi_key]["String_1"]
            str2_lst = umi_dict[umi_key]["String_2"]
            # If any of the keys have an empty list entry, separately calculates the pos_hits and stores it
            if str1_lst:
                if str2_lst:
                    counter += 1
                    base_res["String_1_Hits"] = pos_hits(str1_lst, record_pos)
                    base_res["String_2_Hits"] = pos_hits(str2_lst, record_pos)
                    var_dict = ffpe_finder(base_res, ref_var, ref_base)
                else:
                    str1_s_hits = pos_hits(str1_lst, record_pos)
            elif str2_lst:
                if not str1_lst:
                    str2_s_hits = pos_hits(str2_lst, record_pos)
            sing_dict = {"String_1_Single": str1_s_hits, "String_2_Single": str2_s_hits}
            pos_res[umi_key] = {"Single_Hits": sing_dict, "Mate_Hits": var_dict}
    except KeyError as e:
        print("ERROR: The requested key " + str(e) + " does not exist")
    return pos_res


def pos_hits(inp_lst, record_pos):
    # Loops through each key and its respective reads to extract their variant position data and then counts
    # The no. hits for each respective letter for this position
    n_a = 0
    n_t = 0
    n_g = 0
    n_c = 0
    n_n = 0
    n_d = 0

    base_lst = ["A", "T", "G", "C", "N"]
    # Iterates through every query_name entry within the given UMI-key for the direction
    for query_name, read in inp_lst.items():
        read_base = None
        # If two entries exist for the same query name (ie: forward and reverse strand), calls the base_check function
        # For both reads, then checks if these are identical or if they are different. If they are different, skips the
        # Query-name pair as this  indicates some form of error, as the same molecule should have a identical base
        if len(read) == 2:
            read_base = base_check(read[0], record_pos)
            r2b = base_check(read[1], record_pos)
            if read_base != r2b:
                continue
        elif len(read) == 1:
            read_base = base_check(read[0], record_pos)
            #todo handle umapped mate
        else:
            raise Exception("No entries exceeds 2, cannot handle the no. reads: {}".format(len(read)))

        # Checks the read_base and adds a count if it matches against any of the entries in the base_lst, if so adds
        # To a counter
        if read_base == base_lst[0]:
            n_a += 1
        elif read_base == base_lst[1]:
            n_t += 1
        elif read_base == base_lst[2]:
            n_g += 1
        elif read_base == base_lst[3]:
            n_c += 1
        elif read_base == base_lst[4]:
            n_n += 1
        else:
            n_d += 1

    # Returns a dict with the counts for each base respectively
    pos_dict = {"A": n_a, "T": n_t, "G": n_g, "C": n_c, "N": n_n, "-": n_d}
    return pos_dict


def base_check(read, record_pos):
    try:
        # Gets the positions the sequence maps to in the reference
        # Full length with soft clips is required for the index selection to be correct
        read_pos = read.get_reference_positions(full_length=True)
        ind_pos = read_pos.index(record_pos)
        # Gets the index of the position the sequence maps to
        # Obtains the query sequence (the sequence as it were read)
        read_seq = read.query_sequence
        read_seq = list(read_seq)
        # Gets the base present at the index position of the variant
        read_base = read_seq[ind_pos]
        # Checks the value of the base, adds to a counter based on its value
        return read_base
    except ValueError:
        pass


def ffpe_finder(base_res, ref_var, ref_base):
    ffpe_dict = {}
    mut_dict = {}
    ref_dict = {}
    out_dict = {}
    try:
        # Retrieves the forward and reverse base hits for the key value
        f_hits = base_res["String_1_Hits"]
        r_hits = base_res["String_2_Hits"]
        # Iterates through each variant called by the program
        for var_call in ref_var:
            for r_base in ref_base:
                # First checks if the variant is present in both forward and reverse molecule, if so it is deemed a
                # Mutation
                if f_hits[var_call] > 0 and r_hits[var_call] > 0:
                    mut_dict = {"String_1": f_hits, "String_2": r_hits}
                # If the variant is present in only on of the strings, it is deemed a FFPE
                elif f_hits[var_call] > 0 and r_hits[var_call] == 0:
                    ffpe_dict = {"String_1": f_hits, "String_2": r_hits}
                elif f_hits[var_call] == 0 and r_hits[var_call] > 0:
                    ffpe_dict = {"String_1": f_hits, "String_2": r_hits}
                # If both positions are deemed to be the reference, it is added as a reference to the dict
                elif f_hits[r_base] > 0 and r_hits[r_base] > 0:
                    ref_dict = {"String_1": f_hits, "String_2": r_hits}
                # If none of the above statements are fulfilled, then it is deemed as another form of mutation
                else:
                    out_dict = {"String_1": f_hits, "String_2": r_hits}
    except KeyError as e:
        print("No match for: " + str(e) + " found, comparison not possible")
        # Returns both dictionaries in a list
    var_dict = {"FFPE_Hits": ffpe_dict, "Mutation_Hits": mut_dict, "Reference_Hits": ref_dict,
                "Outlier_Hits": out_dict}
    # print(var_lst)
    return var_dict


def sup_count(input_dict, nuc):
    n_sup = {"String_1_Single": {}, "String_2_Single": {}, "Paired": {}}
    for n in nuc:
        n_sup["String_1_Single"][n] = 0
        n_sup["String_2_Single"][n] = 0
        n_sup["Paired"] = {"String_1": {}, "String_2": {}}
        n_sup["Paired"]["String_1"][n] = 0
        n_sup["Paired"]["String_2"][n] = 0
        for umi_key in input_dict:
            # Iterates through each alternative allele called by the variant caller
            # Counts the no. molecules supporting the variant for each direction in the single hits
            if input_dict[umi_key]["Single_Hits"]:
                for sh_dir in input_dict[umi_key]["Single_Hits"]:
                    if input_dict[umi_key]["Single_Hits"][sh_dir]:
                        if n in input_dict[umi_key]["Single_Hits"][sh_dir]:
                            n_sup[sh_dir][n] += input_dict[umi_key]["Single_Hits"][sh_dir][n]
            # Sees if the dict "Mate Hits" is empty, if not, counts the no molecules supporting the variant for each
            # dict within the "Mate Hits" dictionary
            if input_dict[umi_key]["Mate_Hits"]:
                for v_h in input_dict[umi_key]["Mate_Hits"]:
                    if input_dict[umi_key]["Mate_Hits"][v_h]:
                        if input_dict[umi_key]["Mate_Hits"][v_h]["String_1"]:
                            if n in input_dict[umi_key]["Mate_Hits"][v_h]["String_1"]:
                                n_sup["Paired"]["String_1"][n] += input_dict[umi_key]["Mate_Hits"][v_h][
                                    "String_1"][n]
                        if input_dict[umi_key]["Mate_Hits"][v_h]["String_2"]:
                            if n in input_dict[umi_key]["Mate_Hits"][v_h]["String_2"]:
                                n_sup["Paired"]["String_2"][n] += input_dict[umi_key]["Mate_Hits"][v_h][
                                        "String_2"][n]
    return n_sup


def main():
    # Reads in the bam file
    bam_file = pysam.AlignmentFile("26-sort.REF_chr8.bam", "r")
    vcf_file = pysam.VariantFile("26-ensembl.chr8.vcf", "r")
    # Extracts all reads present in the bam file which have their regions matching to that of a variant call in the VCF
    vcf_extract(vcf_file, bam_file)
    bam_file.close()


if __name__ == "__main__":
    main()
