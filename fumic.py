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
            var_sd = str(pos_sup["Variant"]["Paired"][var]) + ":" + str(pos_sup["Variant"][
                        "Forward Single"][var]) + ":" + str(pos_sup["Variant"]["Reverse Single"][var])
        for ref in n_ref:
            ref_sd = str(pos_sup["Reference"]["Paired"][ref]) + ":" + str(pos_sup["Reference"]["Forward Single"][
                                                ref]) + ":" + str(pos_sup["Reference"]["Reverse Single"][ref])
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
            if ffpe_data[umi_key]["Variant Hits"]:
                if ffpe_data[umi_key]["Variant Hits"]["FFPE Hits"]:
                    n_cop.filter.add("FFPE")
                    break
        # Adds a record to the new VCF-file
        n_vcf.write(n_cop)


def pos_checker(bam_lst, record_pos, ref_var, ref_base):
    umi_dict = {}
    base_res = {}
    var_dict = {}
    pos_res = {}
    sing_dict = {}
    try:
        for read in bam_lst:
            # Gets the name of the sequence
            qn = read.query_name
            # Splits the name based on "_" character
            qn_spl = qn.split("_")
            # Extracts the 1- based leftmost mapping POSition
            # l_pos = read.reference_start
            # r_pos = read.reference_end
            # Selects the last element of the vector
            umi = qn_spl[-1]
            # Obtains the barcode's based on splitting the vector  using the "+" symbol
            brc = umi.split("+")
            umi_l = brc[0]
            umi_r = brc[1]

            # Checks read polarity to see whether or not it is reverse paired
            strand = "Forward_Molecule"
            if read.is_read1:
                if read.is_reverse:
                    strand = "Reverse_Molecule"
                    # umi_id = str(l_pos) + "_" + str(r_pos) + "_" + umi_l + "_" + umi_r
                    umi_id = umi_r + "_" + umi_l
                else:
                    # umi_id = str(l_pos) + "_" + str(r_pos) + "_" + umi_l + "_" + umi_r
                    umi_id = umi_l + "_" + umi_r
            else:
                if read.is_reverse:
                    # umi_id = str(r_pos) + "_" + str(l_pos) + "_" + umi_r + "_" + umi_l
                    umi_id = umi_l + "_" + umi_r
                else:
                    strand = "Reverse_Molecule"
                    # umi_id = str(r_pos) + "_" + str(l_pos) + "_" + umi_r + "_" + umi_l
                    umi_id = umi_r + "_" + umi_l

            # Adds every read to a list corresponding to its polarity, the lists in turn are part of a dict represented
            # By its unique id-tag created from its leftmost position and UMIs
            try:
                umi_dict[umi_id][strand].append(read)
            except KeyError:
                umi_dict[umi_id] = {"Forward_Molecule": [], "Reverse_Molecule": []}
                umi_dict[umi_id][strand].append(read)
        f_s_hits = {}
        r_s_hits = {}
        umi_dict = {k: v for k, v in umi_dict.items() if v}
        counter = 0
        for umi_key in umi_dict.keys():
            f_lst = umi_dict[umi_key]["Forward_Molecule"]
            r_lst = umi_dict[umi_key]["Reverse_Molecule"]
            # If any of the keys have an empty list entry, removes this entry and only calculates the value for the
            # Entry containing the correct value
            if f_lst:
                if r_lst:
                    counter += 1
                    base_res["Forward Reads"] = f_lst
                    base_res["Reverse Reads"] = r_lst
                    base_res["Forward Hits"] = pos_hits(f_lst, record_pos)
                    base_res["Reverse Hits"] = pos_hits(r_lst, record_pos)
                    var_dict = ffpe_finder(base_res, ref_var, ref_base)
            elif f_lst:
                if not r_lst:
                    f_lst = umi_dict[umi_key]["Forward_Molecule"]
                    f_s_hits = pos_hits(f_lst, record_pos)
            elif r_lst:
                if not f_lst:
                    r_lst = umi_dict[umi_key]["Reverse_Molecule"]
                    r_s_hits = pos_hits(r_lst, record_pos)
            sing_dict = {"Forward Single": f_s_hits, "Reverse Single": r_s_hits}
            pos_res[umi_key] = {"Single Hits": sing_dict, "Variant Hits": var_dict}
    except KeyError:
        print("ERROR: The requested key does not exist")
    return pos_res


def pos_hits(inp_lst, record_pos):
    # Loops through each key and its respective reads to extract their variant position data and then counts
    # The no. hits for each respective letter for this position
    n_a_f = 0
    n_t_f = 0
    n_g_f = 0
    n_c_f = 0
    n_n_f = 0
    n_d_f = 0
    n_a_r = 0
    n_t_r = 0
    n_g_r = 0
    n_c_r = 0
    n_n_r = 0
    n_d_r = 0
    # count = 0
    for read in inp_lst:
        # Gets the positions the sequence maps to in the reference
        # Full length with soft clips is required for the index selection to be correct
        read_pos = read.get_reference_positions(full_length=True)
        # read_pos_en = list(enumerate(read_pos, 0))
        # print(read_pos_en)
        try:
            # Gets the index of the position the sequence maps to
            ind_pos = read_pos.index(record_pos)
            # Obtains the query sequence (the sequence as it were read)
            read_seq = read.query_sequence
            read_seq = list(read_seq)
            # Gets the base present at the index position of the variant
            read_base = read_seq[ind_pos]
            # Checks the value of the base, adds to a counter based on its value
            if read.is_reverse:
                if read_base == 'A':
                    n_a_r += 1
                elif read_base == 'T':
                    n_t_r += 1
                elif read_base == 'G':
                    n_g_r += 1
                elif read_base == 'C':
                    n_c_r += 1
                else:
                    n_n_r += 1
            else:
                if read_base == 'A':
                    n_a_f += 1
                elif read_base == 'T':
                    n_t_f += 1
                elif read_base == 'G':
                    n_g_f += 1
                elif read_base == 'C':
                    n_c_f += 1
                else:
                    n_n_f += 1
        # If the reference-base is not found, the mutation is noted as a deletion
        except ValueError:
            if read.is_reverse:
                n_d_r += 1
            else:
                n_d_f += 1
    # Returns a dict with the counts for each base respectively
    pos_dict = {"Forward String": {"A": n_a_f, "T": n_t_f, "G": n_g_f, "C": n_c_f, "N": n_n_f, "-": n_d_f},
                "Reverse String": {"A": n_a_r, "T": n_t_r, "G": n_g_r, "C": n_c_r, "N": n_n_r, "-": n_d_r}}
    return pos_dict


def ffpe_finder(base_res, ref_var, ref_base):
    ffpe_dict = {}
    mut_dict = {}
    ref_dict = {}
    out_dict = {}
    try:
        # Retrieves the forward and reverse base hits for the key value
        f1_hits = base_res["Forward Hits"]["Forward String"]
        r2_hits = base_res["Forward Hits"]["Reverse String"]
        f2_hits = base_res["Reverse Hits"]["Forward String"]
        r1_hits = base_res["Reverse Hits"]["Reverse String"]
        # Iterates through each variant called by the program
        for var_call in ref_var:
            for r_base in ref_base:
                # Sees if the variant exists alongside the reference, if so it is deemed an FFPE
                if f1_hits[var_call] > 0 and r2_hits[var_call] > 0:
                    mut_dict = {"Forward Molecule": base_res["Forward Hits"],
                                "Reverse Molecule": base_res["Reverse Hits"]}
                elif f2_hits[var_call] > 0 and r1_hits[var_call] > 0:
                    mut_dict = {"Forward Molecule": base_res["Forward Hits"],
                                "Reverse Molecule": base_res["Reverse Hits"]}
                elif f1_hits[var_call] > 0 and r2_hits[var_call] == 0:
                    ffpe_dict = {"Forward Molecule": base_res["Forward Hits"],
                                 "Reverse Molecule": base_res["Reverse Hits"]}
                elif f1_hits[var_call] == 0 and r2_hits[var_call] > 0:
                    ffpe_dict = {"Forward Molecule": base_res["Forward Hits"],
                                 "Reverse Molecule": base_res["Reverse Hits"]}
                elif f2_hits[var_call] > 0 and r1_hits[var_call] == 0:
                    ffpe_dict = {"Forward Molecule": base_res["Forward Hits"],
                                 "Reverse Molecule": base_res["Reverse Hits"]}
                elif f2_hits[var_call] == 0 and r2_hits[var_call] > 0:
                    ffpe_dict = {"Forward Molecule": base_res["Forward Hits"],
                                 "Reverse Molecule": base_res["Reverse Hits"]}
                elif f1_hits[r_base] > 0 and r2_hits[r_base] > 0:
                    ref_dict = {"Forward Molecule": base_res["Forward Hits"],
                                "Reverse Molecule": base_res["Reverse Hits"]}
                elif f2_hits[r_base] > 0 and r1_hits[r_base] > 0:
                    ref_dict = {"Forward Molecule": base_res["Forward Hits"],
                                "Reverse Molecule": base_res["Reverse Hits"]}
                # If none of the above statements are fulfilled, then it is deemed as another form of mutation
                else:
                    mut_dict = {"Forward Molecule": base_res["Forward Hits"],
                                "Reverse Molecule": base_res["Reverse Hits"]}
    except KeyError as e:
        print("No match for: " + str(e) + " found, comparison not possible")
        # Returns both dictionaries in a list
    var_dict = {"FFPE Hits": ffpe_dict, "Mutation Hits": mut_dict, "Reference Hits": ref_dict,
                "Other Mutation Hits": out_dict}
    # print(var_lst)
    return var_dict


def sup_count(input_dict, nuc):
    n_sup = {"Forward Single": {}, "Reverse Single": {}, "Paired": {}}

    for n in nuc:
        n_sup["Forward Single"][n] = 0
        n_sup["Reverse Single"][n] = 0
        n_sup["Paired"][n] = 0
        for umi_key in input_dict:
            # Iterates through each alternative allele called by the variant caller
            # Counts the no. molecules supporting the variant in the single forward molecule
            if input_dict[umi_key]["Single Hits"]:
                for dict_entry in input_dict[umi_key]["Single Hits"]:
                    if input_dict[umi_key]["Single Hits"][dict_entry]:
                        if input_dict[umi_key]["Single Hits"][dict_entry]["Forward String"]:
                            if n in input_dict[umi_key]["Single Hits"][dict_entry]["Forward String"]:
                                n_sup[dict_entry][n] += input_dict[umi_key]["Single Hits"][dict_entry][
                                    "Forward String"][n]
                            if n in input_dict[umi_key]["Single Hits"][dict_entry]["Reverse String"]:
                                n_sup[dict_entry][n] += input_dict[umi_key]["Single Hits"]["Forward Single"][
                                    "Reverse String"][n]

            # Sees if the dict variant hits is empty, if not, counts the no molecules supporting the variant as FFPE
            if input_dict[umi_key]["Variant Hits"]:
                for v_h in input_dict[umi_key]["Variant Hits"]:
                    print(v_h)
                    if input_dict[umi_key]["Variant Hits"][v_h]:
                        if input_dict[umi_key]["Variant Hits"][v_h]["Forward Molecule"]:
                            if input_dict[umi_key]["Variant Hits"][v_h]["Forward Molecule"]["Forward String"]:
                                if n in input_dict[umi_key]["Variant Hits"][v_h]["Forward Molecule"]["Forward String"]:
                                    n_sup["Paired"][n] += input_dict[umi_key]["Variant Hits"][v_h][
                                        "Forward Molecule"]["Forward String"][n]
                                if n in input_dict[umi_key]["Variant Hits"][v_h]["Forward Molecule"][
                                        "Reverse String"]:
                                    n_sup["Paired"][n] += input_dict[umi_key]["Variant Hits"][v_h][
                                        "Forward Molecule"]["Reverse String"][n]
                        if input_dict[umi_key]["Variant Hits"][v_h]["Reverse Molecule"]:
                            if input_dict[umi_key]["Variant Hits"][v_h]["Reverse Molecule"]["Forward String"]:
                                if n in input_dict[umi_key]["Variant Hits"][v_h]["Reverse Molecule"]["Forward String"]:
                                    n_sup["Paired"][n] += input_dict[umi_key]["Variant Hits"][v_h][
                                        "Reverse Molecule"]["Forward String"][n]
                                if n in input_dict[umi_key]["Variant Hits"][v_h]["Forward Molecule"]["Reverse String"]:
                                    n_sup["Paired"][n] += input_dict[umi_key]["Variant Hits"][v_h][
                                        "Reverse Molecule"]["Reverse String"][n]
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
