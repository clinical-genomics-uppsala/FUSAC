# First attempt to read in BAM as SAM
# !/usr/bin/python

# FUMIC (FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data)
# By Hugo Swenson, with assistance from Patrik Smeds and Claes Edenvall
# For Klinisk Genomik, Uppsala Akademiska Sjukhus 2018

# Imports modules
import pysam
import argparse
import string


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
    n_vcf = pysam.VariantFile("output_vcf.vcf", mode='w', header=vcf_head)
    # Initialize a null list
    # Retrieve the position from every record in the VCF file
    for record in vcf_file.fetch():
        bam_lst = []
        # Copies the record information
        n_alt = record.alts
        n_pos = record.pos
        n_ref = record.ref
        n_cop = record.copy()
        # The position that is returned to Python is 0 - based, NOT 1 - based as in the VCF file.
        record_pos = (n_pos - 1)

        for alt in n_alt:
            # Removes any entry of the r_base also present in the var_call
            n_ref = list(n_ref)
            n_ref = [r_base for r_base in n_ref if r_base != alt]
        # Use the record position to fetch all reads matching it, then append to the list
        for read in bam_file.fetch('chr8', record_pos, record_pos+1):
            bam_lst.append(read)
        # Calls the pos_checker function to obtain ffpe_data
        ffpe_data = (pos_checker(bam_lst, record_pos, n_alt, n_ref))
        pos_sup = sup_count(ffpe_data, n_ref, n_alt)
        # Checks if any record in the returned dict indicates an FFPE, if so updates the n_fil parameter
        for umi_key in ffpe_data:
            if ffpe_data[umi_key]["Variant Hits"]:
                if ffpe_data[umi_key]["Variant Hits"]["FFPE Hits"]:
                    n_cop.filter.add("FFPE")
                    break
        # Adds a record to the new VCF-file
        n_vcf.write(n_cop)


def sup_count(input_dict, n_ref, n_alt):
    alt_sup = {}
    ref_sup = {}
    pos_sup = {}
    for umi_key in input_dict:
        # Iterates through each alternative allele called by the variant caller
        for alt in n_alt:
            # Counts the no. molecules supporting the variant in the single forward molecule
            if input_dict[umi_key]["Forward Single"]:
                alt_sup[alt] += input_dict[umi_key]["Forward Single"].count(alt)
            # Counts the no. molecules supporting the variant in the single reverse molecule
            if input_dict[umi_key]["Reverse Single"]:
                alt_sup[alt] += input_dict[umi_key]["Reverse Single"].count(alt)
            # Sees if the dict variant hits is empty, if not, counts the no molecules supporting the variant as FFPE
            if input_dict[umi_key]["Variant Hits"]:
                if input_dict[umi_key]["Variant Hits"]["FFPE Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["FFPE Hits"][alt]:
                        alt_sup[alt] += input_dict[umi_key]["Variant Hits"]["FFPE Hits"][alt].count(alt)
                # Sees if the dict variant hits is empty, if not, counts the no molecules supporting the variant as
                # a regular mutation
                if input_dict[umi_key]["Variant Mutation Hits"]["Variant Mutation Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["Variant Mutation Hits"][alt]:
                        alt_sup[alt] += input_dict[umi_key]["Variant Hits"]["Variant Mutation Hits"][alt].count(alt)
                # Sees if the dict Ref Hits is empty, if not, counts the no molecules supporting the variant as
                # a reference (ie: not mutated)
                if input_dict[umi_key]["Variant Hits"]["Reference Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["Reference Hits"][alt]:
                        alt_sup[alt] += input_dict[umi_key]["Variant Hits"]["Reference Hits"][alt].count(alt)
                # Sees if the dict Other Mutation Hits is empty, if not, counts the no molecules supporting the variant
                # as another form of mutation
                if input_dict[umi_key]["Variant Hits"]["Other Mutation Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["Other Mutation Hits"][alt]:
                        alt_sup[alt] += input_dict[umi_key]["Variant Hits"]["Other Mutation Hits"][alt].count(alt)
        # Iterates through each base present in the reference
        for ref in n_ref:
            if input_dict[umi_key]["Forward Single"]:
                ref_sup[ref] += input_dict[umi_key]["Forward Single"].count(ref)
            if input_dict[umi_key]["Reverse Single"]:
                ref_sup[ref] += input_dict[umi_key]["Reverse Single"].count(n_alt)
            if input_dict[umi_key]["Variant Hits"]:
                if input_dict[umi_key]["Variant Hits"]["FFPE Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["FFPE Hits"][ref]:
                        ref_sup[ref] += input_dict[umi_key]["Variant Hits"]["FFPE Hits"][ref].count(ref)
                if input_dict[umi_key]["Variant Hits"]["Variant Mutation Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["Variant Mutation Hits"][ref]:
                        ref_sup[ref] += input_dict[umi_key]["Variant Hits"]["Variant Mutation Hits"][ref].count(ref)
                if input_dict[umi_key]["Variant Hits"]["Reference Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["Reference Hits"][ref]:
                        ref_sup[ref] += input_dict[umi_key]["Variant Hits"]["Reference Hits"][ref].count(ref)
                if input_dict[umi_key]["Variant Hits"]["Other Mutation Hits"]:
                    if input_dict[umi_key]["Variant Hits"]["Other Mutation Hits"][ref]:
                        ref_sup[ref] += input_dict[umi_key]["Variant Hits"]["Other Mutation Hits"][ref].count(ref)
        pos_sup[umi_key] = {"Variant": alt_sup, "Reference": ref_sup}
    return pos_sup


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
                umi_add = {"Forward_Molecule": [], "Reverse_Molecule": []}
                umi_dict[umi_id] = umi_add
                umi_dict[umi_id][strand].append(read)
        f_hits = {}
        r_hits = {}
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
                    # print("Forward_lst")
                    # print(f_lst)
                    # print("Reverse_lst")
                    # print(r_lst)
                    # print(" Variant position : " + str(record_pos))
                    base_res["Forward Reads"] = f_lst
                    base_res["Reverse Reads"] = r_lst
                    f_hits[umi_key] = pos_hits(f_lst, record_pos)
                    r_hits[umi_key] = pos_hits(r_lst, record_pos)
                    base_res["Forward Hits"] = f_hits
                    base_res["Reverse Hits"] = r_hits
                    var_dict = ffpe_finder(base_res, ref_var, ref_base, umi_key)
                    # print("Forward Molecule Hits")
                    # print(f_hits)
                    # print("Reverse Molecule hits")
                    # print(r_hits)
            elif f_lst:
                if not r_lst:
                    f_lst = umi_dict[umi_key]["Forward_Molecule"]
                    f_s_hits = pos_hits(f_lst, record_pos)
            elif r_lst:
                if not f_lst:
                    r_lst = umi_dict[umi_key]["Reverse_Molecule"]
                    r_s_hits = pos_hits(r_lst, record_pos)
            pos_res[umi_key] = {"Forward Single": f_s_hits, "Reverse Single": r_s_hits, "Variant Hits": var_dict}
            # print(pos_res)
        #     if counter == 300:
        #         break
        # exit()
    except KeyError:
        print("ERROR: The requested key does not exist")
    return pos_res


def pos_hits(inp_lst, record_pos):
    # Loops through each key and its respective reads to extract their variant position data and then counts
    # The no. hits for each respective letter for this position
    n_a = 0
    n_t = 0
    n_g = 0
    n_c = 0
    n_n = 0
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
            # # If the read is reverse, get the original read sequence
            # if read.is_reverse:
            #     print("The read aligns to the reverse strand")
            #     print("The query-sequence is : " + str(read_seq))
            #     print("The length of the sequence is : " + str(len(read_seq)))
            #     print("The base position is : " + str(ind_pos))
            #     # Retrieve the base present at the index position where the variant is located
            #     read_base = read_seq[ind_pos]
            #     print("The selected rev base is : " + read_base)
            # else:
            #     print("The read aligns to the forward strand")
            #     print("The forward-sequence is : " + str(read_seq))
            #     print("The base position is : " + str(ind_pos))
            #
            #     read_base = read_seq[ind_pos]
            #     print("The selected forward base is : " + read_base)

            # Checks the value of the base, adds to a counter based on its value
            if read_base == 'A':
                n_a += 1
            elif read_base == 'T':
                n_t += 1
            elif read_base == 'G':
                n_g += 1
            elif read_base == 'C':
                n_c += 1
            else:
                n_n += 1

        except ValueError:
            print("Error: reference position: " + str(record_pos) + " not in the short read")
    # Returns a dict with the counts for each base respectively
    pos_dict = {"A": n_a, "T": n_t, "G": n_g, "C": n_c, "N": n_n}
    return pos_dict


def ffpe_finder(base_res, ref_var, ref_base, umi_key):
    ffpe_dict = {}
    mut_dict = {}
    ref_dict = {}
    out_dict = {}

    try:
        # Checks if the umi-key is present in the reverse molecule
        # Retrieves the forward and reverse base hits for the key value
        f_hits = base_res["Forward Hits"][umi_key]
        r_hits = base_res["Reverse Hits"][umi_key]
        # Iterates through each variant called by the program
        for var_call in ref_var:
            for r_base in ref_base:
                # Sees if the variant is present in both molecules, if so, it is not classified as a FFPE
                if f_hits[var_call] > 0 and r_hits[var_call] > 0:
                    mut_dict[var_call] = {"Forward Molecule": f_hits, "Reverse Molecule": r_hits}
                # Sees if the variant is present in only one molecule, if so it is determined to be a FFPE
                elif f_hits[var_call] > 0 and r_hits[var_call] == 0:
                    ffpe_dict[var_call] = {"Forward Molecule": f_hits, "Reverse Molecule": r_hits}
                elif r_hits[var_call] > 0 and f_hits[var_call] == 0:
                    ffpe_dict[var_call] = {"Forward Molecule": f_hits, "Reverse Molecule": r_hits}
                # If none of the above if statements have been fulfilled, sees if the reference molecule has a count
                # Exceeding 0 in both of the molecules, if so it is deemed not a mutation and added to the ref_dict
                elif f_hits[r_base] > 0 and r_hits[r_base] > 0:
                    ref_dict[r_base] = {"Forward Molecule": f_hits, "Reverse Molecule": r_hits}
                # If none of the above statements are fulfilled, then it is deemed as another form of mutation
                else:
                    mut_dict[var_call] = {"Forward Molecule": f_hits, "Reverse Molecule": r_hits}
    except KeyError as e:
        print("No matching key for key: " + str(e) + " found, comparison not possible")
        # Returns both dictionaries in a list
    var_dict = {"FFPE Hits": ffpe_dict, "Variant Mutation Hits": mut_dict, "Reference Hits": ref_dict,
                "Other Mutation Hits": out_dict}
    # print(var_lst)
    return var_dict


def main():
    # Reads in the bam file
    bam_file = pysam.AlignmentFile("26-sort.REF_chr8.bam", "r")
    vcf_file = pysam.VariantFile("26-ensembl.chr8.vcf", "r")
    # Extracts all reads present in the bam file which have their regions matching to that of a variant call in the VCF
    vcf_extract(vcf_file, bam_file)
    bam_file.close()


main()
