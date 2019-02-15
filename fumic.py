# First attempt to read in BAM as SAM
# !/usr/bin/python

# FUMIC (FFPE-artfact UMI-based Mapper for Imputation in Cancer-sample tissue data)
# By Hugo Swenson, with assistance from Patrik Smeds and Claes Edenvall
# For Klinisk Genomik, Uppsala Akademiska Sjukhus 2018

# Imports modules
import pysam

# def bam_dict(readfile):
#     new_dict = {}
#     for read in readfile:
#         fs = read.get_forward_sequence()
#         # Gets the name of the sequence
#         qn = read.query_name
#         # Splits the name based on "_" character
#         qn_spl = qn.split("_")
#         # Extracts the 1- based leftmost mapping POSition
#         l_pos = read.reference_start
#         r_pos = read.reference_end
#         # Selects the last element of the vector
#         umi = qn_spl[-1]
#         # Obtains the barcode's based on splitting the vector  using the "+" symbol
#         brc = umi.split("+")
#         umi_l = brc[0]
#         umi_r = brc[1]
#         # Checks read polarity to see whether or not it is reverse paired
#         if read.is_reverse:
#             try:
#                 s_dir = "reverse"
#                 seq_d = {"Direction": s_dir, "Left-Pos": r_pos, "Right-Pos": l_pos, "UMI-L": umi_r, "UMI-R": umi_l, "sequence": fs}
#                 new_dict[str(r_pos) + "_" + str(l_pos) + "_" + umi_r + "_" + umi_l] = seq_d
#             except KeyError:
#                 print("The requested key does not exist")
#         else:
#             try:
#                 s_dir = "forward"
#                 seq_d = {"Direction": s_dir, "Left-Pos": l_pos, "Right-Pos": r_pos, "UMI-L": umi_r, "UMI-R": umi_l, "sequence": fs}
#                 new_dict[str(l_pos) + "_" + str(r_pos) + "_" + umi_l + "_" + umi_r] = seq_d
#             except KeyError:
#                 print("The requested key does not exist")
#     return new_dict


def dict_search(in_dict):
    # initialize a null list
    unique_lst = []
    # traverse for all elements
    for key in in_dict.keys():
        # check if exists in unique_list or not
        if key not in unique_lst:
            unique_lst.append(key)
    return unique_lst


def vcf_extract(vcf_file, bam_file):
    # Initialize a null list
    # Retrieve the position from every record in the VCF file
    ffpe_data = {}
    for record in vcf_file.fetch():
        bam_lst = []
        # The position that is returned to Python is 0 - based, NOT 1 - based as in the VCF file.
        record_pos = (record.pos - 1)
        # Use the record position to fetch all reads matching it, then append to the list
        for read in bam_file.fetch('chr8', record_pos, record_pos+1):
            bam_lst.append(read)
        ffpe_data[record_pos] = (pos_checker(bam_lst, record_pos))
    return ffpe_data


def pos_checker(bam_lst, record_pos):
    umi_dict = {}
    base_res = {}
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
                    base_res["Forward hits"] = f_hits[umi_key]
                    base_res["Reverse hits"] = r_hits[umi_key]
                    f_found = ffpe_finder(base_res, record_pos, 0.9)
                    # print("Forward Molecule Hits")
                    # print(f_hits)
                    # print("Reverse Molecule hits")
                    # print(r_hits)
            if counter == 5:
                break
        exit()
    except KeyError:
        print("ERROR: The requested key does not exist")
    return base_res


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
        read_pos_en = list(enumerate(read_pos, 0))
        # print(read_pos_en)
        try:
            # Gets the index of the position the sequence maps to
            ind_pos = read_pos.index(record_pos)
            read_seq = read.query_sequence
            read_seq = list(read_seq)
            # print("The query name is : " + str(read.query_name))
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
    # Returns a dict with the percentile support of each base
    pos_dict = {"A": n_a, "T": n_t, "G": n_g, "C": n_c, "N": n_n}
    return pos_dict


def ffpe_finder(base_res, record_pos, pos_sup):
    ffpe_dict = {}
    artf_dict = {}
    for umi_key in base_res["Forward hits"]:
        try:
            if umi_key in base_res["Reverse hits"]:
                f_hits = base_res["Forward hits"][umi_key]
                r_hits = base_res["Reverse hits"][umi_key]
                f_sum = sum(f_hits.values())
                r_sum = sum(r_hits.values())
                for f_base in f_hits:
                    if (f_hits[f_base]//f_sum) != (r_hits[f_base]/r_sum):
                        if (f_hits[f_base]/f_sum) > pos_sup:
                            ffpe_dict[umi_key]["Forward Molecule"] = f_hits
                            ffpe_dict[umi_key]["Reverse Molecule"] = r_hits
                        # ffpe_dict[umi_key]["Forward Molecule"] = {"Base": f_base, "Positive support": f_hits[f_base],
                        #                                           "Negative support": (f_sum - f_hits[f_base])}
                        # r_base = max(r_hits, key=r_hits)
                        # ffpe_dict[umi_key]["Reverse Molecule"] = {"Base": r_base, "Positive support": r_hits[r_base],
                        #                                           "Negative support": (r_sum - r_hits[r_base])}
                        elif (r_hits[f_base]/r_sum) > pos_sup:
                            ffpe_dict[umi_key]["Forward Molecule"] = f_hits
                            ffpe_dict[umi_key]["Reverse Molecule"] = r_hits
                        elif (f_hits[f_base]/f_sum) < pos_sup and (f_hits[f_base]/f_sum) != 0:
                            artf_dict[umi_key]["Forward Molecule"] = f_hits
                            artf_dict[umi_key]["Reverse Molecule"] = r_hits
                        elif (r_hits[f_base]/r_sum) < pos_sup and (r_hits[f_base]/r_sum) != 0:
                            artf_dict[umi_key]["Forward Molecule"] = f_hits
                            artf_dict[umi_key]["Reverse Molecule"] = r_hits
        except KeyError:
            print(" No matching key found, comparison not possible")
    var_lst = [ffpe_dict, artf_dict]
    return var_lst


# def count_n(unique_k):
#     n_1 = 0
#     n_2 = 0
#     n_3 = 0
#     n_n = 0
#     for letter in unique_k:
#         n_c = letter.count('N')
#         if n_c == 1:
#             n_1 += 1
#         elif n_c == 2:
#             n_2 += 1
#         elif n_c == 3:
#             n_3 += 1
#         elif n_c > 3:
#             n_n += 1
#     t_val = len(unique_k) - (n_1 + n_2 + n_3 + n_n)
#     h_vec = [len(unique_k), t_val, n_1, n_2, n_3, n_n]
#     bars = ('Total', 'True', 'N', 'NN', 'NNN', 'N > 3')
#     y_pos = np.arange(len(bars))
#
#     plt.bar(y_pos, h_vec)
#     plt.ylabel('Counts')
#     plt.xticks(y_pos, h_vec)
#     plt.savefig('Bar_plot.png')


def main():
    # Reads in the bam file
    bam_file = pysam.AlignmentFile("26-sort.REF_chr8.bam", "r")
    vcf_file = pysam.VariantFile("26-ensembl.chr8.vcf", "r")
    # Extracts all reads present in the bam file which have their regions matching to that of a variant call in the VCF
    result_dat = vcf_extract(vcf_file, bam_file)
    new_lst = []
    new_bam = pysam.AlignmentFile("filtered_bam.bam", "wb", template=bam_file)
    for pos in result_dat:
        new_lst.append(result_dat[pos]["Reads"])
    for read in new_lst:
        if read.is_paired:
            new_bam.write(read)
    # Generate new dictionary of all BAM reads in the extracted dict with UMI's as identifier tags
    bam_file.close()

    # Closes the alignment file
    # pp = pprint.PrettyPrinter(indent=4)
    # pp.pprint(c_uniq)


main()
