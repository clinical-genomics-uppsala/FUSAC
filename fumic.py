# !/usr/bin/env python3

# FUMIC - FFPE-artefact UMI-based Mapper for Imputation in Cancer-sample tissue data
# By Hugo Swenson, with assistance from Patrik Smeds and Claes Edenvall
# Made for Klinisk Genetik, Uppsala Akademiska Sjukhus 2019

# Imports modules
import pysam
import time
import argparse
import threading
import queue

parser = argparse.ArgumentParser(description='FUMIC - FFPE-artefact UMI-based Mapper for Imputation in '
                                                 'Cancer-sample tissue data')
parser.add_argument('-b', '--inputBAM', help='Input BAM file (Required)', required=True)
parser.add_argument('-v', '--inputVCF', help='Input VCF file (Required)', required=True)
parser.add_argument('-qs', '--queueSize', help='Input Queue-Size (Optional)', required=False, default=0)
parser.add_argument('-t', '--threads', help='No. threads to run the program (Optional)', required=False, default=1)
args = vars(parser.parse_args())

thr_que = queue.Queue(int(args["queueSize"]))
res_que = queue.Queue()


def vcf_extract(record, bam_file):
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
    # The position that is returned to Python is 0 - based, NOT 1 - based as in the VCF file.
    n_pos = (n_pos - 1)

    # Removes reference bases also present in the variant call
    # for alt in n_alt:
    #     n_ref = [r_base for r_base in n_ref if r_base != alt]

    # Use the record position to fetch all reads matching it, then append to the list
    for read in bam_file.fetch('chr8', n_pos, n_pos+1):
        bam_lst.append(read)

    # Calls the pos_checker function to obtain ffpe_data
    mpd_data = pos_checker(bam_lst, n_pos, n_alt, n_ref)[0]
    unmpd_data = pos_checker(bam_lst, n_pos, n_alt, n_ref)[1]

    mpd_inf = inf_builder(mpd_data, n_alt, n_ref)
    unmpd_inf = inf_builder(unmpd_data, n_alt, n_ref)

    for sample in n_cop.samples:
        n_cop.samples[sample]['UMI'] = str(mpd_inf[0][0]) + ";" + str(mpd_inf[0][1]) + ";" + str(mpd_inf[0][2]) + \
                                           ";" + str(mpd_inf[0][3]) + ";" + str(mpd_inf[0][4]) + ";" + mpd_inf[1] + \
                                           ";" + mpd_inf[2] + ";" + mpd_inf[3] + ";" + mpd_inf[4]
        n_cop.samples[sample]['UUMI'] = str(unmpd_inf[0][0]) + ";" + str(unmpd_inf[0][1]) + ";" + \
                                        str(unmpd_inf[0][2]) + ";" + str(unmpd_inf[0][3]) + ";" + \
                                        str(unmpd_inf[0][4]) + ";" + unmpd_inf[1] + ";" + unmpd_inf[2] \
                                        + ";" + unmpd_inf[3] + ";" + unmpd_inf[4]

        # umi_lst = []
        # for ref in pos_sup["Reference"]:
        #     ref_c = pos_sup["Variant"][ref]
        #     umi_lst.append(ref_c + ":")
        # for var in pos_sup["Variant"]:
        #     umi_lst.append(pos_sup["Variant"][var])

        # Checks if any record in the returned dict indicates an FFPE, if so updates the n_fil parameter
    for umi_key in mpd_data:
        if mpd_data[umi_key]["Mate_Hits"]:
            if mpd_data[umi_key]["Mate_Hits"]["FFPE_Hits"]:
                n_cop.filter.add("FFPE")
                break
    return n_cop


def inf_builder(inp, n_ref, n_alt):
    # Calls the sup_count function to obtain the support for each nucleotide called
    ref_sup = sup_count(inp, n_ref)
    alt_sup = sup_count(inp, n_alt)
    type_sup = mol_count(inp)
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


def pos_checker(bam_lst, record_pos, ref_var, ref_base):
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
                    mpd_b_dict = ffpe_finder(mpd_dict, ref_var, ref_base)
                    unmpd_b_dict = ffpe_finder(unmpd_dict, ref_var, ref_base)
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
    mpd_bas = None
    unmpd_bas = None
    uread_base = None
    read_base = None

    # Iterates through every query_name entry within the given UMI-key for the direction
    for query_name, read in inp_lst.items():
        # If two entries exist for the same query name (ie: forward and reverse strand), calls the base_check function
        # For both reads, then checks if these are identical or if they are different. If they are different, skips the
        # Query-name pair as this  indicates some form of error, as the same molecule should have a identical base
        if len(read) == 2:
            read_base = base_check(read[0], record_pos)
            r2b = base_check(read[1], record_pos)
            # If the read mates do not agree on which base the position represents, the data cannot be used
            if read_base != r2b:
                continue
        elif len(read) == 1:
            # Checks if the read simply is alone as its mate did not align,
            # If it lacks a mate it is categorized as unmapped.
            if read[0].mate_is_unmapped:
                uread_base = base_check(read[0], record_pos)
            else:
                read_base = base_check(read[0], record_pos)
        else:
            raise Exception("No. entries exceeds 2, cannot handle the no. reads: {}".format(len(read)))
    # Selects the most prominent base in the unmapped/mapped dict
    if read_base:
        mpd_dict = pos_count(read_base)
        mpd_bas = max(mpd_dict, key=mpd_dict.get)
    if uread_base:
        unmpd_dict = pos_count(uread_base)
        unmpd_bas = max(unmpd_dict, key=unmpd_dict.get)
    # Returns a list of the mapped and unmapped most prominent base
    bas_lst = [mpd_bas, unmpd_bas]
    return bas_lst


def base_check(read, record_pos):
    try:
        # Gets the positions the sequence maps to in the reference
        # Full length with soft clips is required for the index selection to be correct
        read_pos = read.get_reference_positions(full_length=True)
        # Gets the index of the position the sequence maps to
        ind_pos = read_pos.index(record_pos)
        # Obtains the query sequence (the sequence as it were read)
        read_seq = read.query_sequence
        read_seq = list(read_seq)
        # Gets the base present at the index position of the variant
        read_base = read_seq[ind_pos]
        return read_base
    except ValueError:
        pass


def pos_count(read_base):
    base_lst = ["A", "T", "G", "C", "N", "-"]
    base_cl = [0]*6
    # Checks the read_base against the base_lst and increases the corresponding counter if a match is found
    for base in base_lst:
        if base == read_base:
            base_cl[base_lst.index(base)] += 1
    # Returns a dict with the counts for each entry in the base_lst
    pos_dict = {"A": base_cl[0], "T": base_cl[1], "G": base_cl[2], "C": base_cl[3], "N": base_cl[4], "-": base_cl[5]}
    return pos_dict


def ffpe_finder(base_res, ref_var, ref_base):
    n_sym = "N"
    del_sym = "-"
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

    try:
        # Retrieves the forward and reverse base hits for the key value
        f_hits = base_res["String_1_Hits"]
        r_hits = base_res["String_2_Hits"]
        # Iterates through each variant called by the program
        # First checks if the variant is present in both forward and reverse molecule, if so it is a Mutation
        if f_hits != n_sym and f_hits != del_sym and r_hits != n_sym and r_hits != del_sym and f_hits != r_hits:
            ffpe_dict = {"String_1": f_hits, "String_2": r_hits}
            ffpe_ind += 1
        # If both positions are deemed to be the reference, it is added to the ref_dict
        elif f_hits == r_hits and f_hits == ref_base:
            ref_dict = {"String_1": f_hits, "String_2": r_hits}
            ref_ind += 1
        # If both positions contains the variant, it is added to the mut_dict
        elif f_hits == r_hits and f_hits == ref_var:
            mut_dict = {"String_1": f_hits, "String_2": r_hits}
            mut_ind += 1
        # If either of the positions are equal to an N or a deletion, they are deemed as a N or Del respectively
        elif f_hits == n_sym or r_hits == n_sym:
            n_dict = {"String_1": f_hits, "String_2": r_hits}
            n_ind += 1
        elif f_hits == del_sym or r_hits == del_sym:
            del_dict = {"String_1": f_hits, "String_2": r_hits}
            del_ind += 1
    except KeyError as e:
        print("No match for: " + str(e) + " found, comparison not possible")
    # Returns all dictionaries in a list along with their total counts
    var_dict = {"Reference_Hits": ref_dict, "Mutation_Hits": mut_dict, "FFPE_Hits": ffpe_dict,
                "N_Hits": n_dict, "Del_Hits": del_dict, "Reference_Support": ref_ind, "Mutation_Support": mut_ind,
                "FFPE_Support": ffpe_ind, "N_Support": n_ind, "Del_Support": del_ind}
    return var_dict


def mol_count(input_dict):
    ref_sup = 0
    mut_sup = 0
    ffpe_sup = 0
    n_sup = 0
    del_sup = 0
    # Counts the total support for each variant-option for every UMI-key entry in input_dict
    for umi_key in input_dict:
        if input_dict[umi_key]["Mate_Hits"]:
            ref_sup += input_dict[umi_key]["Mate_Hits"]["Reference_Support"]
            mut_sup += input_dict[umi_key]["Mate_Hits"]["Mutation_Support"]
            ffpe_sup += input_dict[umi_key]["Mate_Hits"]["FFPE_Support"]
            n_sup += input_dict[umi_key]["Mate_Hits"]["N_Support"]
            del_sup += input_dict[umi_key]["Mate_Hits"]["Del_Support"]
        else:
            continue
    return [ref_sup, mut_sup, ffpe_sup, n_sup, del_sup]


def sup_count(input_dict, nuc):
    n_sup = {"Paired": {}, "String_1_Single": {}, "String_2_Single": {}}
    for n in nuc:
        n_sup["Paired"] = {"String_1": {}, "String_2": {}}
        n_sup["Paired"]["String_1"][n] = 0
        n_sup["Paired"]["String_2"][n] = 0
        n_sup["String_1_Single"][n] = 0
        n_sup["String_2_Single"][n] = 0
        for umi_key in input_dict:
            # Iterates through each alternative allele called by the variant caller
            # Counts the no. molecules supporting the variant for each direction in the single hits
            if input_dict[umi_key]["Single_Hits"]:
                for sh_dir in input_dict[umi_key]["Single_Hits"]:
                    if input_dict[umi_key]["Single_Hits"][sh_dir]:
                        if n in input_dict[umi_key]["Single_Hits"][sh_dir]:
                            n_sup[sh_dir][n] += 1
            # Sees if the dict "Mate Hits" is empty, if not, counts the no molecules supporting the variant for each
            # dict within the "Mate Hits" dictionary
            if input_dict[umi_key]["Mate_Hits"]:
                new_lst = list(input_dict[umi_key]["Mate_Hits"].items())
                new_dict = dict(new_lst[0:-7])
                for v_h in new_dict:
                    if new_dict[v_h]:
                        for s_t in new_dict[v_h]:
                            if new_dict[v_h][s_t]:
                                if n in new_dict[v_h][s_t]:
                                    n_sup["Paired"][s_t][n] += 1
    return n_sup


class QueueThread(threading.Thread):
    def __init__(self, vcf_file, target=None, name=None):
        super(QueueThread, self).__init__()
        self.target = target
        self.name = name
        self.vcf_file = vcf_file

    def run(self):
        # Populates the thr_que if not full until every record in the vcf_file has been retrieved
        for record in self.vcf_file:
            thr_que.put(record)


class ResultThread(threading.Thread):
    def __init__(self, bam_path, target=None, name=None):
        super(ResultThread, self).__init__()
        self.target = target
        self.name = name
        self.bam_path = bam_path

    def run(self):
        # Calls upon the function vcf_extract while the queue is not empty, stores the results in res_que if not None
        bam_file = pysam.AlignmentFile(self.bam_path, "r")
        while not thr_que.empty():
            record = thr_que.get()
            n_cop = vcf_extract(record, bam_file)
            if n_cop is not None:
                res_que.put(n_cop)


def main():
    t_start = time.time()
    vcf_file = pysam.VariantFile(args['inputVCF'], "r")
    bam_path = args['inputBAM']
    vcf_head = vcf_file.header
    # Generates a new filter category as well as two new format categories for the generated output
    vcf_head.filters.add('FFPE', None, None, 'FFPE Artefact')
    vcf_head.formats.add("UMI", ".", "String", "UMI information for variant then reference "
                                               "Paired;Single;Paired;Single")
    vcf_head.formats.add("UUMI", ".", "String", "Unmapped UMI information for variant then reference"
                                                " Paired;Single;Paired;Single")
    n_vcf = pysam.VariantFile("fumic_output.vcf", mode='w', header=vcf_head)

    # Starts the producer thread to populate the queue
    p_que = QueueThread(name='producer', vcf_file=vcf_file)
    p_que.start()
    time.sleep(0.5)
    cons = [ResultThread(name='consumer', bam_path=bam_path) for t in range(int(args["threads"]))]
    # Starts the consumer thread to generate output from the queue
    for c in cons:
        c.start()
    p_que.join()
    for c in cons:
        c.join()
    # Writes the consumer output to the vcf-file
    while not res_que.empty():
        n_vcf.write(res_que.get())
    t_end = time.time()
    print("Total runtime: " + str(t_end - t_start) + "s")


if __name__ == "__main__":
    main()
