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

import count_function
import pos_function


def vcf_extract(record, bam_file, ffpe_b, ext_fun, spl_fun, spl_cha):
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
    res_data = pos_function.pos_checker(bam_lst, n_pos, n_alt, n_ref, ffpe_b, ext_fun, spl_fun, spl_cha)
    mpd_data = res_data[0]
    unmpd_data = res_data[1]

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

        # Checks if any record in the returned dict indicates an FFPE, if so updates the n_fil parameter
    for umi_key in mpd_data:
        if mpd_data[umi_key]["Mate_Hits"]:
            if mpd_data[umi_key]["Mate_Hits"]["FFPE_Hits"]:
                n_cop.filter.add("FFPE")
                break
    return n_cop


def inf_builder(inp, n_ref, n_alt):
    # Calls the sup_count function to obtain the support for each nucleotide called
    ref_sup = count_function.sup_count(inp, n_ref)
    alt_sup = count_function.sup_count(inp, n_alt)
    type_sup = count_function.mol_count(inp)
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


class QueueThread(threading.Thread):
    def __init__(self, vcf_file, thr_que, target=None, name=None):
        super(QueueThread, self).__init__()
        self.target = target
        self.name = name
        self.thr_que = thr_que
        self.vcf_file = vcf_file

    def run(self):
        # Populates the thr_que if not full until every record in the vcf_file has been retrieved
        for record in self.vcf_file:
            self.thr_que.put(record)


class ResultThread(threading.Thread):
    def __init__(self, bam_path, thr_que, res_que, ffpe_b, ext_fun, spl_fun, spl_cha, target=None, name=None):
        super(ResultThread, self).__init__()
        self.target = target
        self.name = name
        self.thr_que = thr_que
        self.res_que = res_que
        self.bam_path = bam_path
        self.ffpe_b = ffpe_b
        self.ext_fun = ext_fun
        self.spl_fun = spl_fun
        self.spl_cha = spl_cha

    def run(self):
        # Calls upon the function vcf_extract while the queue is not empty, stores the results in res_que if not None
        bam_file = pysam.AlignmentFile(self.bam_path, "r", check_sq=False)
        while not self.thr_que.empty():
            record = self.thr_que.get()
            n_cop = vcf_extract(record, bam_file, self.ffpe_b, self.ext_fun, self.spl_fun, self.spl_cha)
            if n_cop is not None:
                self.res_que.put(n_cop)


def main():
    t_start = time.time()

    parser = argparse.ArgumentParser(description='FUMIC - FFPE-artefact UMI-based Mapper for Imputation in '
                                                 'Cancer-sample tissue data')
    parser.add_argument('-b', '--inputBAM', help='Input BAM file (Required)', required=True)
    parser.add_argument('-v', '--inputVCF', help='Input VCF file (Required)', required=True)
    parser.add_argument('-t', '--threads', help='No. threads to run the program (Optional)', required=False, default=1)
    parser.add_argument('-qs', '--queueSize', help='Input Queue-Size (Optional)', required=False, default=0)
    parser.add_argument('-fb', '--FFPEBases', help='Choose "all" to include all base transitions in the analysis, '
                                                   '(default: C:G>T:A)', required=False, default="standard")
    parser.add_argument('-dt', '--DataType', help='Data-type: Query-name-based (default, qrn) , Tagged-based (tgg)', required=False, default="rdn")
    parser.add_argument('-sc', '--SplitCharacter', help='Split character for the UMI-tag, qrn default = +, tgg default = "" for splitting in half', required=False, default="+")

    args = vars(parser.parse_args())
    thr_que = queue.Queue(int(args["queueSize"]))
    ffpe_b = str(args["FFPEBases"])
    d_type = str(args["DataType"])
    spl_cha = str(args["SplitCharacter"])

    if d_type == "rdn":
        ext_fun = pos_function.qrn_ext
    else:
        ext_fun = pos_function.tgg_ext
        spl_cha = ""

    if spl_cha == "+":
        spl_fun = pos_function.cha_splt
    else:
        spl_fun = pos_function.hlf_splt

    res_que = queue.Queue()
    vcf_file = pysam.VariantFile(args['inputVCF'], "r")
    bam_path = args['inputBAM']
    vcf_head = vcf_file.header
    # Generates a new filter category as well as two new format categories for the generated output
    vcf_head.filters.add('FFPE', None, None, 'FFPE Artefact')
    vcf_head.formats.add("UMI", ".", "String", "UMI information for variant then reference "
                                               "Paired;Single;Paired;Single")
    vcf_head.formats.add("UUMI", ".", "String", "Unmapped UMI information for variant then reference"
                                                " Paired;Single;Paired;Single")

    n_vcf = pysam.VariantFile('fumic_output.vcf', mode='w', header=vcf_head)

    # Starts the producer thread to populate the queue
    p_que = QueueThread(name='producer', vcf_file=vcf_file, thr_que=thr_que)
    p_que.start()
    cons = [ResultThread(name='consumer', bam_path=bam_path, thr_que=thr_que, res_que=res_que, ffpe_b=ffpe_b, ext_fun=ext_fun, spl_fun=spl_fun, spl_cha=spl_cha)
            for t in range(int(args["threads"]))]
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

