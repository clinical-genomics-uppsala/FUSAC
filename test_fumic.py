# !/usr/bin/env python3


# Import modules
import unittest
import base_function as bf
import build_function as buf
import count_function as cf
import pos_function as pf


class ReadCheck:

    def __init__(self, is_read1, is_read2, is_reverse, query_sequence, query_name):
        self.is_read1 = is_read1
        self.is_read2 = is_read2
        self.is_reverse = is_reverse
        self.query_sequence = query_sequence
        self.query_name = query_name

    def get_reference_positions(self, full_length):
        if full_length:
            return list(range(0, len(self.query_sequence)+1))
        else:
            return list(range(0, len(self.query_sequence)+1))


class TestCase(unittest.TestCase):
    def setUp(self):
        self.rec_pos = 1
        self.var_nuc = "C"
        self.ref_nuc = "T"
        self.unk_sym = "N"
        self.del_sym = "-"
        self.umi_key = "AAATTT_CCCGGG"
        self.sing_dict = {"String_1_Single": {}, "String_2_Single": {}}
        self.ffpe_b_1 = "standard"
        self.ffpe_b_2 = "all"
        self.ext_fun_1 = pf.qrn_ext
        self.ext_fun_2 = pf.rx_ext
        self.spl_fun_1 = pf.cha_splt
        self.spl_fun_2 = pf.hlf_splt

        # Creating clean paired reads
        f1c = ReadCheck(True, False, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_AAATTT+CCCGGG")
        f1c2 = ReadCheck(True, False, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_UMI_AAATTT+CCCGGG")
        f1ch = ReadCheck(True, False, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_UMI_AAATTTCCCGGG")
        f2c = ReadCheck(False, True, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_CCCGGG+AAATTT")
        r1c = ReadCheck(True, False, True,  "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_CCCGGG+AAATTT")
        r2c = ReadCheck(False, True, True,  "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_AAATTT+CCCGGG")
        self.f1c = f1c
        self.f1c2 = f1c2
        self.f1ch = f1ch
        self.f2c = f2c
        self.r1c = r1c
        self.r2c = r2c

        # Creating a Mutation read on position 2 in the forward direction
        f1m = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_AAATTT+CCCGGG")
        f2m = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_CCCGGG+AAATTT")
        r1m = ReadCheck(True, False, True,  "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_CCCGGG+AAATTT")
        r2m = ReadCheck(False, True, True,  "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_AAATTT+CCCGGG")

        # Creating a FFPE artefact read on position 2 in the forward direction
        f1f = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_AAATTT+CCCGGG")
        f2f = ReadCheck(False, True, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_CCCGGG+AAATTT")
        r1f = ReadCheck(True, False, True, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_CCCGGG+AAATTT")
        r2f = ReadCheck(False, True, True, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_AAATTT+CCCGGG")

        # Creating a Unknown read on position 2 in the forward direction
        f1n = ReadCheck(True, False, False, "ANCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_AAATTT+CCCGGG")
        f2n = ReadCheck(False, True, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_CCCGGG+AAATTT")
        r1n = ReadCheck(True, False, True,  "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_CCCGGG+AAATTT")
        r2n = ReadCheck(False, True, True,  "ANCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_AAATTT+CCCGGG")

        # Creating a Deletion read on position 2 in the forward direction
        f1d = ReadCheck(True, False, False, "A-CGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_AAATTT+CCCGGG")
        f2d = ReadCheck(False, True, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_CCCGGG+AAATTT")
        r1d = ReadCheck(True, False, True, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_CCCGGG+AAATTT")
        r2d = ReadCheck(False, True, True, "A-CGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_AAATTT+CCCGGG")

        # Creates fake paired bam-lists
        self.clea_lst = [f1c, f2c, r1c, r2c]
        self.muta_lst = [f1m, f2m, r1m, r2m]
        self.ffpe_lst = [f1f, f2f, r1f, r2f]
        self.n_lst = [f1n, f2n, r1n, r2n]
        self.del_lst = [f1d, f2d, r1d, r2d]

        # Creates single lists
        self.sing_clea_lst = [f1c, r2c]
        self.sing_muta_lst = [f1m, r2m]
        self.sing_ffpe_lst = [f1f, r2f]
        self.sing_n_lst = [f1n, r2n]
        self.sing_del_lst = [f1d, r2d]

        # Creates fake single string_1 lists
        self.str1_clea_lst = {"PairedCleanRead_AAATTT+CCCGGG": [f1c, r2c]}
        self.str1_muta_lst = {"PairedMutaRead_AAATTT+CCCGGG": [f1m, r2m]}
        self.str1_ffpe_lst = {"PairedFFPERead_AAATTT+CCCGGG": [f1f, r2f]}
        self.str1_n_lst = {"PairedUnknRead_AAATTT+CCCGGG": [f1n, r2n]}
        self.str1_del_lst = {"PairedDelRead_AAATTT+CCCGGG": [f1d, r2d]}
        self.str1_nm_lst = {"PairedCleanRead_AAATTT+CCCGGG": [f1f, r2c]}

        # Creates fake single string_2 lists
        self.str2_clea_lst = {"PairedCleanRead_CCCGGG+AAATTT": [f2c, r1c]}
        self.str2_muta_lst = {"PairedMutaRead_CCCGGG+AAATTT+": [f2m, r1m]}
        self.str2_ffpe_lst = {"PairedFFPERead_CCCGGG+AAATTT+": [f2f, r1f]}
        self.str2_n_lst = {"PairedNUnknRead_CCCGGG+AAATTT": [f2n, r1n]}
        self.str2_del_lst = {"PairedDelRead_CCCGGG+AAATTT": [f2d, r1d]}
        self.str2_nm_lst = {"PairedCleanRead_CCCGGG+AAATTT": [f2f, r1c]}

        # Creates fake paired pos_hits output dicts
        self.c_base_res = {"String_1_Hits": "T", "String_2_Hits": "T"}
        self.m_base_res = {"String_1_Hits": "C", "String_2_Hits": "C"}
        self.f_base_res = {"String_1_Hits": "C", "String_2_Hits": "T"}
        self.n_base_res = {"String_1_Hits": "N", "String_2_Hits": "T"}
        self.d_base_res = {"String_1_Hits": "-", "String_2_Hits": "T"}

        self.c_ffpe = {'Reference_Hits': {'String_1': "T", "String_2": "T"}, 'Mutation_Hits': {},
                       'FFPE_Hits': {},  'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 1, 'Mutation_Support': 0,
                       'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 0}

        self.m_ffpe = {'Reference_Hits': {}, 'Mutation_Hits': {'String_1': "C", "String_2": "C"},
                       'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 0, 'Mutation_Support': 1,
                       'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 0}

        self.f_ffpe = {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {
            'String_1': "C", "String_2": "T"}, 'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 0,
                       'Mutation_Support': 0, 'FFPE_Support': 1, 'N_Support': 0, 'Del_Support': 0}

        self.n_ffpe = {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {
                       'String_1': "N", "String_2": "T"}, 'Del_Hits': {}, 'Reference_Support': 0,
                       'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 1, 'Del_Support': 0}

        self.d_ffpe = {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {
                       'String_1': "-", "String_2": "T"}, 'Reference_Support': 0, 'Mutation_Support': 0,
                       'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 1}

        self.c_ffpe_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {'String_1': 'T', 'String_2': 'T'}, 'Mutation_Hits': {}, 'FFPE_Hits': {},
             'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 1, 'Mutation_Support': 0, 'FFPE_Support': 0,
             'N_Support': 0, 'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.f_ffpe_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {'String_1': 'C', 'String_2': 'T'}, 'N_Hits': {},
             'Del_Hits': {}, 'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 1, 'N_Support': 0,
             'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {},'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.m_ffpe_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'Mutation_Hits': {'String_1': 'C', 'String_2': 'C'}, 'FFPE_Hits': {}, 'N_Hits': {},
             'Del_Hits': {}, 'Reference_Support': 0, 'Mutation_Support': 1, 'FFPE_Support': 0, 'N_Support': 0,
             'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits' : {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.n_ffpe_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {'String_1': 'N', 'String_2': 'T'},
             'Del_Hits': {}, 'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 1,
             'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {},'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.d_ffpe_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {'String_1': '-',
                                                                                                    'String_2': 'T'},
             'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 1}}},
                         {self.umi_key: {'Single_Hits': {'String_1_Single': {}, 'String_2_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'Mutation_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'Mutation_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.str1_ffpe_c_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': 'T', 'String_2_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'String_1_Single': None, 'String_2_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_m_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': 'C', 'String_2_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'String_1_Single': None, 'String_2_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_f_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': 'C', 'String_2_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'String_1_Single': None, 'String_2_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_n_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': 'N', 'String_2_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'String_1_Single': None, 'String_2_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_d_d = [{self.umi_key: {'Single_Hits': {'String_1_Single': '-', 'String_2_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'String_1_Single': None, 'String_2_Single': {}},
                                              'Mate_Hits': {}}}]

    def test_qrn_ext(self):
        # Test method for the qrn_ext function with slight variations to the read query-name structure
        self.assertEqual(pf.qrn_ext(self.f1c), "AAATTT+CCCGGG")
        self.assertEqual(pf.qrn_ext(self.f1c2), "AAATTT+CCCGGG")
        self.assertEqual(pf.qrn_ext(self.f1ch), "AAATTTCCCGGG")

    def test_cha_splt(self):
        # Tests the cha_splt function with different chars
        self.assertEqual(pf.cha_splt("AAATTT+CCCGGG", "+"), ["AAATTT", "CCCGGG"])
        self.assertEqual(pf.cha_splt("AAATTT_CCCGGG", "_"), ["AAATTT", "CCCGGG"])
        self.assertEqual(pf.cha_splt("AAATTT-CCCGGG", "-"), ["AAATTT", "CCCGGG"])
        self.assertEqual(pf.cha_splt("AAATTTXCCCGGG", "X"), ["AAATTT", "CCCGGG"])

    def test_hlf_splt(self):
        # Tests the hlf_splt function with different chars
        self.assertEqual(pf.hlf_splt("AAATTTCCCGGG", "+"), ["AAATTT", "CCCGGG"])
        self.assertEqual(pf.hlf_splt("AAATTTCCCGGG", "_"), ["AAATTT", "CCCGGG"])
        self.assertEqual(pf.hlf_splt("CCCGGGAAATTT", "-"), ["CCCGGG", "AAATTT"])
        self.assertEqual(pf.hlf_splt("CCCGGGAAATTT", "X"), ["CCCGGG", "AAATTT"])

    def test_umi_maker_r1f(self):
        # Tests the umi_maker function for a read 1 forward case
        self.assertEqual(pf.umi_maker(self.f1c, pf.cha_splt(pf.qrn_ext(self.f1c), "+")),
                         ["PairedCleanRead_AAATTT+CCCGGG", "String_1", "AAATTT_CCCGGG"])

    def test_umi_maker_r2f(self):
        # Tests the umi_maker function for a read 2 forward case
        self.assertEqual(pf.umi_maker(self.f2c, pf.cha_splt(pf.qrn_ext(self.f2c), "+")),
                         ["PairedCleanRead_CCCGGG+AAATTT", "String_2", "AAATTT_CCCGGG"])

    def test_umi_maker_r1r(self):
        # Tests the umi_maker function for a read 1 reverse case
        self.assertEqual(pf.umi_maker(self.r1c, pf.cha_splt(pf.qrn_ext(self.r1c), "+")),
                         ["PairedCleanRead_CCCGGG+AAATTT", "String_2", "AAATTT_CCCGGG"])

    def test_umi_maker_r2r(self):
        # Tests the umi_maker function for a read 2 reverse case
        self.assertEqual(pf.umi_maker(self.r2c, pf.cha_splt(pf.qrn_ext(self.r2c), "+")),
                         ["PairedCleanRead_AAATTT+CCCGGG", "String_1", "AAATTT_CCCGGG"])

    def test_pos_hits_c(self):
        # Test method for the pos-hits function for "clean" data
        self.assertEqual(pf.pos_hits(self.str1_clea_lst, self.rec_pos), [self.ref_nuc, None])
        self.assertEqual(pf.pos_hits(self.str2_clea_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_f(self):
        # Test method for the pos-hits function for ffpe-artefact data
        self.assertEqual(pf.pos_hits(self.str1_ffpe_lst, self.rec_pos), [self.var_nuc, None])
        self.assertEqual(pf.pos_hits(self.str2_ffpe_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_m(self):
        # Test method for the pos-hits function for SNPs
        self.assertEqual(pf.pos_hits(self.str1_muta_lst, self.rec_pos), [self.var_nuc, None])
        self.assertEqual(pf.pos_hits(self.str2_muta_lst, self.rec_pos), [self.var_nuc, None])

    def test_pos_hits_n(self):
        # Test method for the pos-hits function for n in the mutated position
        self.assertEqual(pf.pos_hits(self.str1_n_lst, self.rec_pos), [self.unk_sym, None])
        self.assertEqual(pf.pos_hits(self.str2_n_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_d(self):
        # Test method for the pos-hits function for a del in the mutated position
        self.assertEqual(pf.pos_hits(self.str1_del_lst, self.rec_pos), [self.del_sym, None])
        self.assertEqual(pf.pos_hits(self.str2_del_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_nm(self):
        # Test method for the pos-hits function for a list where one of the strings differ (case 2)
        self.assertEqual(pf.pos_hits(self.str1_nm_lst, self.rec_pos), [None, None])
        self.assertEqual(pf.pos_hits(self.str2_nm_lst, self.rec_pos), ["T", None])

    def test_ffpe_finder_c(self):
        # Test method for the ffpe_finder function
        # First checks the clean list
        self.assertEqual(bf.ffpe_finder(self.c_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["Reference_Hits"][
                "String_1"], self.ref_nuc)
        self.assertEqual(bf.ffpe_finder(self.c_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["Reference_Hits"][
                "String_2"], self.ref_nuc)

    def test_ffpe_finder_m(self):
        # Test method for the ffpe_finder function, for the mut-list
        self.assertEqual(bf.ffpe_finder(self.m_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["Mutation_Hits"][
                             "String_1"], self.var_nuc)
        self.assertEqual(bf.ffpe_finder(self.m_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["Mutation_Hits"][
                             "String_2"], self.var_nuc)

    def test_ffpe_finder_f(self):
        # Test method for the ffpe_finder function
        # Then checks the ffpe-list
        self.assertEqual(bf.ffpe_finder(self.f_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["FFPE_Hits"][
                             "String_1"], self.var_nuc)
        self.assertEqual(bf.ffpe_finder(self.f_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["FFPE_Hits"][
                             "String_2"], self.ref_nuc)

    def test_ffpe_finder_n(self):
        # Test method for the ffpe_finder function, for the n-list
        self.assertEqual(bf.ffpe_finder(self.n_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["N_Hits"][
                             "String_1"], self.unk_sym)
        self.assertEqual(bf.ffpe_finder(self.n_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["N_Hits"][
                             "String_2"], self.ref_nuc)

    def test_ffpe_finder_d(self):
        # Test method for the ffpe_finder function, for the n-list
        self.assertEqual(bf.ffpe_finder(self.d_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["Del_Hits"][
                             "String_1"], self.del_sym)
        self.assertEqual(bf.ffpe_finder(self.d_base_res, self.var_nuc, self.ref_nuc, self.ffpe_b_1)["Del_Hits"][
                             "String_2"], self.ref_nuc)

    def test_var_extract_c(self):
        self.assertDictEqual(buf.var_extract(self.clea_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                                             self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.c_ffpe_d[0])

    def test_var_extract_m(self):
        self.assertDictEqual(buf.var_extract(self.muta_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                                             self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.m_ffpe_d[0])

    def test_var_extract_f(self):
        self.assertDictEqual(buf.var_extract(self.ffpe_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                                             self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.f_ffpe_d[0])

    def test_var_extract_n(self):
        self.assertDictEqual(buf.var_extract(self.n_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                                             self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.n_ffpe_d[0])

    def test_var_extract_d(self):
        self.assertDictEqual(buf.var_extract(self.del_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                                             self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.d_ffpe_d[0])

    def test_var_extract_sing_c(self):
        self.assertDictEqual(buf.var_extract(self.sing_clea_lst, self.rec_pos, self.var_nuc, self.ref_nuc,
                                             self.ffpe_b_1, self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0],
                             self.str1_ffpe_c_d[0])

    def test_var_extract_sing_m(self):
        self.assertDictEqual(buf.var_extract(self.sing_muta_lst, self.rec_pos, self.var_nuc, self.ref_nuc,
                                             self.ffpe_b_1, self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0],
                             self.str1_ffpe_m_d[0])

    def test_var_extract_sing_f(self):
        self.assertDictEqual(buf.var_extract(self.sing_ffpe_lst, self.rec_pos, self.var_nuc, self.ref_nuc,
                                             self.ffpe_b_1, self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0],
                             self.str1_ffpe_f_d[0])

    def test_var_extract_sing_n(self):
        self.assertDictEqual(
            buf.var_extract(self.sing_n_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                            self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.str1_ffpe_n_d[0])

    def test_var_extract_sing_d(self):
        self.assertDictEqual(
            buf.var_extract(self.sing_del_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_b_1,
                            self.ext_fun_1, self.spl_fun_1, self.spl_cha_1)[0], self.str1_ffpe_d_d[0])


if __name__ == '__main__':
    unittest.main()
