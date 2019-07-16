# !/usr/bin/env python3


# Import modules
import unittest
import nuc_function as nf
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
        self.ref_nuc = "C"
        self.var_nuc = "T"
        self.unk_sym = "N"
        self.del_sym = "-"
        self.umi_key = "AAATTT_CCCGGG"
        self.sing_dict = {"Pos_Str_Single": {}, "Neg_Str_Single": {}}
        self.ffpe_n_1 = "standard"
        self.ffpe_n_2 = "all"
        self.umi_spl_cha = "+"
        self.qrn_spl_cha = "_"
        self.ext_fun_1 = pf.qrn_ext
        self.ext_fun_2 = pf.rx_ext
        self.spl_fun_1 = pf.cha_splt
        self.spl_fun_2 = pf.hlf_splt

        # Creating clean paired reads
        f1c = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_AAATTT+CCCGGG")
        f1c2 = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_UMI_AAATTT+CCCGGG")
        f1ch = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_UMI_AAATTTCCCGGG")
        f2c = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_CCCGGG+AAATTT")
        r1c = ReadCheck(True, False, True,  "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_CCCGGG+AAATTT")
        r2c = ReadCheck(False, True, True,  "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedCleanRead_AAATTT+CCCGGG")
        self.f1c = f1c
        self.f1c2 = f1c2
        self.f1ch = f1ch
        self.f2c = f2c
        self.r1c = r1c
        self.r2c = r2c

        # Creating a variant read on position 2 in the forward direction
        f1v = ReadCheck(True, False, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_AAATTT+CCCGGG")
        f2v = ReadCheck(False, True, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_CCCGGG+AAATTT")
        r1v = ReadCheck(True, False, True, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_CCCGGG+AAATTT")
        r2v = ReadCheck(False, True, True, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedRefRead_AAATTT+CCCGGG")

        # Creating a FFPE artefact read on position 2 in the forward direction
        f1f = ReadCheck(True, False, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_AAATTT+CCCGGG")
        f2f = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_CCCGGG+AAATTT")
        r1f = ReadCheck(True, False, True,  "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_CCCGGG+AAATTT")
        r2f = ReadCheck(False, True, True,  "ATCGATCGAATCGATCGATCGATCGATCGATCG", "PairedFFPERead_AAATTT+CCCGGG")

        # Creating a Unknown read on position 2 in the forward direction
        f1n = ReadCheck(True, False, False, "ANCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_AAATTT+CCCGGG")
        f2n = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_CCCGGG+AAATTT")
        r1n = ReadCheck(True, False, True,  "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_CCCGGG+AAATTT")
        r2n = ReadCheck(False, True, True,  "ANCGATCGAATCGATCGATCGATCGATCGATCG", "PairedUnknRead_AAATTT+CCCGGG")

        # Creating a Deletion read on position 2 in the forward direction
        f1d = ReadCheck(True, False, False, "A-CGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_AAATTT+CCCGGG")
        f2d = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_CCCGGG+AAATTT")
        r1d = ReadCheck(True, False, True, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_CCCGGG+AAATTT")
        r2d = ReadCheck(False, True, True, "A-CGATCGAATCGATCGATCGATCGATCGATCG", "PairedDelRead_AAATTT+CCCGGG")

        # Creates fake paired bam-lists
        self.ref_lst = [f1c, f2c, r1c, r2c]
        self.var_lst = [f1v, f2v, r1v, r2v]
        self.ffpe_lst = [f1f, f2f, r1f, r2f]
        self.n_lst = [f1n, f2n, r1n, r2n]
        self.del_lst = [f1d, f2d, r1d, r2d]

        # Creates single lists
        self.sing_ref_lst = [f1c, r2c]
        self.sing_var_lst = [f1v, r2v]
        self.sing_ffpe_lst = [f1f, r2f]
        self.sing_n_lst = [f1n, r2n]
        self.sing_del_lst = [f1d, r2d]

        # Creates fake single positive strand lists
        self.str1_ref_lst = {"PairedCleanRead_AAATTT+CCCGGG": [f1c, r2c]}
        self.str1_var_lst = {"PairedMutaRead_AAATTT+CCCGGG": [f1v, r2v]}
        self.str1_ffpe_lst = {"PairedFFPERead_AAATTT+CCCGGG": [f1f, r2f]}
        self.str1_n_lst = {"PairedUnknRead_AAATTT+CCCGGG": [f1n, r2n]}
        self.str1_del_lst = {"PairedDelRead_AAATTT+CCCGGG": [f1d, r2d]}
        self.str1_nm_lst = {"PairedCleanRead_AAATTT+CCCGGG": [f1f, r2c]}

        # Creates fake single negative strand lists
        self.str2_ref_lst = {"PairedCleanRead_CCCGGG+AAATTT": [f2c, r1c]}
        self.str2_var_lst = {"PairedMutaRead_CCCGGG+AAATTT+": [f2v, r1v]}
        self.str2_ffpe_lst = {"PairedFFPERead_CCCGGG+AAATTT+": [f2f, r1f]}
        self.str2_n_lst = {"PairedNUnknRead_CCCGGG+AAATTT": [f2n, r1n]}
        self.str2_del_lst = {"PairedDelRead_CCCGGG+AAATTT": [f2d, r1d]}
        self.str2_nm_lst = {"PairedCleanRead_CCCGGG+AAATTT": [f2f, r1c]}

        # Creates fake paired pos_hits output dicts
        self.r_base_res = {"Pos_Str_Hits": "C", "Neg_Str_Hits": "C"}
        self.v_base_res = {"Pos_Str_Hits": "T", "Neg_Str_Hits": "T"}
        self.f_base_res = {"Pos_Str_Hits": "T", "Neg_Str_Hits": "C"}
        self.n_base_res = {"Pos_Str_Hits": "N", "Neg_Str_Hits": "C"}
        self.d_base_res = {"Pos_Str_Hits": "-", "Neg_Str_Hits": "C"}

        self.r_ffpe = {'Reference_Hits': {'Pos_Str': "C", "Neg_Str": "C"}, 'True_Variant_Hits': {},
                       'FFPE_Hits': {},  'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 1, 'True_Variant_Support': 0,
                       'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 0}

        self.v_ffpe = {'Reference_Hits': {}, 'True_Variant_Hits': {'Pos_Str': "T", "Neg_Str": "T"},
                       'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 0, 'True_Variant_Support': 1,
                       'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 0}

        self.f_ffpe = {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {
            'Pos_Str': "T", "Neg_Str": "C"}, 'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 0,
                       'True_Variant_Support': 0, 'FFPE_Support': 1, 'N_Support': 0, 'Del_Support': 0}

        self.n_ffpe = {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {
                       'Pos_Str': "N", "Neg_Str": "C"}, 'Del_Hits': {}, 'Reference_Support': 0,
                       'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 1, 'Del_Support': 0}

        self.d_ffpe = {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {
                       'Pos_Str': "-", "Neg_Str": "C"}, 'Reference_Support': 0, 'True_Variant_Support': 0,
                       'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 1}

        self.r_ffpe_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {'Pos_Str': 'C', 'Neg_Str': 'C'}, 'True_Variant_Hits': {}, 'FFPE_Hits': {},
             'N_Hits': {}, 'Del_Hits': {}, 'Reference_Support': 1, 'True_Variant_Support': 0, 'FFPE_Support': 0,
             'N_Support': 0, 'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg-Str_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.v_ffpe_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'True_Variant_Hits': {'Pos_Str': 'T', 'Neg_Str': 'T'}, 'FFPE_Hits': {}, 'N_Hits': {},
             'Del_Hits': {}, 'Reference_Support': 0, 'True_Variant_Support': 1, 'FFPE_Support': 0, 'N_Support': 0,
             'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {},
                              'Del_Hits': {},
                              'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.f_ffpe_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {'Pos_Str': 'T', 'Neg_Str': 'C'}, 'N_Hits': {},
             'Del_Hits': {}, 'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 1, 'N_Support': 0,
             'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {},'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.n_ffpe_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {'Pos_Str': 'N', 'Neg_Str': 'C'},
             'Del_Hits': {}, 'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 1,
             'Del_Support': 0}}},
                         {self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.d_ffpe_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
            {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {'Pos_Str': '-',
                                                                                                    'Neg_Str': 'C'},
             'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 0, 'Del_Support': 1}}},
                         {self.umi_key: {'Single_Hits': {'Pos_Str_Single': {}, 'Neg_Str_Single': {}}, 'Mate_Hits':
                             {'Reference_Hits': {}, 'True_Variant_Hits': {}, 'FFPE_Hits': {}, 'N_Hits': {}, 'Del_Hits': {},
                              'Reference_Support': 0, 'True_Variant_Support': 0, 'FFPE_Support': 0, 'N_Support': 0,
                              'Del_Support': 0}}}]

        self.str1_ffpe_r_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': 'C', 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'Pos_Str_Single': None, 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_v_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': 'T', 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'Pos_Str_Single': None, 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_f_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': 'T', 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'Pos_Str_Single': None, 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_n_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': 'N', 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'Pos_Str_Single': None, 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}}]

        self.str1_ffpe_d_d = [{self.umi_key: {'Single_Hits': {'Pos_Str_Single': '-', 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}},
                              {self.umi_key: {'Single_Hits': {'Pos_Str_Single': None, 'Neg_Str_Single': {}},
                                              'Mate_Hits': {}}}]

    def test_qrn_ext(self):
        # Test method for the qrn_ext function with slight variations to the read query-name structure
        self.assertEqual(pf.qrn_ext(self.f1c, self.qrn_spl_cha), "AAATTT+CCCGGG")
        self.assertEqual(pf.qrn_ext(self.f1c2, self.qrn_spl_cha), "AAATTT+CCCGGG")
        self.assertEqual(pf.qrn_ext(self.f1ch, self.qrn_spl_cha), "AAATTTCCCGGG")

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
        self.assertEqual(pf.umi_maker(self.f1c, pf.cha_splt(pf.qrn_ext(self.f1c, self.qrn_spl_cha), "+")),
                         ["PairedCleanRead_AAATTT+CCCGGG", "Pos_Str", "AAATTT_CCCGGG"])

    def test_umi_maker_r2f(self):
        # Tests the umi_maker function for a read 2 forward case
        self.assertEqual(pf.umi_maker(self.f2c, pf.cha_splt(pf.qrn_ext(self.f2c, self.qrn_spl_cha), "+")),
                         ["PairedCleanRead_CCCGGG+AAATTT", "Neg_Str", "AAATTT_CCCGGG"])

    def test_umi_maker_r1r(self):
        # Tests the umi_maker function for a read 1 reverse case
        self.assertEqual(pf.umi_maker(self.r1c, pf.cha_splt(pf.qrn_ext(self.r1c, self.qrn_spl_cha), "+")),
                         ["PairedCleanRead_CCCGGG+AAATTT", "Neg_Str", "AAATTT_CCCGGG"])

    def test_umi_maker_r2r(self):
        # Tests the umi_maker function for a read 2 reverse case
        self.assertEqual(pf.umi_maker(self.r2c, pf.cha_splt(pf.qrn_ext(self.r2c, self.qrn_spl_cha), "+")),
                         ["PairedCleanRead_AAATTT+CCCGGG", "Pos_Str", "AAATTT_CCCGGG"])

    def test_pos_hits_c(self):
        # Test method for the pos-hits function for "clean" data
        self.assertEqual(pf.pos_hits(self.str1_ref_lst, self.rec_pos), [self.ref_nuc, None])
        self.assertEqual(pf.pos_hits(self.str2_ref_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_f(self):
        # Test method for the pos-hits function for ffpe-artefact data
        self.assertEqual(pf.pos_hits(self.str1_ffpe_lst, self.rec_pos), [self.var_nuc, None])
        self.assertEqual(pf.pos_hits(self.str2_ffpe_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_m(self):
        # Test method for the pos-hits function for SNPs
        self.assertEqual(pf.pos_hits(self.str1_var_lst, self.rec_pos), [self.var_nuc, None])
        self.assertEqual(pf.pos_hits(self.str2_var_lst, self.rec_pos), [self.var_nuc, None])

    def test_pos_hits_n(self):
        # Test method for the pos-hits function for n in the mutated position
        self.assertEqual(pf.pos_hits(self.str1_n_lst, self.rec_pos), [self.unk_sym, None])
        self.assertEqual(pf.pos_hits(self.str2_n_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_d(self):
        # Test method for the pos-hits function for a del in the mutated position
        self.assertEqual(pf.pos_hits(self.str1_del_lst, self.rec_pos), [self.del_sym, None])
        self.assertEqual(pf.pos_hits(self.str2_del_lst, self.rec_pos), [self.ref_nuc, None])

    def test_pos_hits_nm(self):
        # Test method for the pos-hits function for a list where one of the strands differ (case 2)
        self.assertEqual(pf.pos_hits(self.str1_nm_lst, self.rec_pos), [None, None])
        self.assertEqual(pf.pos_hits(self.str2_nm_lst, self.rec_pos), ["C", None])

    def test_ffpe_finder_c(self):
        # Test method for the ffpe_finder function
        # First checks the clean list
        self.assertEqual(nf.ffpe_finder(self.r_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["Reference_Hits"][
                "Pos_Str"], self.ref_nuc)
        self.assertEqual(nf.ffpe_finder(self.r_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["Reference_Hits"][
                "Neg_Str"], self.ref_nuc)

    def test_ffpe_finder_m(self):
        # Test method for the ffpe_finder function, for the mut-list
        self.assertEqual(nf.ffpe_finder(self.v_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["True_Variant_Hits"][
                             "Pos_Str"], self.var_nuc)
        self.assertEqual(nf.ffpe_finder(self.v_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["True_Variant_Hits"][
                             "Neg_Str"], self.var_nuc)

    def test_ffpe_finder_f(self):
        # Test method for the ffpe_finder function
        self.assertEqual(nf.ffpe_finder(self.f_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)[
                             "FFPE_Hits"]["Pos_Str"], self.var_nuc)
        self.assertEqual(nf.ffpe_finder(self.f_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)[
                             "FFPE_Hits"]["Neg_Str"], self.ref_nuc)

    def test_ffpe_finder_n(self):
        # Test method for the ffpe_finder function, for the n-list
        self.assertEqual(nf.ffpe_finder(self.n_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["N_Hits"][
                             "Pos_Str"], self.unk_sym)
        self.assertEqual(nf.ffpe_finder(self.n_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["N_Hits"][
                             "Neg_Str"], self.ref_nuc)

    def test_ffpe_finder_d(self):
        # Test method for the ffpe_finder function, for the n-list
        self.assertEqual(nf.ffpe_finder(self.d_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["Del_Hits"][
                             "Pos_Str"], self.del_sym)
        self.assertEqual(nf.ffpe_finder(self.d_base_res, self.var_nuc, self.ref_nuc, self.ffpe_n_1)["Del_Hits"][
                             "Neg_Str"], self.ref_nuc)

    def test_var_extract_c(self):
        self.assertDictEqual(buf.var_extract(self.ref_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                                             self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0], self.r_ffpe_d[0])

    def test_var_extract_m(self):
        self.assertDictEqual(buf.var_extract(self.var_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                                             self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0], self.v_ffpe_d[0])

    def test_var_extract_f(self):
        self.assertDictEqual(buf.var_extract(self.ffpe_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                                             self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0], self.f_ffpe_d[0])

    def test_var_extract_n(self):
        self.assertDictEqual(buf.var_extract(self.n_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                                             self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0], self.n_ffpe_d[0])

    def test_var_extract_d(self):
        self.assertDictEqual(buf.var_extract(self.del_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                                             self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0], self.d_ffpe_d[0])

    def test_var_extract_sing_c(self):
        self.assertDictEqual(buf.var_extract(self.sing_ref_lst, self.rec_pos, self.var_nuc, self.ref_nuc,
                                             self.ffpe_n_1, self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0],
                             self.str1_ffpe_r_d[0])

    def test_var_extract_sing_m(self):
        self.assertDictEqual(buf.var_extract(self.sing_var_lst, self.rec_pos, self.var_nuc, self.ref_nuc,
                                             self.ffpe_n_1, self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0],
                             self.str1_ffpe_v_d[0])

    def test_var_extract_sing_f(self):
        self.assertDictEqual(buf.var_extract(self.sing_ffpe_lst, self.rec_pos, self.var_nuc, self.ref_nuc,
                                             self.ffpe_n_1, self.ext_fun_1, self.spl_fun_1,
                                             self.qrn_spl_cha, self.umi_spl_cha)[0],
                             self.str1_ffpe_f_d[0])

    def test_var_extract_sing_n(self):
        self.assertDictEqual(
            buf.var_extract(self.sing_n_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                            self.ext_fun_1, self.spl_fun_1, self.qrn_spl_cha,
                            self.umi_spl_cha)[0],  self.str1_ffpe_n_d[0])

    def test_var_extract_sing_d(self):
        self.assertDictEqual(
            buf.var_extract(self.sing_del_lst, self.rec_pos, self.var_nuc, self.ref_nuc, self.ffpe_n_1,
                            self.ext_fun_1, self.spl_fun_1, self.qrn_spl_cha, self.umi_spl_cha)[0], self.str1_ffpe_d_d[0])


if __name__ == '__main__':
    unittest.main()
