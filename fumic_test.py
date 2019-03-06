# !/usr/bin/env python3


# Import modules
import unittest
import fumic


class ReadCheck:

    def __init__(self, is_read1, is_read2, is_reverse, query_sequence, query_name):
        self.is_read1 = is_read1
        self.is_read2 = is_read2
        self.is_reverse = is_reverse
        self.query_sequence = query_sequence
        self.query_name = query_name

    def get_reference_positions(self, full_length):
            return list(range(0, len(self.query_sequence)+1))


class TestCase(unittest.TestCase):

    def setUp(self):
        # Creating clean reads
        f1c = ReadCheck(True, False, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "Forwardread1_AAATTT+CCCGGG")
        f2c = ReadCheck(False, True, False, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "Forwardread2_CCCGGG+AAATTT")
        r1c = ReadCheck(True, False, True,  "ATCGATCGAATCGATCGATCGATCGATCGATCG", "Reverseread1_CCCGGG+AAATTT")
        r2c = ReadCheck(False, True, True,  "ATCGATCGAATCGATCGATCGATCGATCGATCG", "Reverseread2_AAATTT+CCCGGG")

        # Creating a FFPE artefact read on position 2 in the forward direction
        f1f = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "Forwardread1_" + "AAATTT+CCCGGG")
        f2f = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "Forwardread2_" + "CCCGGG+AAATTT")
        r1f = ReadCheck(True, False, True, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "Reverseread1_" + "CCCGGG+AAATTT")
        r2f = ReadCheck(False, True, True, "ATCGATCGAATCGATCGATCGATCGATCGATCG", "Reverseread1_" + "AAATTT+CCCGGG")

        # Creating a Mutation read on position 2 in the forward direction
        f1m = ReadCheck(True, False, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "Forwardread1_" + "AAATTT+CCCGGG")
        f2m = ReadCheck(False, True, False, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "Forwardread2_" + "CCCGGG+AAATTT")
        r1m = ReadCheck(True, False, True, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "Reverseread1_" + "CCCGGG+AAATTT")
        r2m = ReadCheck(False, True, True, "ACCGATCGAATCGATCGATCGATCGATCGATCG", "Reverseread1_" + "AAATTT+CCCGGG")

        # Creates fake bam-lists
        self.clea_lst = [f1c, f2c, r1c, r2c]
        self.ffpe_lst = [f1f, f2f, r1f, r2f]
        self.muta_lst = [f1m, f2m, r1m, r2m]

        # Creates fake forward lists
        self.f_clea_lst = [f1c, r2c]
        self.f_ffpe_lst = [f1f, r2f]
        self.f_muta_lst = [f1m, r2m]

        # Creates fake reverse-lists
        self.r_clea_lst = [f2c, r1c]
        self.r_ffpe_lst = [f2f, r1f]
        self.r_muta_lst = [f2m, r1m]

        fwd_ch = {"Forward String": {"A": 0, "T": 1, "G": 0, "C": 0, "N": 0, "-": 0},
                  "Reverse String": {"A": 0, "T": 1, "G": 0, "C": 0, "N": 0, "-": 0}}
        rev_ch = {"Forward String": {"A": 0, "T": 1, "G": 0, "C": 0, "N": 0, "-": 0},
                  "Reverse String": {"A": 0, "T": 1, "G": 0, "C": 0, "N": 0, "-": 0}}

        fwd_fh = {"Forward String": {"A": 0, "T": 0, "G": 0, "C": 1, "N": 0, "-": 0},
                  "Reverse String": {"A": 0, "T": 1, "G": 0, "C": 0, "N": 0, "-": 0}}
        rev_fh = {"Forward String": {"A": 0, "T": 0, "G": 0, "C": 1, "N": 0, "-": 0},
                  "Reverse String": {"A": 0, "T": 1, "G": 0, "C": 0, "N": 0, "-": 0}}

        fwd_mh = {"Forward String": {"A": 0, "T": 0, "G": 0, "C": 1, "N": 0, "-": 0},
                  "Reverse String": {"A": 0, "T": 0, "G": 0, "C": 1, "N": 0, "-": 0}}
        rev_mh = {"Forward String": {"A": 0, "T": 0, "G": 0, "C": 1, "N": 0, "-": 0},
                  "Reverse String": {"A": 0, "T": 0, "G": 0, "C": 1, "N": 0, "-": 0}}

        self.c_base_res = {"Forward Reads": self.f_clea_lst, "Reverse Reads": self.r_clea_lst,
                           "Forward Hits": fwd_ch, "Reverse Hits": rev_ch}
        self.f_base_res = {"Forward Reads": self.f_ffpe_lst, "Reverse Reads": self.r_ffpe_lst,
                           "Forward Hits": fwd_fh, "Reverse Hits": rev_fh}
        self.m_base_res = {"Forward Reads": self.f_muta_lst, "Reverse Reads": self.r_muta_lst,
                           "Forward Hits": fwd_mh, "Reverse Hits": rev_mh}

    def test_pos_hits(self):
        # Test method for the pos-hits function
        rec_pos = 1
        ref_var = "C"
        ref_bas = "T"

        # First checks the clean list
        self.assertEqual(fumic.pos_hits(self.f_clea_lst, rec_pos)["Forward String"][ref_bas], 1)
        self.assertEqual(fumic.pos_hits(self.f_clea_lst, rec_pos)["Reverse String"][ref_bas], 1)
        self.assertEqual(fumic.pos_hits(self.r_clea_lst, rec_pos)["Forward String"][ref_bas], 1)
        self.assertEqual(fumic.pos_hits(self.r_clea_lst, rec_pos)["Reverse String"][ref_bas], 1)

        # Then checks the ffpe-lists
        self.assertEqual(fumic.pos_hits(self.f_ffpe_lst, rec_pos)["Forward String"][ref_var], 1)
        self.assertEqual(fumic.pos_hits(self.r_ffpe_lst, rec_pos)["Reverse String"][ref_bas], 1)
        self.assertEqual(fumic.pos_hits(self.f_ffpe_lst, rec_pos)["Forward String"][ref_var], 1)
        self.assertEqual(fumic.pos_hits(self.r_ffpe_lst, rec_pos)["Reverse String"][ref_bas], 1)

    def test_ffpe_finder(self):
        # Test method for the ffpe_finder function
        ref_var = "C"
        ref_bas = "T"

        # First checks the clean list
        self.assertEqual(fumic.ffpe_finder(self.c_base_res, ref_var, ref_bas)["Reference Hits"][
                "Forward Molecule"]["Forward String"][ref_bas], 1)
        self.assertEqual(fumic.ffpe_finder(self.c_base_res, ref_var, ref_bas)["Reference Hits"][
                "Forward Molecule"]["Reverse String"][ref_bas], 1)

        self.assertEqual(fumic.ffpe_finder(self.c_base_res, ref_var, ref_bas)["Reference Hits"][
                "Reverse Molecule"]["Forward String"][ref_bas], 1)
        self.assertEqual(fumic.ffpe_finder(self.c_base_res, ref_var, ref_bas)["Reference Hits"][
                "Reverse Molecule"]["Reverse String"][ref_bas], 1)

        # Then checks the ffpe-list
        self.assertEqual(fumic.ffpe_finder(self.f_base_res, ref_var, ref_bas)["FFPE Hits"][
                             "Forward Molecule"]["Forward String"][ref_var], 1)
        self.assertEqual(fumic.ffpe_finder(self.f_base_res, ref_var, ref_bas)["FFPE Hits"][
                             "Forward Molecule"]["Reverse String"][ref_bas], 1)
        self.assertEqual(fumic.ffpe_finder(self.f_base_res, ref_var, ref_bas)["FFPE Hits"][
                             "Reverse Molecule"]["Forward String"][ref_var], 1)
        self.assertEqual(fumic.ffpe_finder(self.f_base_res, ref_var, ref_bas)["FFPE Hits"][
                             "Reverse Molecule"]["Reverse String"][ref_bas], 1)

        # Then checks the mut-list
        self.assertEqual(fumic.ffpe_finder(self.m_base_res, ref_var, ref_bas)["Mutation Hits"][
                             "Forward Molecule"]["Forward String"][ref_var], 1)
        self.assertEqual(fumic.ffpe_finder(self.m_base_res, ref_var, ref_bas)["Mutation Hits"][
                             "Forward Molecule"]["Reverse String"][ref_var], 1)
        self.assertEqual(fumic.ffpe_finder(self.m_base_res, ref_var, ref_bas)["Mutation Hits"][
                             "Reverse Molecule"]["Forward String"][ref_var], 1)
        self.assertEqual(fumic.ffpe_finder(self.m_base_res, ref_var, ref_bas)["Mutation Hits"][
                             "Reverse Molecule"]["Reverse String"][ref_var], 1)

    def test_sup_count(self):
        ref_var = "C"
        ref_bas = "T"
        umi_key = "AAATTT_CCCGGG"

        c_ffpe_d = {}
        f_ffpe_d = {}
        m_ffpe_d = {}

        c_ffpe = fumic.ffpe_finder(self.c_base_res, ref_var, ref_bas)
        f_ffpe = fumic.ffpe_finder(self.f_base_res, ref_var, ref_bas)
        m_ffpe = fumic.ffpe_finder(self.m_base_res, ref_var, ref_bas)

        sing_dict = {"Forward Single": {}, "Reverse Single": {}}

        c_ffpe_d[umi_key] = {"Single Hits": sing_dict, "Variant Hits": c_ffpe}
        f_ffpe_d[umi_key] = {"Single Hits": sing_dict, "Variant Hits": f_ffpe}
        m_ffpe_d[umi_key] = {"Single Hits": sing_dict, "Variant Hits": m_ffpe}

        self.assertEqual(fumic.sup_count(c_ffpe_d, ref_bas)["Paired"]["Forward Molecule"][ref_bas], 2)
        self.assertEqual(fumic.sup_count(c_ffpe_d, ref_var)["Paired"]["Forward Molecule"][ref_var], 0)
        self.assertEqual(fumic.sup_count(c_ffpe_d, ref_bas)["Paired"]["Reverse Molecule"][ref_bas], 2)
        self.assertEqual(fumic.sup_count(c_ffpe_d, ref_var)["Paired"]["Reverse Molecule"][ref_var], 0)

        self.assertEqual(fumic.sup_count(f_ffpe_d, ref_bas)["Paired"]["Forward Molecule"][ref_bas], 1)
        self.assertEqual(fumic.sup_count(f_ffpe_d, ref_var)["Paired"]["Forward Molecule"][ref_var], 1)
        self.assertEqual(fumic.sup_count(f_ffpe_d, ref_bas)["Paired"]["Reverse Molecule"][ref_bas], 1)
        self.assertEqual(fumic.sup_count(f_ffpe_d, ref_var)["Paired"]["Reverse Molecule"][ref_var], 1)

        self.assertEqual(fumic.sup_count(m_ffpe_d, ref_bas)["Paired"]["Forward Molecule"][ref_bas], 0)
        self.assertEqual(fumic.sup_count(m_ffpe_d, ref_var)["Paired"]["Forward Molecule"][ref_var], 2)
        self.assertEqual(fumic.sup_count(m_ffpe_d, ref_bas)["Paired"]["Reverse Molecule"][ref_bas], 0)
        self.assertEqual(fumic.sup_count(m_ffpe_d, ref_var)["Paired"]["Reverse Molecule"][ref_var], 2)


    def test_pos_checker(self):
        rec_pos = 1
        ref_var = "C"
        ref_bas = "T"
        umi_key = "AAATTT_CCCGGG"

        c_ffpe_d = {}
        f_ffpe_d = {}
        m_ffpe_d = {}

        c_ffpe = fumic.ffpe_finder(self.c_base_res, ref_var, ref_bas)
        f_ffpe = fumic.ffpe_finder(self.f_base_res, ref_var, ref_bas)
        m_ffpe = fumic.ffpe_finder(self.m_base_res, ref_var, ref_bas)

        sing_dict = {"Forward Single": {}, "Reverse Single": {}}

        c_ffpe_d[umi_key] = {"Single Hits": sing_dict, "Variant Hits": c_ffpe}
        f_ffpe_d[umi_key] = {"Single Hits": sing_dict, "Variant Hits": f_ffpe}
        m_ffpe_d[umi_key] = {"Single Hits": sing_dict, "Variant Hits": m_ffpe}

        self.assertDictEqual(fumic.pos_checker(self.clea_lst, rec_pos, ref_var, ref_bas), c_ffpe_d)
        self.assertDictEqual(fumic.pos_checker(self.ffpe_lst, rec_pos, ref_var, ref_bas), f_ffpe_d)
        self.assertDictEqual(fumic.pos_checker(self.muta_lst, rec_pos, ref_var, ref_bas), m_ffpe_d)


if __name__ == '__main__':
    unittest.main()
