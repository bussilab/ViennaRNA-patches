import unittest
import RNA


SEQ = "GCGCUUCGCCGC"


def pf_with_up(e1, e2=None):
    fc = RNA.fold_compound(SEQ)
    fc.mfe()
    fc.sc_add_up(3, e1)
    if e2 is not None:
        fc.sc_add_up(3, e2)
    return fc.pf()[1]


def pf_with_bp(e1, e2=None):
    fc = RNA.fold_compound(SEQ)
    fc.mfe()
    fc.sc_add_bp(2, 11, e1)
    if e2 is not None:
        fc.sc_add_bp(2, 11, e2)
    return fc.pf()[1]


class TestContinuousSoftConstraints(unittest.TestCase):

    def test_unpaired_constraints_accumulate_continuously(self):
        f_split = pf_with_up(0.003, 0.003)
        f_once = pf_with_up(0.006)

        self.assertAlmostEqual(f_split, f_once, places=10)

    def test_basepair_constraints_accumulate_continuously(self):
        f_split = pf_with_bp(0.003, 0.003)
        f_once = pf_with_bp(0.006)

        self.assertAlmostEqual(f_split, f_once, places=10)

    def test_unpatched_behavior_would_differ(self):
        f_small = pf_with_up(0.004)
        f_zero = pf_with_up(0.0)

        self.assertNotAlmostEqual(f_small, f_zero, places=10)

    def test_mfe_rounds_accumulated_unpaired_contributions(self):
        fc = RNA.fold_compound("AAAAAA")
        fc.sc_add_up(3, 0.003)
        fc.sc_add_up(3, 0.003)

        self.assertAlmostEqual(fc.mfe()[1], 0.01, places=5)
        self.assertAlmostEqual(fc.pf()[1], 0.006, places=5)

    def test_mfe_rounds_accumulated_base_pair_contributions(self):
        structure = "((((....))))"

        def energies(contributions):
            fc = RNA.fold_compound(SEQ)
            for value in contributions:
                fc.sc_add_bp(2, 11, value)
            return fc.eval_structure(structure), fc.pf()[1]

        mfe_zero, pf_zero = energies([])
        mfe_split, pf_split = energies([0.003, 0.003])
        mfe_once, pf_once = energies([0.006])

        self.assertAlmostEqual(mfe_split - mfe_zero, 0.01, places=5)
        self.assertAlmostEqual(mfe_split, mfe_once, places=5)
        self.assertAlmostEqual(pf_split, pf_once, places=5)
        self.assertNotAlmostEqual(pf_split, pf_zero, places=5)

    def test_comparative_base_pair_constraints(self):
        alignment = [SEQ, SEQ]

        def energies(contributions):
            fc = RNA.fold_compound(alignment)

            if contributions:
                # Exercise comparative initialization on a fresh fold compound.
                self.assertEqual(
                    fc.sc_add_bp_comparative(
                        [2, 2], [8, 8], [contributions[0], contributions[0]],
                    ),
                    2,
                )

            fc.mfe()
            fc.pf()

            for value in contributions[1:]:
                self.assertEqual(
                    fc.sc_add_bp_comparative(
                        [2, 2], [8, 8], [value, value],
                        RNA.OPTION_MFE | RNA.OPTION_PF,
                    ),
                    2,
                )
            structure = "(((...)))..."
            mfe = fc.eval_structure(structure)
            pf = fc.pf()[1]
            return mfe, pf

        mfe_zero, pf_zero = energies([])
        mfe_split, pf_split = energies([0.003, 0.003])
        mfe_once, pf_once = energies([0.006])

        self.assertAlmostEqual(mfe_split - mfe_zero, 0.01, places=5)
        self.assertAlmostEqual(mfe_split, mfe_once, places=5)
        self.assertAlmostEqual(pf_split, pf_once, places=5)
        self.assertNotAlmostEqual(pf_split, pf_zero, places=5)


if __name__ == "__main__":
    unittest.main()
