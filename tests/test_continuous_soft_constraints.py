import unittest
import RNA


SEQ = "GCGCUUCGCCGC"


def pf_with_up(e1, e2=None):
    fc = RNA.fold_compound(SEQ)
    fc.mfe()
    fc.exp_params_rescale(0.01 * fc.mfe()[1])
    fc.sc_add_up(3, e1)
    if e2 is not None:
        fc.sc_add_up(3, e2)
    return fc.pf()[1]


def pf_with_bp(e1, e2=None):
    fc = RNA.fold_compound(SEQ)
    fc.mfe()
    fc.exp_params_rescale(0.01 * fc.mfe()[1])
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


if __name__ == "__main__":
    unittest.main()
