import unittest
import RNA


SEQ = "CGACGUACCGUUUUGCAAAGGCGUGGCGGCCCCCAUGAACAUUGACCGUCACUGUUUCCACGUAUGUUCU"


class TestSuboptSoftConstraints(unittest.TestCase):

    def test_sc_add_up_subopt_is_filtered_by_evaluated_energy(self):
        baseline = RNA.fold_compound(SEQ).subopt(130)

        fc = RNA.fold_compound(SEQ)
        fc.sc_add_up(12, -0.01)

        # sc_add_up() contributes to unpaired nucleotides.  This is equivalent
        # to the wrapper's paired-nucleotide convention up to a structure-
        # independent energy shift, which does not affect the energy range.
        expected_energies = {
            s.structure: s.energy - (0.01 if s.structure[11] == "." else 0.0)
            for s in baseline
        }
        best = min(expected_energies.values())
        expected = {
            structure
            for structure, energy in expected_energies.items()
            if energy - best <= 1.0001
        }

        subopt = fc.subopt(100)
        structures = [s.structure for s in subopt]

        self.assertEqual(len(structures), len(set(structures)))
        self.assertEqual(set(structures), expected)
        for s in subopt:
            self.assertEqual("%6.2f" % s.energy,
                             "%6.2f" % fc.eval_structure(s.structure))

    def test_sc_add_up_subopt_has_no_duplicates(self):
        seq = "GCCAUCGACAUUGUAUACCUGUGCCCAAGA"
        lambdas = [
            0.61, 0.60, 0.79, -2.40, 1.17,
            -0.20, 0.95, -0.51, -0.37, 0.10,
            2.60, -0.77, -0.92, 0.19, -1.34,
            -0.46, -1.01, 0.07, 1.36, 0.90,
            -0.10, 0.66, -0.64, 0.41, -0.11,
            0.21, 0.50, 0.27, -1.42, 0.17,
        ]
        md = RNA.md()
        md.uniq_ML = 1
        fc = RNA.fold_compound(seq, md)
        for i, value in enumerate(lambdas):
            fc.sc_add_up(i + 1, -value)

        structures = [s.structure for s in fc.subopt(600)]
        self.assertEqual(len(structures), len(set(structures)))


if __name__ == "__main__":
    unittest.main()
