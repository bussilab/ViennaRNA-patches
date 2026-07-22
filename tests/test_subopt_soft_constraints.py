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


if __name__ == "__main__":
    unittest.main()
