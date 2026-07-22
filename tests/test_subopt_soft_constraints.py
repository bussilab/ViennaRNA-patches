import unittest
import RNA


class TestSuboptSoftConstraints(unittest.TestCase):

    def test_boundary_completeness(self):
        seq = "CGACGUACCGUUUUGCAAAGGCGUGGCGGCCCCCAUGAACAUUGACCGUCACUGUUUCCACGUAUGUUCU"
        baseline = RNA.fold_compound(seq).subopt(130)

        fc = RNA.fold_compound(seq)
        fc.sc_add_up(12, -0.01)

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

        solution = fc.subopt(100)
        structures = [s.structure for s in solution]

        self.assertEqual(set(structures), expected)
        self.assertEqual(len(structures), len(set(structures)))
        for s in solution:
            self.assertEqual("%6.2f" % s.energy,
                             "%6.2f" % fc.eval_structure(s.structure))

    def test_no_duplicates(self):
        seq = "GCCAUCGACAUUGUAUACCUGUGCCCAAGA"
        lambdas = [
            0.61, 0.60, 0.79, -2.40, 1.17, -0.20, 0.95, -0.51, -0.37, 0.10,
            2.60, -0.77, -0.92, 0.19, -1.34, -0.46, -1.01, 0.07, 1.36, 0.90,
            -0.10, 0.66, -0.64, 0.41, -0.11, 0.21, 0.50, 0.27, -1.42, 0.17,
        ]

        md = RNA.md()
        md.uniq_ML = 1
        fc = RNA.fold_compound(seq, md)
        for i, value in enumerate(lambdas):
            fc.sc_add_up(i + 1, -value)

        structures = [s.structure for s in fc.subopt(600)]
        self.assertEqual(len(structures), len(set(structures)))

    def test_multiloop_backtracking(self):
        seq = "AGAAAAAUGCUUUUGGUCUCCGGUAGACAAUGCGUCUUACCUGGUGUCACUCGUUGACAUGUCGUCAGCUCGAUUUAACCUAAUCCAAUGAUCCGCGUU"
        lambdas = [
            0.59814674, 0.63414093, -0.93836241, -1.47296183, 0.7413845,
            -1.12585136, 0.84930008, -0.3000255, -1.11879364, 0.1682144,
            -1.65509097, 0.24366851, 1.5856395, 1.14631001, -1.42703653,
            0.06955476, -0.82710223, -0.12502622, 0.46963834, 0.11080732,
            -1.78100577, -0.40946286, -0.20849948, 0.23935086, -0.22252288,
            -0.76962024, -0.60420856, -1.17886947, -1.34621151, 0.33801942,
            -0.87748157, 0.45155889, 0.14420838, 0.12269876, 0.3448331,
            0.0842953, 1.45964308, -1.04130463, 0.5860351, -0.97677242,
            1.62478009, -0.52916401, -1.28500793, -0.57730234, 2.47363993,
            -2.13843494, -1.19901456, 0.30513793, 0.57730671, 0.41269904,
            0.52871422, -0.485297, 1.3887759, -0.11510145, 2.84733117,
            1.68232046, 0.00299492, -0.30354263, 0.47890643, -0.58661871,
            1.09265141, -1.60499188, 0.68584452, 1.97162116, 1.0474096,
            -0.60308955, 0.18457419, -0.58659373, -0.56666659, 1.16731899,
            0.4722467, 0.26563025, 0.72977425, -0.04399878, 1.35543012,
            1.19618532, 1.00137383, -1.46020754, 0.02119784, -1.32222561,
            1.27981108, 0.6496825, -0.47264512, -0.78978465, -1.26327039,
            -0.18946186, 0.12064246, -0.08047584, -0.40224871, 0.04757904,
            -1.50177435, -1.1083075, 1.24885175, -1.16564824, 1.2691256,
            -0.09189248, -0.65950989, 0.16997237, 0.77994649,
        ]
        expected = {
            ".....(((((....(((..(((((((((.....))).)))).)).(((.......)))..((.....))........)))..(((....)))..)))))",
            ".....(((((....(((..(((((((((.....))).))).))).(((.......)))..((.....))........)))..(((....)))..)))))",
        }

        fc = RNA.fold_compound(seq)
        for i, value in enumerate(lambdas):
            fc.sc_add_up(i + 1, -value)

        mfe_structure, mfe = fc.mfe()
        solution = fc.subopt(47)
        structures = [s.structure for s in solution]

        self.assertEqual(set(structures), expected)
        self.assertEqual(len(structures), len(expected))
        self.assertIn(mfe_structure, expected)
        for s in solution:
            self.assertLessEqual(s.energy, mfe + 0.47 + 0.0001)
            self.assertAlmostEqual(s.energy,
                                   fc.eval_structure(s.structure),
                                   places=4)


if __name__ == "__main__":
    unittest.main()
