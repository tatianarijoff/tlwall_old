import unittest
import numpy as np
from matplotlib import pyplot as plt
import tlwall


class TestCfgIo(unittest.TestCase):
    def test_onelayer(self):
        print('\nTesting one layer cfg')
        read_cfg = tlwall.CfgIo()
        chamber = read_cfg.read_chamber('input/one_layer.cfg')
        self.assertEqual('prova', chamber.component_name)
        self.assertEqual(2e-2, chamber.pipe_hor_m)
        self.assertEqual('CW', chamber.layers[0].layer_type)
        self.assertEqual('PEC', chamber.layers[1].layer_type)

    def test_frequency(self):
        print('\nTesting frequency')
        read_cfg = tlwall.CfgIo()
        freq = read_cfg.read_freq('input/one_layer.cfg')
        self.assertEqual(2, freq.fstep)

    def test_beam(self):
        print('\nTesting beam')
        read_cfg = tlwall.CfgIo()
        beam = read_cfg.read_beam('input/one_layer.cfg')
        self.assertEqual(0.01, beam.test_beam_shift)
        self.assertEqual(1.006, round(beam.gammarel, 3))
        self.assertEqual(5.3, round(beam.Ekin_MeV, 1))
        self.assertEqual(100.0, round(beam.p_MeV_c, 1))

    def test_tlwall(self):
        print('\nTesting tlwall')
        read_cfg = tlwall.CfgIo('input/one_layer.cfg')
        mywall = read_cfg.read_tlwall()
        self.assertEqual(0.01, mywall.beam.test_beam_shift)
        self.assertEqual(1.006, round(mywall.beam.gammarel, 3))
        self.assertEqual(5.3, round(mywall.beam.Ekin_MeV, 1))
        self.assertEqual(100.0, round(mywall.beam.p_MeV_c, 1))
        fig = plt.figure
        plt.plot(mywall.f, mywall.ZTrans.real)
        plt.yscale('symlog')
        plt.xscale('symlog')
        plt.show()


if __name__ == '__main__':
    unittest.main()
