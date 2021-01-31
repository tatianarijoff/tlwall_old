import unittest
import numpy as np
import tlwall


class TestChamber(unittest.TestCase):
    def test_default(self):
        print('\nTesting default values')
        chamber = tlwall.Chamber()
        self.assertEqual(1, chamber.pipe_len_m)
        self.assertEqual('CIRCULAR', chamber.chamber_shape)

    def test_pipe_radius(self):
        print('\nTesting pipe radius')
        chamber = tlwall.Chamber(pipe_rad_m=0.02)
        self.assertEqual(0.02, chamber.pipe_rad_m)
        self.assertEqual(0.02, chamber.pipe_hor_m)
        self.assertEqual(0.02, chamber.pipe_ver_m)

    def test_horizontal_vertical(self):
        print('\nTesting pipe horizontal and vertical dimension')
        chamber = tlwall.Chamber(pipe_hor_m=0.02, pipe_ver_m=0.01)
        self.assertEqual(0.02, chamber.pipe_hor_m)
        self.assertEqual(0.01, chamber.pipe_ver_m)
        self.assertEqual(0.01, chamber.pipe_rad_m)


if __name__ == '__main__':
    unittest.main()
