import unittest
import numpy as np
import flap
import flap_jet_dataapi as dataapi

# --------------INDEPENDENT FUNCTIONS FROM FLAP STORAGE------------------------


def test_basic():
    ky6 = dataapi.getsignal(95531, 'JPF/DH/Y6-EMITER<VLT', options={"Check Time Equidistant": True})
    return int(np.mean(ky6.data)), ky6.get_coordinate_object("Time").mode.equidistant
    

class StandaloneTest(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(test_basic(), (773, False))

# -------------------------------------COMPATIBILITY WITH FLAP----------------------------------------------------------


def test_register():
    flap.register_data_source('JET_API',
                              get_data_func=dataapi.get_data,
                              add_coord_func=dataapi.add_coordinate)
    return 'JET_API' in flap.list_data_sources()

def test_reading():
    flap.register_data_source('JET_API',
                              get_data_func=dataapi.get_data,
                              add_coord_func=dataapi.add_coordinate)
    ky6 = flap.get_data('JET_API', name='JPF/DH/Y6-EMITER<VLT',
                      exp_id=95531,
                      object_name='JPF Voltage', options={})
    return int(np.mean(ky6.data))

class FLAPTest(unittest.TestCase):
    def test_register(self):
        self.assertTrue(test_register())

    def test_reading(self):
        self.assertEqual(test_reading(), 773)

if __name__ == '__main__':
    unittest.main()
