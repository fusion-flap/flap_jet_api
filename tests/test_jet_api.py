import os
import unittest
import numpy as np
import flap
import flap_jet_api as jetapi

# --------------INDEPENDENT FUNCTIONS FROM FLAP STORAGE------------------------


def test_ppf():
    ky6 = jetapi.getsignal(95531, 'JPF/DH/Y6-EMITER<VLT', options={"Check Time Equidistant": True,
                                                                   "Datapath": "PPF",
                                                                   "Cache Data": False})
    return int(np.mean(ky6.data)), ky6.get_coordinate_object("Time").mode.equidistant

def test_sal():
    ky6 = jetapi.getsignal(95531, 'JPF/DH/Y6-EMITER<VLT', options={"Check Time Equidistant": True,
                                                                   "Datapath": "SAL",
                                                                   "Cache Data": False})
    return int(np.mean(ky6.data)), ky6.get_coordinate_object("Time").mode.equidistant

def test_cache():
    cached_file = os.sep.join(jetapi.__file__.split(os.sep)[:-1] + ['cached','ppf_ky6i_cali-mvecsei-95531.hdf5'])
    existed = os.path.exists(cached_file)
    if existed is True:
        print("Can't properly test data caching, file already exists")
    try:
        jetapi.getsignal(95531, 'ppf/ky6i/cali', options={"Check Time Equidistant": True,
                                                          "Datapath": "SAL",
                                                          "Cache Data": True,
                                                          "UID": "mvecsei"})
    except ImportError:
        jetapi.getsignal(95531, 'ppf/ky6i/cali', options={"Check Time Equidistant": True,
                                                          "Datapath": "PPF",
                                                          "Cache Data": True,
                                                          "UID": "mvecsei"})
    exists = os.path.exists(cached_file)
    if not existed:
        os.remove(cached_file)
    return exists

def test_team():
    try:
        jetapi.getsignal(95531, 'ppf/ky6i/cali', options={"Check Time Equidistant": True,
                                                          "Datapath": "SAL",
                                                          "Cache Data": False,
                                                          "UID": "KY6-team"})
        return True

    except ImportError:
        jetapi.getsignal(95531, 'ppf/ky6i/cali', options={"Check Time Equidistant": True,
                                                          "Datapath": "PPF",
                                                          "Cache Data": False,
                                                          "UID": "KY6-team"})
        return True

def test_signal_list():
    try:
        jetapi.getsignal(95531, 'KY6-CrossCalib', options={"Check Time Equidistant": True,
                                                          "Datapath": "SAL",
                                                          "Cache Data": False,
                                                          "UID": "mvecsei"})
        return True

    except ImportError:
        jetapi.getsignal(95531, 'KY6-CrossCalib', options={"Check Time Equidistant": True,
                                                          "Datapath": "PPF",
                                                          "Cache Data": False,
                                                          "UID": "mvecsei"})
        return True
    

class StandaloneTest(unittest.TestCase):
    def test_ppf(self):
        self.assertEqual(test_ppf(), (773, False))
    def test_sal(self):
        self.assertEqual(test_sal(), (773, False))
    def test_cache(self):
        self.assertTrue(test_cache())
    def test_team(self):
        self.assertTrue(test_team())
    def test_signal_list(self):
        self.assertTrue(test_signal_list())

# -------------------------------------COMPATIBILITY WITH FLAP----------------------------------------------------------


def test_register():
    flap.register_data_source('JET_API',
                              get_data_func=jetapi.get_data,
                              add_coord_func=jetapi.add_coordinate)
    return 'JET_API' in flap.list_data_sources()

def test_reading():
    flap.register_data_source('JET_API',
                              get_data_func=jetapi.get_data,
                              add_coord_func=jetapi.add_coordinate)
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
