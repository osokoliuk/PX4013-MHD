"""
Module      : cvpy.get_cvs.py
Description : Parser for cataclysmic variable lightcurves
Copyright   : (c) Oleksii Sokoliuk, 2026
License     : MIT
Maintainer  : oleksii.sokoliuk@mao.kiev.ua
Stability   : experimental
Portability : portable

This script finds and downloads all of the available lightcurves from Kepler,
TESS and K2 missions for a given source coordinates.
"""

import matplotlib.pyplot as plt
import numpy as np
import lightkurve as lk
import ctypes

lib = ctypes.CDLL("./cvpy/libhsffi.so")


CDouble = ctypes.c_double
CInt = ctypes.c_int
CDoublePtr = ctypes.POINTER(CDouble)

lib.hs_start.argtypes = []
lib.hs_start.restype = None

lib.hs_stop.argtypes = []
lib.hs_stop.restype = None

lib.polyspectraPython.argtypes = [
    ctypes.c_char_p,
    CDouble,
    CDouble,
    CDouble,
    CDouble,
    CDoublePtr,
    CInt,
    CDoublePtr,
    CDoublePtr,
    CDoublePtr,
    CDoublePtr,
]
lib.polyspectraPython.restype = None


class CV:
    ########################################################################
    # Initialize a class HMF (Halo Mass Function)
    # float a - scale factor value (related to redshift via a = 1/(1+z))
    # array of floats k - wavenumber in units of 1/Mpc
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # float par1, par2 - corresponding MG parameters
    # Masses - array of CDM halo masses
    ########################################################################

    def __init__(self, source, telescope, nobs):
        self.source = source
        self.telescope = telescope
        self.nobs = nobs

    def extract_lightcurves_from_tpf(self):
        lc_arr = []
        search_result = lk.search_targetpixelfile(self.source, author=self.telescope)
        search_download = search_result.download_all()[: self.nobs]
        for dataset in search_download:
            lc_arr.append(dataset.to_lightcurve())
        return lc_arr

    def plot_lightcurves(self, lc_arr, filename):
        fig, ax = plt.subplots()
        for lc in lc_arr:
            lc.plot(c="k", ax=ax)
        plt.savefig(filename)

    def calculate_polyspectrum(self, lc):

        sigmaT = 0.14
        m = 10
        n = 512
        t = n * 0.1
        length = int(t)

        a = np.empty(length, dtype=np.float64)
        b = np.empty((length, length), dtype=np.float64)
        c = np.empty((length, length, length), dtype=np.float64)
        d = np.empty((length, length), dtype=np.float64)

        lib.hs_start()

        try:
            lib.polyspectraPython(
                b"DiscreteFourier",
                m,
                n,
                t,
                sigmaT,
                lc.ctypes.data_as(CDoublePtr),
                CInt(lc.size),
                a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                d.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )
        finally:
            lib.hs_stop()

        print(a, b, c, d)
