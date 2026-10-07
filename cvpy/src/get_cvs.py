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
