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

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import numpy as np
import lightkurve as lk
import ctypes
import smplotlib

plt.rcParams["font.family"] = "AVHershey Complex"
plt.rcParams["font.sans-serif"] = ["AVHershey Complex:heavy"]
plt.rcParams["font.serif"] = ["AVHershey Complex:heavy"]
mpl.rcParams["axes.linewidth"] = 2
plt.rcParams.update({"text.usetex": False})
colors = ["#4169e1", "#dc143c", "#2e8b58", "#d1691e", "#9d3cce"]


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
    ctypes.c_char_p,
    ctypes.c_char_p,
    CDouble,
    CDouble,
    CDouble,
    CDouble,
    CDouble,
    CDouble,
    CDoublePtr,
    CDoublePtr,
    CInt,
    CDoublePtr,
    CDoublePtr,
    CDoublePtr,
]
lib.polyspectraPython.restype = None


class CV:
    ########################################################################
    # Initialize a class CV (Cataclysmic Variable)
    # string source - cataclysmic variable name / coordinates
    # string telescope - telescope name (Kepler/K2/TESS)
    # int nobs - number of observations to be retrieved
    ########################################################################

    def __init__(self, source, telescope, nobs):
        self.source = source
        self.telescope = telescope
        self.nobs = nobs

    def extract_lightcurves_from_tpf(self):
        lc_arr = []
        if self.telescope == "":
            search_result = lk.search_targetpixelfile(self.source)
        else:
            search_result = lk.search_targetpixelfile(
                self.source, author=self.telescope
            )
        search_download = search_result.download_all()[: self.nobs]
        for dataset in search_download:
            lc_arr.append(dataset.to_lightcurve())
        return lc_arr

    def plot_lightcurves(self, lc_arr, filename):

        fig, ax = plt.subplots()
        for lc in lc_arr:
            lc.plot(c="k", ax=ax)
        plt.savefig(filename)

    def calculate_polyspectrum(
        self, lc, times, m, n, t, dt, sigmaT, forb, spectraType, estimatorType
    ):
        if spectraType == "Full":
            s2 = np.empty(n, dtype=np.float64)
            s3 = np.empty((n, n), dtype=np.float64)
            compactS4 = np.empty((n, n), dtype=np.float64)
        else:
            s2 = np.empty(n, dtype=np.float64)
            s3 = np.empty(n, dtype=np.float64)
            compactS4 = np.empty(n, dtype=np.float64)

        lib.hs_start()

        try:
            lib.polyspectraPython(
                b"DiscreteFourier",
                spectraType.encode("ascii"),
                estimatorType.encode("ascii"),
                m,
                n,
                t,
                dt,
                sigmaT,
                forb,
                lc.ctypes.data_as(CDoublePtr),
                times.ctypes.data_as(CDoublePtr),
                CInt(lc.size),
                s2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                s3.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                compactS4.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            )
        finally:
            lib.hs_stop()

        return s2, s3, compactS4

    def plot_polyspectrum(self, s2, s3, compactS4, n, t, vmin, spectraType):
        if spectraType == "Full":
            plt.cla()
            fig = plt.figure(figsize=(8, 9))
            ax = plt.subplot(321)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                labelsize=14,
                length=6,
                width=1.7,
                zorder=10,
            )
            ax.tick_params(
                axis="both",
                which="minor",
                direction="out",
                labelsize=11,
                length=4,
                width=1.1,
                zorder=10,
            )

            for spine in ax.spines.values():
                spine.set_zorder(11)

            x = np.linspace(0, int(round(t)), int(round(t)) + 1) / t

            plt.plot(x, s2, c=colors[2])

            plt.xscale("log")
            plt.yscale("log")
            plt.title(r"$\rm Power\;spectrum$")
            plt.xlabel(r"$\rm f\;(Hz)$")
            plt.ylabel(r"$\rm P(f)$")

            ax = plt.subplot(323)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())

            ax.xaxis.set_major_locator(plt.MaxNLocator(6))

            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                labelsize=14,
                length=6,
                width=1.7,
                zorder=10,
            )
            ax.tick_params(
                axis="both",
                which="minor",
                direction="out",
                labelsize=11,
                length=4,
                width=1.1,
                zorder=10,
            )

            for spine in ax.spines.values():
                spine.set_zorder(11)

            s3 = s3
            s3[s3 <= vmin] = vmin
            plt.pcolormesh(x, x, np.log10(s3), cmap="cubehelix", rasterized=True)

            plt.title(r"$\rm Bispectrum$")
            plt.xlabel(r"$\rm f_1\;(Hz)$")
            plt.ylabel(r"$\rm f_2\;(Hz)$")
            plt.colorbar()

            ax = plt.subplot(325)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())

            ax.xaxis.set_major_locator(plt.MaxNLocator(6))

            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                labelsize=14,
                length=6,
                width=1.7,
                zorder=10,
            )
            ax.tick_params(
                axis="both",
                which="minor",
                direction="out",
                labelsize=11,
                length=4,
                width=1.1,
                zorder=10,
            )

            for spine in ax.spines.values():
                spine.set_zorder(11)

            compactS4 = compactS4
            compactS4[compactS4 <= vmin] = vmin

            plt.pcolormesh(x, x, np.log10(compactS4), cmap="cubehelix", rasterized=True)

            plt.title(r"$\rm Trispectrum$")
            plt.xlabel(r"$\rm f_1\;(Hz)$")
            plt.ylabel(r"$\rm f_2\;(Hz)$")
            plt.colorbar()
            plt.tight_layout()
            plt.savefig("spectra.png", bbox_inches="tight")
        elif spectraType == "Reduced":
            plt.cla()
            fig = plt.figure(figsize=(8, 9))
            ax = plt.subplot(321)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                labelsize=14,
                length=6,
                width=1.7,
                zorder=10,
            )
            ax.tick_params(
                axis="both",
                which="minor",
                direction="out",
                labelsize=11,
                length=4,
                width=1.1,
                zorder=10,
            )

            for spine in ax.spines.values():
                spine.set_zorder(11)

            x = np.linspace(0, int(round(t)), int(round(t)) + 1) / t

            plt.plot(x, s2, c=colors[2])

            plt.xscale("log")
            plt.yscale("log")
            plt.title(r"$\rm Power\;spectrum$")
            plt.xlabel(r"$\rm f\;(Hz)$")
            plt.ylabel(r"$\rm P(f)$")

            ax = plt.subplot(323)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                labelsize=14,
                length=6,
                width=1.7,
                zorder=10,
            )
            ax.tick_params(
                axis="both",
                which="minor",
                direction="out",
                labelsize=11,
                length=4,
                width=1.1,
                zorder=10,
            )

            for spine in ax.spines.values():
                spine.set_zorder(11)

            x = np.linspace(0, int(round(t)), int(round(t)) + 1) / t

            plt.plot(x, np.log10(s3), c=colors[2])

            plt.xscale("log")
            plt.yscale("log")
            plt.title(r"$\rm Power\;spectrum$")
            plt.xlabel(r"$\rm f\;(Hz)$")
            plt.ylabel(r"$\rm P(f)$")

            ax = plt.subplot(325)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.xaxis.set_major_locator(plt.MaxNLocator(6))
            ax.yaxis.set_major_locator(plt.MaxNLocator(6))

            ax.tick_params(
                axis="both",
                which="major",
                direction="out",
                labelsize=14,
                length=6,
                width=1.7,
                zorder=10,
            )
            ax.tick_params(
                axis="both",
                which="minor",
                direction="out",
                labelsize=11,
                length=4,
                width=1.1,
                zorder=10,
            )

            for spine in ax.spines.values():
                spine.set_zorder(11)

            x = np.linspace(0, int(round(t)), int(round(t)) + 1) / t

            plt.plot(x, np.log10(compactS4), c=colors[2])

            plt.xscale("log")
            plt.yscale("log")
            plt.title(r"$\rm Power\;spectrum$")
            plt.xlabel(r"$\rm f\;(Hz)$")
            plt.ylabel(r"$\rm P(f)$")

            plt.savefig("spectra_reduced.png", bbox_inches="tight")
