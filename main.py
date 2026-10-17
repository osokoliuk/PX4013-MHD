import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt
import numpy as np


def main():
    filename = "test.png"

    CV_library = CV(source="RR Lyr", telescope="Kepler", nobs=1)
    # lc_array = CV_library.extract_lightcurves_from_tpf()
    # CV_library.plot_lightcurves(lc_array, filename)
    x = np.linspace(0, 0.1 * 1e5, 100000)
    lc_array = np.sin(x)
    CV_library.calculate_polyspectrum(lc_array)


if __name__ == "__main__":
    main()
