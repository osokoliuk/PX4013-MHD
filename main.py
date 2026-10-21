import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt
import numpy as np


def main():
    filename = "test.png"
    CV_library = CV(source="V452 Lyr", telescope="Kepler", nobs=1)
    lc_array = np.array(CV_library.extract_lightcurves_from_tpf()[0].flux)
    # CV_library.plot_lightcurves(lc_array, filename)
    # print(min(lc_array), max(lc_array))
    x = np.linspace(0, 1e4 - 0.1, 100000)
    lc_array = np.sin(x) + np.cos(x) ** 2
    m = 10
    n = 1024
    t = n / 10
    sigmaT = 0.14
    forb = 0.5
    spectraType = "full"
    s2, s3, s4, compactS4 = CV_library.calculate_polyspectrum(
        lc_array, m, n, t, sigmaT, forb, spectraType
    )
    print("Done")
    CV_library.plot_polyspectrum(
        s2, s3, compactS4, n, t, vmin=0, spectraType=spectraType
    )


if __name__ == "__main__":
    main()
