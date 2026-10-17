import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt
import numpy as np


def main():
    filename = "test.png"
    print("test")
    CV_library = CV(source="RR Lyr", telescope="Kepler", nobs=1)
    lc_array = np.array(CV_library.extract_lightcurves_from_tpf().flux)
    print("test")
    # CV_library.plot_lightcurves(lc_array, filename)
    m = 10
    n = 128
    t = n / 10
    sigmaT = 0.14
    print(lc_array)
    s2, s3, s4, compactS4 = CV_library.calculate_polyspectrum(lc_array, m, n, t, sigmaT)
    CV_library.plot_polyspectrum(s2, s3, compactS4, n, t)


if __name__ == "__main__":
    main()
