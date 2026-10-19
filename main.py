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
    print(min(lc_array), max(lc_array))
    m = 5
    n = 1024
    t = n / 10
    sigmaT = 0.14
    print(lc_array)
    s2, s3, s4, compactS4 = CV_library.calculate_polyspectrum(lc_array, m, n, t, sigmaT)
    CV_library.plot_polyspectrum(s2, s3, compactS4, n, t, vmin=np.sqrt(9.2 / n))


if __name__ == "__main__":
    main()
