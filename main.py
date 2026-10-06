import os
import sys
from get_cvs import CV
import matplotlib.pyplot as plt


def main():
    ax = plt.subplots(111)
    filename = "test.png"

    CV_library = CV(source="V344 Lyr", telescope="Kepler", nobs=5)
    lc_array = CV_library.extract_lightcurves_from_tpf()
    CV_library.plot_lightcurves(lc_array, ax, filename)


if __name__ == "__main__":
    main()
