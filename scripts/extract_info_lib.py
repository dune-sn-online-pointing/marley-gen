from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import uproot
import os
import pandas as pd

import seaborn as sns
# Apply the default theme
sns.set_theme(palette='dark')
#sns.set_style("whitegrid")

from ipywidgets import *
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = REPO_ROOT / "plots"
DEFAULT_OUT_DIR.mkdir(parents=True, exist_ok=True)


# ------------------ Take data from root ----------------------
def get_data(file_path):
    root_path = Path(file_path)
    if not root_path.is_absolute():
        root_path = (REPO_ROOT / root_path).resolve()

    marleydata = uproot.open(root_path)["mst"]
    marleyframe = marleydata.arrays(library='pd')
    marleyarr = marleydata.arrays(library='np')
    return marleyframe, marleyarr

# ----------------------------------------

# ------------------ Plot Generated Energy Spectrum ----------------------
def plot_energy_spectrum(marleyarr, outfolder=DEFAULT_OUT_DIR, plot_name="generated_energy_spectrum.png"):
    out_path = Path(outfolder)
    out_path.mkdir(parents=True, exist_ok=True)

    plt.figure()
    plt.hist(marleyarr["Ev"], bins=80, range=(0, 80), label='E')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Counts')
    plt.title('Generated Energy Spectrum')
    plt.legend()
    plt.savefig(out_path / plot_name)
    plt.close()

# ----------------------------------------
    
# ------------------ Extract Fluxes from Snowglobes ----------------------
def get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "CC", model = "gvkm"):
    if mode == "CC":
        sg_livermore = "livermore_nue_Ar40_ar40kt_events.dat"
        sg_gvkm = "gvkm_nue_Ar40_ar40kt_events.dat"
        sg_pinched = "pinched_0_paperoriginal_nue_Ar40_ar40kt_events.dat"
        if model == "gvkm":
            sg = np.genfromtxt(snowglobes_path + sg_gvkm, delimiter=' ', skip_footer=2)
            sg_energies = sg[:, 0] * 1000
            sg_spectrum = sg[:, 1]
            sg_norm = np.sum(sg_spectrum) * np.diff(sg_energies)[0]

        elif model == "livermore":
            sg = np.genfromtxt(snowglobes_path + sg_livermore, skip_footer=2)
            sg_energies = sg[:, 0] * 1000
            sg_spectrum = sg[:, 1]
            sg_norm = np.sum(sg_spectrum) * np.diff(sg_energies)[0]

        elif model == "pinched":
            sg = np.genfromtxt(snowglobes_path + sg_pinched, delimiter=' ', skip_footer=2)
            sg_energies = sg[:, 0] * 1000
            sg_spectrum = sg[:, 1]
            sg_norm = np.sum(sg_spectrum) * np.diff(sg_energies)[0]
        return sg_energies, sg_spectrum, sg_norm

    elif mode == "ES":
        neutrino_types = ["nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
        sg_livermore = ["livermore_{}_e_ar40kt_events.dat".format(nu) for nu in neutrino_types]
        sg_gvkm = ["gvkm_{}_e_ar40kt_events.dat".format(nu) for nu in neutrino_types]
        sg_pinched = ["pinched_0_paperoriginal_{}_e_ar40kt_events.dat".format(nu) for nu in neutrino_types]

        if model == "gvkm":
            sg_data = [np.genfromtxt(snowglobes_path + s, delimiter=' ', skip_footer=2) for s in sg_gvkm]
            sg_energies = sg_data[0][:, 0] * 1000
            sg_spectrum = np.sum([s[:, 1] for s in sg_data], axis=0)
            sg_norm = np.sum(sg_spectrum) * np.diff(sg_energies)[0]
        elif model == "livermore":
            sg_data = [np.genfromtxt(snowglobes_path + s, skip_footer=2) for s in sg_livermore]
            sg_energies = sg_data[0][:, 0] * 1000
            sg_spectrum = np.sum([s[:, 1] for s in sg_data], axis=0)
            sg_norm = np.sum(sg_spectrum) * np.diff(sg_energies)[0]
        elif model == "pinched":
            sg_data = [np.genfromtxt(snowglobes_path + s, delimiter=' ', skip_footer=2) for s in sg_pinched]
            sg_energies = sg_data[0][:, 0] * 1000
            sg_spectrum = np.sum([s[:, 1] for s in sg_data], axis=0)
            sg_norm = np.sum(sg_spectrum) * np.diff(sg_energies)[0]
            
        return sg_energies, sg_spectrum, sg_norm
    else:
        print("Mode not recognized")
        return 0, 0, 0

# ----------------------------------------
    
# ------------------ Create Weights ----------------------
def create_weights(marleyarr, energy, spectrum, norm):
    weights_mc = np.zeros_like(marleyarr["Ev"])
    bins = np.arange(0, 80, 0.5)
    for i in range(len(bins)-1):
        # Get the events in the bin
        bin_events = marleyarr["Ev"][(marleyarr["Ev"] > bins[i]) & (marleyarr["Ev"] < bins[i+1])]
        if bin_events.size == 0:
            continue
        # Weigh the events by the inverse bin size
        weights_mc[(marleyarr["Ev"] > bins[i]) & (marleyarr["Ev"] < bins[i+1])] = 1/len(bin_events)
    spl = CubicSpline(energy, spectrum/norm)
    weights = spl(marleyarr["Ev"]) * weights_mc
    weights = weights / np.sum(weights)
    return weights

# ----------------------------------------