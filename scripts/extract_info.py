from pathlib import Path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import uproot
import os
import pandas as pd
mpl.rcParams['agg.path.chunksize'] = 10000
from extract_info_lib import *
import seaborn as sns
# Apply the default theme
sns.set_theme(palette='dark')
#sns.set_style("whitegrid")

from ipywidgets import *
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline


REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

out_folder = str(OUT_DIR)

PDG_CODES = {22: 'gamma', 11: 'e', 13: 'mu', 211: 'pi',
              321: 'K', 2212: 'p', 2112: 'n', 
              1000010020: 'deuteron', 1000010030: 'triton', 1000020040: 'alpha', 1000020030: 'He3'}

# ------------------ Take data from root ----------------------

#file_path = "crpa/events45.root"
file_path_CC = REPO_ROOT / "data/flat.root"
file_path_ES = REPO_ROOT / "data/flat_ES.root"
    
_, marleyarr_CC = get_data(file_path_CC)

_, marleyarr_ES = get_data(file_path_ES)

# ----------------------------------------

# ------------------ Plot Generated Energy Spectrum ----------------------

plot_energy_spectrum(marleyarr_CC, outfolder = out_folder, plot_name = "generated_energy_spectrum_CC.png")
plot_energy_spectrum(marleyarr_ES, outfolder = out_folder, plot_name = "generated_energy_spectrum_ES.png")

# ----------------------------------------




# ------------------ Plot Fluxes from Snowglobes ----------------------
liv_energies_CC, liv_spectrum_CC, liv_norm_CC = get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "CC", model = "livermore")
gvkm_energies_CC, gvkm_spectrum_CC, gvkm_norm_CC = get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "CC", model = "gvkm")
pinched_energies_CC, pinched_spectrum_CC, pinched_norm_CC = get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "CC", model = "pinched")

liv_energies_ES, liv_spectrum_ES, liv_norm_ES = get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "ES", model = "livermore")
gvkm_energies_ES, gvkm_spectrum_ES, gvkm_norm_ES = get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "ES", model = "gvkm")
pinched_energies_ES, pinched_spectrum_ES, pinched_norm_ES = get_fluxes(snowglobes_path = "interacted-fluxes/", mode = "ES", model = "pinched")

plt.figure()

# plt.plot(liv_energies_CC, liv_spectrum_CC/liv_norm_CC, label='Livermore')
plt.plot(gvkm_energies_CC, gvkm_spectrum_CC/gvkm_norm_CC, label='GKVM')
# plt.plot(pinched_energies_CC, pinched_spectrum_CC/pinched_norm_CC, label='Pinched')
plt.xlabel('Energy [MeV]')
plt.ylabel('Flux [arb. units]')
plt.title('Interacted Fluxes - CC')
plt.legend()
plt.savefig(out_folder + "interacted_fluxes_CC_norm.png")
plt.clf()

# plt.plot(liv_energies_ES, liv_spectrum_ES/liv_norm_ES, label='Livermore')
plt.plot(gvkm_energies_ES, gvkm_spectrum_ES/gvkm_norm_ES, label='GKVM')
# plt.plot(pinched_energies_ES, pinched_spectrum_ES/pinched_norm_ES, label='Pinched')
plt.xlabel('Energy [MeV]')
plt.ylabel('Flux [arb. units]')
plt.title('Interacted Fluxes - ES')
plt.legend()
plt.savefig(out_folder + "interacted_fluxes_ES_norm.png")
plt.clf()

plt.plot(gvkm_energies_ES, gvkm_spectrum_ES/gvkm_norm_ES, label="ES")
plt.plot(gvkm_energies_CC, gvkm_spectrum_CC/gvkm_norm_CC, label='CC')
plt.xlabel('Energy [MeV]')
plt.ylabel('Flux [arb. units]')
plt.title('Interacted Fluxes - GKVM')
plt.legend(fontsize=18)
plt.savefig(out_folder + "interacted_fluxes_both_norm.png")
plt.clf()

plt.plot(liv_energies_CC, liv_spectrum_CC, label='Livermore')
plt.plot(gvkm_energies_CC, gvkm_spectrum_CC, label='GKVM')
plt.plot(pinched_energies_CC, pinched_spectrum_CC, label='Pinched')
plt.xlabel('Energy [MeV]')
plt.ylabel('Flux [arb. units]')
plt.title('Interacted Fluxes - CC')
plt.legend()
plt.savefig(out_folder + "interacted_fluxes_CC.png")
plt.clf()

plt.plot(liv_energies_ES, liv_spectrum_ES, label='Livermore')
plt.plot(gvkm_energies_ES, gvkm_spectrum_ES, label='GKVM')
plt.plot(pinched_energies_ES, pinched_spectrum_ES, label='Pinched')
plt.xlabel('Energy [MeV]')
plt.ylabel('Flux [arb. units]')
plt.title('Interacted Fluxes - ES')
plt.legend()
plt.savefig(out_folder + "interacted_fluxes_ES.png")
plt.clf()

# ----------------------------------------

# ------------------ Create GKVM weights ----------------------

gkvm_weights_CC = create_weights(marleyarr_CC, gvkm_energies_CC, gvkm_spectrum_CC, gvkm_norm_CC)
gkvm_weights_ES = create_weights(marleyarr_ES, gvkm_energies_ES, gvkm_spectrum_ES, gvkm_norm_ES)
print(gkvm_weights_CC.shape)
print(gkvm_weights_ES.shape)
plt.plot(marleyarr_CC["Ev"], gkvm_weights_CC, label='CC', linewidth=0., marker='o', markersize=0.5)
plt.plot(marleyarr_ES["Ev"], gkvm_weights_ES, label='ES', linewidth=0., marker='o', markersize=0.5)
plt.xlabel('Energy [MeV]')
plt.ylabel('Weights')
plt.title('GKVM Weights')
plt.legend()
plt.savefig(out_folder + "gkvm_weights.png")
plt.clf()

# ----------------------------------------

# ------------------ Create COS ----------------------

mod_pl_CC = np.sqrt(marleyarr_CC["pxl"]**2 + marleyarr_CC["pyl"]**2 + marleyarr_CC["pzl"]**2)
cos_thetal_CC = marleyarr_CC["pzl"]/mod_pl_CC
mod_pl_ES = np.sqrt(marleyarr_ES["pxr"]**2 + marleyarr_ES["pyr"]**2 + marleyarr_ES["pzr"]**2)
cos_thetal_ES = marleyarr_ES["pzr"]/mod_pl_ES

plt.hist(cos_thetal_CC, bins=20, range=(-1, 1), label='CC', weights=gkvm_weights_CC, alpha=0.5)
plt.hist(cos_thetal_ES, bins=20, range=(-1, 1), label='ES', weights=gkvm_weights_ES, alpha=0.5)
plt.xlabel('cos($\\theta$)')
plt.ylabel('Counts')
plt.title('cos($\\theta$) - GKVM')
plt.yscale('log')
plt.legend(fontsize=18)
plt.savefig(out_folder + "cos_theta_gkvm_weights.png")
plt.clf()


plt.hist(cos_thetal_CC, bins=20, range=(-1, 1), label='CC', alpha=0.5)
plt.hist(cos_thetal_ES, bins=20, range=(-1, 1), label='ES', alpha=0.5)
plt.xlabel('cos($\\theta$)')
plt.ylabel('Counts')
plt.title('cos($\\theta$) - GKVM Weights')
plt.yscale('log')
plt.legend()
plt.savefig(out_folder + "cos_theta_gkvm.png")
plt.clf()
# ----------------------------------------

# ------------------ Extract cos < 0 ----------------------

marleyarr_CC_neg_Ev = marleyarr_CC["Ev"][cos_thetal_CC < 0]
marleyarr_ES_neg_Ev = marleyarr_ES["Ev"][cos_thetal_ES < 0]

gkvm_weights_CC_neg = gkvm_weights_CC[cos_thetal_CC < 0]
gkvm_weights_ES_neg = gkvm_weights_ES[cos_thetal_ES < 0]

plt.hist(marleyarr_CC_neg_Ev, bins=80, range=(0, 80), label='CC', weights=gkvm_weights_CC_neg, alpha=0.5)
plt.hist(marleyarr_ES_neg_Ev, bins=80, range=(0, 80), label='ES', weights=gkvm_weights_ES_neg, alpha=0.5)
plt.xlabel('Energy [MeV]')
plt.ylabel('Counts')
plt.title('cos($\\theta$) < 0 - GKVM Weights')
plt.yscale('log')
plt.legend()
plt.savefig(out_folder + "cos_theta_neg_gkvm_weights.png")
plt.clf()

# ----------------------------------------

# ------------------ 2D plot E vs cos ----------------------
plt.figure(figsize=(8, 5))
plt.hist2d(marleyarr_CC["Ev"], cos_thetal_CC, bins=100, range=((0, 70), (-1, 1)), weights=gkvm_weights_CC, cmap='viridis', density=True, norm=mpl.colors.LogNorm())
# plt.hist2d(marleyarr_CC["Ev"], cos_thetal_CC, bins=100, range=((0, 70), (-1, 1)), weights=gkvm_weights_CC, cmap='viridis', cmin=0.0005, density=True, norm=mpl.colors.LogNorm())
# plt.hist2d(marleyarr_CC["Ev"], cos_thetal_CC, bins=100, range=((0, 70), (-1, 1)), cmap='viridis', cmin=0.0005, density=True, norm=mpl.colors.LogNorm())

plt.colorbar(label='Counts')
plt.xlabel('Energy [MeV]')
plt.ylabel('cos($\\theta$)')
plt.title('Energy vs cos($\\theta$) - CC - GKVM Weights')
plt.savefig(out_folder + "E_vs_cos_CC_gkvm_weights.png")
plt.clf()

# plt.hist2d(marleyarr_ES["Ev"], cos_thetal_ES, bins=100, range=((5, 70), (-1, 1)), weights=gkvm_weights_ES, cmap='viridis', cmin=0.005, density=True, norm=mpl.colors.LogNorm())
plt.hist2d(marleyarr_ES["Ev"], cos_thetal_ES, bins=100, range=((0, 70), (0, 1)), weights=gkvm_weights_ES, cmap='viridis', density=True, norm=mpl.colors.LogNorm())
# plt.hist2d(marleyarr_ES["Ev"], cos_thetal_ES, bins=100, range=((0, 70), (-1, 1)), cmap='viridis', cmin=0.0005, density=True, norm=mpl.colors.LogNorm())
plt.colorbar(label='Counts')
plt.xlabel('Energy [MeV]')
plt.ylabel('cos($\\theta$)')
plt.title('Energy vs cos($\\theta$) - ES - GKVM Weights')
plt.savefig(out_folder + "E_vs_cos_ES_gkvm_weights.png")
plt.clf()

# ----------------------------------------

# ------------------ Plot number of secondary particles ----------------------

plt.figure()
plt.hist(marleyarr_CC["np"], bins=20, range=(0, 20), label='CC', alpha=0.5, weights=gkvm_weights_CC)
plt.hist(marleyarr_ES["np"], bins=20, range=(0, 20), label='ES', alpha=0.5, weights=gkvm_weights_ES)
plt.xlabel('Number of secondary particles')
plt.ylabel('Counts')
plt.title('Number of secondary particles - GKVM Weights')
plt.legend()
plt.savefig(out_folder + "n_sec_particles_gkvm_weights.png")
plt.clf()

# ----------------------------------------

# ------------------ Plot spectrum of secondary particles ----------------------

particles_CC = {}
particle_weights_CC = {}
print(marleyarr_CC["np"].shape)
for i in range(marleyarr_CC["np"].shape[0]):
    for j in range(marleyarr_CC["np"][i]):
        if not marleyarr_CC["pdgp"][i][j] in particles_CC:
            particles_CC[marleyarr_CC["pdgp"][i][j]] = [marleyarr_CC["KEp"][i][j]]
            particle_weights_CC[marleyarr_CC["pdgp"][i][j]] = [gkvm_weights_CC[i]]
        else:
            particles_CC[marleyarr_CC["pdgp"][i][j]].append(marleyarr_CC["KEp"][i][j])
            particle_weights_CC[marleyarr_CC["pdgp"][i][j]].append(gkvm_weights_CC[i])

for key in particles_CC:
    print(key, len(particles_CC[key]), np.min(particles_CC[key]), np.max(particles_CC[key]))
    plt.hist(particles_CC[key], bins=50, range=(0, 25), label=PDG_CODES[key], alpha=0.3, weights=particle_weights_CC[key])

plt.xlabel('Energy [MeV]')
plt.ylabel('Counts')
plt.title('Energy spectrum of secondary particles - CC - GKVM Weights')
plt.legend()
plt.yscale('log')
plt.savefig(out_folder + "sec_particles_spectrum_CC_gkvm_weights.png")
plt.clf()
# create a stack plot with the particles  
sorted_keys = sorted(particles_CC.keys(), key=lambda x: -len(particles_CC[x]))
sorted_keys.reverse()

plt.hist([particles_CC[key] for key in sorted_keys], bins=50, range=(0, 25), label=[PDG_CODES[key] for key in sorted_keys], stacked=True, alpha=0.7, weights=[particle_weights_CC[key] for key in sorted_keys])
plt.xlabel('Energy [MeV]')
plt.ylabel('Counts')
plt.title('Energy spectrum of de-excitation products - CC - GKVM Weights - Stacked')
plt.legend()
plt.yscale('log')
plt.savefig(out_folder + "sec_particles_spectrum_CC_gkvm_weights_stack.png")
plt.clf()

plt.hist([particles_CC[key] for key in sorted_keys], bins=50, range=(0, 25), label=[PDG_CODES[key] for key in sorted_keys], stacked=True, alpha=0.7)
plt.xlabel('Energy [MeV]')
plt.ylabel('Counts')
plt.title('Energy spectrum of de-excitation products - CC  - Stacked')
plt.legend()
plt.yscale('log')
plt.savefig(out_folder + "sec_particles_spectrum_CC_stack.png")
plt.clf()

plt.hist(particles_CC[22], bins=50, range=(0, 25), label=PDG_CODES[22], alpha=0.3, weights=particle_weights_CC[22])
plt.xlabel('Energy [MeV]')
plt.ylabel('Counts')
plt.title('Energy spectrum of gammas - CC - GKVM Weights')
plt.legend()
plt.yscale('log')
plt.savefig(out_folder + "sec_particles_spectrum_CC_gamma_gkvm_weights.png")

# ----------------------------------------

