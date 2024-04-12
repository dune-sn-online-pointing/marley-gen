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


out_folder = "plots/"
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

PDG_CODES = {22: 'gamma', 11: 'e', 13: 'mu', 211: 'pi',
              321: 'K', 2212: 'p', 2112: 'n', 
              1000010020: 'deuteron', 1000010030: 'triton', 1000020040: 'alpha', 1000020030: 'He3'}


#file_path = "crpa/events45.root"
file_path = "marley-old/flat.root"
# file_path = "marley-old/flat_ES.root"
    
marleydata = uproot.open(file_path)["mst"]
marleyframe = marleydata.arrays(library='pd')
marleyarr = marleydata.arrays(library='np')
marleyarr["pdgp"][0]

for key in marleyarr.keys():
    continue
    print(key)

# Plot the energy distribution
fig, ax = plt.subplots()
ax.hist(marleyarr["Ev"], bins=100, range=(0, 100), label='E')

print(marleyarr["Ev"])

plt.savefig(out_folder + "marley_energy.png")


print(type(marleyframe))
print(type(marleyframe["pdgp"]))

# Count the number of events with a given pdg code, without using the value_counts() method
#marleyframe.pdgp.nunique()
print(type(marleyframe["pdgp"][0]))

# Each entry in marleyframe["pdgp"] is a list of pdg codes for the particles in the event
# Merge all these lists into a single list
pdglist = np.concatenate(marleyframe["pdgp"])
print(len(pdglist))

# Count the number of times each pdg code appears in the list
pdgcounts = np.unique(pdglist, return_counts=True)

print(pdgcounts)

# Plot the number of times each pdg code appears
fig, ax = plt.subplots()
particle_types = [PDG_CODES[pdg] for pdg in pdgcounts[0]]
#ax.bar(pdgcounts[0], pdgcounts[1])
ax.bar(particle_types, pdgcounts[1])
ax.set_yscale('log')
# ax.set_xlabel('PDG code')
# ax.set_ylabel('Number of occurrences')
plt.savefig(out_folder + "marley_pdg.png")


gamma, proton, neutron, deuteron, triton, alpha, he3 = 0, 0, 0, 0, 0, 0, 0

for i, arr in enumerate(marleyarr["pdgp"]):
    if 22 in arr:
        gamma += 1
    if 2112 in arr:
        neutron += 1
    if 2212 in arr:
        proton += 1
    if 1000010020 in arr:
        deuteron += 1
    if 1000010030 in arr:
        triton += 1
    if 1000020040 in arr:
        alpha += 1
    if 1000020030 in arr:
        he3 += 1

event_num = len(marleyarr["pdgp"])
print(gamma, proton, neutron, deuteron, triton, alpha, he3)
print(gamma/event_num * 100, proton/event_num * 100, neutron/event_num * 100, deuteron/event_num * 100, triton/event_num * 100, alpha/event_num * 100, he3/event_num * 100)

mod_pl = np.sqrt(marleyarr["pxl"]**2 + marleyarr["pyl"]**2 + marleyarr["pzl"]**2)
cos_thetal = marleyarr["pzl"]/mod_pl

# Momentum transfer
q_sq = (marleyarr["El"] - marleyarr["Ev"])**2 - ((marleyarr["pxl"] - marleyarr["pxv"])**2 + (marleyarr["pyl"] - marleyarr["pyv"])**2 + (marleyarr["pzl"] - marleyarr["pzv"])**2)
kinectic_l = marleyarr["KEl"]

q0 = marleyarr["El"] - marleyarr["Ev"]
q3 = marleyarr["pzl"] - marleyarr["pzv"]

plt.figure()
plt.hist2d(-q_sq, cos_thetal, bins=100, density=False, norm=mpl.colors.LogNorm(), cmap='gnuplot');
plt.xlabel(r'$Q^2$ (MeV$^2$)')
plt.ylabel(r'$\cos{\theta_l}$')
plt.colorbar()

plt.savefig(out_folder + "marley_q2_cos.png")
print(np.sqrt(-q_sq))
print(kinectic_l)


plt.figure()
plt.hist2d(marleyarr["Ex"], -q0, bins=100, density=True, norm=mpl.colors.LogNorm());
plt.colorbar()
plt.savefig(out_folder + "marley_ex_q0.png")
# plt.figure()
# plt.scatter(-q3, -q0, s=0.1)

cut = 10
q0_cont = q0[q0 < -cut]
prop = len(q0_cont) / len(q0)

print("MARLEY flux-averaged cross section:", marleyarr["xsec"][0])
print("MARLEY flux-averaged cross section above {} MeV: {} ".format(cut, marleyarr["xsec"][0] * prop))

print("Proportion:", prop)

unique_e = np.unique(np.round(marleyarr["Ex"], 4))
print(unique_e)

flux_averaged = marleyarr["xsec"][0]
_, mbins = np.histogram(-q0, bins=100)
bin_width = np.diff(mbins)[0]

weights = flux_averaged / event_num * 1/bin_width * np.ones(event_num);
# plt.figure(1)
# my, _, _ = plt.hist(-q0, bins=mbins, weights=weights);
my, _ = np.histogram(-q0, bins=mbins, weights=weights)

cut = 0
plt.figure()
plt.hist(cos_thetal[-q0 < cut], bins=100);
plt.xlabel(r'$\cos{\theta_l}$')
plt.figure()
plt.hist(cos_thetal[-q0 > cut], bins=100);
plt.xlabel(r'$\cos{\theta_l}$')

plt.savefig(out_folder + "marley_cos.png")

ar0_abs = -35039.9000
k0_abs = -33535.50

k0 = 0
k1 = 29.8299
k2 = 800.1431
k3 = 891.394
k4 = 1643.638
k5 = 1959.071

delta = k0_abs - ar0_abs

cl0_abs = -27560

delta_cl = cl0_abs - ar0_abs

delta + k2


# mode = "ES" # "ES" or "CC"
mode = "CC" # "ES" or "CC"

#snowglobes_path = "/Users/pbarham/OneDrive/workspace/snowglobes/out/"
snowglobes_path = "interacted-fluxes/"

if mode == "CC":
    sg_livermore = "livermore_nue_Ar40_ar40kt_events.dat"
    sg_gvkm = "gvkm_nue_Ar40_ar40kt_events.dat"
    sg_pinched = "pinched_0_paperoriginal_nue_Ar40_ar40kt_events.dat"

    sg_livermore_data = np.genfromtxt(snowglobes_path + sg_livermore, skip_footer=2)
    sg_gvkm_data = np.genfromtxt(snowglobes_path + sg_gvkm, delimiter=' ', skip_footer=2)
    sg_pinched_data = np.genfromtxt(snowglobes_path + sg_pinched, delimiter=' ', skip_footer=2)

    sg_livermore_energies = sg_livermore_data[:, 0] * 1000
    sg_livermore_spectrum = sg_livermore_data[:, 1] 
    sg_livermore_norm = np.sum(sg_livermore_spectrum) * np.diff(sg_livermore_energies)[0]

    sg_gvkm_energies = sg_gvkm_data[:, 0] * 1000
    sg_gvkm_spectrum = sg_gvkm_data[:, 1] 
    sg_gkvm_norm = np.sum(sg_gvkm_spectrum) * np.diff(sg_gvkm_energies)[0]

    sg_pinched_energies = sg_pinched_data[:, 0] * 1000
    sg_pinched_spectrum = sg_pinched_data[:, 1]
    sg_pinched_norm = np.sum(sg_pinched_spectrum) * np.diff(sg_pinched_energies)[0]

elif mode == "ES":
    neutrino_types = ["nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]
    sg_livermore = ["livermore_{}_e_ar40kt_events.dat".format(nu) for nu in neutrino_types]
    sg_gvkm = ["gvkm_{}_e_ar40kt_events.dat".format(nu) for nu in neutrino_types]
    sg_pinched = ["pinched_0_paperoriginal_{}_e_ar40kt_events.dat".format(nu) for nu in neutrino_types]

    sg_livermore_data = [np.genfromtxt(snowglobes_path + s, skip_footer=2) for s in sg_livermore]
    sg_gvkm_data = [np.genfromtxt(snowglobes_path + s, delimiter=' ', skip_footer=2) for s in sg_gvkm]
    sg_pinched_data = [np.genfromtxt(snowglobes_path + s, delimiter=' ', skip_footer=2) for s in sg_pinched]

    sg_livermore_energies = sg_livermore_data[0][:, 0] * 1000
    sg_livermore_spectrum = np.sum([s[:, 1] for s in sg_livermore_data], axis=0)
    sg_livermore_norm = np.sum(sg_livermore_spectrum) * np.diff(sg_livermore_energies)[0]

    sg_gvkm_energies = sg_gvkm_data[0][:, 0] * 1000
    sg_gvkm_spectrum = np.sum([s[:, 1] for s in sg_gvkm_data], axis=0)
    sg_gkvm_norm = np.sum(sg_gvkm_spectrum) * np.diff(sg_gvkm_energies)[0]

    sg_pinched_energies = sg_pinched_data[0][:, 0] * 1000
    sg_pinched_spectrum = np.sum([s[:, 1] for s in sg_pinched_data], axis=0)
    sg_pinched_norm = np.sum(sg_pinched_spectrum) * np.diff(sg_pinched_energies)[0]


plt.figure()
plt.plot(sg_livermore_energies, sg_livermore_spectrum/sg_livermore_norm, label="Livermore")
plt.plot(sg_gvkm_energies, sg_gvkm_spectrum/sg_gkvm_norm, label="GKVM", c='orange')
plt.plot(sg_pinched_energies, sg_pinched_spectrum/sg_pinched_norm, c='red', label="Garching")

#plt.hist(marleyarr["Ev"], bins=100, density=True, color='green', alpha=0.5, label="Livermore MARLEY");
plt.xlabel('Energy (MeV)')
plt.ylabel('Event density')
plt.legend()
plt.savefig(out_folder + "fluxes.png")
# What percentage of these interacted fluxes is above 10 MeV?
# sg_livermore_above_10 = np.sum(sg_livermore_spectrum[sg_livermore_energies > 10]) / sg_livermore_norm * np.diff(sg_livermore_energies)[0]
# sg_gvkm_above_10 = np.sum(sg_gvkm_spectrum[sg_gvkm_energies > 10]) / sg_gkvm_norm * np.diff(sg_gvkm_energies)[0]
# sg_pinched_above_10 = np.sum(sg_pinched_spectrum[sg_pinched_energies > 10]) / sg_pinched_norm * np.diff(sg_pinched_energies)[0]

# print(sg_livermore_above_10, sg_gvkm_above_10, sg_pinched_above_10)

# Plot all GKVM fluxes
if mode == "ES":
    plt.figure()
    for s in sg_gvkm_data:
        plt.plot(s[:, 0] * 1000, s[:, 1])

# For a flat spectrum of energies, let's weigh the events by a given interacted flux
bins = np.arange(0, 80, 0.5)
fig, ax = plt.subplots()
h, _, _ = ax.hist(marleyarr["Ev"], bins=bins, label='E')

# Step 0: as we know the spectrum should be identically flat, let's weigh the events by their (inverse) bin size (a 0.5 MeV bin size should work)
weights_mc = np.zeros_like(marleyarr["Ev"])
for i in range(len(h)):
    # Get the events in the bin
    bin_events = marleyarr["Ev"][(marleyarr["Ev"] > bins[i]) & (marleyarr["Ev"] < bins[i+1])]
    if bin_events.size == 0:
        continue
    # Weigh the events by the inverse bin size
    weights_mc[(marleyarr["Ev"] > bins[i]) & (marleyarr["Ev"] < bins[i+1])] = 1/len(bin_events)

fig, ax = plt.subplots()
ax.hist(marleyarr["Ev"], bins=bins, weights=weights_mc, label='E')

# Now, let's weigh the events by the interacted flux
# Step 1: Interpolate the flux so we can get the value of the flux at the energy of each event
spl = CubicSpline(sg_gvkm_energies, sg_gvkm_spectrum/sg_gkvm_norm)
weights = spl(marleyarr["Ev"]) * weights_mc

# Step 2: As the flat energy range doesn't cover the full flux (and we have the weights_mc term), normalise the weights so the add up to 1
weights = weights / np.sum(weights)

# Step 3: Weigh the events by the flux, plot the histogram
plt.figure()
plt.hist(marleyarr["Ev"], bins=bins, weights=weights, density=True, color='green', alpha=0.5, label="GKVM MARLEY");
plt.savefig(out_folder + "flux_weighted.png")
# Some histograms

fig, ax = plt.subplots()
h = ax.hist2d(marleyarr["Ev"], cos_thetal, bins=100, density=True, weights=weights, cmap='viridis', norm=mpl.colors.LogNorm());
cbar = plt.colorbar(h[3], ax=ax)

ax.set_xlabel('Energy (MeV)')
ax.set_ylabel(r'$\cos{\theta_l}$')
plt.savefig(out_folder + "flux_weighted_2d.png")

