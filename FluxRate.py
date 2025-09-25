# ANNIE flux histograms, produced using gsimple files
# Author: Steven Doran
# based on James Minock's scripts

print('\nLoading...\n'); 
import os          
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl
font = {'family' : 'serif', 'size' : 12 }
mpl.rc('font', **font)
mpl.rcParams['mathtext.fontset'] = 'cm' # Set the math font to Computer Modern
mpl.rcParams['legend.fontsize'] = 1
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

# ***************************************************** #
# # # # # # # # # # # Configuration # # # # # # # # # # #
# ----------------------------------------------------- #

plot_name = 'ANNIE BNB flux.png'


generate_rootfile = True           # generate root file with neutrino fluxes (if False, just plot)
rootfile_name = 'ANNIE_FLUX.root'  # will contain neutrino fluxes for the 4 relevant flavors


Full_Volume = True  # will use the entire ANNIE detector for the flux calculation (radius = 1.524m, half-height = 1.98m)

                    # If Full_Volume = False, specify FV (centered at (0,0,0), geometry assumed to be cylindrical)
FV_radius = 1       # radius [m]
FV_half_z = 1.5     # half-height [m]


TEST_RUN = False    # will only run over 20 gsimple files (used for fast debugging and testing)


# Histogram details
bin_size = 0.05        # [GeV] -> default = 0.05, corresponds to the 50 MeV bins we want for comparison to GENIE XS spline files
upper_hist_limit = 10  # largest bin [GeV] -> default = 10, corresponds to maximum value in the GENIE XS spline files

# currently, default in code for flux units set to (keep commented, won't work up here):
# scaling_factor = 1 / total_POT / surface_area    # proper unit conversion for flux
                                                   # default: v / 50MeV / POT / cm^2
# this can be changed in the script (you'll also want to change the label below to reflect the correct units)
hist_y_label = r'$\Phi (\nu)$' + ' / 50MeV / cm' + r'$^2$' + ' / POT'   # default: r'$\Phi (\nu)$' + ' / 50MeV / cm' + r'$^2$' + ' / POT'

upper_plot_limit = 5   # upper energy limit displayed in plot image (default = 5)


gsimple_directory = '/pnfs/annie/persistent/flux/annie_gsimple/gsimple_may2006_baseline_root/'
file_names = [os.path.join(gsimple_directory, f) for f in os.listdir(gsimple_directory)]

# ----------------------------------------------------- #
# ***************************************************** #


# detector parameters
# ---- full detector ----
if Full_Volume:
    tank_radius = 1.524  # [m]
    tank_height = 1.98   # [m] (half-height)

# -------- FV -----------
else:
    tank_radius = FV_radius
    tank_height = FV_half_z
    
p_c = np.array([0, -0.1446, 1.681])  # center position for simulated geometry [m]
surface_area = (2 * tank_radius * 2 * tank_height) * 1e4  # cyclinder cross section = rectangle, convert to cm^2

print(f'\nVolume selected: {tank_radius}m radius, {tank_height}m half-height\n')

# ----------------------------------------------------- #

nu_e_energy = []
nu_e_bar_energy = []
nu_mu_energy = []
nu_mu_bar_energy = []

# determine if the neutrino has entered the detector volume
def inDetector(x,y,z):

    z_c = p_c[2]
    y_c = p_c[1]
    x_c = p_c[0]

    y_min = y_c - tank_height    # if FV not centered at (0,0,0), you can modify these lines for y
    y_max = y_c + tank_height

                                 # and this line for z / x
    if (y_min < y < y_max) and (np.sqrt( (z - z_c)**2 + (x - x_c)**2 ) < tank_radius):
        return True
    else:
        return False
    

# determine if the neutrino ever hits the detector
def HitDetector(x,y,z,px,py,pz):

    z_c = p_c[2]
    
    if pz <= 0.:                     # if neutrino is going the opposite way (upstream), it'll never hit the tank
        return False
    
    # create unit vectors for direction
    p = np.sqrt(px*px + py*py + pz*pz)
    ux = px/p
    uy = py/p
    uz = pz/p

    # step up to front of detector
    z_front = z_c - tank_radius      # position of the front of the tank (in z)
    first_step = (z_front - z)/uz
    x = x + ux*first_step
    y = y + uy*first_step
    z = z_front
    
    # iterate until middle of tank (cylinder = it won't enter if it didnt already pass through the front half) or return true if in tank
    N_steps = 50      # somewhat arbitrary but more steps doesn't seem to make a signficiant difference
    step = z_c/(N_steps*uz)
    while z <= z_c:
        z = z + uz*step
        x = x + ux*step
        y = y + uy*step
        if inDetector(x, y, z):      # if during propagation the neutrino passes through the detector volume
            return True
    return False

# ----------------------------------------------------- #

total_count = 0; passed = 0; count = 1; total_POT = 0
for file in range(len(file_names)):

    print(count, '/', len(file_names), '  (', round(100*count/len(file_names),2), '% )')
    count += 1

    if TEST_RUN == True:
        if count > 10:
            continue

    with uproot.open(file_names[file]) as root_file:
        
        flux = root_file["flux/entry"]
        E = flux["E"].array(library="np")            # neutrino energy [GeV]
        PDG = flux["pdg"].array(library="np")        # particle ID: 12/-12: nu e / nu e bar, 14/-14: nu mu / nu mu bar
        wgt = flux["wgt"].array(library="np")        # corresponding neutrino weight (all gsimple files have weights stripped = 1.0)
        vtxx = flux["vtxx"].array(library="np")      # neutrino ray origin [m]
        vtxy = flux["vtxy"].array(library="np")
        vtxz = flux["vtxz"].array(library="np")
        px = flux["px"].array(library="np")          # neutrino 3 momentum components [GeV/c]
        py = flux["py"].array(library="np")
        pz = flux["pz"].array(library="np")

        meta = root_file["meta/meta/"]
        protons = meta["protons"].array(library="np")
        total_POT += protons[0]

        for i in range(len(E)):
            if PDG[i] != 12 and PDG[i] != -12 and PDG[i] != 14 and PDG[i] != -14:
                print('NU TAU DETECTED!!!!!!!!!')    # we don't expect any of these
                continue
            else:
                total_count += 1
                if not HitDetector(vtxx[i],vtxy[i],vtxz[i],px[i],py[i],pz[i]):
                    continue
                else:
                    passed += 1
                    if PDG[i] == 12:      # nu e
                        nu_e_energy.append(E[i])
                    elif PDG[i] == -12:   # nu e bar
                        nu_e_bar_energy.append(E[i])
                    elif PDG[i] == 14:    # nu mu
                        nu_mu_energy.append(E[i])
                    elif PDG[i] == -14:   # nu mu bar
                        nu_mu_bar_energy.append(E[i])


# calculate flux mean and median by flavor (and totals)
avg_nu_mu = round(1000*np.average(nu_mu_energy),2)      # x1000 to get from GeV -> MeV
med_nu_mu = round(1000*np.median(nu_mu_energy),2)
avg_nu_e = round(1000*np.average(nu_e_energy),2)
med_nu_e = round(1000*np.median(nu_e_energy),2)
avg_nu_mu_bar = round(1000*np.average(nu_mu_bar_energy),2)
med_nu_mu_bar = round(1000*np.median(nu_mu_bar_energy),2)
avg_nu_e_bar = round(1000*np.average(nu_e_bar_energy),2)
med_nu_e_bar = round(1000*np.median(nu_e_bar_energy),2)

nu_avg_energy = round(1000*np.average(nu_mu_energy + nu_e_energy),2)
nu_med_energy = round(1000*np.median(nu_mu_energy + nu_e_energy),2)
nu_bar_avg_energy = round(1000*np.average(nu_mu_bar_energy + nu_e_bar_energy),2)
nu_bar_med_energy = round(1000*np.median(nu_mu_bar_energy + nu_e_bar_energy),2)

# ----------------------------------------------------- #

plt.figure(figsize=(5, 4))
binning = np.arange(0,upper_hist_limit+bin_size,bin_size)
bin_centers = (binning[:-1] + binning[1:]) / 2

# raw counts before applying units
counts_nu_mu, _ = np.histogram(nu_mu_energy, bins=binning)
counts_nu_mu_bar, _ = np.histogram(nu_mu_bar_energy, bins=binning)
counts_nu_e, _ = np.histogram(nu_e_energy, bins=binning)
counts_nu_e_bar, _ = np.histogram(nu_e_bar_energy, bins=binning)

# errors
errors_nu_mu = np.sqrt(counts_nu_mu)
errors_nu_mu_bar = np.sqrt(counts_nu_mu_bar)
errors_nu_e = np.sqrt(counts_nu_e)
errors_nu_e_bar = np.sqrt(counts_nu_e_bar)

# unit conversion (* CHANGE IF NEEDED *)
scaling_factor = 1 / total_POT / surface_area        # v / 50 MeV / POT / cm^2

# scale counts for appropriate units
counts_nu_mu = counts_nu_mu * scaling_factor
counts_nu_mu_bar = counts_nu_mu_bar * scaling_factor
counts_nu_e = counts_nu_e * scaling_factor
counts_nu_e_bar = counts_nu_e_bar * scaling_factor

# also need to scale the errors
errors_nu_mu = errors_nu_mu * scaling_factor
errors_nu_mu_bar = errors_nu_mu_bar * scaling_factor
errors_nu_e = errors_nu_e * scaling_factor
errors_nu_e_bar = errors_nu_e_bar * scaling_factor


plt.hist(binning[:-1], binning, weights = counts_nu_mu, histtype = 'step',
    label= r'$\nu_{\mu}$', color = 'red', linewidth = 1.5)
plt.hist(binning[:-1], binning, weights = counts_nu_mu_bar, histtype = 'step',
    label= r'$\bar{\nu}_{\mu}$', color = 'blue', linewidth = 1.5)
plt.hist(binning[:-1], binning, weights = counts_nu_e, histtype = 'step', linestyle='--',
    label= r'$\nu_{e}$', color = 'red', linewidth = 1)
plt.hist(binning[:-1], binning, weights = counts_nu_e_bar, histtype = 'step', linestyle='--',
    label= r'$\bar{\nu}_{e}$', color = 'blue', linewidth = 1)

# error bands
plt.fill_between(bin_centers, counts_nu_mu - errors_nu_mu, counts_nu_mu + errors_nu_mu,
                 color='red', alpha=0.4, step='mid')
plt.fill_between(bin_centers, counts_nu_mu_bar - errors_nu_mu_bar, counts_nu_mu_bar + errors_nu_mu_bar,
                 color='blue', alpha=0.4, step='mid')
plt.fill_between(bin_centers, counts_nu_e - errors_nu_e, counts_nu_e + errors_nu_e,
                 color='red', alpha=0.2, step='mid')
plt.fill_between(bin_centers, counts_nu_e_bar - errors_nu_e_bar, counts_nu_e_bar + errors_nu_e_bar,
                 color='blue', alpha=0.2, step='mid')

plt.yscale('log')
plt.xlabel('Energy (GeV)', loc = 'right')
plt.ylabel(hist_y_label)  # see configurations
plt.legend(fontsize = 11, frameon = False, loc = 'upper right')     # legend text may be too big depending on xlim
plt.xlim(0, upper_plot_limit)   # adjust if you want in configuration
ax = plt.gca()
ax.tick_params(axis='x', which = 'both', direction= 'in', top = True)
ax.tick_params(axis='y', which = 'both', direction= 'in', right = True)
ax.xaxis.set_major_locator(MultipleLocator(1))   # 0.5 for 3 GeV xlim
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.text(0.6, 0.9, 'ANNIE Preliminary', transform=ax.transAxes, fontsize=9, color='grey', ha='center')
plt.tight_layout()
plt.savefig(plot_name,dpi=300,bbox_inches='tight',pad_inches=.3,facecolor = 'w')
plt.close()

# ----------------------------------------------------- #

# output flux to ROOT
print(root_file)
with uproot.recreate(rootfile_name) as file:
    file["numu_cv"] = (counts_nu_mu, binning)
    file["numubar_cv"] = (counts_nu_mu_bar, binning)
    file["nue_cv"] = (counts_nu_e, binning)
    file["nuebar_cv"] = (counts_nu_e_bar, binning)

# ----------------------------------------------------- #

print('\n')
print('****************************************')
print(total_count, 'total neutrinos produced')
print(passed, 'neutrinos that hit the FV (', round(100*passed/total_count,3), '% )')
print(round(total_POT/(1e6),4), 'e6 POT')
print('\n')
print('Flux breakdown')
print('----------------------------------------------------------------------')
print('muon neutrino:         ', str(round(100*len(nu_mu_energy)/passed,2)), '%, <E> =', avg_nu_mu, 'MeV, E_med =', med_nu_mu, 'MeV')
print('muon antineutrino:     ', str(round(100*len(nu_mu_bar_energy)/passed,2)), '%,  <E> =', avg_nu_mu_bar, 'MeV, E_med =', med_nu_mu_bar, 'MeV')
print('electron neutrino:     ', str(round(100*len(nu_e_energy)/passed,2)), '%,  <E> =', avg_nu_e, 'MeV, E_med =', med_nu_e, 'MeV')
print('electron antineutrino: ', str(round(100*len(nu_e_bar_energy)/passed,2)), '%,  <E> =', avg_nu_e_bar, 'MeV, E_med =', med_nu_e_bar, 'MeV')
print('**********************************************************************')
print('Neutrino Avg Energy      =', nu_avg_energy, 'MeV,  Median Energy = ', nu_med_energy, 'MeV')
print('Anti-Neutrino Avg Energy =', nu_bar_avg_energy, 'MeV,  Median Energy = ', nu_bar_med_energy, 'MeV')

print('\n\ndone\n')
