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

# ----------------------------------------------------- #

plot_name = 'ANNIE BNB flux _ 10.png'

gsimple_directory = '/pnfs/annie/persistent/flux/annie_gsimple/gsimple_may2006_baseline_root/'
#gsimple_directory = "gsimple_flux/"
file_names = []
for file_name in os.listdir(gsimple_directory):
    file_names.append(gsimple_directory + file_name)

# ----------------------------------------------------- #

nu_e_energy = [[], []]; nu_e_bar_energy = [[], []]     # energy, weight
nu_mu_energy = [[], []]; nu_mu_bar_energy = [[], []]

tank_radius = 1.524  # [m]
tank_height = 1.98   # [m]
p_c = [0, -0.1446, 1.681]  # center position for simulated geometry [m]
surface_area = (2*tank_radius*2*tank_height)  # [m], 2*radius * height (tank_height is only half height so we must multiply by 2)

# ----------------------------------------------------- #

# determine if the neutrino has entered the detector volume
def inDetector(x,y,z):

    z_c = p_c[2]
    y_c = p_c[1]
    x_c = p_c[0]

    y_min = y_c - tank_height
    y_max = y_c + tank_height

    if (y_min < y < y_max) and (np.sqrt( (z - z_c)**2 + (x - x_c)**2 ) < tank_radius):
        return True
    else:
        return False
    

# determine if the neutrino ever hits the detector
def HitDetector(x,y,z,px,py,pz):

    z_c = p_c[2]
    
    if pz <= 0.:        # if neutrino is going the opposite way (upstream), it'll never hit the tank
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
    N_steps = 50
    step = z_c/(N_steps*uz)
    while z <= z_c:
        z = z + uz*step
        x = x + ux*step
        y = y + uy*step
        if inDetector(x, y, z):    # if during propagation the neutrino passes through the detector volume
            return True
    return False

# ----------------------------------------------------- #

total_count = 0; passed = 0; count = 1; total_POT = 0
for file in range(len(file_names)):

    print(count, '/', len(file_names), '  (', round(100*count/len(file_names),2), '% )')
    count += 1

    with uproot.open(file_names[file]) as root_file:
        
        flux = root_file["flux/entry"]
        E = flux["E"].array(library="np")            # neutrino energy [GeV]
        PDG = flux["pdg"].array(library="np")        # particle ID: 12/-12: nu e / nu e bar, 14/-14: nu mu / nu mu bar
        wgt = flux["wgt"].array(library="np")        # corresponding neutrino weight (default = 1.0)
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
            if wgt[i] != 1.0:
                print('Non 1.0 weight found!!!!')    # This assumption is critical for calculating the normalization. Has been tested for ~50 files and they're all 1.0
            else:
                total_count += 1
                if not HitDetector(vtxx[i],vtxy[i],vtxz[i],px[i],py[i],pz[i]):
                    continue
                else:
                    passed += 1
                    if PDG[i] == 12:      # nu e
                        nu_e_energy[0].append(E[i])
                    elif PDG[i] == -12:   # nu e bar
                        nu_e_bar_energy[0].append(E[i])
                    elif PDG[i] == 14:    # nu mu
                        nu_mu_energy[0].append(E[i])
                    elif PDG[i] == -14:   # nu mu bar
                        nu_mu_bar_energy[0].append(E[i])


# ----------------------------------------------------- #

plt.figure(figsize=(5, 4))
bin_size = .05   # GeV --> corresponds to the 50 MeV bins we want
binning = np.arange(0,5+bin_size+bin_size,bin_size)   # extend it slightly so we don't see the end of the line plot

# apply proper units
counts_nu_mu, _ = np.histogram(nu_mu_energy[0], bins=binning)
counts_nu_mu_bar, _ = np.histogram(nu_mu_bar_energy[0], bins=binning)
counts_nu_e, _ = np.histogram(nu_e_energy[0], bins=binning)
counts_nu_e_bar, _ = np.histogram(nu_e_bar_energy[0], bins=binning)

#'''
# for units of 10^12 POT per cm^2
surface_area = surface_area * 10000    # convert from m^2 to cm^2
scaling_factor = 1 * (1e12) / total_POT / surface_area      # v / 50MeV / 10^12 POT / cm^2
#'''

#scaling_factor = 1 * (1e6) / total_POT / surface_area       # v / 50MeV / 10^6 POT / m^2

counts_nu_mu = counts_nu_mu * scaling_factor
counts_nu_mu_bar = counts_nu_mu_bar * scaling_factor
counts_nu_e = counts_nu_e * scaling_factor
counts_nu_e_bar = counts_nu_e_bar * scaling_factor


plt.hist(binning[:-1], binning, weights = counts_nu_mu, histtype = 'step',
    label= r'$\nu_{\mu}$', color = 'red', linewidth = 1.5)
plt.hist(binning[:-1], binning, weights = counts_nu_mu_bar, histtype = 'step',
    label= r'$\bar{\nu}_{\mu}$', color = 'blue', linewidth = 1.5)
plt.hist(binning[:-1], binning, weights = counts_nu_e, histtype = 'step', linestyle='--',
    label= r'$\nu_{e}$', color = 'red', linewidth = 1)
plt.hist(binning[:-1], binning, weights = counts_nu_e_bar, histtype = 'step', linestyle='--',
    label= r'$\bar{\nu}_{e}$', color = 'blue', linewidth = 1)

plt.yscale('log')
plt.xlabel('Energy (GeV)', loc = 'right')
#plt.ylabel(r'$\Phi (\nu)$' + ' / 50MeV / m' + r'$^2$' + ' / 10' + r'$^6$' + ' POT')
plt.ylabel(r'$\Phi (\nu)$' + ' / 50MeV / cm' + r'$^2$' + ' / 10' + r'$^{12}$' + ' POT')
plt.legend(fontsize = 11, frameon = False, loc = 'upper right')
plt.xlim(0, 5)
ax = plt.gca()
ax.tick_params(axis='x', which = 'both', direction= 'in', top = True)
ax.tick_params(axis='y', which = 'both', direction= 'in', right = True)
ax.xaxis.set_major_locator(MultipleLocator(1))   # 0.5 for 3 GeV scale
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.text(0.6, 0.9, 'ANNIE Preliminary', transform=ax.transAxes, fontsize=9, color='grey', ha='center')
plt.tight_layout()
plt.savefig(plot_name,dpi=300,bbox_inches='tight',pad_inches=.3,facecolor = 'w')
plt.close()

# ----------------------------------------------------- #

print('\n')
print('****************************************')
print(total_count, 'total neutrinos produced')
print(passed, 'neutrinos that hit the detector (', round(100*passed/total_count,3), '% )')
print(round(total_POT/(1e6),4), 'e6 POT')
print('\n')
print('Flux breakdown')
print('-------------------------------------')
print('muon neutrino:         ', str(round(100*len(nu_mu_energy[0])/passed,2)), '%')
print('muon antineutrino:     ', str(round(100*len(nu_mu_bar_energy[0])/passed,2)), '%')
print('electron neutrino:     ', str(round(100*len(nu_e_energy[0])/passed,2)), '%')
print('electron antineutrino: ', str(round(100*len(nu_e_bar_energy[0])/passed,2)), '%')
print('****************************************')

print('\n\ndone\n')
