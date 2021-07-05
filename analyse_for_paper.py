from astropy import coordinates as coord
from astropy import units as u
import argparse
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit # mod ly
import sys
from matplotlib.ticker import FormatStrFormatter
import numpy.ma as ma
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid.inset_locator import inset_axes

parser = argparse.ArgumentParser(description='Plot routine for fitting over a glitch.')
parser.add_argument('-p', '--parfile', help='Path to ephemeris', required=True)
parser.add_argument('-s', '--stride', help='Stride data text file', required=True)
args = parser.parse_args()
par = args.parfile
strdat = args.stride

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

#plt.rcParams["figure.figsize"] = [7.5,10.6]

pglep=np.zeros(100)
pglf0=np.zeros(100)
pglf1=np.zeros(100)
pglf2=np.zeros(100)
pglf0d=np.zeros(100)
pgltd=np.ones(100)
pglf0d2=np.zeros(100)
pgltd2=np.ones(100)
max_glitch=0

#Open par file and extract parameters
with open(par) as f:
    for line in f:
        line = line.strip()
        e=line.split()
        if e[0] == "PSRJ":
            psrn=e[1]
        if e[0].startswith("GLEP_"):
            i=int(e[0][5:])
            pglep[i-1] = float(e[1])
            max_glitch = max(i,max_glitch)
        if e[0].startswith("GLF0_"):
            i=int(e[0][5:])
            pglf0[i-1] = float(e[1])
        if e[0].startswith("GLF1_"):
            i=int(e[0][5:])
            pglf1[i-1] = float(e[1])
        if e[0].startswith("GLF2_"):
            i=int(e[0][5:])
            pglf2[i-1] = float(e[1])
        if e[0].startswith("GLTD_"):
            i=int(e[0][5:])
            pgltd[i-1] = float(e[1])
        if e[0].startswith("GLF0D_"):
            i=int(e[0][6:])
            pglf0d[i-1] = float(e[1])
# mod ly
        if e[0].startswith("GLTD2_"):
            i=int(e[0][6:])
            pgltd2[i-1] = float(e[1])
        if e[0].startswith("GLF0D2_"):
            i=int(e[0][7:])
            pglf0d2[i-1] = float(e[1])
# mod ly
        if e[0] == "F0":
            F0=float(e[1])
        if e[0] == "PB":
            PB=float(e[1])
        if e[0] == "F1":
            F1=float(e[1])
        if e[0] == "F2":
            F2=float(e[1])
        if e[0] == "START":
            start=float(e[1])
        if e[0] == "FINISH":
            finish=float(e[1])
        if e[0] == "PEPOCH":
            pepoch=float(e[1])

print("")
print("Parameters in par file")
print("F0:", F0)
print("F1:", F1)
print("F2:", F2)
print("")

# for i in range(max_glitch):

print("Glitch epoch:", pglep[0])
print("GLF0:", pglf0[0])
print("GLF1:", pglf1[0])
print("GLF2:", pglf2[0])
print("")
print("GLF0D_1:", pglf0d[0], " - GLTD_1", pgltd[0], "GLF0D2_1:", pglf0d2[0], " - GLTD2_1", pgltd2[0])
print("Initial jump:", pglf0d[0]+pglf0d2[0]+pglf0[0])

#print("GLF0D_2:", pglf0d[1], " - T1", pgltd[1])
#print("GLF0D_3:", pglf0d[2], " - T1", pgltd[2])

#print("Initial jump:", np.sum(pglf0d)+pglf0[0])

dat="deltanu_{}.asc".format(psrn)
#par="test_final.par"

t, nu = np.loadtxt(dat,unpack=True)
t, no_nudot, nudot, nudot_model = np.loadtxt("nudot_{}.asc".format(psrn), unpack=True)
#this_mjd, this_dnu, this_dnu_err = np.loadtxt("mjd_dnu_dnuerr_nof2.txt", unpack=True)
panel1_mjd, panel1_nu, panel1_err = np.loadtxt("panel1_{}.txt".format(psrn), unpack=True)
panel2_mjd, panel2_nu, panel2_err = np.loadtxt("panel2_{}.txt".format(psrn), unpack=True)
panel3_mjd, panel3_nu, panel3_err = np.loadtxt("panel3_{}.txt".format(psrn), unpack=True)
str_mjd, str_nudot, str_err , str_nu2dot, str_2err = np.loadtxt(strdat, unpack=True, usecols=[0,3,4,5,6])
#print(np.max(str_err)) ; sys.exit(9)


def glexp(xx,td,f0d):
    '''
    xx = time since glitch epoch
    td = decay timescale
    f0d = decay amplitude
    '''

    ee = np.zeros_like(xx)
    tau1=td*86400.0
    ee[xx>0] = f0d * np.exp(-xx[xx>0]/tau1)
    return ee


#Time since period epoch
x = (t-pepoch)*86400.0

#First derivative term of taylor series
tf1 = F1 * x

#second derivative term of taylor series
#f = f0 + f1 t + 0.5 f2 t^2
tf2 = 0.5 * x * x * F2

glf0 = np.zeros_like(t)
glf1 = np.zeros_like(t)
glf2 = np.zeros_like(t)
exp1 = np.zeros_like(t)
exp2 = np.zeros_like(t)
mdglf1 = np.zeros_like(t)
mdglf2 = np.zeros_like(t)
mdcor = np.zeros_like(t)

#str2nu = np.zeros_like(str_mjd)
strglf1 = np.zeros_like(str_mjd)
strglf2 = np.zeros_like(str_mjd)
strexp1 = np.zeros_like(str_mjd)
strexp2 = np.zeros_like(str_mjd)
strcor = np.zeros_like(str_mjd)

gleps = []
for gi in range(len(pglep)):
    if float(pglep[gi]) != 0:
        glep = pglep[gi]
        gleps.append(pglep[gi])
        print("The {} glitch at {}".format(gi+1, glep))
        #Time since glitch epoch
        xx = (t-glep)*86400.0 # add

        #Permanent change term (constant)
        glf0[xx>0] += pglf0[gi]

        #GLF1 term
        glf1[xx>0] += xx[xx>0] * pglf1[gi]

        #GLF2 term
        glf2[xx>0] += 0.5* (xx[xx>0]**2) * pglf2[gi]

        #Change in spin-down rate
        mdglf1[xx>0] += pglf1[gi]  

        #Change in spin-down rate derivative
        mdglf2[xx>0] += xx[xx>0] * pglf2[gi]  

        #transient terms
        exp1 += glexp(xx,pgltd[gi],pglf0d[gi])
        if pglf0d2[gi] != 0:
            exp2 += glexp(xx,pgltd2[gi],pglf0d2[gi])

        #subtract exp
        mdcor += glexp(xx,pgltd[gi],pglf0d[gi])/(pgltd[gi]*86400) + glexp(xx,pgltd2[gi],pglf0d2[gi])/(pgltd2[gi]*86400)

        #stide data terms
        strx = (str_mjd-glep)*86400
        #str2nu[strx>0] = strx[strx>0] * str_nu2dot

        #stride GLF1 term
        strglf1[strx>0] += pglf1[gi]

        #stride GLF2 term
        strglf2[strx>0] += strx[strx>0] * pglf2[gi]

        #transient terms
        strexp1 += glexp(strx,pgltd[gi],pglf0d[gi])
        if pglf0d2[gi] != 0:
            strexp2 += glexp(strx,pgltd2[gi],pglf0d2[gi])

        #subtract exp
        strcor += glexp(strx,pgltd[gi],pglf0d[gi])/(pgltd[gi]*86400) + glexp(strx,pgltd2[gi],pglf0d2[gi])/(pgltd2[gi]*86400)

# mod ly
glep = pglep[0]

mask_len = []
for i in range(len(t)):
    for gi in range(len(gleps)):
        if t[i] <= pglep[gi] < t[i+1]:
            mask_len.append(i) # or i+1?
# mask length is the length of pre-glitch data

#frequency evolution with second derivative and spin-down change subtracted
numod = nu - tf2 # why only f2? deltanu=nu-F0-f1?
mc = ma.array(numod)
mc[mask_len] = ma.masked  # mc: mask data at glep

# mod ly
#sf2 = x * F2
#nudot_mod = nudot_model - sf2
# mod ly
md = ma.array(nudot_model) 
md[mask_len] = ma.masked # md: mask data at glep


if all(f0d==0 for f0d in pglf0d): # remove the 3rd panel
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(7.5, 10.6)) #10.6
    #fig = plt.figure(figsize=(7.5, 10.6))
    #gs = gridspec.GridSpec(4, 1)
    #gs.update(wspace=0.025, hspace=0.0001)

    plt.subplot(411)
    plt.plot(t-glep,1e6*mc, 'k-', zorder=2)
    plt.errorbar(panel1_mjd - glep, panel1_nu, yerr=panel1_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1)
    #plt.errorbar(this_mjd - glep, this_dnu, yerr=this_dnu_err, marker='.', color='k', ecolor='k', linestyle='None', alpha=0.2, markersize=4)
    plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15)
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    #plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])

    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])

    plt.subplot(412)
    plt.plot(t-glep,1e6*(mc-glf0-glf1-glf2), 'k-', zorder=2)
    plt.errorbar(panel2_mjd - glep, panel2_nu, yerr=panel2_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([-30, 30])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    #plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15, labelpad=15)

    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])

    plt.subplot(413)
    plt.plot(t-glep, md/1e5, 'k-', zorder=2)
    plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.ylim([-3.725, -3.62])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    plt.ylabel(r'$\dot{\nu}$ ($10^{-10}$ Hz s$^{-1}$)', fontsize=15)
    plt.xlabel("Days since glitch epoch", fontsize=15)
    plt.subplots_adjust(wspace=0, hspace=0.002)

    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])

    #INSET
    #inset_axes2 = inset_axes(ax[3], 
    #                    width="40%", # width = 30% of parent_bbox
    #                    height=0.6, # height : 1 inch
    #                    loc=1,
    #                    borderpad=1)
    #plt.plot(t-glep, md/1e5, 'k-', zorder=2)
    #plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([0,50])
    #frame = plt.gca()
    ##frame.axes.xaxis.set_ticklabels([])
    #frame.axes.yaxis.set_ticklabels([])

    #plt.plot(t-glep, md/1e5, 'k-', zorder=2)
    #plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([-0.001, 5]) 
    #plt.xticks([])
    #plt.yticks([])
    #END INSET

    plt.subplot(414)
    plt.plot(t-glep, md/1e5-(mdglf1+mdglf2-mdcor)/1e-10, 'k-', zorder=2)
    plt.errorbar(str_mjd - glep, (str_nudot-strglf1-strglf2+strcor)/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.ylim([-3.725, -3.62])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    plt.ylabel(r'$\dot{\nu}_{\mathrm{ng}}$ ($10^{-10}$ Hz s$^{-1}$)', fontsize=15)
    plt.xlabel("Days since glitch epoch", fontsize=15)
    plt.subplots_adjust(wspace=0, hspace=0.002)

    #frame = plt.gca()
    #frame.axes.xaxis.set_ticklabels([])

    plt.tight_layout()
    plt.savefig("nu_nudot_gp_{}.pdf".format(psrn), format='pdf', dpi=400)
    plt.show()


else:
    fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(7.5, 12.8)) #10.6
    #fig = plt.figure(figsize=(7.5, 10.6))
    #gs = gridspec.GridSpec(4, 1)
    #gs.update(wspace=0.025, hspace=0.0001)

    plt.subplot(511)
    plt.plot(t-glep,1e6*mc, 'k-', zorder=2)
    plt.errorbar(panel1_mjd - glep, panel1_nu, yerr=panel1_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1)
    #plt.errorbar(this_mjd - glep, this_dnu, yerr=this_dnu_err, marker='.', color='k', ecolor='k', linestyle='None', alpha=0.2, markersize=4)
    plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15)
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    #plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])

    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])

    plt.subplot(512)
    plt.plot(t-glep,1e6*(mc-glf0-glf1), 'k-', zorder=2)
    plt.errorbar(panel2_mjd - glep, panel2_nu, yerr=panel2_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([-30, 30])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    #plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15, labelpad=15)

    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])

    #If exists recovery
    plt.subplot(513)
    plt.ylabel(r'$\delta \nu$ ($\mu$Hz)', fontsize=15)
    #plt.ylim([-0.3, 0.08])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    plt.plot(t-glep,1e6*(mc-glf0-glf1-exp1-exp2), 'k-', zorder=2)
    plt.errorbar(panel3_mjd - glep, panel3_nu, yerr=panel3_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    #plt.axvline(0, color='k', linestyle='dashed', alpha=0.5)
    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])
    #frame.axes.yaxis.set_ticklabels([])

    #INSET
    #inset_axes1 = inset_axes(ax[2], 
    #                    width="40%", # width = 30% of parent_bbox
    #                    height=0.8, # height : 1 inch
    #                    loc=4, 
    #                    borderpad=1
    #                    )
    #plt.plot(t-glep,1e6*(mc-glf0-glf1-exp1), 'k-', zorder=2)
    #plt.errorbar(panel3_mjd - glep, panel3_nu, yerr=panel3_err, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([0, 50]) 
    #frame = plt.gca()
    #   #frame.axes.xaxis.set_ticklabels([])
    #frame.axes.yaxis.set_ticklabels([])
    #END INSET

    plt.subplot(514)
    plt.plot(t-glep, md/1e5, 'k-', zorder=2)
    plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.ylim([-3.725, -3.62])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    plt.ylabel(r'$\dot{\nu}$ ($10^{-10}$ Hz s$^{-1}$)', fontsize=15)
    plt.xlabel("Days since glitch epoch", fontsize=15)
    plt.subplots_adjust(wspace=0, hspace=0.002)

    frame = plt.gca()
    frame.axes.xaxis.set_ticklabels([])

    #INSET
    #inset_axes2 = inset_axes(ax[3], 
    #                    width="40%", # width = 30% of parent_bbox
    #                    height=0.6, # height : 1 inch
    #                    loc=1,
    #                    borderpad=1)
    #plt.plot(t-glep, md/1e5, 'k-', zorder=2)
    #plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([0,50])
    #frame = plt.gca()
    ##frame.axes.xaxis.set_ticklabels([])
    #frame.axes.yaxis.set_ticklabels([])

    #plt.plot(t-glep, md/1e5, 'k-', zorder=2)
    #plt.errorbar(str_mjd - glep, str_nudot/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.xlim([-0.001, 5]) 
    #plt.xticks([])
    #plt.yticks([])
    #END INSET

    plt.subplot(515)
    plt.plot(t-glep, md/1e5-(mdglf1+mdglf2-mdcor)/1e-10, 'k-', zorder=2)
    plt.errorbar(str_mjd - glep, (str_nudot-strglf1-strglf2+strcor)/1e-10, yerr=str_err/1e-10, marker='.', color='r', ecolor='r', linestyle='None', alpha=1, zorder=1, markersize=4)
    #plt.ylim([-3.725, -3.62])
    #plt.xlim([-30, 50])
    #plt.xlim([-10, 10])
    for gls in gleps:
        plt.axvline(gls-glep, color='k', linestyle='dotted', alpha=0.3, linewidth=2)
    plt.ylabel(r'$\dot{\nu}_{\mathrm{ng}}$ ($10^{-10}$ Hz s$^{-1}$)', fontsize=15)
    plt.xlabel("Days since glitch epoch", fontsize=15)
    plt.subplots_adjust(wspace=0, hspace=0.002)

    #frame = plt.gca()
    #frame.axes.xaxis.set_ticklabels([])

    plt.tight_layout()
    plt.savefig("nu_nudot_gp_{}.pdf".format(psrn), format='pdf', dpi=400)
    plt.show()
