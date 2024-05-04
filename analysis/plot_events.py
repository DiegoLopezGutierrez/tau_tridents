import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

fig1 = plt.figure(figsize=(20,16))
gs1 = fig1.add_gridspec(3, hspace=0)
ax11,ax21,ax31 = gs1.subplots(sharex=True, sharey=False)

fig2 = plt.figure(figsize=(20,16))
gs2 = fig2.add_gridspec(3, hspace=0)
ax12,ax22,ax32 = gs2.subplots(sharex=True, sharey=False)

LEFT_LIMIT_COH =  1e-8
LEFT_LIMIT_INCOH = 1e-5

#parameters = ['1t, 1yr, DUNE', '147t, 3yr, DUNE', '1t, 1yr, Alt.', '147t, 3yr, Alt.', '1t, 1yr, Alt. 120', '147t, 3yr, Alt. 120', '1t, 1yr, DUNE tau-opt', '147t, 3yr, DUNE tau-opt']
parameters = ['1t, 1yr, DUNE', '147t, 3yr, DUNE', '1t, 1yr, Alt. 120', '147t, 3yr, Alt. 120', '1t, 1yr, DUNE tau-opt', '147t, 3yr, DUNE tau-opt']
N_PARAMETERS = len(parameters) - 1

#mu2_coh = np.array([8.696e-1, 3.835e2, 8.269e-1, 3.647e2, 1.092, 4.816e2, 5.922, 2.612e3])
mu2_coh = np.array([1.177, 5.190e2, 1.170, 5.160e2, 5.922, 2.612e3])
#mu2_incoh_p = np.array([3.249e-1, 1.433e2, 3.100e-1, 1.367e2, 3.838e-1, 1.693e2, 1.922, 8.478e2])
mu2_incoh = np.array([4.926e-1, 2.172e2, 4.902e-1, 2.162e2, 2.296, 1.012e3])

#tau1_coh = np.array([3.929e-5, 1.733e-2, 5.213e-5, 2.299e-2, 2.510e-3, 1.107, 1.424e-2, 6.280])
tau1_coh = np.array([3.123e-3, 1.377, 2.689e-3, 1.186, 1.424e-2, 6.280])
#tau1_incoh_p = np.array([2.823e-3, 1.245, 2.866e-3, 1.264, 1.172e-2, 5.167, 6.266e-2, 2.763e1])
tau1_incoh = np.array([1.734e-2, 7.646, 1.675e-2, 7.385, 8.356e-2, 3.685e1])

#tau2_coh = np.array([3.244e-11, 1.430e-8, 7.024e-11, 3.097e-8, 1.148e-7, 5.061e-5, 8.890e-7, 3.920e-4])
tau2_coh = np.array([2.107e-7, 9.291e-5, 1.230e-7, 5.423e-5, 8.890e-7, 3.920e-4])
#tau2_incoh_p = np.array([7.324e-7, 3.230e-4, 1.295e-6, 5.709e-4, 9.068e-5, 3.999e-2, 5.145e-4, 2.269e-1])
tau2_incoh = np.array([1.611e-4, 7.104e-2, 1.424e-4, 6.278e-2, 7.535e-4, 3.323e-1])

cols = ['rebeccapurple','teal','goldenrod']

X_TEXT_COH = 2.3e3
X_TEXT_INCOH = 1.5e3

# Note: the y-axis labels (i.e. parameters) get assigned an integer y-value indexed at 0. So, e.g., 1t, 1yr, DUNE will have y=3, 147t, 3yr, DUNE will have y=2 and so on ending at y=0

for ax, data, col in zip([ax11, ax21, ax31],[mu2_coh, tau1_coh, tau2_coh], cols):
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_COH, align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5)
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

for ax, data, col in zip([ax12, ax22, ax32],[mu2_incoh, tau1_incoh, tau2_incoh], cols):
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_INCOH, align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5)
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

ax31.set_xscale('log')
ax31.set_xlim([LEFT_LIMIT_COH,3e3])

ax32.set_xscale('log')
ax32.set_xlim([LEFT_LIMIT_INCOH,2e3])

# Adding text labels
t_shift = [0.01,0.01,0.01,0.0,0.0,0.0, 0.0, 0.0] # Manual vertical shift for now
N_COH = len(mu2_coh)
N_INCOH = len(mu2_incoh)

def val_to_str(value):
    if value >= 0.1:
        string = r"{\bf %.2f}" % (value)
        return string
    else:
        string = r"{\bf %.2e}" % (value)
        return string

for i in range(0, N_COH):
    mu2_coh_val = mu2_coh[i]
    tau1_coh_val = tau1_coh[i]
    tau2_coh_val = tau2_coh[i]

    mu2_coh_str = val_to_str(mu2_coh_val)
    tau1_coh_str = val_to_str(tau1_coh_val)
    tau2_coh_str = val_to_str(tau2_coh_val)

    ax11.text(mu2_coh_val - 2*mu2_coh_val/5, (N_COH-1)-i-t_shift[i], mu2_coh_str,va='center',ha='right',fontsize=22)
    ax21.text(tau1_coh_val - 2*tau1_coh_val/5, (N_COH-1)-i-t_shift[i], tau1_coh_str,va='center',ha='right',fontsize=22)
    ax31.text(tau2_coh_val - 2*tau2_coh_val/5, (N_COH-1)-i-t_shift[i], tau2_coh_str,va='center',ha='right',fontsize=22)

for i in range(0, N_INCOH):
    mu2_incoh_val = mu2_incoh[i]
    tau1_incoh_val = tau1_incoh[i]
    tau2_incoh_val = tau2_incoh[i]

    mu2_incoh_str = val_to_str(mu2_incoh_val)
    tau1_incoh_str = val_to_str(tau1_incoh_val)
    tau2_incoh_str = val_to_str(tau2_incoh_val)

    ax12.text(mu2_incoh_val - 2*mu2_incoh_val/5, (N_INCOH-1)-i-t_shift[i], mu2_incoh_str,va='center',ha='right',fontsize=22)
    ax22.text(tau1_incoh_val - 2*tau1_incoh_val/5, (N_INCOH-1)-i-t_shift[i], tau1_incoh_str,va='center',ha='right',fontsize=22)
    ax32.text(tau2_incoh_val - 2*tau2_incoh_val/5, (N_INCOH-1)-i-t_shift[i], tau2_incoh_str,va='center',ha='right',fontsize=22)

# Fixing ticks
locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=50)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=100)
ax31.xaxis.set_major_locator(locmaj)
ax31.xaxis.set_minor_locator(locmin)
ax32.xaxis.set_major_locator(locmaj)
ax32.xaxis.set_minor_locator(locmin)

# Grid lines
ax11.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
ax21.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
ax31.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)

ax12.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
ax22.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
ax32.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)

# Tick parameters
ax11.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
ax21.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
ax31.tick_params(which='both',right=False,tickdir='out',top=False)

ax12.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
ax22.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
ax32.tick_params(which='both',right=False,tickdir='out',top=False)

# Labels
y_label_loc = N_PARAMETERS - 0.05
ax11.text(X_TEXT_COH,y_label_loc,r'$\nu_\mu \mu^+ \mu^-$',fontsize=30,alpha=0.7,color=cols[0],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')
ax21.text(X_TEXT_COH,y_label_loc,r'$\nu_\tau \tau^+ \mu^-$',fontsize=30,alpha=0.7,color=cols[1],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')
ax31.text(X_TEXT_COH,y_label_loc,r'$\nu_\mu \tau^+ \tau^-$',fontsize=30,alpha=0.7,color=cols[2],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

ax12.text(X_TEXT_INCOH,y_label_loc,r'$\nu_\mu \mu^+ \mu^-$',fontsize=30,alpha=0.7,color=cols[0],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')
ax22.text(X_TEXT_INCOH,y_label_loc,r'$\nu_\tau \tau^+ \mu^-$',fontsize=30,alpha=0.7,color=cols[1],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')
ax32.text(X_TEXT_INCOH,y_label_loc,r'$\nu_\mu \tau^+ \tau^-$',fontsize=30,alpha=0.7,color=cols[2],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

ax31.set_xlabel(r'Number of Events')
ax11.set_title(r'Coherent scattering off $^{40}$Ar')

ax32.set_xlabel(r'Number of Events')
ax12.set_title(r'Incoherent scattering off $^{40}$Ar')

fig1.savefig('../plots/events_coh.png',transparent=False,bbox_inches='tight')
fig2.savefig('../plots/events_incoh.png',transparent=False,bbox_inches='tight')
