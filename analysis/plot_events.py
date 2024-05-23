import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

fig1 = plt.figure(figsize=(20,16))
#fig1 = plt.figure()
gs1 = fig1.add_gridspec(3, hspace=0)
ax11,ax21,ax31 = gs1.subplots(sharex=True, sharey=False)

fig2 = plt.figure(figsize=(20,16))
#fig2 = plt.figure()
gs2 = fig2.add_gridspec(3, hspace=0)
ax12,ax22,ax32 = gs2.subplots(sharex=True, sharey=False)

#LEFT_LIMIT_COH =  1e-11
#LEFT_LIMIT_INCOH = 1e-11

LEFT_LIMIT_COH =  1e-9
LEFT_LIMIT_INCOH = 1e-7

#parameters = ['1t, 1yr, DUNE', '67t, 3yr, DUNE', '1t, 1yr, DUNE tau-opt', '67t, 3yr, DUNE tau-opt', '1t, 1yr, FASER', '1.2t, 3yr, FASER']
parameters = [r'\textit{1t, 1yr}, \textbf{DUNE Std.}', r'\textit{67t, 3yr}, \textbf{DUNE Std.}', r'\textit{1t, 1yr}, \textbf{DUNE $\tau$-opt}', r'\textit{67t, 3yr}, \textbf{DUNE $\tau$-opt}']
N_PARAMETERS = len(parameters) - 1

#mu2_coh = np.array([1.177, 2.366e2, 5.922, 1.190e3, 1, 1])
#mu2_incoh = np.array([4.926e-1, 9.901e1, 2.296, 4.614e2, 2.287e-10, 8.235e-10])
#
#tau1_coh = np.array([1.425e-3, 2.865e-1, 6.504e-3, 1.307, 3.162e-8, 1.138e-7])
#tau1_incoh = np.array([1.734e-2, 3.485, 8.356e-2, 1.679e1, 1.119e-8, 4.030e-8])
#
#tau2_coh = np.array([2.107e-7, 4.235e-5, 8.890e-7, 1.787e-4, 7.892e-10, 2.841e-9])
#tau2_incoh = np.array([1.611e-4, 3.238e-2, 7.535e-4, 1.515e-1, 9.122e-10, 3.284e-9])

# Data #
mu2_coh = np.array([1.177, 2.366e2, 5.922, 1.190e3])
mu2_incoh = np.array([4.926e-1, 9.901e1, 2.296, 4.614e2])

tau1_coh = np.array([1.425e-3, 2.865e-1, 6.504e-3, 1.307])
tau1_incoh = np.array([1.789e-2, 3.596, 8.626e-2, 1.734e1])

tau2_coh = np.array([2.107e-7, 4.235e-5, 8.890e-7, 1.787e-4])
tau2_incoh = np.array([1.771e-4, 3.560e-2, 8.316e-4, 1.672e-1])

# Data uncertainty #
mu2_coh_err = np.array([7.021e-2, 1.411e1, 3.532e-1, 7.100e1])
mu2_incoh_err = np.array([1.515e-1, 3.045e1, 7.059e-1, 1.419e2])

tau1_coh_err = np.array([8.526e-5, 1.714e-2, 3.891e-4, 7.821e-2])
tau1_incoh_err = np.array([5.501e-3, 1.106, 2.652e-2, 5.331])

tau2_coh_err = np.array([1.276e-8, 2.565e-6, 5.378e-8, 1.081e-5])
tau2_incoh_err = np.array([5.445e-5, 1.094e-2, 2.556e-4, 5.138e-2])


cols = ['rebeccapurple','teal','goldenrod']

X_TEXT_COH = 2.3e3
X_TEXT_INCOH = 1.5e3

# Note: the y-axis labels (i.e. parameters) get assigned an integer y-value indexed at 0. So, e.g., 1t, 1yr, DUNE will have y=3, 147t, 3yr, DUNE will have y=2 and so on ending at y=0

for ax, data, data_err, col in zip([ax11, ax21, ax31],[mu2_coh, tau1_coh, tau2_coh], [mu2_coh_err, tau1_coh_err, tau2_coh_err], cols):
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_COH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

for ax, data, data_err, col in zip([ax12, ax22, ax32],[mu2_incoh, tau1_incoh, tau2_incoh], [mu2_incoh_err, tau1_incoh_err, tau2_incoh_err], cols):
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(parameters), np.flipud(data), left=LEFT_LIMIT_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
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

def val_to_str(value, error):
    if value >= 0.1:
        string = r"{\bf %.2f} $\pm$ {\bf %.2f}" % (value, error)
        return string
    else:
        string = r"{\bf %.1e} $\pm$ {\bf %.1e}" % (value, error)
        return string

for i in range(0, N_COH):
    mu2_coh_val = mu2_coh[i]
    tau1_coh_val = tau1_coh[i]
    tau2_coh_val = tau2_coh[i]

    mu2_coh_error = mu2_coh_err[i]
    tau1_coh_error = tau1_coh_err[i]
    tau2_coh_error = tau2_coh_err[i]

    mu2_coh_str = val_to_str(mu2_coh_val, mu2_coh_error)
    tau1_coh_str = val_to_str(tau1_coh_val, tau1_coh_error)
    tau2_coh_str = val_to_str(tau2_coh_val, tau2_coh_error)

    ax11.text(mu2_coh_val - 2*mu2_coh_val/5, (N_COH-1)-i-t_shift[i], mu2_coh_str,va='center',ha='right',fontsize=22)
    ax21.text(tau1_coh_val - 2*tau1_coh_val/5, (N_COH-1)-i-t_shift[i], tau1_coh_str,va='center',ha='right',fontsize=22)
    ax31.text(tau2_coh_val - 2*tau2_coh_val/5, (N_COH-1)-i-t_shift[i], tau2_coh_str,va='center',ha='right',fontsize=22)

for i in range(0, N_INCOH):
    mu2_incoh_val = mu2_incoh[i]
    tau1_incoh_val = tau1_incoh[i]
    tau2_incoh_val = tau2_incoh[i]

    mu2_incoh_error = mu2_incoh_err[i]
    tau1_incoh_error = tau1_incoh_err[i]
    tau2_incoh_error = tau2_incoh_err[i]

    mu2_incoh_str = val_to_str(mu2_incoh_val, mu2_incoh_error)
    tau1_incoh_str = val_to_str(tau1_incoh_val, tau1_incoh_error)
    tau2_incoh_str = val_to_str(tau2_incoh_val, tau2_incoh_error)

    ax12.text(mu2_incoh_val - 2*mu2_incoh_val/5, (N_INCOH-1)-i-t_shift[i], mu2_incoh_str,va='center',ha='right',fontsize=22)
    ax22.text(tau1_incoh_val - 2*tau1_incoh_val/5, (N_INCOH-1)-i-t_shift[i], tau1_incoh_str,va='center',ha='right',fontsize=22)
    ax32.text(tau2_incoh_val - 2*tau2_incoh_val/5, (N_INCOH-1)-i-t_shift[i], tau2_incoh_str,va='center',ha='right',fontsize=22)

# Fixing ticks
locmaj = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=1000)
locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=2000)
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

ax31.set_xlabel(r'{\bf Number of Events}')
ax11.set_title(r'{\bf Coherent scattering off $^{40}$Ar}', fontsize=60)

ax32.set_xlabel(r'{\bf Number of Events}')
ax12.set_title(r'{\bf Incoherent scattering off $^{40}$Ar}', fontsize=60)

fig1.savefig('../plots/events_coh.png',transparent=False,bbox_inches='tight')
fig2.savefig('../plots/events_incoh.png',transparent=False,bbox_inches='tight')
