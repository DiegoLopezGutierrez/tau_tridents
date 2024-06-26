import numpy as np
import matplotlib
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

STYLE_DIR = '../plots/styles/'
plt.style.use(STYLE_DIR+'sty.mplstyle')

## Neutrino Mode -- Coherent -- DUNE ##
fig1 = plt.figure(figsize=(30,20))
gs1 = fig1.add_gridspec(8, hspace=0)
ax11,ax21,ax31,ax41,ax51,ax61,ax71,ax81 = gs1.subplots(sharex=True, sharey=False)

## Neutrino Mode -- Incoherent -- DUNE ##
fig2 = plt.figure(figsize=(30,20))
gs2 = fig2.add_gridspec(8, hspace=0)
ax12,ax22,ax32,ax42,ax52,ax62,ax72,ax82 = gs2.subplots(sharex=True, sharey=False)

## Antineutrino Mode -- Coherent -- DUNE ##
fig3 = plt.figure(figsize=(30,20))
gs3 = fig3.add_gridspec(8, hspace=0)
ax13,ax23,ax33,ax43,ax53,ax63,ax73,ax83 = gs3.subplots(sharex=True, sharey=False)

## Antineutrino Mode -- Incoherent -- DUNE ##
fig4 = plt.figure(figsize=(30,20))
gs4 = fig4.add_gridspec(8, hspace=0)
ax14,ax24,ax34,ax44,ax54,ax64,ax74,ax84 = gs4.subplots(sharex=True, sharey=False)

## Neutrino Mode -- Coherent -- MINOS and MINOS+##
fig5 = plt.figure(figsize=(20,10))
gs5 = fig5.add_gridspec(4, hspace=0)
ax15,ax25,ax35,ax45 = gs5.subplots(sharex=True, sharey=False)

## Neutrino Mode -- Incoherent -- MINOS and MINOS+##
fig6 = plt.figure(figsize=(20,10))
gs6 = fig6.add_gridspec(4, hspace=0)
ax16,ax26,ax36,ax46 = gs6.subplots(sharex=True, sharey=False)

## Antineutrino Mode -- Coherent -- MINOS and MINOS+##
fig7 = plt.figure(figsize=(22,8))
gs7 = fig7.add_gridspec(4, hspace=0)
ax17,ax27,ax37,ax47 = gs7.subplots(sharex=True, sharey=False)

## Antineutrino Mode -- Incoherent -- MINOS and MINOS+##
fig8 = plt.figure(figsize=(22,8))
gs8 = fig8.add_gridspec(4, hspace=0)
ax18,ax28,ax38,ax48 = gs8.subplots(sharex=True, sharey=False)

## Coherent -- T2K INGRID Phase 1 and Phase 2 ##
fig9 = plt.figure(figsize=(25,8))
gs9 = fig9.add_gridspec(4, hspace=0)
ax19,ax29,ax39,ax49 = gs9.subplots(sharex=True, sharey=False)

## Incoherent -- T2K INGRID Phase 1 and Phase 2 ##
fig10 = plt.figure(figsize=(25,8))
gs10 = fig10.add_gridspec(4, hspace=0)
ax110,ax210,ax310,ax410 = gs10.subplots(sharex=True, sharey=False)

## Coherent -- FASER ##
fig11 = plt.figure(figsize=(20,6))
gs11 = fig11.add_gridspec(3, hspace=0)
ax111,ax211,ax311 = gs11.subplots(sharex=True, sharey=False)

## Incoherent -- FASER ##
fig12 = plt.figure(figsize=(20,6))
gs12 = fig12.add_gridspec(3, hspace=0)
ax112,ax212,ax312 = gs12.subplots(sharex=True, sharey=False)

## DUNE Neutrino Mode -- Events ##
# vmu #
DUNE_neutrino_vmu_2mu_coh = np.array([2.366e2, 1.190e3])
DUNE_neutrino_vmu_2mu_incoh = np.array([1.007e2, 4.710e2])

DUNE_neutrino_vmu_1tau_coh = np.array([2.865e-1, 1.307])
DUNE_neutrino_vmu_1tau_incoh = np.array([3.819, 1.843e1])

DUNE_neutrino_vmu_2tau_coh = np.array([4.235e-5, 1.787e-4])
DUNE_neutrino_vmu_2tau_incoh = np.array([4.193e-2, 1.979e-1])

# vmubar #
DUNE_neutrino_vmubar_2mu_coh = np.array([4.718e1, 5.047e1])
DUNE_neutrino_vmubar_2mu_incoh = np.array([1.834e1, 1.913e1])

DUNE_neutrino_vmubar_1tau_coh = np.array([4.267e-2, 5.655e-2])
DUNE_neutrino_vmubar_1tau_incoh = np.array([9.262e-1, 1.090])

DUNE_neutrino_vmubar_2tau_coh = np.array([2.815e-6, 4.679e-6])
DUNE_neutrino_vmubar_2tau_incoh = np.array([7.050e-3, 9.102e-03])

# ve #
DUNE_neutrino_ve_1tau_coh = np.array([1.188e-02, 3.245e-02])
DUNE_neutrino_ve_1tau_incoh = np.array([1.159e-01, 4.085e-01])

# vebar #
DUNE_neutrino_vebar_1tau_coh = np.array([3.376e-03, 7.815e-03])
DUNE_neutrino_vebar_1tau_incoh = np.array([5.295e-02, 8.499e-02])

## DUNE Antineutrino Mode -- Events ##
# vmu #
DUNE_antineutrino_vmu_2mu_coh = np.array([1.035e+02, 1.070e+02])
DUNE_antineutrino_vmu_2mu_incoh = np.array([3.822e+01, 3.875e+01])

DUNE_antineutrino_vmu_1tau_coh = np.array([2.999e-01, 4.227e-01])
DUNE_antineutrino_vmu_1tau_incoh = np.array([3.205, 3.749])

DUNE_antineutrino_vmu_2tau_coh = np.array([6.591e-05, 9.933e-05])
DUNE_antineutrino_vmu_2tau_incoh = np.array([4.236e-02, 5.767e-02])

# vmubar #
DUNE_antineutrino_vmubar_2mu_coh = np.array([1.453e+02, 7.535e+02])
DUNE_antineutrino_vmubar_2mu_incoh = np.array([6.619e+01, 3.080e+02])

DUNE_antineutrino_vmubar_1tau_coh = np.array([5.725e-02, 2.707e-01])
DUNE_antineutrino_vmubar_1tau_incoh = np.array([1.299, 6.541])

DUNE_antineutrino_vmubar_2tau_coh = np.array([3.934e-06, 2.240e-05])
DUNE_antineutrino_vmubar_2tau_incoh = np.array([9.316e-03, 4.474e-02])

# ve #
DUNE_antineutrino_ve_1tau_coh = np.array([8.844e-03, 1.712e-02])
DUNE_antineutrino_ve_1tau_incoh = np.array([1.044e-01, 1.447e-01])

# vebar #
DUNE_antineutrino_vebar_1tau_coh = np.array([8.378e-03, 1.059e-02])
DUNE_antineutrino_vebar_1tau_incoh = np.array([1.004e-01, 1.749e-01])

## DUNE Neutrino Mode -- Event Uncertainty ##
# vmu #
DUNE_neutrino_vmu_2mu_coh_err = np.array([1.411e+01, 7.100e+01])
DUNE_neutrino_vmu_2mu_incoh_err = np.array([3.097e+01, 1.449e+02])

DUNE_neutrino_vmu_1tau_coh_err = np.array([1.714e-02, 7.821e-02])
DUNE_neutrino_vmu_1tau_incoh_err = np.array([1.174, 5.668])

DUNE_neutrino_vmu_2tau_coh_err = np.array([2.565e-06, 1.081e-05])
DUNE_neutrino_vmu_2tau_incoh_err = np.array([1.289e-02, 6.085e-02])

# vmubar #
DUNE_neutrino_vmubar_2mu_coh_err = np.array([2.815, 3.011])
DUNE_neutrino_vmubar_2mu_incoh_err = np.array([5.640, 5.884])

DUNE_neutrino_vmubar_1tau_coh_err = np.array([2.555e-03, 3.385e-03])
DUNE_neutrino_vmubar_1tau_incoh_err = np.array([2.848e-01, 3.351e-01])

DUNE_neutrino_vmubar_2tau_coh_err = np.array([1.709e-07, 2.836e-07])
DUNE_neutrino_vmubar_2tau_incoh_err = np.array([2.167e-03, 2.798e-03])

# ve #
DUNE_neutrino_ve_1tau_coh_err = np.array([7.138e-04, 1.951e-03])
DUNE_neutrino_ve_1tau_incoh_err = np.array([3.568e-02, 1.257e-01])

# vebar #
DUNE_neutrino_vebar_1tau_coh_err = np.array([2.031e-04, 4.697e-04])
DUNE_neutrino_vebar_1tau_incoh_err = np.array([1.630e-02, 2.616e-02])

## DUNE Antineutrino Mode -- Event Uncertainty ##
# vmu #
DUNE_antineutrino_vmu_2mu_coh_err = np.array([6.174, 6.384])
DUNE_antineutrino_vmu_2mu_incoh_err = np.array([1.175e+01, 1.192e+01])

DUNE_antineutrino_vmu_1tau_coh_err = np.array([1.794e-02, 2.527e-02])
DUNE_antineutrino_vmu_1tau_incoh_err = np.array([9.856e-01, 1.153])

DUNE_antineutrino_vmu_2tau_coh_err = np.array([3.980e-06, 6.000e-06])
DUNE_antineutrino_vmu_2tau_incoh_err = np.array([1.302e-02, 1.773e-02])

# vmubar #
DUNE_antineutrino_vmubar_2mu_coh_err = np.array([8.670, 4.494e+01])
DUNE_antineutrino_vmubar_2mu_incoh_err = np.array([2.036e+01, 9.471e+01])

DUNE_antineutrino_vmubar_1tau_coh_err = np.array([3.428e-03, 1.621e-02])
DUNE_antineutrino_vmubar_1tau_incoh_err = np.array([3.993e-01, 2.011])

DUNE_antineutrino_vmubar_2tau_coh_err = np.array([2.388e-07, 1.355e-06])
DUNE_antineutrino_vmubar_2tau_incoh_err = np.array([2.864e-03, 1.375e-02])

# ve #
DUNE_antineutrino_ve_1tau_coh_err = np.array([5.315e-04, 1.028e-03])
DUNE_antineutrino_ve_1tau_incoh_err = np.array([3.212e-02, 4.454e-02])

# vebar #
DUNE_antineutrino_vebar_1tau_coh_err = np.array([5.037e-04, 6.370e-04])
DUNE_antineutrino_vebar_1tau_incoh_err = np.array([3.089e-02, 5.383e-02])

## MINOS and MINOS+ Neutrino Mode -- Events ##
# vmu #
MINOS_neutrino_vmu_2mu_coh = np.array([1.979e+01, 9.448e+01])
MINOS_neutrino_vmu_2mu_incoh = np.array([6.939, 3.218e+01])

MINOS_neutrino_vmu_1tau_coh = np.array([1.106e-04, 1.814e-04])
MINOS_neutrino_vmu_1tau_incoh = np.array([4.691e-02, 1.156e-01])

# vmubar #
MINOS_neutrino_vmubar_2mu_coh = np.array([5.225, 4.553])
MINOS_neutrino_vmubar_2mu_incoh = np.array([1.690, 1.468])

MINOS_neutrino_vmubar_1tau_coh = np.array([4.246e-05, 3.931e-05])
MINOS_neutrino_vmubar_1tau_incoh = np.array([1.844e-02, 1.673e-02])

## MINOS and MINOS+ Antineutrino Mode -- Events ##
# vmu #
MINOS_antineutrino_vmu_2mu_coh = np.array([2.973])
MINOS_antineutrino_vmu_2mu_incoh = np.array([9.515e-01])

MINOS_antineutrino_vmu_1tau_coh = np.array([2.618e-05])
MINOS_antineutrino_vmu_1tau_incoh = np.array([1.128e-02])

# vmubar #
MINOS_antineutrino_vmubar_2mu_coh = np.array([3.451])
MINOS_antineutrino_vmubar_2mu_incoh = np.array([1.267])

MINOS_antineutrino_vmubar_1tau_coh = np.array([1.180e-05])
MINOS_antineutrino_vmubar_1tau_incoh = np.array([5.237e-03])

## MINOS and MINOS+ Neutrino Mode -- Event Uncertainty ##
# vmu #
MINOS_neutrino_vmu_2mu_coh_err = np.array([1.609, 7.680])
MINOS_neutrino_vmu_2mu_incoh_err = np.array([2.142, 9.935])

MINOS_neutrino_vmu_1tau_coh_err = np.array([6.969e-06, 1.142e-05])
MINOS_neutrino_vmu_1tau_incoh_err = np.array([1.445e-02, 3.559e-02])

# vmubar #
MINOS_neutrino_vmubar_2mu_coh_err = np.array([4.247e-01, 3.701e-01])
MINOS_neutrino_vmubar_2mu_incoh_err = np.array([5.219e-01, 4.534e-01])

MINOS_neutrino_vmubar_1tau_coh_err = np.array([2.674e-06, 2.476e-06])
MINOS_neutrino_vmubar_1tau_incoh_err = np.array([5.679e-03, 5.154e-03])

## MINOS and MINOS+ Antineutrino Mode -- Event Uncertainty ##
# vmu #
MINOS_antineutrino_vmu_2mu_coh_err = np.array([2.417e-01])
MINOS_antineutrino_vmu_2mu_incoh_err = np.array([2.937e-01])

MINOS_antineutrino_vmu_1tau_coh_err = np.array([1.649e-06])
MINOS_antineutrino_vmu_1tau_incoh_err = np.array([3.475e-03])

# vmubar #
MINOS_antineutrino_vmubar_2mu_coh_err = np.array([2.805e-01])
MINOS_antineutrino_vmubar_2mu_incoh_err = np.array([3.911e-01])

MINOS_antineutrino_vmubar_1tau_coh_err = np.array([7.429e-07])
MINOS_antineutrino_vmubar_1tau_incoh_err = np.array([1.613e-03])

## T2K INGRID Phase 1 and Phase 2 -- Events ##
# vmu #
INGRID_neutrino_vmu_2mu_coh = np.array([2.649e+01, 6.792e+01])
INGRID_neutrino_vmu_2mu_incoh = np.array([1.288e+01, 3.302e+01])

INGRID_neutrino_vmu_1tau_coh = np.array([5.901e-11, 1.513e-10])
INGRID_neutrino_vmu_1tau_incoh = np.array([2.256e-10, 5.784e-10])

# vmubar #
INGRID_neutrino_vmubar_2mu_coh = np.array([1.661, 4.259])
INGRID_neutrino_vmubar_2mu_incoh = np.array([8.074e-01, 2.070])

INGRID_neutrino_vmubar_1tau_coh = np.array([3.700e-12, 9.487e-12])
INGRID_neutrino_vmubar_1tau_incoh = np.array([1.414e-11, 3.626e-11])

## T2K INGRID Phase 1 and Phase 2 -- Event Uncertainty ##
# vmu #
INGRID_neutrino_vmu_2mu_coh_err = np.array([2.154, 5.522])
INGRID_neutrino_vmu_2mu_incoh_err = np.array([3.975, 1.019e+01])

INGRID_neutrino_vmu_1tau_coh_err = np.array([3.752e-12, 9.619e-12])
INGRID_neutrino_vmu_1tau_incoh_err = np.array([6.950e-11, 1.782e-10])

# vmubar #
INGRID_neutrino_vmubar_2mu_coh_err = np.array([1.350e-01, 3.462e-01])
INGRID_neutrino_vmubar_2mu_incoh_err = np.array([2.493e-01, 6.391e-01])

INGRID_neutrino_vmubar_1tau_coh_err = np.array([2.352e-13, 6.032e-13])
INGRID_neutrino_vmubar_1tau_incoh_err = np.array([4.358e-12, 1.117e-11])

## FASER -- Events ##
# vmu + vmubar #
FASER_2mu_coh = np.array([1.859e-01, 6.762e+01])
FASER_2mu_incoh = np.array([2.226e-02, 8.094])

FASER_1tau_coh = np.array([5.497e-02, 1.999e+01])
FASER_1tau_incoh = np.array([1.712e-02, 6.226])

FASER_2tau_coh = np.array([1.977e-03, 7.188e-01])
FASER_2tau_incoh = np.array([1.475e-03, 5.362e-01])

## FASER -- Event Uncertainty ##
# vmu + vmubar #
FASER_2mu_coh_err = np.array([1.514e-02, 5.507])
FASER_2mu_incoh_err = np.array([6.957e-03, 2.530])

FASER_1tau_coh_err = np.array([4.474e-03, 1.627])
FASER_1tau_incoh_err = np.array([5.350e-03, 1.946])

FASER_2tau_coh_err = np.array([1.606e-04, 5.838e-02])
FASER_2tau_incoh_err = np.array([4.605e-04, 1.674e-01])

## X (Events) Limits ##
# DUNE #
LEFT_LIMIT_DUNE_NEUTRINO_COH =  1e-8
LEFT_LIMIT_DUNE_NEUTRINO_INCOH = 4e-4

LEFT_LIMIT_DUNE_ANTINEUTRINO_COH =  1e-8
LEFT_LIMIT_DUNE_ANTINEUTRINO_INCOH = 4e-4

RIGHT_LIMIT_DUNE_NEUTRINO_COH =  3e3
RIGHT_LIMIT_DUNE_NEUTRINO_INCOH = 2e3

RIGHT_LIMIT_DUNE_ANTINEUTRINO_COH = 2e3
RIGHT_LIMIT_DUNE_ANTINEUTRINO_INCOH = 2e3

X_TEXT_DUNE_NEUTRINO_COH = 2.3e3
X_TEXT_DUNE_NEUTRINO_INCOH = 1.5e3
X_TEXT_DUNE_ANTINEUTRINO_COH = 1.5e3
X_TEXT_DUNE_ANTINEUTRINO_INCOH = 1.5e3

# MINOS #
LEFT_LIMIT_MINOS_NEUTRINO_COH =  1e-6
LEFT_LIMIT_MINOS_NEUTRINO_INCOH = 1e-3

LEFT_LIMIT_MINOS_ANTINEUTRINO_COH =  3e-7
LEFT_LIMIT_MINOS_ANTINEUTRINO_INCOH = 6e-4

RIGHT_LIMIT_MINOS_NEUTRINO_COH =  3e2
RIGHT_LIMIT_MINOS_NEUTRINO_INCOH = 1e2

RIGHT_LIMIT_MINOS_ANTINEUTRINO_COH = 1e2
RIGHT_LIMIT_MINOS_ANTINEUTRINO_INCOH = 7e0

X_TEXT_MINOS_NEUTRINO_COH = 2.5e2
X_TEXT_MINOS_NEUTRINO_INCOH = 8e1
X_TEXT_MINOS_ANTINEUTRINO_COH = 8e1
X_TEXT_MINOS_ANTINEUTRINO_INCOH = 6e0

# T2K INGRID #
LEFT_LIMIT_INGRID_NEUTRINO_COH =  2e-14
LEFT_LIMIT_INGRID_NEUTRINO_INCOH = 1e-13

RIGHT_LIMIT_INGRID_NEUTRINO_COH =  1e3
RIGHT_LIMIT_INGRID_NEUTRINO_INCOH = 5e2

X_TEXT_INGRID_NEUTRINO_COH = 9e2
X_TEXT_INGRID_NEUTRINO_INCOH = 4e2

# FASER #
LEFT_LIMIT_FASER_COH =  1e-4
LEFT_LIMIT_FASER_INCOH = 1e-4

RIGHT_LIMIT_FASER_COH =  3e2
RIGHT_LIMIT_FASER_INCOH = 2e1

X_TEXT_FASER_COH = 2.5e2
X_TEXT_FASER_INCOH = 1.75e1


### Plotting Setup ###
## DUNE ##
# Colors #
# Extracted from a color-blind-friendly palette: https://davidmathlogic.com/colorblind/#%23663399-%23008080-%23DAA520-%23A71342-%231E3282-%23c6b276-%2325B5BF-%23C52FC5
DUNE_cols = ['#663399',
             '#008080',
             '#DAA520',
             '#A71342',
             '#1E3282',
             '#C6B276',
             '#25B5BF',
             '#C52FC5']

# Coherent #
DUNE_neutrino_coherent_events = [DUNE_neutrino_vmu_2mu_coh,
                                 DUNE_neutrino_vmubar_2mu_coh,
                                 DUNE_neutrino_vmu_1tau_coh,     
                                 DUNE_neutrino_vmubar_1tau_coh,
                                 DUNE_neutrino_ve_1tau_coh,
                                 DUNE_neutrino_vebar_1tau_coh,
                                 DUNE_neutrino_vmu_2tau_coh,
                                 DUNE_neutrino_vmubar_2tau_coh]

DUNE_antineutrino_coherent_events = [DUNE_antineutrino_vmu_2mu_coh,
                                     DUNE_antineutrino_vmubar_2mu_coh,
                                     DUNE_antineutrino_vmu_1tau_coh,
                                     DUNE_antineutrino_vmubar_1tau_coh,
                                     DUNE_antineutrino_ve_1tau_coh,
                                     DUNE_antineutrino_vebar_1tau_coh,
                                     DUNE_antineutrino_vmu_2tau_coh,
                                     DUNE_antineutrino_vmubar_2tau_coh]

DUNE_neutrino_coherent_err = [DUNE_neutrino_vmu_2mu_coh_err,
                              DUNE_neutrino_vmubar_2mu_coh_err,
                              DUNE_neutrino_vmu_1tau_coh_err,
                              DUNE_neutrino_vmubar_1tau_coh_err,
                              DUNE_neutrino_ve_1tau_coh_err,
                              DUNE_neutrino_vebar_1tau_coh_err,
                              DUNE_neutrino_vmu_2tau_coh_err,
                              DUNE_neutrino_vmubar_2tau_coh_err]

DUNE_antineutrino_coherent_err = [DUNE_antineutrino_vmu_2mu_coh_err,
                                  DUNE_antineutrino_vmubar_2mu_coh_err,
                                  DUNE_antineutrino_vmu_1tau_coh_err,
                                  DUNE_antineutrino_vmubar_1tau_coh_err,
                                  DUNE_antineutrino_ve_1tau_coh_err,
                                  DUNE_antineutrino_vebar_1tau_coh_err,
                                  DUNE_antineutrino_vmu_2tau_coh_err,
                                  DUNE_antineutrino_vmubar_2tau_coh_err]

# Incoherent #
DUNE_neutrino_incoherent_events = [DUNE_neutrino_vmu_2mu_incoh,
                                   DUNE_neutrino_vmubar_2mu_incoh,
                                   DUNE_neutrino_vmu_1tau_incoh,
                                   DUNE_neutrino_vmubar_1tau_incoh,
                                   DUNE_neutrino_ve_1tau_incoh,
                                   DUNE_neutrino_vebar_1tau_incoh,
                                   DUNE_neutrino_vmu_2tau_incoh,
                                   DUNE_neutrino_vmubar_2tau_incoh]

DUNE_antineutrino_incoherent_events = [DUNE_antineutrino_vmu_2mu_incoh,
                                       DUNE_antineutrino_vmubar_2mu_incoh,
                                       DUNE_antineutrino_vmu_1tau_incoh,
                                       DUNE_antineutrino_vmubar_1tau_incoh,
                                       DUNE_antineutrino_ve_1tau_incoh,
                                       DUNE_antineutrino_vebar_1tau_incoh,
                                       DUNE_antineutrino_vmu_2tau_incoh,
                                       DUNE_antineutrino_vmubar_2tau_incoh]

DUNE_neutrino_incoherent_err = [DUNE_neutrino_vmu_2mu_incoh_err,
                                DUNE_neutrino_vmubar_2mu_incoh_err,
                                DUNE_neutrino_vmu_1tau_incoh_err,
                                DUNE_neutrino_vmubar_1tau_incoh_err,
                                DUNE_neutrino_ve_1tau_incoh_err,
                                DUNE_neutrino_vebar_1tau_incoh_err,
                                DUNE_neutrino_vmu_2tau_incoh_err,
                                DUNE_neutrino_vmubar_2tau_incoh_err]

DUNE_antineutrino_incoherent_err = [DUNE_antineutrino_vmu_2mu_incoh_err,
                                    DUNE_antineutrino_vmubar_2mu_incoh_err,
                                    DUNE_antineutrino_vmu_1tau_incoh_err,
                                    DUNE_antineutrino_vmubar_1tau_incoh_err,
                                    DUNE_antineutrino_ve_1tau_incoh_err,
                                    DUNE_antineutrino_vebar_1tau_incoh_err,
                                    DUNE_antineutrino_vmu_2tau_incoh_err,
                                    DUNE_antineutrino_vmubar_2tau_incoh_err]

# Axes #
DUNE_neutrino_coherent_axes = [ax11, ax21, ax31, ax41, ax51, ax61, ax71, ax81]
DUNE_neutrino_incoherent_axes = [ax12, ax22, ax32, ax42, ax52, ax62, ax72, ax82]
DUNE_antineutrino_coherent_axes = [ax13, ax23, ax33, ax43, ax53, ax63, ax73, ax83]
DUNE_antineutrino_incoherent_axes = [ax14, ax24, ax34, ax44, ax54, ax64, ax74, ax84]

# All Data #
all_data = [[DUNE_neutrino_coherent_events, DUNE_neutrino_coherent_err, DUNE_neutrino_coherent_axes],
            [DUNE_neutrino_incoherent_events, DUNE_neutrino_incoherent_err, DUNE_neutrino_incoherent_axes],
            [DUNE_antineutrino_coherent_events, DUNE_antineutrino_coherent_err, DUNE_antineutrino_coherent_axes],
            [DUNE_antineutrino_incoherent_events, DUNE_antineutrino_incoherent_err, DUNE_antineutrino_incoherent_axes]]

all_axes = [DUNE_neutrino_coherent_axes, DUNE_neutrino_incoherent_axes, DUNE_antineutrino_coherent_axes, DUNE_antineutrino_incoherent_axes]

all_processes = [r'$\nu_\mu \mu^+ \mu^-$',
                 r'$\bar{\nu}_\mu \mu^+ \mu^-$',
                 r'$\nu_\tau \tau^+ \mu^-$',
                 r'$\bar{\nu}_\tau \mu^+ \tau^-$',
                 r'$\nu_\tau \tau^+ e^-$',
                 r'$\bar{\nu}_\tau e^+ \tau^-$',
                 r'$\nu_\mu \tau^+ \tau^-$',
                 r'$\bar{\nu}_\mu \tau^+ \tau^-$']

all_label_pos = [X_TEXT_DUNE_NEUTRINO_COH, X_TEXT_DUNE_NEUTRINO_INCOH, X_TEXT_DUNE_ANTINEUTRINO_COH, X_TEXT_DUNE_ANTINEUTRINO_INCOH]


DUNE_parameters = [r'\textbf{Standard}',
                   r'\textbf{$\tau$-Optimized}']

## MINOS ##
# Colors #
# Extracted from a color-blind-friendly palette: https://davidmathlogic.com/colorblind/#%23663399-%23008080-%23DAA520-%23A71342-%231E3282-%23c6b276-%2325B5BF-%23C52FC5
MINOS_cols = ['#663399',
             '#008080',
             '#DAA520',
             '#A71342']

# Coherent #
MINOS_neutrino_coherent_events = [MINOS_neutrino_vmu_2mu_coh,
                                 MINOS_neutrino_vmubar_2mu_coh,
                                 MINOS_neutrino_vmu_1tau_coh,     
                                 MINOS_neutrino_vmubar_1tau_coh]

MINOS_antineutrino_coherent_events = [MINOS_antineutrino_vmu_2mu_coh,
                                     MINOS_antineutrino_vmubar_2mu_coh,
                                     MINOS_antineutrino_vmu_1tau_coh,
                                     MINOS_antineutrino_vmubar_1tau_coh]

MINOS_neutrino_coherent_err = [MINOS_neutrino_vmu_2mu_coh_err,
                              MINOS_neutrino_vmubar_2mu_coh_err,
                              MINOS_neutrino_vmu_1tau_coh_err,
                              MINOS_neutrino_vmubar_1tau_coh_err]

MINOS_antineutrino_coherent_err = [MINOS_antineutrino_vmu_2mu_coh_err,
                                  MINOS_antineutrino_vmubar_2mu_coh_err,
                                  MINOS_antineutrino_vmu_1tau_coh_err,
                                  MINOS_antineutrino_vmubar_1tau_coh_err]

# Incoherent #
MINOS_neutrino_incoherent_events = [MINOS_neutrino_vmu_2mu_incoh,
                                   MINOS_neutrino_vmubar_2mu_incoh,
                                   MINOS_neutrino_vmu_1tau_incoh,
                                   MINOS_neutrino_vmubar_1tau_incoh]

MINOS_antineutrino_incoherent_events = [MINOS_antineutrino_vmu_2mu_incoh,
                                       MINOS_antineutrino_vmubar_2mu_incoh,
                                       MINOS_antineutrino_vmu_1tau_incoh,
                                       MINOS_antineutrino_vmubar_1tau_incoh]

MINOS_neutrino_incoherent_err = [MINOS_neutrino_vmu_2mu_incoh_err,
                                MINOS_neutrino_vmubar_2mu_incoh_err,
                                MINOS_neutrino_vmu_1tau_incoh_err,
                                MINOS_neutrino_vmubar_1tau_incoh_err]

MINOS_antineutrino_incoherent_err = [MINOS_antineutrino_vmu_2mu_incoh_err,
                                    MINOS_antineutrino_vmubar_2mu_incoh_err,
                                    MINOS_antineutrino_vmu_1tau_incoh_err,
                                    MINOS_antineutrino_vmubar_1tau_incoh_err]

# Axes #
MINOS_neutrino_coherent_axes = [ax15, ax25, ax35, ax45]
MINOS_neutrino_incoherent_axes = [ax16, ax26, ax36, ax46]
MINOS_antineutrino_coherent_axes = [ax17, ax27, ax37, ax47]
MINOS_antineutrino_incoherent_axes = [ax18, ax28, ax38, ax48]

# All Data #
all_data_MINOS_neutrino = [[MINOS_neutrino_coherent_events, MINOS_neutrino_coherent_err, MINOS_neutrino_coherent_axes],
                           [MINOS_neutrino_incoherent_events, MINOS_neutrino_incoherent_err, MINOS_neutrino_incoherent_axes]]
all_data_MINOS_antineutrino = [[MINOS_antineutrino_coherent_events, MINOS_antineutrino_coherent_err, MINOS_antineutrino_coherent_axes],
                               [MINOS_antineutrino_incoherent_events, MINOS_antineutrino_incoherent_err, MINOS_antineutrino_incoherent_axes]]

all_axes_MINOS_neutrino = [MINOS_neutrino_coherent_axes, MINOS_neutrino_incoherent_axes] 
all_axes_MINOS_antineutrino = [MINOS_antineutrino_coherent_axes, MINOS_antineutrino_incoherent_axes]

all_processes_MINOS = [r'$\nu_\mu \mu^+ \mu^-$',
                       r'$\bar{\nu}_\mu \mu^+ \mu^-$',
                       r'$\nu_\tau \tau^+ \mu^-$',
                       r'$\bar{\nu}_\tau \mu^+ \tau^-$']

all_label_pos_MINOS_neutrino = [X_TEXT_MINOS_NEUTRINO_COH, X_TEXT_MINOS_NEUTRINO_INCOH]
all_label_pos_MINOS_antineutrino = [X_TEXT_MINOS_ANTINEUTRINO_COH, X_TEXT_MINOS_ANTINEUTRINO_INCOH]


MINOS_neutrino_parameters = [r'\textbf{MINOS}',
                             r'\textbf{MINOS+}']

MINOS_antineutrino_parameters = [r'\textbf{MINOS}']

## T2K INGRID ##
# Colors #
# Extracted from a color-blind-friendly palette: https://davidmathlogic.com/colorblind/#%23663399-%23008080-%23DAA520-%23A71342-%231E3282-%23c6b276-%2325B5BF-%23C52FC5
INGRID_cols = ['#663399',
               '#008080',
               '#DAA520',
               '#A71342']

# Coherent #
INGRID_neutrino_coherent_events = [INGRID_neutrino_vmu_2mu_coh,
                                   INGRID_neutrino_vmubar_2mu_coh,
                                   INGRID_neutrino_vmu_1tau_coh,     
                                   INGRID_neutrino_vmubar_1tau_coh]

INGRID_neutrino_coherent_err = [INGRID_neutrino_vmu_2mu_coh_err,
                                INGRID_neutrino_vmubar_2mu_coh_err,
                                INGRID_neutrino_vmu_1tau_coh_err,
                                INGRID_neutrino_vmubar_1tau_coh_err]

# Incoherent #
INGRID_neutrino_incoherent_events = [INGRID_neutrino_vmu_2mu_incoh,
                                     INGRID_neutrino_vmubar_2mu_incoh,
                                     INGRID_neutrino_vmu_1tau_incoh,
                                     INGRID_neutrino_vmubar_1tau_incoh]

INGRID_neutrino_incoherent_err = [INGRID_neutrino_vmu_2mu_incoh_err,
                                  INGRID_neutrino_vmubar_2mu_incoh_err,
                                  INGRID_neutrino_vmu_1tau_incoh_err,
                                  INGRID_neutrino_vmubar_1tau_incoh_err]

# Axes #
INGRID_neutrino_coherent_axes = [ax19, ax29, ax39, ax49]
INGRID_neutrino_incoherent_axes = [ax110, ax210, ax310, ax410]

# All Data #
all_data_INGRID_neutrino = [[INGRID_neutrino_coherent_events, INGRID_neutrino_coherent_err, INGRID_neutrino_coherent_axes],
                            [INGRID_neutrino_incoherent_events, INGRID_neutrino_incoherent_err, INGRID_neutrino_incoherent_axes]]

all_axes_INGRID_neutrino = [INGRID_neutrino_coherent_axes, INGRID_neutrino_incoherent_axes] 

all_processes_INGRID = [r'$\nu_\mu \mu^+ \mu^-$',
                        r'$\bar{\nu}_\mu \mu^+ \mu^-$',
                        r'$\nu_\tau \tau^+ \mu^-$',
                        r'$\bar{\nu}_\tau \mu^+ \tau^-$']

all_label_pos_INGRID_neutrino = [X_TEXT_INGRID_NEUTRINO_COH, X_TEXT_INGRID_NEUTRINO_INCOH]

INGRID_neutrino_parameters = [r'\textbf{Phase 1}',
                             r'\textbf{Phase 2}']

## FASER ##
# Colors #
# Extracted from a color-blind-friendly palette: https://davidmathlogic.com/colorblind/#%23663399-%23008080-%23DAA520-%23A71342-%231E3282-%23c6b276-%2325B5BF-%23C52FC5
FASER_cols = ['#663399',
               '#008080',
               '#DAA520']

# Coherent #
FASER_coherent_events = [FASER_2mu_coh,
                         FASER_1tau_coh,
                         FASER_2tau_coh]

FASER_coherent_err = [FASER_2mu_coh_err,
                      FASER_1tau_coh_err,
                      FASER_2tau_coh_err]

# Incoherent #
FASER_incoherent_events = [FASER_2mu_incoh,
                           FASER_1tau_incoh,
                           FASER_2tau_incoh]

FASER_incoherent_err = [FASER_2mu_incoh_err,
                        FASER_1tau_incoh_err,
                        FASER_2tau_incoh_err]

# Axes #
FASER_coherent_axes = [ax111, ax211, ax311]
FASER_incoherent_axes = [ax112, ax212, ax312]

# All Data #
all_data_FASER = [[FASER_coherent_events, FASER_coherent_err, FASER_coherent_axes],
                  [FASER_incoherent_events, FASER_incoherent_err, FASER_incoherent_axes]]

all_axes_FASER = [FASER_coherent_axes, FASER_incoherent_axes] 

all_processes_FASER = [r'$\nu_\mu \mu^+ \mu^-$',
                       r'$\nu_\tau \tau^+ \mu^-$',
                       r'$\nu_\mu \tau^+ \tau^-$']

all_label_pos_FASER = [X_TEXT_FASER_COH, X_TEXT_FASER_INCOH]

FASER_parameters = [r'\textbf{FASER$\nu$}', r'\textbf{FASER$\nu$2}']

N_DUNE_PARAMETERS = len(DUNE_parameters) - 1
N_MINOS_NEUTRINO_PARAMETERS = len(MINOS_neutrino_parameters) - 1
N_MINOS_ANTINEUTRINO_PARAMETERS = len(MINOS_antineutrino_parameters) - 1
N_INGRID_NEUTRINO_PARAMETERS = len(INGRID_neutrino_parameters) - 1
N_FASER_PARAMETERS = len(FASER_parameters) - 1

### Plotting ###
# Note: the y-axis labels (i.e. parameters) get assigned an integer y-value indexed at 0.
# So, e.g., 1t, 1yr, DUNE will have y=3, 147t, 3yr, DUNE will have y=2 and so on ending at y=0

# DUNE Neutrino Coherent #
for ax, data, data_err, col in zip(DUNE_neutrino_coherent_axes, DUNE_neutrino_coherent_events, DUNE_neutrino_coherent_err, DUNE_cols):
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_NEUTRINO_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_NEUTRINO_COH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# DUNE Neutrino Incoherent #
for ax, data, data_err, col in zip(DUNE_neutrino_incoherent_axes, DUNE_neutrino_incoherent_events, DUNE_neutrino_incoherent_err, DUNE_cols):
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_NEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_NEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# DUNE Antineutrino Coherent #
for ax, data, data_err, col in zip(DUNE_antineutrino_coherent_axes, DUNE_antineutrino_coherent_events, DUNE_antineutrino_coherent_err, DUNE_cols):
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_ANTINEUTRINO_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_ANTINEUTRINO_COH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# DUNE Antineutrino Incoherent #
for ax, data, data_err, col in zip(DUNE_antineutrino_incoherent_axes, DUNE_antineutrino_incoherent_events, DUNE_antineutrino_incoherent_err, DUNE_cols):
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_ANTINEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(DUNE_parameters), np.flipud(data), left=LEFT_LIMIT_DUNE_ANTINEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# MINOS Neutrino Coherent #
for ax, data, data_err, col in zip(MINOS_neutrino_coherent_axes, MINOS_neutrino_coherent_events, MINOS_neutrino_coherent_err, MINOS_cols):
    ax.barh(np.flipud(MINOS_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_NEUTRINO_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(MINOS_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_NEUTRINO_COH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# MINOS Neutrino Incoherent #
for ax, data, data_err, col in zip(MINOS_neutrino_incoherent_axes, MINOS_neutrino_incoherent_events, MINOS_neutrino_incoherent_err, MINOS_cols):
    ax.barh(np.flipud(MINOS_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_NEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(MINOS_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_NEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# MINOS Antineutrino Coherent #
for ax, data, data_err, col in zip(MINOS_antineutrino_coherent_axes, MINOS_antineutrino_coherent_events, MINOS_antineutrino_coherent_err, MINOS_cols):
    ax.barh(np.flipud(MINOS_antineutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_ANTINEUTRINO_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(MINOS_antineutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_ANTINEUTRINO_COH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# MINOS Antineutrino Incoherent #
for ax, data, data_err, col in zip(MINOS_antineutrino_incoherent_axes, MINOS_antineutrino_incoherent_events, MINOS_antineutrino_incoherent_err, MINOS_cols):
    ax.barh(np.flipud(MINOS_antineutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_ANTINEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(MINOS_antineutrino_parameters), np.flipud(data), left=LEFT_LIMIT_MINOS_ANTINEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# INGRID Neutrino Coherent #
for ax, data, data_err, col in zip(INGRID_neutrino_coherent_axes, INGRID_neutrino_coherent_events, INGRID_neutrino_coherent_err, INGRID_cols):
    ax.barh(np.flipud(INGRID_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_COH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(INGRID_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_COH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# INGRID Neutrino Incoherent #
for ax, data, data_err, col in zip(INGRID_neutrino_incoherent_axes, INGRID_neutrino_incoherent_events, INGRID_neutrino_incoherent_err, INGRID_cols):
    ax.barh(np.flipud(INGRID_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(INGRID_neutrino_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# FASER Coherent #
for ax, data, data_err, col in zip(FASER_coherent_axes, FASER_coherent_events, FASER_coherent_err, FASER_cols):
    ax.barh(np.flipud(FASER_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(FASER_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# FASER Incoherent #
for ax, data, data_err, col in zip(FASER_incoherent_axes, FASER_incoherent_events, FASER_incoherent_err, FASER_cols):
    ax.barh(np.flipud(FASER_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_INCOH, align='center', height=1, edgecolor='k', facecolor='w', lw=1.5, zorder=-10, alpha=1)
    ax.barh(np.flipud(FASER_parameters), np.flipud(data), left=LEFT_LIMIT_INGRID_NEUTRINO_INCOH, xerr=np.flipud(data_err), align='center', height=0.98, facecolor=col, lw=0, zorder=-10, alpha=0.5, error_kw=dict(markeredgewidth=3,markeredgecolor='k',solid_joinstyle='round',capsize=8,elinewidth=3, ecolor=col))
    ax.patch.set_facecolor(col)
    ax.patch.set_alpha(0.1)

# Adding text labels
t_shift = [0.01,0.01,0.01,0.0,0.0,0.0,0.0,0.0] # Manual vertical shift for now
N_VALUES_DUNE = N_DUNE_PARAMETERS + 1
N_VALUES_MINOS_NEUTRINO = N_MINOS_NEUTRINO_PARAMETERS + 1
N_VALUES_MINOS_ANTINEUTRINO = N_MINOS_ANTINEUTRINO_PARAMETERS + 1
N_VALUES_INGRID_NEUTRINO = N_INGRID_NEUTRINO_PARAMETERS + 1
N_VALUES_FASER = N_FASER_PARAMETERS + 1

def val_to_str(value, error):
    if value >= 0.1:
        string = r"{\bf %.2f} $\pm$ {\bf %.2f}" % (value, error)
        return string
    else:
        string = r"{\bf %.1e} $\pm$ {\bf %.1e}" % (value, error)
        return string

for collection in all_data:
    event_array = collection[0]
    error_array = collection[1]
    axes = collection[2]
    for j in range(len(axes)):
        for i in range(0, N_VALUES_DUNE):
            event_val = event_array[j][i]
            error_val = error_array[j][i]
            string = val_to_str(event_val, error_val)
            axes[j].text(event_val - 2*event_val/5, (N_VALUES_DUNE-1)-i-t_shift[i], string, va='center', ha='right', fontsize=30)

for collection in all_data_MINOS_neutrino:
    event_array = collection[0]
    error_array = collection[1]
    axes = collection[2]
    for j in range(len(axes)):
        for i in range(0, N_VALUES_MINOS_NEUTRINO):
            event_val = event_array[j][i]
            error_val = error_array[j][i]
            string = val_to_str(event_val, error_val)
            axes[j].text(event_val - 2*event_val/5, (N_VALUES_MINOS_NEUTRINO-1)-i-t_shift[i], string, va='center', ha='right', fontsize=22)

for collection in all_data_MINOS_antineutrino:
    event_array = collection[0]
    error_array = collection[1]
    axes = collection[2]
    for j in range(len(axes)):
        for i in range(0, N_VALUES_MINOS_ANTINEUTRINO):
            event_val = event_array[j][i]
            error_val = error_array[j][i]
            string = val_to_str(event_val, error_val)
            axes[j].text(event_val - 2*event_val/5, (N_VALUES_MINOS_ANTINEUTRINO-1)-i-t_shift[i], string, va='center', ha='right', fontsize=22)

for collection in all_data_INGRID_neutrino:
    event_array = collection[0]
    error_array = collection[1]
    axes = collection[2]
    for j in range(len(axes)):
        for i in range(0, N_VALUES_INGRID_NEUTRINO):
            event_val = event_array[j][i]
            error_val = error_array[j][i]
            string = val_to_str(event_val, error_val)
            axes[j].text(event_val - 2*event_val/5, (N_VALUES_INGRID_NEUTRINO-1)-i-t_shift[i], string, va='center', ha='right', fontsize=22)

for collection in all_data_FASER:
    event_array = collection[0]
    error_array = collection[1]
    axes = collection[2]
    for j in range(len(axes)):
        for i in range(0, N_VALUES_FASER):
            event_val = event_array[j][i]
            error_val = error_array[j][i]
            string = val_to_str(event_val, error_val)
            axes[j].text(event_val - 2*event_val/5, (N_VALUES_FASER-1)-i-t_shift[i], string, va='center', ha='right', fontsize=22)

## Limits ##
# DUNE #
ax81.set_xscale('log')
ax81.set_xlim((LEFT_LIMIT_DUNE_NEUTRINO_COH,RIGHT_LIMIT_DUNE_NEUTRINO_COH))

ax82.set_xscale('log')
ax82.set_xlim([LEFT_LIMIT_DUNE_NEUTRINO_INCOH,RIGHT_LIMIT_DUNE_NEUTRINO_INCOH])

ax83.set_xscale('log')
ax83.set_xlim([LEFT_LIMIT_DUNE_ANTINEUTRINO_COH,RIGHT_LIMIT_DUNE_ANTINEUTRINO_COH])

ax84.set_xscale('log')
ax84.set_xlim([LEFT_LIMIT_DUNE_ANTINEUTRINO_INCOH,RIGHT_LIMIT_DUNE_ANTINEUTRINO_INCOH])

# MINOS
ax45.set_xscale('log')
ax45.set_xlim((LEFT_LIMIT_MINOS_NEUTRINO_COH,RIGHT_LIMIT_MINOS_NEUTRINO_COH))

ax46.set_xscale('log')
ax46.set_xlim([LEFT_LIMIT_MINOS_NEUTRINO_INCOH,RIGHT_LIMIT_MINOS_NEUTRINO_INCOH])

ax47.set_xscale('log')
ax47.set_xlim([LEFT_LIMIT_MINOS_ANTINEUTRINO_COH,RIGHT_LIMIT_MINOS_ANTINEUTRINO_COH])

ax48.set_xscale('log')
ax48.set_xlim([LEFT_LIMIT_MINOS_ANTINEUTRINO_INCOH,RIGHT_LIMIT_MINOS_ANTINEUTRINO_INCOH])

# INGRID
ax49.set_xscale('log')
ax49.set_xlim((LEFT_LIMIT_INGRID_NEUTRINO_COH,RIGHT_LIMIT_INGRID_NEUTRINO_COH))

ax410.set_xscale('log')
ax410.set_xlim([LEFT_LIMIT_INGRID_NEUTRINO_INCOH,RIGHT_LIMIT_INGRID_NEUTRINO_INCOH])

# FASER
ax311.set_xscale('log')
ax311.set_xlim((LEFT_LIMIT_FASER_COH,RIGHT_LIMIT_FASER_COH))

ax312.set_xscale('log')
ax312.set_xlim([LEFT_LIMIT_FASER_INCOH,RIGHT_LIMIT_FASER_INCOH])

## Fixing ticks ##
# VERY IMPORTANT NOTE ON FORMATTING: Assign a separate major and minor locator to each Axis. The locator stores references to the Axis data and view limits so assigning the same
# major or minor locator to multiple axes will give strange behavior, such as not giving the correct ticks for the given range.

# DUNE
locmaj1 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin1 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj2 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin2 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj3 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin3 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj4 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin4 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)

ax81.xaxis.set_major_locator(locmaj1)
ax81.xaxis.set_minor_locator(locmin1)

ax82.xaxis.set_major_locator(locmaj2)
ax82.xaxis.set_minor_locator(locmin2)

ax83.xaxis.set_major_locator(locmaj3)
ax83.xaxis.set_minor_locator(locmin3)

ax84.xaxis.set_major_locator(locmaj4)
ax84.xaxis.set_minor_locator(locmin4)

# MINOS
locmaj5 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin5 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj6 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin6 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj7 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin7 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj8 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin8 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)

ax45.xaxis.set_major_locator(locmaj5)
ax45.xaxis.set_minor_locator(locmin5)

ax46.xaxis.set_major_locator(locmaj6)
ax46.xaxis.set_minor_locator(locmin6)

ax47.xaxis.set_major_locator(locmaj7)
ax47.xaxis.set_minor_locator(locmin7)

ax48.xaxis.set_major_locator(locmaj8)
ax48.xaxis.set_minor_locator(locmin8)

# INGRID
locmaj9 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin9 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj10 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin10 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)

ax49.xaxis.set_major_locator(locmaj9)
ax49.xaxis.set_minor_locator(locmin9)

ax410.xaxis.set_major_locator(locmaj10)
ax410.xaxis.set_minor_locator(locmin10)

# FASER
locmaj11 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin11 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)
locmaj12 = matplotlib.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=500)
locmin12 = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=500)

ax311.xaxis.set_major_locator(locmaj11)
ax311.xaxis.set_minor_locator(locmin11)

ax312.xaxis.set_major_locator(locmaj12)
ax312.xaxis.set_minor_locator(locmin12)

## Grid Lines, Tick Parameters and Labels ##

# DUNE
y_label_loc_DUNE = N_DUNE_PARAMETERS - 0.05

for l, axes in enumerate(all_axes):
    for k, ax in enumerate(axes):
        # Grid lines
        ax.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
        # Tick parameters
        ax.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
        ax.tick_params(axis='x',labelsize=40)
        # Labels
        ax.text(all_label_pos[l],y_label_loc_DUNE,all_processes[k],fontsize=35,alpha=0.7,color=DUNE_cols[k],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

ax81.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax82.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax83.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax84.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)

#ax81.tick_params(axis='x',labelsize=40)

# MINOS
y_label_loc_MINOS_neutrino = N_MINOS_NEUTRINO_PARAMETERS - 0.05
y_label_loc_MINOS_antineutrino = N_MINOS_ANTINEUTRINO_PARAMETERS - 0.05

for l, axes in enumerate(all_axes_MINOS_neutrino):
    for k, ax in enumerate(axes):
        # Grid lines
        ax.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
        # Tick parameters
        ax.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
        # Labels
        ax.text(all_label_pos_MINOS_neutrino[l],y_label_loc_MINOS_neutrino,all_processes_MINOS[k],fontsize=35,alpha=0.7,color=MINOS_cols[k],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

for l, axes in enumerate(all_axes_MINOS_antineutrino):
    for k, ax in enumerate(axes):
        # Grid lines
        ax.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
        # Tick parameters
        ax.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
        # Labels
        ax.text(all_label_pos_MINOS_antineutrino[l],y_label_loc_MINOS_antineutrino,all_processes_MINOS[k],fontsize=35,alpha=0.7,color=MINOS_cols[k],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

ax45.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax46.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax47.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax48.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)

# INGRID
y_label_loc_INGRID_neutrino = N_INGRID_NEUTRINO_PARAMETERS - 0.05

for l, axes in enumerate(all_axes_INGRID_neutrino):
    for k, ax in enumerate(axes):
        # Grid lines
        ax.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
        # Tick parameters
        ax.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
        # Labels
        ax.text(all_label_pos_INGRID_neutrino[l],y_label_loc_INGRID_neutrino,all_processes_INGRID[k],fontsize=35,alpha=0.7,color=INGRID_cols[k],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

ax49.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax410.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)

# FASER
y_label_loc_FASER = N_FASER_PARAMETERS - 0.05

for l, axes in enumerate(all_axes_FASER):
    for k, ax in enumerate(axes):
        # Grid lines
        ax.xaxis.grid(True, linestyle='--', which='major',color='grey', alpha=.45)
        # Tick parameters
        ax.tick_params(which='both',right=False,tickdir='out',top=False,bottom=False)
        # Labels
        ax.text(all_label_pos_FASER[l],y_label_loc_FASER,all_processes_FASER[k],fontsize=35,alpha=0.7,color=FASER_cols[k],path_effects=[pe.Stroke(linewidth=1.5, foreground='k',alpha=0.5), pe.Normal()],ha='right')

ax311.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)
ax312.tick_params(which='both',right=False,tickdir='out',top=False,bottom=True)

## Axes Labels ##
# DUNE
ax81.set_xlabel(r'{\bf Number of Events}')
ax11.set_title(r'{\bf Neutrino Mode - Coherent scattering}', fontsize=60)

ax82.set_xlabel(r'{\bf Number of Events}')
ax12.set_title(r'{\bf Neutrino Mode - Incoherent scattering}', fontsize=60)

ax83.set_xlabel(r'{\bf Number of Events}')
ax13.set_title(r'{\bf Antineutrino Mode - Coherent scattering}', fontsize=60)

ax84.set_xlabel(r'{\bf Number of Events}')
ax14.set_title(r'{\bf Antineutrino Mode - Incoherent scattering}', fontsize=60)

# MINOS
ax45.set_xlabel(r'{\bf Number of Events}')
ax15.set_title(r'{\bf Neutrino Mode - Coherent scattering}', fontsize=60)

ax46.set_xlabel(r'{\bf Number of Events}')
ax16.set_title(r'{\bf Neutrino Mode - Incoherent scattering}', fontsize=60)

ax47.set_xlabel(r'{\bf Number of Events}')
ax17.set_title(r'{\bf Antineutrino Mode - Coherent scattering}', fontsize=60)

ax48.set_xlabel(r'{\bf Number of Events}')
ax18.set_title(r'{\bf Antineutrino Mode - Incoherent scattering}', fontsize=60)

# INGRID
ax49.set_xlabel(r'{\bf Number of Events}')
ax19.set_title(r'{\bf Neutrino Mode - Coherent scattering}', fontsize=60)

ax410.set_xlabel(r'{\bf Number of Events}')
ax110.set_title(r'{\bf Neutrino Mode - Incoherent scattering}', fontsize=60)

# INGRID
ax311.set_xlabel(r'{\bf Number of Events}')
ax111.set_title(r'{\bf Coherent scattering}', fontsize=60)

ax312.set_xlabel(r'{\bf Number of Events}')
ax112.set_title(r'{\bf Incoherent scattering}', fontsize=60)

## Save Figures ##
# DUNE
fig1.savefig('../plots/events_DUNE_neutrino_coh.png',transparent=False,bbox_inches='tight',dpi=100)
fig2.savefig('../plots/events_DUNE_neutrino_incoh.png',transparent=False,bbox_inches='tight',dpi=100)
fig3.savefig('../plots/events_DUNE_antineutrino_coh.png',transparent=False,bbox_inches='tight',dpi=100)
fig4.savefig('../plots/events_DUNE_antineutrino_incoh.png',transparent=False,bbox_inches='tight',dpi=100)

# MINOS
fig5.savefig('../plots/events_MINOS_neutrino_coh.png',transparent=False,bbox_inches='tight',dpi=100)
fig6.savefig('../plots/events_MINOS_neutrino_incoh.png',transparent=False,bbox_inches='tight',dpi=100)
fig7.savefig('../plots/events_MINOS_antineutrino_coh.png',transparent=False,bbox_inches='tight',dpi=100)
fig8.savefig('../plots/events_MINOS_antineutrino_incoh.png',transparent=False,bbox_inches='tight',dpi=100)

# INGRID
fig9.savefig('../plots/events_INGRID_neutrino_coh.png',transparent=False,bbox_inches='tight',dpi=100)
fig10.savefig('../plots/events_INGRID_neutrino_incoh.png',transparent=False,bbox_inches='tight',dpi=100)

# FASER
fig11.savefig('../plots/events_FASER_coh.png',transparent=False,bbox_inches='tight',dpi=100)
fig12.savefig('../plots/events_FASER_incoh.png',transparent=False,bbox_inches='tight',dpi=100)
