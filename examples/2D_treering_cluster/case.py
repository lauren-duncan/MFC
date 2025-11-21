#!/usr/bin/env python3
import math
import json
import numpy as np
# Bubble screen
# Description: A planar acoustic wave interacts with a bubble cloud
# in water. The background field is modeled in using an Eulerian framework,
# while the bubbles are tracked using a Lagrangian framework.

# Reference values for nondimensionalization
x0 = 1.0e-3  # length - m
rho0 = 700.0  # density - kg/m3
patm = 101325.0  # Atmospheric pressure - Pa
pi_inf_host = 7e9  # Stiffness - Pa
G_modulus = 7.E7
gamma_host = 2.  # Specific heat ratio
c0 = np.sqrt(gamma_host*(patm+pi_inf_host)/rho0)
p0 = rho0 * c0 * c0  # pressure - Pa
T0 = 298  # temperature - K
#print(c0)
# Host properties (water)
mu_host = 1e-3  # Dynamic viscosity - Pa.s
c_host = c0  # speed of sound - m/s
T_host = 298  # temperature K

# air properties
rho_outer = rho0  # kg/m^3
rho_mid = rho0  # kg/m^3
rho_inner = rho0  # kg/m^3
rho_gas = 1.225  # kg/m^3

# Lagrangian bubbles' properties
R_uni = 8314  # Universal gas constant - J/kmol/K
MW_g = 28.0  # Molar weight of the gas - kg/kmol
MW_v = 18.0  # Molar weight of the vapor - kg/kmol
gamma_g = 1.4  # Specific heat ratio of the gas
gamma_v = 1.333  # Specific heat ratio of the vapor
pv = 2350  # Vapor pressure of the host - Pa
cp_g = 1.0e3  # Specific heat of the gas - J/kg/K
cp_v = 2.1e3  # Specific heat of the vapor - J/kg/K
k_g = 0.025  # Thermal conductivity of the gas - W/m/K
k_v = 0.02  # Thermal conductivity of the vapor - W/m/K
diffVapor = 2.5e-5  # Diffusivity coefficient of the vapor - m2/s
sigBubble = 0.069  # Surface tension of the bubble - N/m
mu_g = 1.48e-5
Re_inv_host = 1.0 / (mu_host / (rho0 * c0 * x0))
Re_inv_gas = 1.0 / (mu_g / (rho0 * c0 * x0))
Cau_inv = G_modulus / p0
# Acoustic source properties
pamp = 1.0e5  # Amplitude of the acoustic source - Pa
freq = 300e03  # Source frequency - Hz
wlen = c_host / freq  # Wavelength - m

# Domain and time set up

L = 2e-3
xb = -L  # Domain boundaries - m (x direction)
xe = L
yb = -L  # Domain boundaries - m (y direction)
ye = L
outer_radius = 1.4

Nx = 399  # number of elements into x direction
Ny = 399  # number of elements into y direction

# uniform grid
dx = 2.*L/Nx 
z_virtual = L  # Virtual depth (z direction)
#CFL = 0.8
dt = 0.0005 #0.009 # CFL * dx / c0
tfinal = 1500 * dt
tstop = int(tfinal/dt)

nframes = 200.
tframes = int(tstop/nframes)

# Configuring case dictionary
print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": xb / x0,
            "x_domain%end": xe / x0,
            "y_domain%beg": yb / x0,
            "y_domain%end": ye / x0,
            "stretch_y": "F",
            "stretch_x": "F",
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": tstop,
            "t_step_save": tframes,
            # Simulation Algorithm Parameters
            "model_eqns": 2,
            "hypoelasticity": "T",
            "fd_order": 4,
            "time_stepper": 3,
            "num_fluids": 2,
            "num_patches": 2,
            "viscous": "T",
            "mpp_lim": "T",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -6,
            "bc_y%end": -6,
            # Acoustic source
            "acoustic_source": "F",
            "num_source": 1,
            "acoustic(1)%support": 2,
            "acoustic(1)%pulse": 1,
            "acoustic(1)%npulse": 1,
            "acoustic(1)%mag": pamp / p0,
            "acoustic(1)%wavelength": wlen / x0,
            "acoustic(1)%length": 2 * (ye - yb) / x0,
            "acoustic(1)%loc(1)": -7.0e-03 / x0,
            "acoustic(1)%loc(2)": 0.0,
            "acoustic(1)%dir": 0.0,
            "acoustic(1)%delay": 0.0,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "lag_db_wrt": "T",
            # Patch 1: Outside air
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": (xe - xb) / x0,
            "patch_icpp(1)%length_y": (ye - yb) / x0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": patm / p0,
            "patch_icpp(1)%alpha_rho(1)": 0.0,
            "patch_icpp(1)%alpha_rho(2)": rho_gas / rho0,
            "patch_icpp(1)%alpha(1)": 0.0,
            "patch_icpp(1)%alpha(2)": 1.0,
            # Patch 2: Outer tree
            "patch_icpp(2)%geometry": 2,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": 0.0,
            "patch_icpp(2)%y_centroid": 0.0,
            "patch_icpp(2)%radius": outer_radius,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": patm / p0,
            "patch_icpp(2)%alpha_rho(1)": rho_outer / rho0,
            "patch_icpp(2)%alpha_rho(2)": 0.0,
            "patch_icpp(2)%alpha(1)": 1.0,
            "patch_icpp(2)%alpha(2)": 0.0,
            # Lagrangian Bubbles
            "bubbles_lagrange": "T",
            "bubble_model": 2,  # Keller-Miksis model
            "Cau_inv": Cau_inv,
            "lag_params%nBubs_glb": 200,  # Number of bubbles
            "lag_params%solver_approach": 2,
            "lag_params%cluster_type": 2,
            "lag_params%pressure_corrector": "T",
            "lag_params%smooth_type": 1,
            "lag_params%heatTransfer_model": "T",
            "lag_params%massTransfer_model": "T",
            "lag_params%epsilonb": 1.0,
            "lag_params%valmaxvoid": 0.9,
            "lag_params%write_bubbles": "F",
            "lag_params%write_bubbles_stats": "F",
            "lag_params%c0": c0,
            "lag_params%rho0": rho0,
            "lag_params%T0": T0,
            "lag_params%x0": x0,
            "lag_params%diffcoefvap": diffVapor,
            "lag_params%Thost": T_host,
            "lag_params%charwidth": z_virtual / x0,
            # Fluids Physical Parameters
            # Host medium
            "fluid_pp(1)%gamma": 1.0 / (gamma_host - 1.0),
            "fluid_pp(1)%pi_inf": gamma_host * ((pi_inf_host) / p0) / (gamma_host - 1.0),
            "fluid_pp(1)%Re(1)": Re_inv_host,
            "fluid_pp(1)%mul0": mu_host,
            "fluid_pp(1)%ss": sigBubble,
            "fluid_pp(1)%pv": pv,
            "fluid_pp(1)%gamma_v": gamma_v,
            "fluid_pp(1)%M_v": MW_v,
            "fluid_pp(1)%k_v": k_v,
            "fluid_pp(1)%cp_v": cp_v,
            "fluid_pp(1)%G": Cau_inv,
            # Bubble gas state
            "fluid_pp(2)%gamma": 1.0 / (gamma_g - 1.0),
            "fluid_pp(2)%pi_inf": 0.0e00,
            "fluid_pp(2)%Re(1)": Re_inv_gas,
            "fluid_pp(2)%gamma_v": gamma_g,
            "fluid_pp(2)%M_v": MW_g,
            "fluid_pp(2)%k_v": k_g,
            "fluid_pp(2)%cp_v": cp_g,
            "fluid_pp(2)%G": 0,
        }
    )
)

