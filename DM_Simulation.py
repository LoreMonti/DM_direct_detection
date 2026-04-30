# ==========================================================
# Dark Matter Direct Detection Simulation
#
# Simulation of the expected event rate in direct detection
# experiments, comparing Xenon and NaI targets for a 70 GeV WIMP.
#
# Includes:
# - Helm form factor
# - Standard Halo Model (SHM)
# - Annual modulation of the signal
#
# Author: Lorenzo Monti
# ==========================================================


import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.special import erf


# =============================================================================
# Constants and parameters for conversion
# =============================================================================
# Constants for converting from fm to GeV^-1 (in natural units)
fm_to_GeV = 5.068
GeVinv_to_fm = 1.0 / fm_to_GeV

# Energy conversion: 1 keV = 1e-6 GeV
keV_to_GeV = 1.0e-6
# Conversion from GeV to kg (1 GeV = 1.783e-27 kg)
GeV_to_kg = 1.783e-27

# Nucleon mass (approximate proton/neutron mass ~ 0.9315 GeV)
m_n = 0.9315

# Local dark matter density (typical value)
rho_chi = 0.3  # GeV/cm^3

# Days in a year
days_per_year = 365.0
days_to_seconds = 86400

# Speed of light in km/s
c_km_s = 3.0e5

# Mass of detector in kg
M_det = 1.0


# =============================================================================
# Problem Data
# =============================================================================
# WIMP mass (GeV)
m_chi = 70.0

# Target parameters (mass numbers)
A_Xe = 132  # Xenon
A_NaI = 150  # NaI

# Spin-independent cross section (cm^2)
sigma_SI = 1e-41

# Set Energy
E_range = np.linspace(0.1, 60, 100)    # (avoid zero energy to prevent division by zero)
E_range2 = np.linspace(10, 60, 100)


# =============================================================================
# Function: Form Factor
# Function to compute the Helm form factor F_SI^2(q)
# where q is the momentum transfer in GeV (natural units)
# =============================================================================
def form_factor(E_keV, A):
    """
    Calculate the Helm form factor F_SI^2(q) using Helm's parametrization.

    Parameters:
      E_keV : Nuclear recoil energy in keV
      A     : Mass number of the nucleus

    Returns:
      F_SI^2 : The square of the spin-independent form factor
    """
    # Convert energy from keV to GeV
    E_GeV = E_keV * keV_to_GeV

    # Nuclear mass (approximated as A * nucleon mass)
    m_N = A * m_n

    # If energy is zero or negative, return the conventional limit (q -> 0)
    if E_GeV <= 0:
        return 1.0

    # Compute the momentum transfer: q^2 = 2 * m_N * E (with c=1)
    q = np.sqrt(2.0 * m_N * E_GeV)

    # Helm parameters (in fm)
    s = 1.0   # Surface thickness in fm
    Rn = 1.2 * A**(1.0/3.0)  # Nuclear radius in fm

    # Convert s and Rn to GeV^-1
    s_GeV = s * fm_to_GeV
    Rn_GeV = Rn * fm_to_GeV

    # Calculate R1
    R1 = np.sqrt(Rn_GeV**2 - 5*(s_GeV**2))

    # Define q * R1
    qR1 = q * R1

    # Define the spherical Bessel function of order 1: j1(x)
    def j1(x):
        return np.sin(x) / x**2 - np.cos(x) / x

    # Avoid numerical issues for very small qR1
    if qR1 < 1e-8:
        prefactor = 1.0  # Limit as q -> 0
    else:
        prefactor = 3.0 * j1(qR1) / qR1

    # Compute the form factor squared including the exponential term
    FF2 = (prefactor**2) * np.exp(- (q**2) * (s_GeV**2))
    return FF2


# Calculate the form factor for each energy in the defined range
FF2_Xe = np.array([form_factor(E, A_Xe) for E in E_range])
FF2_NaI = np.array([form_factor(E, A_NaI) for E in E_range])

# Plotting the form factor vs recoil energy
plt.figure(figsize=(8, 6))
plt.plot(E_range, FF2_Xe, label="Xenon", color="black")
plt.plot(E_range, FF2_NaI, label="NaI", color="red", linestyle="--")
plt.xlabel("Recoil Energy (keV)")
plt.ylabel(r"$F_{SI}^2$")
plt.yscale("log")
plt.title("Helm Form Factor for Xe and NaI")
plt.legend()
plt.grid(True)

filename = "Images/Helm_Form_Factor.pdf"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()

print("----------------------------------------------------------------------")


# =============================================================================
# Earth's Orbital Speed Function
# =============================================================================
def v_earth(t):
    """
    Calculate the Earth's orbital speed using a complex model in galactic coordinates.

    Parameter:
      t : Time (in days)

    Returns:
      Earth's speed (km/s)
    """
    v_orb = 29.8  # Earth's orbital speed in km/s
    omega = 2 * np.pi / 365  # Angular frequency in days^-1 (approximately 0.0172 d^-1)

    # Components of Earth's velocity (in km/s)
    ve_x = 11.1 + v_orb * (0.9941 * np.cos(omega * t) - 0.00504 * np.sin(omega * t))
    ve_y = 238 + 12.2 + v_orb * (0.1088 * np.cos(omega * t) + 0.4946 * np.sin(omega * t))
    ve_z = 7.3 + v_orb * (0.0042 * np.cos(omega * t) - 0.8677 * np.sin(omega * t))
    return np.sqrt(ve_x**2 + ve_y**2 + ve_z**2)


# Define time grid and calculate Earth's speed
t_values = np.linspace(0, 365, 300)
v_e = np.array([v_earth(t) for t in t_values])

# Determine the time and value at which the Earth's speed is maximum
t_max = t_values[np.argmax(v_e)]
max_v_e = np.max(v_e)
print("The maximum Earth's speed occurs at t =", t_max, "days, with v_e =", max_v_e, "km/s")

# Plot Earth's orbital speed using the complex model
plt.figure(figsize=(8,6))
plt.plot(t_values, v_e, color="blue", label=r"$v_e$")
plt.xlabel("Time (days)")
plt.ylabel(r"$v_e$ [km/s]")
plt.title("Earth's Orbital Speed (Complex Model)")
plt.legend()
plt.grid(True)

filename = "Images/Earth_Orbital_Speed.pdf"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()

print("----------------------------------------------------------------------")


# =============================================================================
# Analytic η(E, t, A, m_chi) for the Standard Halo Model
# =============================================================================
def eta(E_keV, t, A, m_chi):
    """
    Calculate the mean inverse speed η(E, t, A, m_chi) for the Standard Halo Model.

    Parameters:
      E_keV : Recoil energy (keV)
      t     : Time (days)
      A     : Mass number of the target nucleus
      m_chi : WIMP mass (GeV)

    Returns:
      η(E, t) in (km/s)^-1
    """
    # SHM model parameters
    v0 = 238.0      # km/s
    v_esc = 544.0   # km/s

    m_N = A * m_n
    E_GeV = E_keV * keV_to_GeV
    mu = (m_chi * m_N) / (m_chi + m_N)

    # Calculate v_min (in km/s)
    v_min = np.sqrt((E_GeV * m_N) / (2 * mu**2)) * c_km_s
    x = v_min / v0
    y = v_earth(t) / v0
    z = v_esc / v0

    N_val = erf(z) - (2 * z / np.sqrt(np.pi)) * np.exp(-z**2)
    if x < 0:
        return 0.0
    elif x <= (z - y):
        val = erf(x + y) - erf(x - y) - (4 * y / np.sqrt(np.pi)) * np.exp(-z**2)
        return (1 / (2 * y * v0)) * (1 / N_val) * val
    elif (z - y) < x <= (z + y):
        val = erf(z) - erf(x - y) - (2 / np.sqrt(np.pi)) * (z + y - x) * np.exp(-z**2)
        return (1 / (2 * y * v0)) * (1 / N_val) * val
    else:
        return 0.0


# Plot η as a function of recoil energy (for a fixed time)
t_fixed = t_max
eta_values_Xe = np.array([eta(E, t_fixed, A_Xe, m_chi) for E in E_range])
eta_values_NaI = np.array([eta(E, t_fixed, A_NaI, m_chi) for E in E_range])

eta_max_Xe = np.max(eta_values_Xe)
eta_max_NaI = np.max(eta_values_NaI)
print("The maximum η for Xe is =", eta_max_Xe, "(km/s)^-1")
print("The maximum η for NaI is =", eta_max_NaI, "(km/s)^-1")

plt.figure(figsize=(8,6))
plt.plot(E_range, eta_values_Xe, label="Xenon", color="black")
plt.plot(E_range, eta_values_NaI, label="NaI", color="red", linestyle="--")
plt.xlabel("Recoil Energy (keV)")
plt.ylabel(r"$\eta(E, t_{max})$ [$(km/s)^{-1}$]")
plt.title("Mean Inverse Speed vs Recoil Energy")
plt.legend()
plt.grid(True)

filename = "Images/Mean_Inverse_Speed_vs_Recoil_Energy.pdf"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()

print("----------------------------------------------------------------------")


# =============================================================================
# Differential Rate and Total Rate Calculation
# =============================================================================
def differential_rate(E_keV, t, A, m_chi, sigma_SI):
    """
    Calculate the differential rate dR/dE (in counts/keV/kg/d) for a given nucleus.

    Parameters:
      E_keV    : Recoil energy (keV)
      t        : Time (days)
      A        : Mass number of the nucleus
      m_chi    : WIMP mass (GeV)
      sigma_SI : Spin-independent cross section (cm^2)

    Returns:
      Differential rate dR/dE
    """
    # Compute the nucleus mass in GeV and kg
    m_N_GeV = A * m_n
    m_N_kg = m_N_GeV * GeV_to_kg

    # Convert energy from keV to GeV
    E_GeV = E_keV * keV_to_GeV
    mu = (m_chi * m_n) / (m_chi + m_n)

    # v_min: minimum velocity (km/s)
    v_min = np.sqrt((E_GeV * m_N_GeV) / (2 * mu**2))

    FF2 = form_factor(E_keV, A)

    # Number of nuclei per kg: N_T = 1 / (m_N in kg)
    N_T = 1.0 / m_N_kg

    # Conversion factor:
    # - dE: from GeV to keV --> 1/keV_to_GeV = 1e6
    # - Time: from seconds to days --> days_to_seconds = 86400
    conv_factor = days_to_seconds / keV_to_GeV  # 86400/1e-6 = 8.64e10

    # Corrected global constant
    C_new = N_T * rho_chi / (2 * m_chi * mu**2) * conv_factor * 1.e6

    return C_new * (A**2) * sigma_SI * FF2 * eta(E_keV, t, A, m_chi)


def total_rate(E_min, E_max, t, A, m_chi):
    """
    Integrate the differential rate over the energy range [E_min, E_max] to obtain the total rate R(t).

    Parameters:
      E_min : Minimum energy (keV)
      E_max : Maximum energy (keV)
      t     : Time (days)
      A     : Mass number of the nucleus
      m_chi : WIMP mass (GeV)

    Returns:
      Total rate R(t) in counts/kg/day
    """
    # Create an energy grid for integration
    energies = np.linspace(E_min, E_max, 100)
    dR = np.array([differential_rate(E, t, A, m_chi, sigma_SI) for E in energies])
    return np.trapezoid(dR, energies)


# Compute total rate R(t) for the four target combinations and energy ranges
R_Xe_0_60  = np.array([total_rate(0, 60, t, A_Xe, m_chi) for t in t_values])
R_Xe_10_60 = np.array([total_rate(10, 60, t, A_Xe, m_chi) for t in t_values])
R_NaI_0_60 = np.array([total_rate(0, 60, t, A_NaI, m_chi) for t in t_values])
R_NaI_10_60 = np.array([total_rate(10, 60, t, A_NaI, m_chi) for t in t_values])


# =============================================================================
# Analysis: Compute R_max, R_min, R_avg, R_mod, and t_max for each case
# =============================================================================
def analyze_rate(R, t_values):
    R_max = np.max(R)
    R_min = np.min(R)
    R_avg = np.mean(R)
    R_mod = 0.5 * (R_max - R_min)
    t_max = t_values[np.argmax(R)]
    return R_max, R_min, R_avg, R_mod, t_max

results_Xe_0_60  = analyze_rate(R_Xe_0_60, t_values)
results_Xe_10_60 = analyze_rate(R_Xe_10_60, t_values)
results_NaI_0_60 = analyze_rate(R_NaI_0_60, t_values)
results_NaI_10_60 = analyze_rate(R_NaI_10_60, t_values)

# Print the results (including R_mod)
print("Xenon [0-60] keV: R_max = {:.3e}, R_min = {:.3e}, R_avg = {:.3e}, R_mod = {:.3e}, t_max = {:.1f} days"
      .format(*results_Xe_0_60))
print("Xenon [10-60] keV: R_max = {:.3e}, R_min = {:.3e}, R_avg = {:.3e}, R_mod = {:.3e}, t_max = {:.1f} days"
      .format(*results_Xe_10_60))
print()
print("NaI [0-60] keV: R_max = {:.3e}, R_min = {:.3e}, R_avg = {:.3e}, R_mod = {:.3e}, t_max = {:.1f} days"
      .format(*results_NaI_0_60))
print("NaI [10-60] keV: R_max = {:.3e}, R_min = {:.3e}, R_avg = {:.3e}, R_mod = {:.3e}, t_max = {:.1f} days"
      .format(*results_NaI_10_60))
print()


# =============================================================================
# Approximation of the rate: R_approx(t) = R0 + R_mod*cos(omega*(t-t_max))
# with omega = 2*pi/365
# =============================================================================
def R_approx(t, R0, R_mod, t_max):
    return R0 + R_mod * np.cos((2*np.pi/365) * (t - t_max))

# --- Approximation for Xenon [0-60] keV ---
R0_Xe, R_mod_Xe, t_max_Xe = results_Xe_0_60[2], results_Xe_0_60[3], results_Xe_0_60[4]
R_approx_vals_Xe = R_approx(t_values, R0_Xe, R_mod_Xe, t_max_Xe)

R_tmax_Xe = total_rate(0, 60, t_max_Xe, A_Xe, m_chi)
R_approx_tmax_Xe = R_approx(t_max_Xe, R0_Xe, R_mod_Xe, t_max_Xe)
error_percent_Xe = abs(R_tmax_Xe - R_approx_tmax_Xe) / R_tmax_Xe * 100
print("Percentage error in the approximation for Xe at t_max: {:.5f}%".format(error_percent_Xe))

# --- Approximation for NaI [0-60] keV ---
R0_NaI, R_mod_NaI, t_max_NaI = results_NaI_0_60[2], results_NaI_0_60[3], results_NaI_0_60[4]
R_approx_vals_NaI = R_approx(t_values, R0_NaI, R_mod_NaI, t_max_NaI)

R_tmax_NaI = total_rate(0, 60, t_max_NaI, A_NaI, m_chi)
R_approx_tmax_NaI = R_approx(t_max_NaI, R0_NaI, R_mod_NaI, t_max_NaI)
error_percent_NaI = abs(R_tmax_NaI - R_approx_tmax_NaI) / R_tmax_NaI * 100
print("Percentage error in the approximation for NaI at t_max: {:.5f}%".format(error_percent_NaI))
print()


# =============================================================================
# Plot: Total Rate vs Time for the different target combinations and sinusoidal approximations
# =============================================================================
plt.figure(figsize=(10, 6))
plt.plot(t_values, R_Xe_0_60, label='Xe [0-60] keV', color='black')
plt.plot(t_values, R_NaI_0_60, label='NaI [0-60] keV', color='red')
plt.plot(t_values, R_Xe_10_60, label='Xe [10-60] keV', color='black')
plt.plot(t_values, R_NaI_10_60, label='NaI [10-60] keV', color='red')

plt.xlabel("Time (days)")
plt.ylabel("Total Rate (counts/kg/d)")
plt.title("Total Rate vs Time for Different Targets and Energy Ranges")
plt.legend()
plt.grid(True)

filename = "Images/Total_Rate.pdf"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()

#-----------------------------------------------------------------------
# With approx. and energy [0-60]
plt.figure(figsize=(10, 6))
plt.plot(t_values, R_Xe_0_60, label='Xe [0-60] keV', color='black')
plt.plot(t_values, R_NaI_0_60, label='NaI [0-60] keV', color='red')
plt.plot(t_values, R_approx_vals_Xe, label='Approx. (Xe [0-60] keV)', color='magenta', linestyle="--")
plt.plot(t_values, R_approx_vals_NaI, label='Approx. (NaI [0-60] keV)', color='green', linestyle="--")

plt.xlabel("Time (days)")
plt.ylabel("Total Rate (counts/kg/d)")
plt.title("Total Rate vs Time for Different Targets and Energy [0-60] with approx.")
plt.legend()
plt.grid(True)

filename = "Images/Total_Rate_[0-60].pdf"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()

#-----------------------------------------------------------------------
# With energy [10-60]
plt.figure(figsize=(10, 6))
plt.plot(t_values, R_Xe_10_60, label='Xe [0-60] keV', color='black')
plt.plot(t_values, R_NaI_10_60, label='NaI [0-60] keV', color='red')

plt.xlabel("Time (days)")
plt.ylabel("Total Rate (counts/kg/d)")
plt.title("Total Rate vs Time for Different Targets and Energy [10-60]")
plt.legend()
plt.grid(True)

filename = "Images/Total_Rate_[10-60].pdf"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.close()