# particle_in_waves
This is dataset from "Phase-resolved settling velocity change of inertial particles in surface gravity waves", which has been submitted to the Journal of Fluid Mechanics (JFM).

## Data file descriptions
### 1. Input CSV files
Particle and wave information used in the experiments are stored in `particle_stag_settling.CSV` and `wavInfo.CSV`.

#### 1-1. `particle_stag_settling.CSV`
- Column 1 (`Case`): particle case ID
- Column 2 (`v_s`): terminal settling velocity measured in quiescent water [mm/s]
- Column 3 (`std`): standard deviation of `v_s` estimated via bootstrapping [mm/s]
- Column 4 (`d`): particle diameter, \(d_p\) [μm]
- Column 5 (`Re`): particle Reynolds number based on the terminal settling velocity, \(Re_{p,t}\)
- Column 6 (`rho_p`): particle density [g/cm³]

#### 1-2. `waveInfo.CSV`
Wave information is provided for each wave condition:
- `H`: wave height [m]
- `T`: wave period [s]
- `omega`: angular frequency [rad/s]
- `k`: wavenumber [rad/m]
- `h`: water depth [m]
- `a`: wave amplitude [m]
- `c`: phase speed [m/s]

---

### 2. Trial-wise slip velocity results (`slip_vs_results_*.mat`)
In each wave–particle case folder, slip-velocity measurement results for each trial are saved as `slip_vs_results_*.mat`. Variables:
- `meta`: logs from the original PTV files and the quiescent-water settling velocity
- `time`: time indexing consistent with PIV
- `x_raw_pix`, `y_raw_pix`: particle positions in image coordinates [px]
- `x_fit`, `y_fit`: spline-fitting results [m]
- `u_p`, `v_p`: particle velocities in the x and z directions [m/s]
- `u_f`, `v_f`: fluid velocities at the particle location in the x and z directions [m/s]
- `u_slip`, `v_slip`: slip (relative) velocities in the x and z directions \((u_p-u_f,\ v_p-v_f)\) [m/s]
- `u_phi_f`, `v_phi_f`: wave phase estimated from the fluid velocity in the x and z directions
- `u_phi_p`, `v_phi_p`: particle-motion phase estimated from the particle velocity in the x and z directions

---

### 3. Ensemble-averaged results (`*_vs_avg_unit.mat`)
In each case folder, ensemble-averaged results are stored in `*_vs_avg_unit.mat`:
- `E_cosh_tracks`, `E_sinh_tracks`: depth factors as a function of \(z_{\mathrm{eff}}\) for each trial
- `E_cosh_ens_out`, `E_sinh_ens_out`: ensemble averages of the depth factors
- `T`: wave period for the corresponding wave case [s]
- `h`: water depth for the corresponding wave case [m]
- `k`: wavenumber for the corresponding wave case [rad/m]
- `k_y_centers`: reference depths for depth-wise binning [m]
- `kz_w_bin_z_mean`, `kz_w_bin_z_std`: wave-averaged particle position \((\overline{z_p})\) [-] and its standard deviation from ensemble averaging [-]
- `kz_w_bin_v_mean`, `kz_w_bin_v_std`: wave-averaged vertical slip velocity [m/s] and its standard deviation from ensemble averaging [m/s]
- `r_avg_u`, `r_se_u`: phase-averaged horizontal slip velocity [m/s] and its standard deviation from ensemble averaging [m/s]
- `r_avg_v`, `r_se_v`: phase-averaged vertical slip velocity [m/s] and its standard deviation from ensemble averaging [m/s]
- `thetas_std`: standardized wave phase \(\phi\)
- `v_s`: measured terminal settling velocity in quiescent water [m/s]

---

### 4. `PiW_figures.m`
Figure-drawing code for Figures 4, 5, 6, 8, and 9.

### 5. `PiW_phase_comp.m`
Code for comparing theory and measurements (Figures 10 and 11).
```
