clear all; close all; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};
parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};

wN = numel(waves); 
pN = numel(parts);

% Wave-averaged (representative track) values
kz_bin_z_mean = cell(wN,pN);   % <kz> for each bin
kz_bin_z_std  = cell(wN,pN);   % std(kz)
kz_bin_v_mean = cell(wN,pN);   % <v>  [m/s]
kz_bin_v_std  = cell(wN,pN);   % std(v) [m/s]

phi_std_tracks     = cell(wN,pN);
phi_std_u_tracks   = cell(wN,pN);
phi_std_v_tracks   = cell(wN,pN);
phi_std_orb_tracks = cell(wN,pN);
dphi_uv_tracks     = cell(wN,pN);
dphi_orb_u_tracks  = cell(wN,pN);

% Std-phase arrays
r_avg_u = cell(wN,pN); 
r_avg_v = cell(wN,pN);
r_se_u  = cell(wN,pN); 
r_se_v  = cell(wN,pN);
thetas_std = cell(wN,pN);      % phi_std bins [rad]

% Scalars per case
v_s = zeros(wN,pN); 
k   = v_s; 
T   = v_s; 
h   = v_s; 
Ec  = v_s; 
Es  = v_s;

nu    = waterNu(24.5);
rho_w = waterRho(24.5);
g     = 9.81;

pInfo = readtable('particle_stag_settling.CSV');

wave_cond = readtable('waveInfo.CSV');
H0 = wave_cond.H;

% Wave height matrix (wN x pN)
H = repmat(H0(:), 1, pN);

% --- Load per-case results generated with common orbital phase ---
avg_suffix = 'vh_avg_unit';

for wI = 1:wN
  wave = waves{wI};

  for pI = 1:pN
    partinfo = parts{pI};

    avgfile = fullfile([wave partinfo], ...
        sprintf('%s%s_%s.mat', wave, partinfo, avg_suffix));

    if ~isfile(avgfile)
        error('Missing orbital-phase averaged file: %s', avgfile);
    end

    S = load(avgfile);

    if ~isfield(S,'phi_std_orb_tracks')
        error('File %s does not contain phi_std_orb_tracks. Re-run avg generation with orb_phi_f.', avgfile);
    end

    % Phase bins and slip statistics
    thetas_std{wI,pI} = S.thetas_std;     % common orbital phase bins [rad]
    r_avg_u{wI,pI}    = S.r_avg_u;        % [m/s]
    r_avg_v{wI,pI}    = S.r_avg_v;        % [m/s]
    r_se_u{wI,pI}     = S.r_se_u;         % [m/s]
    r_se_v{wI,pI}     = S.r_se_v;         % [m/s]

    % Optional diagnostics
    phi_std_tracks{wI,pI}     = S.phi_std_tracks;
    phi_std_u_tracks{wI,pI}   = S.phi_std_u_tracks;
    phi_std_v_tracks{wI,pI}   = S.phi_std_v_tracks;
    phi_std_orb_tracks{wI,pI} = S.phi_std_orb_tracks;
    dphi_uv_tracks{wI,pI}     = S.dphi_uv_tracks;
    dphi_orb_u_tracks{wI,pI}  = S.dphi_orb_u_tracks;

    % Parameters
    v_s(wI,pI) = S.v_s;
    k(wI,pI)   = S.k;
    h(wI,pI)   = S.h;
    T(wI,pI)   = S.T;
    Ec(wI,pI)  = S.E_cosh_ens_out;
    Es(wI,pI)  = S.E_sinh_ens_out;

    % Wave-averaged profiles
    kz_bin_z_mean{wI,pI} = S.kz_w_bin_z_mean;
    kz_bin_z_std{wI,pI}  = S.kz_w_bin_z_std;
    kz_bin_v_mean{wI,pI} = S.kz_w_bin_v_mean;
    kz_bin_v_std{wI,pI}  = S.kz_w_bin_v_std;
  end
end

%% Particle parameters

d_p   = zeros(wN,pN); 
rho_p = zeros(wN,pN);

for pI = 1:pN
  case_mask = strcmp(pInfo.Case, parts{pI});
  d_p(:,pI)   = pInfo.d(case_mask)*1e-6;         % [m]
  rho_p(:,pI) = pInfo.rho_p(case_mask)*1000;     % [kg/m3]
end

omega  = 2*pi./T;
c      = omega./k;

Re_0   = v_s.*d_p/nu;
SchNau = 1 + 0.15*Re_0.^0.687;

tau_p  = d_p.^2.*(rho_p/rho_w + 0.5)./18./nu./SchNau;
St     = omega.*tau_p;
beta   = 3./(2*rho_p/rho_w + 1);

% Theoretical terminal slip scale (Schiller-Naumann-based)
v_t    = ( (omega./k).*St.*(1-beta) ./ tanh(k.*h) );

K = 0.687*(1-SchNau)./SchNau./(v_s./c)./St;

v_g    = sqrt((rho_p/rho_w-1).*g.*d_p);          % [m/s]
Ga     = v_g.*d_p./nu;                           % Galileo number [-]

epsilon = (k.*H/2)./tanh(k.*h);
A_eff   = Ec./cosh(k.*h);
B_eff   = Es./cosh(k.*h);

%% theory vs. measurements (figures 11, 12)
close all;

% Preallocate outputs used in plots/tables
Rx_theo = nan(wN,pN);  Rz_theo = nan(wN,pN);
Rx_exp  = nan(wN,pN);  Rz_exp  = nan(wN,pN);

phix_peak_exp  = nan(wN,pN); phiz_peak_exp  = nan(wN,pN);
phix_peak_theo = nan(wN,pN); phiz_peak_theo = nan(wN,pN);

theta_x1_exp = nan(wN,pN); theta_z1_exp = nan(wN,pN);
theta_x2_exp = nan(wN,pN); theta_z2_exp = nan(wN,pN);

theta_x1_theo = nan(wN,pN); theta_z1_theo = nan(wN,pN);
theta_x2_theo = nan(wN,pN); theta_z2_theo = nan(wN,pN);

offx_exp  = nan(wN,pN);
offx_theo = nan(wN,pN);
offz_exp  = nan(wN,pN);
offz_theo = nan(wN,pN);

Rx1_exp = nan(wN,pN);  Rx2_exp = nan(wN,pN);
Rz1_exp = nan(wN,pN);  Rz2_exp = nan(wN,pN);

Rx1_theo = nan(wN,pN); Rx2_theo = nan(wN,pN);
Rz1_theo = nan(wN,pN); Rz2_theo = nan(wN,pN);

% full-wave O(eps^2) alignment metrics
delta_x_star = nan(wN,pN);   % lag maximizing waveform similarity
delta_z_star = nan(wN,pN);
rho_x_star   = nan(wN,pN);   % max normalized similarity
rho_z_star   = nan(wN,pN);

% Coefficients saved for table export
pxM  = nan(wN,pN); qxM  = nan(wN,pN); pzM  = nan(wN,pN); qzM  = nan(wN,pN);
a0LM = nan(wN,pN); acLM = nan(wN,pN); asLM = nan(wN,pN);
b0LM = nan(wN,pN); bcLM = nan(wN,pN); bsLM = nan(wN,pN);

pxE = nan(wN,pN); qxE = nan(wN,pN);
pzE = nan(wN,pN); qzE = nan(wN,pN);
acE = nan(wN,pN); asE = nan(wN,pN);
bcE = nan(wN,pN); bsE = nan(wN,pN);

% Plot styling
c_wI = { ...
    [0 0 0]                % W1: black
    [0.698 0.1333 0.1333]  % W2: dark reddish
    [0 0.5451 0.5451]      % W3: teal-ish
    };
ms_wI= {'o','^','s'};
p_name = {'M1','M2','','G1','G2','G3','E1','E2','E3'};
label_name = 'ab-cdefgh';

% fig = figure(1); 
% fig.Position = [80 80 1200 500];

for wI = 1:wN
  wave = waves{wI};
  for pI = 1:pN
    partinfo = parts{pI};

    % fprintf('%s %s\n', wave, partinfo);

    % --- Build perturbation/asymptotic solution (extract required variables) ---
    beta_ = beta(wI,pI);
    epsv  = epsilon(wI,pI); 
    Aef   = A_eff(wI,pI); 
    Bef   = B_eff(wI,pI);

    % Measured settling velocity and nondimensionalization by celerity
    vs_meas = v_s(wI,pI);     % [m/s]
    c0      = c(wI,pI);       % [m/s]
    vs0_    = vs_meas / c0;   % [-]

    D0drag  = SchNau(wI,pI);
    St_SN   = St(wI,pI);
    St_     = St_SN;
    St2     = St_^2;
    mSN     = 0.687;

    K_ = mSN*(1-D0drag)/D0drag/vs0_/St_;

    % --- Data: phi_std, r = v_slip/c and u_slip/c with SE ---
    phi_std = thetas_std{wI,pI}(:);                 % [rad], length = nbin

    % Nondimensionalize by celerity
    rv_mean = r_avg_v{wI,pI}(:) ./ c0;              % mean(-v_slip/c)
    ru_mean = r_avg_u{wI,pI}(:) ./ c0;              % mean(u_slip/c)
    rv_se   = r_se_v{wI,pI}(:)  ./ c0;              % SE
    ru_se   = r_se_u{wI,pI}(:)  ./ c0;              % SE

    % Sort angles
    [phi_plot, ord] = sort(mod(phi_std, 2*pi));
    rv_mean = rv_mean(ord); ru_mean = ru_mean(ord);
    rv_se   = rv_se(ord);   ru_se   = ru_se(ord);

    % --- O(eps^2) solution (both w_x and w_z included) ---
    % 1st-order coefficients:
    % w_x1 = p_x cos(phi) + q_x sin(phi), w_z1 = p_z cos(phi) + q_z sin(phi)
    px = ( -St2*Aef*(1-beta_) + St_*Bef*vs0_ )/(St2+1);
    qx = ( -St_*Aef*(1-beta_) - St2*Bef*vs0_ )/(St2+1);
    pz = (  St2*Aef*vs0_ + St_*(1-beta_)*Bef )/(St2+1);
    qz = (  St_*Aef*vs0_ - St2*(1-beta_)*Bef )/(St2+1);   % corrected sign in A22b

    % 2nd-order forcing coefficients (A26, corrected signs)
    % Forcing convention:
    %   d w_x2/dphi = (1/St_)*w_x2 + Gamma0 + Gamma1*cos(2phi) + Gamma2*sin(2phi)
    %   d w_z2/dphi = (1/St_)*w_z2 + Delta0 + Delta1*cos(2phi) + Delta2*sin(2phi)
    % 2nd-order forcing coefficients (corrected Appendix A)
    Gamma0 = 0.5*( ...
          Bef*pz ...
        + (Aef*px + px^2 + qx^2)/St_ ...
        - beta_*Aef*qx ...
        - vs0_*Bef*(Aef + px) ...
        + K_*(px*pz + qx*qz) );
    
    Gamma1 = 0.5*( ...
          Bef*pz ...
        + (Aef*px + px^2 - qx^2)/St_ ...
        + beta_*Aef*qx ...
        - vs0_*Bef*(Aef + px) ...
        + K_*(px*pz - qx*qz) );
    
    Gamma2 = 0.5*( ...
          Bef*qz ...
        + (Aef*qx + 2*px*qx)/St_ ...
        - beta_*Aef*px ...
        + (1-beta_)*Bef^2 ...
        - vs0_*Bef*qx ...
        + K_*(px*qz + qx*pz) );
    
    Delta0 = 0.5*( ...
          Aef*qz ...
        + (Aef*pz + px*pz + qx*qz)/St_ ...
        + (1-beta_)*Aef*Bef ...
        + beta_*Bef*px ...
        - vs0_*Aef*qx ...
        + K_*(pz^2 + qz^2) );
    
    Delta1 = 0.5*( ...
        - Aef*qz ...
        + (Aef*pz + px*pz - qx*qz)/St_ ...
        - (1-beta_)*Aef*Bef ...
        + beta_*Bef*px ...
        + vs0_*Aef*qx ...
        + K_*(pz^2 - qz^2) );
    
    Delta2 = 0.5*( ...
          Aef*pz ...
        + (Aef*qz + px*qz + qx*pz)/St_ ...
        + beta_*Bef*qx ...
        - vs0_*Aef*(Aef + px) ...
        + 2*K_*pz*qz );
    
    a0L = -St_ * Gamma0;
    acL = (-Gamma1*St_ - 2*Gamma2*St2) / (1 + 4*St2);
    asL = (-Gamma2*St_ + 2*Gamma1*St2) / (1 + 4*St2);
    
    b0L = -St_ * Delta0;
    bcL = (-Delta1*St_ - 2*Delta2*St2) / (1 + 4*St2);
    bsL = (-Delta2*St_ + 2*Delta1*St2) / (1 + 4*St2);

    wx_lin2 =  epsv*( px*cos(phi_plot) + qx*sin(phi_plot) ) ...
            + epsv^2*( a0L + acL*cos(2*phi_plot) + asL*sin(2*phi_plot) );

    % Full w_z ≈ -vs_ + eps*w_z1 + eps^2*w_z2  -> r = -w_z
    wz_lin2 = -vs0_ + epsv*( pz*cos(phi_plot) + qz*sin(phi_plot) ) ...
                    + epsv^2*( b0L + bcL*cos(2*phi_plot) + bsL*sin(2*phi_plot) );
    r_lin2  = -wz_lin2;

    % Valid mask
    mexp = isfinite(phi_plot) & isfinite(rv_mean);
    x  = phi_plot(mexp); 
    yv = rv_mean(mexp); 
    yu = ru_mean(mexp);
    
    % Experimental/theoretical mean offsets for Fig. 11
    offx_exp(wI,pI)  = mean(yu, 'omitnan');
    offx_theo(wI,pI) = mean(wx_lin2, 'omitnan');
    
    offz_exp(wI,pI)  = mean(yv, 'omitnan');
    offz_theo(wI,pI) = mean(r_lin2, 'omitnan');

    % =========================================================
    % Vertical fit: y = Q + eps*R1*sin(phi+th1) + eps^2*R2*sin(2phi+th2)
    % Fit unknowns: [R1, th1, R2, th2] (Q is fixed = mean(y))
    % =========================================================
    xv   = x(:);
    yv_c = yv(:);
    
    sev  = rv_se(mexp);
    wv   = 1./max(sev(:).^2, eps);
    W12v = sqrt(wv(:));
    
    Qexp_mean     = mean(yv_c);
    
    % detrend for 4-parameter fit
    yv0 = yv_c - Qexp_mean;
    
    % initial guess from linear 1st/2nd harmonic fit (for stability)
    Xlin = [cos(xv), sin(xv), cos(2*xv), sin(2*xv), ones(numel(xv),1)];
    coef0 = (Xlin.*W12v) \ (yv_c.*W12v);
    A1_0 = coef0(1); B1_0 = coef0(2);
    A2_0 = coef0(3); B2_0 = coef0(4);
    
    R1_0 = hypot(A1_0,B1_0) / max(epsv, eps);
    th1_0 = atan2(A1_0, B1_0);
    R2_0 = hypot(A2_0,B2_0) / max(epsv^2, eps);
    th2_0 = atan2(A2_0, B2_0);
    
    % --- reparameterization: R2 = alpha * R1, with 0<=alpha<=1 ---
    alpha0 = R2_0 / max(R1_0, eps);        % initial alpha
    alpha0 = min(max(alpha0, 0), 1);
    
    % decision variable q = [R1; th1; alpha; th2]
    q0 = [R1_0; th1_0; alpha0; th2_0];
    lb = [0;   -pi; 0;   -pi];
    ub = [Inf;  pi; 1;    pi];
    
    resV = @(q) ( ...
        epsv*q(1)*sin(xv + q(2)) + ...
        (epsv^2)*(q(3)*q(1))*sin(2*xv + q(4)) - ...
        yv0 ) .* W12v;
    
    opts = optimoptions('lsqnonlin','Display','off', ...
        'MaxIterations',300,'MaxFunctionEvaluations',5000);
    
    qV = lsqnonlin(resV, q0, lb, ub, opts);
    
    Rz1  = qV(1);
    thz1 = qV(2);
    alpha= qV(3);
    thz2 = qV(4);
    
    Rz2 = alpha * Rz1;   % <-- enforced: 0 <= Rz2 <= Rz1
    
    Rz1_exp(wI,pI) = Rz1;
    Rz2_exp(wI,pI) = Rz2;

    % convert back to A/B (keep downstream compatible)
    A1v = epsv * Rz1 * sin(thz1);
    B1v = epsv * Rz1 * cos(thz1);
    A2v = epsv^2 * Rz2 * sin(thz2);
    B2v = epsv^2 * Rz2 * cos(thz2);
    
    rv_fit_exp = Qexp_mean ...
               + epsv  * Rz1 * sin(phi_plot + thz1) ...
               + epsv^2* Rz2 * sin(2*phi_plot + thz2);

    fprintf('%s%s: epsilon=%.4f / (R1,R2)=(%.4g, %.4g) / (th1,th2)=(%.3f, %.3f)\n', ...
            wave, partinfo, epsv, Rz1, Rz2, thz1, thz2);
    
    
    % =========================================================
    % Horizontal fit: y = Q + eps*R1*sin(phi+th1) + eps^2*R2*sin(2phi+th2)
    % Fit unknowns: [R1, th1, R2, th2] (Q fixed = mean(y))
    % =========================================================
    xu   = x(:);
    yu_c = yu(:);
    
    seu  = ru_se(mexp);
    wu   = 1./max(seu(:).^2, eps);
    W12u = sqrt(wu(:));
    
    Qu_mean = mean(yu_c);
    yu0 = yu_c - Qu_mean;
    
    Xlin = [cos(xu), sin(xu), cos(2*xu), sin(2*xu), ones(numel(xu),1)];
    coef0 = (Xlin.*W12u) \ (yu_c.*W12u);
    A1_0 = coef0(1); B1_0 = coef0(2);
    A2_0 = coef0(3); B2_0 = coef0(4);
    
    R1_0 = hypot(A1_0,B1_0) / max(epsv, eps);
    th1_0 = atan2(A1_0, B1_0);
    R2_0 = hypot(A2_0,B2_0) / max(epsv^2, eps);
    th2_0 = atan2(A2_0, B2_0);
    
    p0 = [R1_0; th1_0; R2_0; th2_0];
    
    resU = @(p) ( epsv*p(1)*sin(xu + p(2)) + (epsv^2)*p(3)*sin(2*xu + p(4)) - yu0 ) .* W12u;
    
    pU = lsqnonlin(resU, p0, lb, ub, opts);
    
    Rx1 = pU(1);  thx1 = pU(2);
    Rx2 = pU(3);  thx2 = pU(4);
    
    Rx1_exp(wI,pI) = Rx1;
    Rx2_exp(wI,pI) = Rx2;
    
    A1u = epsv * Rx1 * sin(thx1);
    B1u = epsv * Rx1 * cos(thx1);
    A2u = epsv^2 * Rx2 * sin(thx2);
    B2u = epsv^2 * Rx2 * cos(thx2);
    
    Qu = Qu_mean;
    
    ru_fit_exp = Qu_mean ...
               + epsv * Rx1 * sin(phi_plot + thx1) ...
               + epsv^2 * Rx2 * sin(2*phi_plot + thx2);
    

    % --- Experimental peaks/amplitudes ---
    [xmax_exp,imax_x] = max(ru_fit_exp);
    [xmin_exp,~]      = min(ru_fit_exp);
    [ymax_exp,imax_y] = max(rv_fit_exp);
    [ymin_exp,~]      = min(rv_fit_exp);

    Rx_exp(wI,pI) = xmax_exp - xmin_exp;
    Rz_exp(wI,pI) = ymax_exp - ymin_exp;

    Rx1_exp(wI,pI) = Rx1;
    Rx2_exp(wI,pI) = Rx2;
    Rz1_exp(wI,pI) = Rz1;
    Rz2_exp(wI,pI) = Rz2;

    theta_x1_exp(wI,pI) = atan2(A1u,B1u);
    theta_z1_exp(wI,pI) = atan2(A1v,B1v);
    theta_x2_exp(wI,pI) = atan2(A2u,B2u);
    theta_z2_exp(wI,pI) = atan2(A2v,B2v);

    phix_peak_exp(wI,pI) = phi_plot(imax_x);
    phiz_peak_exp(wI,pI) = phi_plot(imax_y);

    offx_exp(wI,pI) = Qu_mean;
    offz_exp(wI,pI) = Qexp_mean;

    % --- Theoretical peaks/amplitudes from perturbation solution ---
    [xmax_theo,imax_x] = max(wx_lin2);
    [xmin_theo,~]      = min(wx_lin2);

    [ymax_theo,imax_y] = max(r_lin2);
    [ymin_theo,~]      = min(r_lin2);

    Rx_theo(wI,pI) = xmax_theo - xmin_theo;
    Rz_theo(wI,pI) = ymax_theo - ymin_theo;

    Rx1_theo(wI,pI) = hypot(px,qx);
    Rx2_theo(wI,pI) = hypot(acL,asL);
    Rz1_theo(wI,pI) = hypot(pz,qz);
    Rz2_theo(wI,pI) = hypot(bcL,bsL);

    % These phases are for w_x and w_z as defined in the theory
    theta_x1_theo(wI,pI) = atan2(px,qx);
    theta_z1_theo(wI,pI) = atan2(pz,qz);
    theta_x2_theo(wI,pI) = atan2(acL,asL);
    theta_z2_theo(wI,pI) = atan2(bcL,bsL);

    phix_peak_theo(wI,pI) = phi_plot(imax_x);
    phiz_peak_theo(wI,pI) = phi_plot(imax_y);

    offx_theo(wI,pI) = mean(wx_lin2,'omitnan');
    offz_theo(wI,pI) = mean(r_lin2,'omitnan');

    % --- full-wave O(eps^2) alignment metrics ---
    % x-direction: direct comparison in w_x convention
    [delta_x_star(wI,pI), rho_x_star(wI,pI)] = wave_align_metric( ...
        epsv, ...
        Rx1_exp(wI,pI), theta_x1_exp(wI,pI), Rx2_exp(wI,pI), theta_x2_exp(wI,pI), ...
        Rx1_theo(wI,pI), theta_x1_theo(wI,pI), Rx2_theo(wI,pI), theta_x2_theo(wI,pI) );

    % z-direction: experimental comparison is based on r = -w_z
    % so the theoretical phases must be converted from w_z to -w_z convention
    theta_z1_theo_r = atan2(-pz,  -qz);
    theta_z2_theo_r = atan2(-bcL, -bsL);

    [delta_z_star(wI,pI), rho_z_star(wI,pI)] = wave_align_metric( ...
        epsv, ...
        Rz1_exp(wI,pI), theta_z1_exp(wI,pI), Rz2_exp(wI,pI), theta_z2_exp(wI,pI), ...
        Rz1_theo(wI,pI), theta_z1_theo_r,    Rz2_theo(wI,pI), theta_z2_theo_r );

    % Save coefficients for table export
    pxM(wI,pI)  = px;   qxM(wI,pI)  = qx;   pzM(wI,pI)  = pz;   qzM(wI,pI)  = qz;
    a0LM(wI,pI) = a0L;  acLM(wI,pI) = acL;  asLM(wI,pI) = asL;
    b0LM(wI,pI) = b0L;  bcLM(wI,pI) = bcL;  bsLM(wI,pI) = bsL;

    pxE(wI,pI) = A1u; qxE(wI,pI) = B1u;
    pzE(wI,pI) = A1v; qzE(wI,pI) = B1v;
    acE(wI,pI) = A2u; asE(wI,pI) = B2u;
    bcE(wI,pI) = A2v; bsE(wI,pI) = B2v;

    % --- Plot update: data ± SE, experimental fit, theory (perturbation) ---
    % In 3x3, slot (3,3,3) is reserved for legend dummy, so map pI=1..8 to 1,2,4,5,6,7,8,9
    spI = pI + (pI>=3);

    % Figure 1: rv_mean and theory curve
    % figure(1)
    % subplot(3,3,spI), hold on; grid on
    % 
    % xfill = [phi_plot; flipud(phi_plot)];
    % yfill = [rv_mean - rv_se; flipud(rv_mean + rv_se)];
    % if all(isfinite(yfill))
    %     f1 = fill(xfill, yfill, [0.88 0.88 0.88], 'EdgeColor','none');
    %     uistack(f1,'bottom');
    % end
    % 
    % plot(phi_plot, rv_mean, '.', 'Color', c_wI{wI}, 'LineWidth',1.6, ...
    %      'DisplayName',sprintf('%s meas.',wave));
    % plot(phi_plot, rv_fit_exp,'--', 'Color', c_wI{wI}, 'LineWidth',1.6, ...
    %      'DisplayName',sprintf('%s fitted',wave));
    % plot(phi_plot, r_lin2, '-.', 'Color', c_wI{wI}, 'LineWidth',1.2, ...
    %      'DisplayName',sprintf('%s asym.',wave));
    % 
    % xlim([0 2*pi]); 
    % xlabel('$\phi$','Interpreter','latex'); 
    % ylabel('$w_z$','Interpreter','latex');
    % ylim([0 0.12])
    % 
    % set(gca,'XTick',[0 pi/2, pi , 1.5*pi, 2*pi], ...
    %         'XTickLabel',{'0','0.5\pi','\pi','1.5\pi','2\pi'})
    % title(partinfo); 
    % set(gca,'FontSize',14)


    % Figure 2: detrended w_z comparison (exp vs theory)
    figure(2)
    subplot(3,3,spI), hold on; grid on

    plot(phi_plot, -rv_mean + Qexp_mean, '-', ...
        'Color',1-(1-c_wI{wI})*0.3,'LineWidth',1.3, ...
        'DisplayName',sprintf('$\\langle w_z \\rangle_{[\\phi]} - \\overline{\\langle w_z \\rangle_{[\\phi]}} \\mathrm{(%s)}$',wave));

    xfill = [phi_plot; flipud(phi_plot)];
    yfill = [-rv_mean + Qexp_mean - rv_se; flipud(-rv_mean + Qexp_mean + rv_se)];
    f1 = fill(xfill, yfill, 1-(1-c_wI{wI})*0.5, 'EdgeColor','none','FaceAlpha',0.5);
    uistack(f1,'bottom')

    plot(phi_plot, -rv_fit_exp + Qexp_mean, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\hat{w_z} - \\overline{\\hat{w_z}} \\, \\mathrm{(%s)}$',wave));
    wz_theo_mod = wz_lin2 - mean(wz_lin2,'omitnan');
    
    plot(phi_plot, wz_theo_mod,'--','Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$w_z + v_s \\, \\mathrm{(Eq.(2.15),\\,%s)}$',wave));

    xlim([0 2*pi]);
    xlabel('$\phi\,\,[\mathrm{rad}]$','Interpreter','latex'); 
    ylabel('$w_z - \overline{w_z}$','Interpreter','latex');
    ylim([-0.005 0.005])
    title(p_name{spI}); 
    set(gca,'FontSize',14)
    set(gca,'XTick',[0 pi/2, pi , 1.5*pi, 2*pi], ...
            'XTickLabel',{'0','0.5\pi','\pi','1.5\pi','2\pi'},'FontSize',16)

    text(-0.18,1.02,sprintf('(%c)',label_name(spI)), 'Units','normalized', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Clipping','off')

    % Figure 3: detrended w_x comparison (exp vs theory)
    figure(3)
    subplot(3,3,spI), hold on; grid on

    xfill_u = [phi_plot; flipud(phi_plot)];
    yfill_u = [ru_mean - mean(yu) - ru_se; flipud(ru_mean - mean(yu) + ru_se)];
    if all(isfinite(yfill_u))
        f2 = fill(xfill_u, yfill_u, 1-(1-c_wI{wI})*0.5, 'EdgeColor','none','FaceAlpha',0.5);
        uistack(f2,'bottom');
    end

    plot(phi_plot, ru_mean - mean(yu), '-', ...
         'Color', 1-(1-c_wI{wI})*0.3, 'LineWidth',1.3, ...
         'DisplayName',sprintf('$\\langle w_x \\rangle_{[\\phi]} - \\overline{\\langle w_x \\rangle_{[\\phi]}} \\mathrm{(%s)}$',wave));
    plot(phi_plot, ru_fit_exp - mean(yu),'-.', 'Color', c_wI{wI}, 'LineWidth',2, ...
         'DisplayName',sprintf('$\\hat{w_x} - \\overline{\\hat{w_x}} \\, \\mathrm{(%s)}$',wave));
    plot(phi_plot, wx_lin2, '--',  'Color',c_wI{wI}, 'LineWidth',2, ...
         'DisplayName',sprintf('$w_x \\, \\mathrm{(Eq.(2.15),\\,%s)}$',wave));

    xlim([0 2*pi]); 
    xlabel('$\phi\,\,[rad]$','Interpreter','latex'); 
    ylabel('$w_x$','Interpreter','latex');
    ylim([-0.005 0.005])
    set(gca,'XTick',[0 pi/2, pi , 1.5*pi, 2*pi], ...
            'XTickLabel',{'0','0.5\pi','\pi','1.5\pi','2\pi'})
    title(p_name{spI}); 
    set(gca,'FontSize',14)

    text(-0.18,1.02,sprintf('(%c)',label_name(spI)), 'Units','normalized', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Clipping','off')

    % Figure 4: nondimensionalize with v_s | w_z comparison (exp vs theory)
    % figure(4)
    % subplot(3,3,spI), hold on; grid on
    % 
    % plot(phi_plot, (-rv_mean + Qexp_mean)*c0/vs0_, '-', ...
    %     'Color',1-(1-c_wI{wI})*0.3,'LineWidth',1.3, ...
    %     'DisplayName',sprintf('$\\langle w_z \\rangle_\\phi - \\overline{\\langle w_z \\rangle_\\phi} \\mathrm{(%s)}$',wave));
    % 
    % xfill = [phi_plot; flipud(phi_plot)];
    % yfill = [(-rv_mean + Qexp_mean - rv_se)*c0/vs0_; flipud((-rv_mean + Qexp_mean + rv_se)*c0/vs0_)];
    % f1 = fill(xfill, yfill, 1-(1-c_wI{wI})*0.5, 'EdgeColor','none','FaceAlpha',0.5);
    % uistack(f1,'bottom')
    % 
    % plot(phi_plot, (-rv_fit_exp + Qexp_mean)*c0/vs0_, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
    %     'DisplayName',sprintf('$\\hat{w_z} - \\overline{\\hat{w_z}} \\, \\mathrm{(%s)}$',wave));
    % wz_theo_mod = wz_lin2 - mean(wz_lin2,'omitnan');
    % 
    % plot(phi_plot, wz_theo_mod*c0/vs0_,'--','Color',c_wI{wI}, 'LineWidth',2, ...
    % 'DisplayName',sprintf('$w_z + v_s \\, \\mathrm{(Eq.(2.15),\\,%s)}$',wave));
    % 
    % xlim([0 2*pi]);
    % xlabel('$\phi\,\,[\mathrm{rad}]$','Interpreter','latex'); 
    % ylabel('$(w_z - \overline{w_z})/v_s$','Interpreter','latex');
    % % ylim([-0.005 0.005])
    % ylim([-0.2 0.2])
    % title(p_name{spI}); 
    % set(gca,'FontSize',14)
    % set(gca,'XTick',[0 pi/2, pi , 1.5*pi, 2*pi], ...
    %         'XTickLabel',{'0','0.5\pi','\pi','1.5\pi','2\pi'},'FontSize',16)
    % 
    % text(-0.18,1.02,sprintf('(%c)',label_name(spI)), 'Units','normalized', ...
    %     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Clipping','off')
  end
end

% Build table and legends after loop

wave_col = cell(wN*pN,1);
part_col = cell(wN*pN,1);
ki = 0;
for wI = 1:wN
  for pI = 1:pN
    ki = ki + 1;
    wave_col{ki} = waves{wI};
    part_col{ki} = parts{pI};
  end
end

% Legend dummy slot: subplot(3,3,3)

% % Figure 1 legend (meas + asym)
% figure(1)
% axL = subplot(3,3,3); cla(axL); axis(axL,'off'); set(axL,'Visible','off'); hold(axL,'on')
% h1 = gobjects(wN,1); h4 = gobjects(wN,1);
% for wI = 1:wN
%     wave = waves{wI};
%     h1(wI) = plot(axL, nan, nan, '.', 'Color', c_wI{wI}, 'LineWidth',1.6, ...
%         'DisplayName',sprintf('%s meas.',wave));
% end
% for wI = 1:wN
%     wave = waves{wI};
%     h4(wI) = plot(axL, nan, nan, '-.', 'Color', c_wI{wI}, 'LineWidth',1.2, ...
%         'DisplayName',sprintf('%s asym.',wave));
% end
% legend(axL, [h1; h4], 'Location','best','Interpreter','latex','FontSize',20);

% Figure 2 legend (pp1 + pp2 + pp3)
fig2 = figure(2); fig2.Position = [-500 100 1000 1000];
axL = subplot(3,3,3); cla(axL); axis(axL,'off'); set(axL,'Visible','off'); hold(axL,'on')
hpp1 = gobjects(wN,1); hpp2 = gobjects(wN,1); hpp3 = gobjects(wN,1);
for wI = 1:wN
    wave = waves{wI};
    hpp1(wI) = plot(axL, nan, nan, '-', 'Color',1-(1-c_wI{wI})*0.3, 'LineWidth',1.3, ...
        'DisplayName',sprintf('$\\langle w_z \\rangle_{[\\phi]} - \\overline{\\langle w_z \\rangle_{[\\phi]}} \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hpp2(wI) = plot(axL, nan, nan, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,\\hat{w_z} - \\overline{\\hat{w_z}} \\, \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hpp3(wI) = plot(axL, nan, nan, '--', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,w_z + v_s \\, \\mathrm{(Eq.(2.18),\\,%s)}$',wave));
end
legend(axL, [hpp1; hpp2; hpp3], 'Location','best','Interpreter','latex','FontSize',16);

% Figure 3 legend (q1 + q3 + q4)
fig3 = figure(3); fig3.Position = [-500 100 1000 1000];
axL = subplot(3,3,3); cla(axL); axis(axL,'off'); set(axL,'Visible','off'); hold(axL,'on')
hq1 = gobjects(wN,1); hq3 = gobjects(wN,1); hq4 = gobjects(wN,1);
for wI = 1:wN
    wave = waves{wI};
    hq1(wI) = plot(axL, nan, nan, '-', 'Color',1-(1-c_wI{wI})*0.3, 'LineWidth',1.3, ...
        'DisplayName',sprintf('$\\langle w_x \\rangle_{[\\phi]} - \\overline{\\langle w_x \\rangle_{[\\phi]}} \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hq3(wI) = plot(axL, nan, nan, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,\\hat{w_x} - \\overline{\\hat{w_x}} \\, \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hq4(wI) = plot(axL, nan, nan, '--', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,w_x \\, \\mathrm{(Eq.(2.18),\\,%s)}$',wave));
end
legend(axL, [hq1; hq3; hq4], 'Location','best','Interpreter','latex','FontSize',16);

% Export table of coefficients
T_coeff_eps2 = table( ...
  wave_col, part_col, ...
  pxM(:), qxM(:), pzM(:), qzM(:), ...
  a0LM(:), acLM(:), asLM(:), b0LM(:), bcLM(:), bsLM(:), ...
  pzE(:), qzE(:), bcE(:), bsE(:), ...
  'VariableNames',{ ...
    'wave','part', ...
    'px','qx','pz','qz', ...
    'a0L','acL','asL','b0L','bcL','bsL','pzE','qzE','bcE','bsE'});

% exportgraphics(fig2,'E:\2.ParticlesInWave\Manuscript\figures\z_slip_phase.pdf')
% exportgraphics(fig3,'E:\2.ParticlesInWave\Manuscript\figures\x_slip_phase.pdf')

%% theory vs. measurements (figure 13): 2x3 with full-wave O(eps^2) lag metric

% Quantitative comparison mask:
% retain all cases, but mark W1-G1/G2/G3 and W1-E1/E2/E3 separately
use_case = true(wN,pN);

% parts = {'M1','M2','G1','G2','G3','E1','E2','E3'}
mark_W1_effect = false(wN,pN);
mark_W1_effect(1,3:8) = true;   % W1-G1/G2/G3 and W1-E1/E2/E3

markers = {'o','^','s'};        % W1, W2, W3
specialMarker = 'd';            % W1-affected cases
cmap_ = turbo(400);

% marker size from rho*
msz0   = 10;
mszAmp = 200;

fig15 = figure(); clf;
fig15.Position = [80 80 1300 780];
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% -----------------------------
% (a) x-direction phase mean (diagnostic)
% -----------------------------
ax1 = nexttile(tl,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
colormap(fig15, cmap_);

for wI = 1:wN
    for pI = 1:pN
        if ~use_case(wI,pI), continue; end
        if isfinite(offx_exp(wI,pI)) && isfinite(offx_theo(wI,pI)) && isfinite(Ga(wI,pI))

            if mark_W1_effect(wI,pI)
                mk = specialMarker;
                ms = 75;
                lw = 1.2;
            else
                mk = markers{wI};
                ms = 55;
                lw = 0.5;
            end

            scatter(ax1, offx_exp(wI,pI), offx_theo(wI,pI), ms, Ga(wI,pI), ...
                'filled', 'Marker', mk, 'MarkerEdgeColor','k', 'LineWidth',lw);
        end
    end
end

xline(ax1,0,'k--','LineWidth',0.8);
yline(ax1,0,'k--','LineWidth',0.8);

vals = [offx_exp(use_case); offx_theo(use_case)];
vals = vals(isfinite(vals));
if ~isempty(vals)
    lim0 = [min(vals), max(vals)];
    if diff(lim0) == 0
        lim0 = lim0 + [-1 1]*1e-4;
    end
    pad = 0.15 * diff(lim0);
    lim1 = lim0 + [-pad pad];
    xlim(ax1, lim1); ylim(ax1, lim1);
end

xlabel(ax1, '$\overline{w}_{x,\mathrm{exp.}}$', 'Interpreter','latex');
ylabel(ax1, '$\overline{w}_{x,\mathrm{theo.}}$', 'Interpreter','latex');
text(ax1, -0.20, 1.04, '(a)', 'Units','normalized', 'FontSize',18);

HL = gobjects(1,wN+1);
for wI = 1:wN
    HL(wI) = plot(ax1, NaN, NaN, ...
        'LineStyle','none', ...
        'Marker', markers{wI}, ...
        'MarkerEdgeColor','k', ...
        'MarkerFaceColor','w', ...
        'DisplayName', waves{wI});
end
HL(wN+1) = plot(ax1, NaN, NaN, ...
    'LineStyle','none', ...
    'Marker', specialMarker, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',1.2, ...
    'DisplayName','W1-G1-3, E1-3');

legend(ax1, HL, 'Location','northwest', 'NumColumns',2);

% -----------------------------
% (b) x-direction peak-to-peak
% -----------------------------
ax2 = nexttile(tl,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

for wI = 1:wN
    for pI = 1:pN
        if ~use_case(wI,pI), continue; end
        if isfinite(Rx_exp(wI,pI)) && isfinite(Rx_theo(wI,pI)) && isfinite(Ga(wI,pI))

            if mark_W1_effect(wI,pI)
                mk = specialMarker;
                ms = 75;
                lw = 1.2;
            else
                mk = markers{wI};
                ms = 55;
                lw = 0.5;
            end

            scatter(ax2, Rx_exp(wI,pI), Rx_theo(wI,pI), ms, Ga(wI,pI), ...
                'filled', 'Marker', mk, 'MarkerEdgeColor','k', 'LineWidth',lw);
        end
    end
end

vals = [Rx_exp(use_case); Rx_theo(use_case)];
vals = vals(isfinite(vals));
if ~isempty(vals)
    lim0 = [min(vals), max(vals)];
    if diff(lim0) == 0
        lim0 = lim0 + [-1 1]*1e-4;
    end
    pad = 0.10 * diff(lim0);
    lim1 = lim0 + [-pad pad];
    plot(ax2, lim1, lim1, 'k--', 'LineWidth',0.8);
    xlim(ax2, lim1); ylim(ax2, lim1);
end

xlabel(ax2, '$w_{x,\mathrm{exp.}}^{\mathrm{p\mbox{-}p}}$', 'Interpreter','latex');
ylabel(ax2, '$w_{x,\mathrm{theo.}}^{\mathrm{p\mbox{-}p}}$', 'Interpreter','latex');
text(ax2, -0.20, 1.04, '(b)', 'Units','normalized', 'FontSize',18);

% -----------------------------
% (c) x-direction full-wave lag delta*_x
% marker size proportional to rho*_x
% -----------------------------
ax3 = nexttile(tl,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');

for wI = 1:wN
    for pI = 1:pN
        if ~use_case(wI,pI), continue; end
        if isfinite(St(wI,pI)) && isfinite(delta_x_star(wI,pI)) && ...
           isfinite(rho_x_star(wI,pI)) && isfinite(Ga(wI,pI))

            msz_case = msz0 + mszAmp * max(0, min(1, rho_x_star(wI,pI)));

            if mark_W1_effect(wI,pI)
                mk = specialMarker;
                ms = 1.25 * msz_case;
                lw = 1.2;
            else
                mk = markers{wI};
                ms = msz_case;
                lw = 0.4;
            end

            scatter(ax3, St(wI,pI), delta_x_star(wI,pI), ms, Ga(wI,pI), ...
                'filled', 'Marker', mk, 'MarkerEdgeColor','k', 'LineWidth',lw);
        end
    end
end

yline(ax3,0,'k--','LineWidth',0.8);
xlabel(ax3, '$St$', 'Interpreter','latex');
ylabel(ax3, '$\delta_x^\ast$', 'Interpreter','latex');
ylim(ax3, [-pi pi]);
set(ax3, 'YTick', [-pi -0.5*pi 0 0.5*pi pi], ...
         'YTickLabel', {'-\pi','-0.5\pi','0','0.5\pi','\pi'});
text(ax3, 0.75, 0.95, '$\mathrm{size}\propto \rho_x^\ast$', ...
    'Units','normalized','Interpreter','latex','FontSize',12);
text(ax3, -0.20, 1.04, '(c)', 'Units','normalized', 'FontSize',18);

cb = colorbar(ax3);
cb.Label.String = '$Ga$';
cb.Label.Interpreter = 'latex';
cb.FontSize = 14;

% -----------------------------
% (d) z-direction phase mean
% -----------------------------
ax4 = nexttile(tl,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');

for wI = 1:wN
    for pI = 1:pN
        if ~use_case(wI,pI), continue; end
        if isfinite(offz_exp(wI,pI)) && isfinite(offz_theo(wI,pI)) && isfinite(Ga(wI,pI))

            if mark_W1_effect(wI,pI)
                mk = specialMarker;
                ms = 75;
                lw = 1.2;
            else
                mk = markers{wI};
                ms = 55;
                lw = 0.5;
            end

            scatter(ax4, offz_exp(wI,pI), offz_theo(wI,pI), ms, Ga(wI,pI), ...
                'filled', 'Marker', mk, 'MarkerEdgeColor','k', 'LineWidth',lw);
        end
    end
end

vals = [offz_exp(use_case); offz_theo(use_case)];
vals = vals(isfinite(vals));
if ~isempty(vals)
    lim0 = [min(vals), max(vals)];
    if diff(lim0) == 0
        lim0 = lim0 + [-1 1]*1e-4;
    end
    pad = 0.10 * diff(lim0);
    lim1 = lim0 + [-pad pad];
    plot(ax4, lim1, lim1, 'k--', 'LineWidth',0.8);
    xlim(ax4, lim1); ylim(ax4, lim1);
end

xlabel(ax4, '$-\overline{w}_{z,\mathrm{exp.}}$', 'Interpreter','latex');
ylabel(ax4, '$-\overline{w}_{z,\mathrm{theo.}}$', 'Interpreter','latex');
text(ax4, -0.20, 1.04, '(d)', 'Units','normalized', 'FontSize',18);

% -----------------------------
% (e) z-direction peak-to-peak
% -----------------------------
ax5 = nexttile(tl,5); hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');

for wI = 1:wN
    for pI = 1:pN
        if ~use_case(wI,pI), continue; end
        if isfinite(Rz_exp(wI,pI)) && isfinite(Rz_theo(wI,pI)) && isfinite(Ga(wI,pI))

            if mark_W1_effect(wI,pI)
                mk = specialMarker;
                ms = 75;
                lw = 1.2;
            else
                mk = markers{wI};
                ms = 55;
                lw = 0.5;
            end

            scatter(ax5, Rz_exp(wI,pI), Rz_theo(wI,pI), ms, Ga(wI,pI), ...
                'filled', 'Marker', mk, 'MarkerEdgeColor','k', 'LineWidth',lw);
        end
    end
end

vals = [Rz_exp(use_case); Rz_theo(use_case)];
vals = vals(isfinite(vals));
if ~isempty(vals)
    lim0 = [min(vals), max(vals)];
    if diff(lim0) == 0
        lim0 = lim0 + [-1 1]*1e-4;
    end
    pad = 0.10 * diff(lim0);
    lim1 = lim0 + [-pad pad];
    plot(ax5, lim1, lim1, 'k--', 'LineWidth',0.8);
    xlim(ax5, lim1); ylim(ax5, lim1);
end

xlabel(ax5, '$w_{z,\mathrm{exp.}}^{\mathrm{p\mbox{-}p}}$', 'Interpreter','latex');
ylabel(ax5, '$w_{z,\mathrm{theo.}}^{\mathrm{p\mbox{-}p}}$', 'Interpreter','latex');
text(ax5, -0.20, 1.04, '(e)', 'Units','normalized', 'FontSize',18);

% -----------------------------
% (f) z-direction full-wave lag delta*_z
% marker size proportional to rho*_z
% -----------------------------
ax6 = nexttile(tl,6); hold(ax6,'on'); box(ax6,'on'); grid(ax6,'on');

for wI = 1:wN
    for pI = 1:pN
        if ~use_case(wI,pI), continue; end
        if isfinite(St(wI,pI)) && isfinite(delta_z_star(wI,pI)) && ...
           isfinite(rho_z_star(wI,pI)) && isfinite(Ga(wI,pI))

            msz_case = msz0 + mszAmp * max(0, min(1, rho_z_star(wI,pI)));

            if mark_W1_effect(wI,pI)
                mk = specialMarker;
                ms = 1.25 * msz_case;
                lw = 1.2;
            else
                mk = markers{wI};
                ms = msz_case;
                lw = 0.4;
            end

            scatter(ax6, St(wI,pI), delta_z_star(wI,pI), ms, Ga(wI,pI), ...
                'filled', 'Marker', mk, 'MarkerEdgeColor','k', 'LineWidth',lw);
        end
    end
end

yline(ax6,0,'k--','LineWidth',0.8);
xlabel(ax6, '$St$', 'Interpreter','latex');
ylabel(ax6, '$\delta_z^\ast$', 'Interpreter','latex');
ylim(ax6, [-pi pi]);
set(ax6, 'YTick', [-pi -0.5*pi 0 0.5*pi pi], ...
         'YTickLabel', {'-\pi','-0.5\pi','0','0.5\pi','\pi'});
text(ax6, 0.75, 0.95, '$\mathrm{size}\propto \rho_z^\ast$', ...
    'Units','normalized','Interpreter','latex','FontSize',12);
text(ax6, -0.20, 1.04, '(f)', 'Units','normalized', 'FontSize',18);

% common style
set([ax1 ax2 ax3 ax4 ax5 ax6], 'FontSize',16, 'LineWidth',0.8);

for ax = [ax1 ax2 ax3 ax4 ax5 ax6]
    caxis(ax, [0 400]);
end

colormap(fig15, cmap_);
cb = colorbar(ax6);
cb.Label.String = '$Ga$';
cb.Label.Interpreter = 'latex';
cb.FontSize = 14;

% exportgraphics(fig15,'E:\2.ParticlesInWave\Manuscript\figures\sol_comparison.pdf')


%%
function [delta_star, rho_star] = wave_align_metric( ...
    epsv, R1_exp, th1_exp, R2_exp, th2_exp, ...
          R1_theo, th1_theo, R2_theo, th2_theo)
% FULL-WAVE O(eps^2) alignment metric based on the modulation part only:
%
%   w_tilde(phi) = eps*R1*sin(phi+th1) + eps^2*R2*sin(2phi+th2)
%
% delta_star : lag d that maximizes normalized waveform similarity
% rho_star   : maximum normalized similarity
%
% NOTE:
% This compares the entire O(eps^2) modulation waveform, not the peak only.

    delta_star = NaN;
    rho_star   = NaN;

    vals = [epsv, R1_exp, th1_exp, R2_exp, th2_exp, ...
                  R1_theo, th1_theo, R2_theo, th2_theo];
    if any(~isfinite(vals))
        return
    end

    phi = linspace(0, 2*pi, 721);
    phi(end) = [];   % remove duplicated endpoint

    w_exp = epsv * R1_exp * sin(phi + th1_exp) ...
          + epsv^2 * R2_exp * sin(2*phi + th2_exp);

    Eexp = trapz(phi, w_exp.^2);
    if ~(isfinite(Eexp) && Eexp > eps)
        return
    end

    delta_grid = linspace(-pi, pi, 1441);
    rho_all = nan(size(delta_grid));

    for ii = 1:numel(delta_grid)
        d = delta_grid(ii);

        w_theo = epsv * R1_theo * sin(phi + d + th1_theo) ...
               + epsv^2 * R2_theo * sin(2*(phi + d) + th2_theo);

        Etheo = trapz(phi, w_theo.^2);
        if ~(isfinite(Etheo) && Etheo > eps)
            continue
        end

        rho_all(ii) = trapz(phi, w_exp .* w_theo) / sqrt(Eexp * Etheo);
    end

    if all(~isfinite(rho_all))
        return
    end

    rho_max = max(rho_all, [], 'omitnan');
    idx_all = find(abs(rho_all - rho_max) < 1e-10);

    if isempty(idx_all)
        return
    end

    % If multiple maxima exist, choose the one with smallest |delta|
    [~, jj] = min(abs(delta_grid(idx_all)));
    idx = idx_all(jj);

    rho_star   = rho_all(idx);
    delta_star = delta_grid(idx);
    delta_star = mod(delta_star + pi, 2*pi) - pi;
end
