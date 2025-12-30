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

% --- Load per-case results ---
for wI = 1:wN
  wave = waves{wI};
  for pI = 1:pN
    partinfo = parts{pI};
    
    S = load(fullfile([wave partinfo], sprintf('%s%s_vs_avg_unit.mat', wave, partinfo)));

    % Phase bins and slip statistics (already binned in saved .mat)
    thetas_std{wI,pI} = S.thetas_std;     % phi_std [rad]
    r_avg_u{wI,pI}    = S.r_avg_u;        % [m/s]
    r_avg_v{wI,pI}    = S.r_avg_v;        % [m/s]
    r_se_u{wI,pI}     = S.r_se_u;         % [m/s]
    r_se_v{wI,pI}     = S.r_se_v;         % [m/s]

    % Wave/particle parameters from saved results
    v_s(wI,pI)        = S.v_s;
    k(wI,pI)          = S.k;
    h(wI,pI)          = S.h;
    T(wI,pI)          = S.T;
    Ec(wI,pI)         = S.E_cosh_ens_out;
    Es(wI,pI)         = S.E_sinh_ens_out;

    % Wave-averaged profiles (kz-binned)
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

%% theory vs. measurements (figure 10)
close all

% Preallocate outputs used in plots/tables
Rx_theo = nan(wN,pN);  Rz_theo = nan(wN,pN);
Rx_exp  = nan(wN,pN);  Rz_exp  = nan(wN,pN);

phix_peak_exp  = nan(wN,pN); phiz_peak_exp  = nan(wN,pN);
phix_peak_theo = nan(wN,pN); phiz_peak_theo = nan(wN,pN);

theta_x1_exp = nan(wN,pN); theta_z1_exp = nan(wN,pN);
theta_x2_exp = nan(wN,pN); theta_z2_exp = nan(wN,pN);

theta_x1_theo = nan(wN,pN); theta_z1_theo = nan(wN,pN);
theta_x2_theo = nan(wN,pN); theta_z2_theo = nan(wN,pN);

off_exp  = nan(wN,pN);
off_theo = nan(wN,pN);

% Coefficients saved for table export
pxM  = nan(wN,pN); qxM  = nan(wN,pN); pzM  = nan(wN,pN); qzM  = nan(wN,pN);
a0LM = nan(wN,pN); acLM = nan(wN,pN); asLM = nan(wN,pN);
b0LM = nan(wN,pN); bcLM = nan(wN,pN); bsLM = nan(wN,pN);

pzE = nan(wN,pN); qzE = nan(wN,pN);
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

fig = figure(1); 
fig.Position = [80 80 1200 500];

for wI = 1:wN
  wave = waves{wI};
  for pI = 1:pN
    partinfo = parts{pI};

    fprintf('%s %s\n', wave, partinfo);

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
    qz = (  St_*Aef*vs0_ + St2*(1-beta_)*Bef )/(St2+1);

    % 2nd-order harmonic coefficients (0th and 2nd harmonics)
    C0 = 0.5*( Aef*qx + Bef*pz ) ...
        + 0.5/St_*( Aef*px + px^2 + qx^2 ) ...
        + 0.5*(1-beta_)*Aef*qx ...
        - 0.5*vs0_*Bef*( Aef + px );

    Cc = 0.5*( -Aef*qx + Bef*pz ) ...
        + 0.5/St_*( Aef*px + px^2 - qx^2 ) ...
        - 0.5*(1-beta_)*Aef*qx ...
        - 0.5*vs0_*Bef*( Aef + px );

    Cs = 0.5*( Aef*px + Bef*qz ) ...
        + 0.5/St_*( Aef*qx + 2*px*qx ) ...
        + 0.5*(1-beta_)*Aef*(Aef + px) ...
        - 0.5*vs0_*Bef*qx + 0.5*(1-beta_)*(Bef^2 - Aef^2);

    E0 = 0.5*( Bef*px + Aef*qz ) ...
       + 0.5/St_*( Aef*pz + px*pz + qx*qz ) ...
       - 0.5*(1-beta_)*Bef*(Aef + px) ...
       - 0.5*vs0_*Aef*qx;

    Ec2 = 0.5*( Bef*px - Aef*qz ) ...
       + 0.5/St_*( Aef*pz + px*pz - qx*qz ) ...
       - 0.5*(1-beta_)*Bef*(Aef + px) ...
       + 0.5*vs0_*Aef*qx ...
       + (1-beta_)*Aef*Bef;

    Es2 = 0.5*( Bef*qx + Aef*pz ) ...
       + 0.5/St_*( Aef*qz + px*qz + qx*pz ) ...
       - 0.5*(1-beta_)*Bef*qx ...
       - 0.5*vs0_*(Aef^2 + Aef*px) ...
       + (1-beta_)*Aef*Bef;

    Gamma0 = C0 + 0.5*K_*(px*pz + qx*qx);
    Gamma1 = Cc + 0.5*K_*(px*pz - qx*qx);
    Gamma2 = Cs + 0.5*K_*(px*qz + qx*px);

    Delta0 = E0 + 0.5*K_*(pz^2 + qz^2);
    Delta1 = Ec2 + 0.5*K_*(pz^2 - qz^2);
    Delta2 = Es2 +     K_*(pz * qz);

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

    % --- Vertical data fit (1st and 2nd harmonics) ---
    Xv = [cos(x), sin(x), cos(2*x), sin(2*x), ones(numel(x),1)];
    sev = rv_se(mexp);
    wv  = 1./max(sev.^2, eps);
    W12v = sqrt(wv(:));
    coef_v = (Xv.*W12v) \ (yv.*W12v);

    A1v=coef_v(1); B1v=coef_v(2); A2v=coef_v(3); B2v=coef_v(4); Qexp=coef_v(5);
    Qexp_mean     = mean(yv);
    off_exp(wI,pI)= Qexp_mean;

    rv_fit_exp = A1v*cos(phi_plot) + B1v*sin(phi_plot) + ...
                 A2v*cos(2*phi_plot) + B2v*sin(2*phi_plot) + Qexp;

    % --- Horizontal data fit (1st and 2nd harmonics) ---
    Xu = [cos(x), sin(x), cos(2*x), sin(2*x), ones(numel(x),1)];
    seu = ru_se(mexp);
    wu  = 1./max(seu.^2, eps);
    W12u = sqrt(wu(:));
    coef_u = (Xu.*W12u) \ (yu.*W12u);

    A1u = coef_u(1); B1u = coef_u(2);
    A2u = coef_u(3); B2u = coef_u(4);
    Qu  = coef_u(5);

    ru_fit_exp = A1u*cos(phi_plot) + B1u*sin(phi_plot) + ...
                 A2u*cos(2*phi_plot) + B2u*sin(2*phi_plot) + Qu;

    % --- Experimental peaks/amplitudes ---
    [xmax_exp,imax_x] = max(ru_fit_exp);
    [xmin_exp,~]      = min(ru_fit_exp);
    [ymax_exp,imax_y] = max(rv_fit_exp);
    [ymin_exp,~]      = min(rv_fit_exp);

    Rx_exp(wI,pI) = xmax_exp - xmin_exp;
    Rz_exp(wI,pI) = ymax_exp - ymin_exp;

    theta_x1_exp(wI,pI) = atan2(A1u,B1u); 
    theta_z1_exp(wI,pI) = atan2(A1v,B1v);
    theta_x2_exp(wI,pI) = atan2(A2u,B2u); 
    theta_z2_exp(wI,pI) = atan2(A2v,B2v);

    phix_peak_exp(wI,pI) = phi_plot(imax_x);
    phiz_peak_exp(wI,pI) = phi_plot(imax_y);

    % --- Theoretical peaks/amplitudes from perturbation solution ---
    [xmax_theo,imax_x] = max(wx_lin2);
    [xmin_theo,~]      = min(wx_lin2);

    [ymax_theo,imax_y] = max(r_lin2);
    [ymin_theo,~]      = min(r_lin2);

    Rx_theo(wI,pI) = xmax_theo - xmin_theo;
    Rz_theo(wI,pI) = ymax_theo - ymin_theo;

    theta_x1_theo(wI,pI) = atan2(px,qx);   
    theta_z1_theo(wI,pI) = atan2(pz,qz);
    theta_x2_theo(wI,pI) = atan2(acL,asL); 
    theta_z2_theo(wI,pI) = atan2(bcL,bsL);

    phix_peak_theo(wI,pI) = phi_plot(imax_x);
    phiz_peak_theo(wI,pI) = phi_plot(imax_y);

    % Theoretical offset (kept for later scatter plot)
    Qthe2 = mean(wz_lin2);
    off_theo(wI,pI) = Qthe2;

    % Save coefficients for table export
    pxM(wI,pI)  = px;   qxM(wI,pI)  = qx;   pzM(wI,pI)  = pz;   qzM(wI,pI)  = qz;
    a0LM(wI,pI) = a0L;  acLM(wI,pI) = acL;  asLM(wI,pI) = asL;
    b0LM(wI,pI) = b0L;  bcLM(wI,pI) = bcL;  bsLM(wI,pI) = bsL;

    pzE(wI,pI) = A1v; qzE(wI,pI) = B1v;
    bcE(wI,pI) = A2v; bsE(wI,pI) = B2v;

    % --- Plot update: data ± SE, experimental fit, theory (perturbation) ---
    % In 3x3, slot (3,3,3) is reserved for legend dummy, so map pI=1..8 to 1,2,4,5,6,7,8,9
    spI = pI + (pI>=3);

    % Figure 1: rv_mean and theory curve
    figure(1)
    subplot(3,3,spI), hold on; grid on

    xfill = [phi_plot; flipud(phi_plot)];
    yfill = [rv_mean - rv_se; flipud(rv_mean + rv_se)];
    if all(isfinite(yfill))
        f1 = fill(xfill, yfill, [0.88 0.88 0.88], 'EdgeColor','none');
        uistack(f1,'bottom');
    end

    plot(phi_plot, rv_mean, '.', 'Color', c_wI{wI}, 'LineWidth',1.6, ...
         'DisplayName',sprintf('%s meas.',wave));
    plot(phi_plot, rv_fit_exp,'--', 'Color', c_wI{wI}, 'LineWidth',1.6, ...
         'DisplayName',sprintf('%s fitted',wave));
    plot(phi_plot, r_lin2, '-.', 'Color', c_wI{wI}, 'LineWidth',1.2, ...
         'DisplayName',sprintf('%s asym.',wave));

    xlim([0 2*pi]); 
    xlabel('$\phi$','Interpreter','latex'); 
    ylabel('$w_z$','Interpreter','latex');
    ylim([0 0.12])

    set(gca,'XTick',[0 pi/2, pi , 1.5*pi, 2*pi], ...
            'XTickLabel',{'0','0.5\pi','\pi','1.5\pi','2\pi'})
    title(partinfo); 
    set(gca,'FontSize',14)

    % Skip some panels for W1 and large particles (as in original)
    if wI == 1 && pI >= 7; 
        continue; 
    end

    % Figure 2: detrended w_z comparison (exp vs theory)
    figure(2)
    subplot(3,3,spI), hold on; grid on

    plot(phi_plot, -rv_mean + Qexp_mean, '-', ...
        'Color',1-(1-c_wI{wI})*0.3,'LineWidth',1.3, ...
        'DisplayName',sprintf('$\\langle w_z \\rangle_\\phi - \\overline{\\langle w_z \\rangle_\\phi} \\mathrm{(%s)}$',wave));

    xfill = [phi_plot; flipud(phi_plot)];
    yfill = [-rv_mean + Qexp_mean - rv_se; flipud(-rv_mean + Qexp_mean + rv_se)];
    f1 = fill(xfill, yfill, 1-(1-c_wI{wI})*0.5, 'EdgeColor','none','FaceAlpha',0.5);
    uistack(f1,'bottom')

    plot(phi_plot, -rv_fit_exp + Qexp, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\hat{w_z} - \\overline{\\hat{w_z}} \\, \\mathrm{(%s)}$',wave));
    plot(phi_plot, wz_lin2 + vs0_,'--','Color',c_wI{wI}, 'LineWidth',2, ...
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
         'DisplayName',sprintf('$\\langle w_x \\rangle_\\phi - \\overline{\\langle w_x \\rangle_\\phi} \\mathrm{(%s)}$',wave));
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

% Figure 1 legend (meas + asym)
figure(1)
axL = subplot(3,3,3); cla(axL); axis(axL,'off'); set(axL,'Visible','off'); hold(axL,'on')
h1 = gobjects(wN,1); h4 = gobjects(wN,1);
for wI = 1:wN
    wave = waves{wI};
    h1(wI) = plot(axL, nan, nan, '.', 'Color', c_wI{wI}, 'LineWidth',1.6, ...
        'DisplayName',sprintf('%s meas.',wave));
end
for wI = 1:wN
    wave = waves{wI};
    h4(wI) = plot(axL, nan, nan, '-.', 'Color', c_wI{wI}, 'LineWidth',1.2, ...
        'DisplayName',sprintf('%s asym.',wave));
end
legend(axL, [h1; h4], 'Location','best','Interpreter','latex','FontSize',20);

% Figure 2 legend (pp1 + pp2 + pp3)
fig2 = figure(2); fig2.Position = [-500 100 1000 1000];
axL = subplot(3,3,3); cla(axL); axis(axL,'off'); set(axL,'Visible','off'); hold(axL,'on')
hpp1 = gobjects(wN,1); hpp2 = gobjects(wN,1); hpp3 = gobjects(wN,1);
for wI = 1:wN
    wave = waves{wI};
    hpp1(wI) = plot(axL, nan, nan, '-', 'Color',1-(1-c_wI{wI})*0.3, 'LineWidth',1.3, ...
        'DisplayName',sprintf('$\\langle w_z \\rangle_\\phi - \\overline{\\langle w_z \\rangle_\\phi} \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hpp2(wI) = plot(axL, nan, nan, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,\\hat{w_z} - \\overline{\\hat{w_z}} \\, \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hpp3(wI) = plot(axL, nan, nan, '--', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,w_z + v_s \\, \\mathrm{(Eq.(2.17),\\,%s)}$',wave));
end
legend(axL, [hpp1; hpp2; hpp3], 'Location','best','Interpreter','latex','FontSize',16);

% Figure 3 legend (q1 + q3 + q4)
fig3 = figure(3); fig3.Position = [-500 100 1000 1000];
axL = subplot(3,3,3); cla(axL); axis(axL,'off'); set(axL,'Visible','off'); hold(axL,'on')
hq1 = gobjects(wN,1); hq3 = gobjects(wN,1); hq4 = gobjects(wN,1);
for wI = 1:wN
    wave = waves{wI};
    hq1(wI) = plot(axL, nan, nan, '-', 'Color',1-(1-c_wI{wI})*0.3, 'LineWidth',1.3, ...
        'DisplayName',sprintf('$\\langle w_x \\rangle_\\phi - \\overline{\\langle w_x \\rangle_\\phi} \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hq3(wI) = plot(axL, nan, nan, '-.', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,\\hat{w_x} - \\overline{\\hat{w_x}} \\, \\mathrm{(%s)}$',wave));
end
for wI = 1:wN
    wave = waves{wI};
    hq4(wI) = plot(axL, nan, nan, '--', 'Color',c_wI{wI}, 'LineWidth',2, ...
        'DisplayName',sprintf('$\\,w_x \\, \\mathrm{(Eq.(2.17),\\,%s)}$',wave));
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

%% Summary scatter plots (figure 11)

markers = {'o','^','s'};    % W1, W2, W3
msz = 60;

phi_diff_z = phiz_peak_theo - phiz_peak_exp;

fig = figure(); 
fig.Position = [300 01 1800 500];

ax_figs(1) = subplot(1,3,1); hold on; grid on;
cmap_ = turbo(400);
colormap(cmap_);
caxis([0 400])

HH = gobjects(1,wN);
for wI = 1:wN
    if wI == 1
        HH(wI) = scatter(off_exp(wI,1:end-2), -off_theo(wI,1:end-2), msz, Ga(wI,1:end-2), ...
                            'filled', 'Marker', markers{wI}, 'DisplayName', waves{wI});
    else
        HH(wI) = scatter(off_exp(wI,:), -off_theo(wI,:), msz, Ga(wI,:), ...
                        'filled', 'Marker', markers{wI}, 'DisplayName', waves{wI});
    end
end
plot([0 1], [0 1],'k--')
xlim([0 0.15]), ylim([0 0.15])
xlabel('$-\overline{w_{z,exp.}}$','Interpreter','latex')
ylabel('$-\overline{w_{z,asym.}}$','Interpreter','latex')
set(gca,'FontSize',20)

% Dummy handles for legend (black markers)
HL = gobjects(1,wN);
for wI = 1:wN
    HL(wI) = plot(NaN, NaN, ...
                  'LineStyle','none', ...
                  'Marker', markers{wI}, ...
                  'MarkerEdgeColor','k', ...
                  'MarkerFaceColor','w', ...
                  'DisplayName', waves{wI});
end
lgd = legend(HL);
set(lgd, 'Location','northeast','NumColumns',3);

ax_figs(2) = subplot(1,3,2); hold on; grid on;
colormap(cmap_);
caxis([0 400])

HH = gobjects(1,wN);
for wI = 1:wN
    if wI == 1
        HH(wI) = scatter(Rz_exp(wI,1:end-2), Rz_theo(wI,1:end-2), msz, Ga(wI,1:end-2), ...
                        'filled', 'Marker', markers{wI}, 'DisplayName', waves{wI});
    else
        HH(wI) = scatter(Rz_exp(wI,:), Rz_theo(wI,:), msz, Ga(wI,:), ...
                        'filled', 'Marker', markers{wI}, 'DisplayName', waves{wI});
    end
end
plot([0 1], [0 1],'k--')

xlabel('$-w_{z,\mathrm{exp.}}^{\mathrm{p\mbox{-}p}}$','Interpreter','latex')
ylabel('$-w_{z,\mathrm{asym.}}^{\mathrm{p\mbox{-}p}}$','Interpreter','latex')
xlim([0 0.005])
ylim([0 0.005])
set(gca,'FontSize',20)

ax_figs(3) = subplot(1,3,3); hold on; grid on;
colormap(cmap_);
caxis([0 400])

HH = gobjects(1,wN);
for wI = 1:wN
    if wI==1
        HH(wI) = scatter(St(wI,1:end-2), phi_diff_z(wI,1:end-2), msz, Ga(wI,1:end-2), ...
                    'filled', 'Marker', markers{wI}, 'DisplayName', waves{wI});
    else
        HH(wI) = scatter(St(wI,:), phi_diff_z(wI,:), msz, Ga(wI,:), ...
                    'filled', 'Marker', markers{wI}, 'DisplayName', waves{wI});
    end
end
yline(0,'k--','LineWidth',1)

cb = colorbar();
cb.Label.String = '$Ga$';
cb.Label.Interpreter = 'latex';

xlabel('$St$','Interpreter','latex')
ylabel('$\Delta \phi_{\mathrm{peak}}$','Interpreter','latex')
ylim([-pi pi])
set(gca,'FontSize',20, ...
       'YTick',[-pi -0.5*pi 0 0.5*pi pi], ...
       'YTickLabel',{'-\pi','-0.5\pi','0','0.5\pi','\pi'})

labels = 'abc';
for sp = 1:3
    ax = ax_figs(sp);
    text(ax,-0.15,1.02,sprintf('(%c)',labels(sp)), 'Units','normalized', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Clipping','off')
end

% exportgraphics(fig,'E:\2.ParticlesInWave\Manuscript\figures\sol_comparison.pdf')