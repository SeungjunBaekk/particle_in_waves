%% ----------------- slip profile (figure 4) ------

clear; close all; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};
parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};

wN = numel(waves);
pN = numel(parts);

% PIV y is positive downward; linear-wave z is positive upward with z=0 at SWL.
y_sign = -1;

num_bins_ky  = 10;     % vertical bins in k*y
B_boot       = 10000;  % bootstrap iterations
alpha_ci     = 0.05;   % e.g. 95% CI

% Storage for binning based on wave-averaged location (track-mean ky)
kz_bin_z_mean = cell(wN,pN);   % <k y_p> in each bin
kz_bin_z_std  = cell(wN,pN);   % std(k y_p) in each bin
kz_bin_v_mean = cell(wN,pN);   % <w_slip> in each bin [m/s]
kz_bin_v_std  = cell(wN,pN);   % std(w_slip) in each bin [m/s]

v_s_mat     = NaN(wN,pN);      % v_s for each (wave, part)
ky_tau5_mat = NaN(wN,pN);      % mean k*y_p at t'=5*tau_p

% Wave colors and markers
c_wI = { ...
    [0 0 0]                % W1: black
    [0.698 0.1333 0.1333]  % W2: red brick
    [0 0.5451 0.5451]      % W3: teal-ish
    };
m_wI = {'o','^','s'};      % W1, W2, W3

% Helper to lighten a color (alpha=0: original color, alpha=1: white)
lighten = @(c, alpha) 1 - alpha*(1 - c);  % 0<=alpha<=1

% Lightness level for trials
alpha_trial = 0.4;

% ----------------- Particle relaxation time from still-water data --------
pInfo = readtable('particle_stag_settling.CSV');

nu    = waterNu(24.5);
rho_w = waterRho(24.5);

tau_p_part = NaN(1,pN);

for pI = 1:pN
    partinfo = parts{pI};
    case_mask = strcmp(pInfo.Case, partinfo);

    if nnz(case_mask) ~= 1
        error('particle_stag_settling.CSV must contain exactly one row for %s.', partinfo);
    end

    d_p   = pInfo.d(case_mask)     * 1e-6;   % [m]
    rho_p = pInfo.rho_p(case_mask) * 1000;   % [kg/m3]
    Re_0  = pInfo.Re(case_mask);             % [-], based on still-water settling

    SchNau = 1 + 0.15*Re_0.^0.687;

    tau_p_part(pI) = d_p^2*(rho_p/rho_w + 0.5) ...
                   / (18*nu*SchNau);         % [s]
end

% ----------------- Prepare 3x3 subplot (figure 4) ------------------------
fig_profile_all = figure('Color','w');
ax_profile = gobjects(1,9);

for sp = 1:9
    ax_profile(sp) = subplot(3,3,sp);

    if sp == 3
        % 3rd slot is reserved for legend
        axis(ax_profile(sp),'off');
    else
        hold(ax_profile(sp),'on');
        grid(ax_profile(sp),'on');
        set(ax_profile(sp),'YDir','reverse','FontSize',16);
        xlabel(ax_profile(sp),'$-w_z/v_s$','Interpreter','latex');
        ylabel(ax_profile(sp),'$z_p$','Interpreter','latex');
        xlim(ax_profile(sp),[0.5 1.5]);
    end
end

% Place 8 particles into 3x3 subplots, leaving the 3rd slot empty
slotForPart = [1 2 4 5 6 7 8 9];

for pI = 1:pN
    ax = ax_profile(slotForPart(pI));
    title(ax, parts{pI}, 'Interpreter','none');
end

% Read wave parameters
wave_cond = readtable('waveInfo.CSV');

% ----------------- Main loop: waves × particles --------------------------
for wI = 1:wN
    wave = waves{wI};

    % --- Extract parameters for this wave ---
    H = wave_cond.H(strcmp(wave_cond.Case,wave));
    T = wave_cond.T(strcmp(wave_cond.Case,wave));
    h = wave_cond.h(strcmp(wave_cond.Case,wave));

    if isempty(H) || isempty(T) || isempty(h)
        warning('Wave %s not found in waveInfo.CSV. Skip.', wave);
        continue;
    end

    omega = 2*pi/T;
    k     = dispersion(h,T);
    a     = H/2;

    fprintf('=== Wave %s: H=%.3f, T=%.3f, h=%.3f, k=%.5f ===\n', ...
            wave, H, T, h, k);

    for pI = 1:pN
        partinfo = parts{pI};
        fprintf('  -> Particle %s\n', partinfo);

        % -----------------------------------------------------------------
        % 1) Load data: slip_vs_results_*.mat
        % -----------------------------------------------------------------
        filepath = [wave, partinfo];
        flist    = dir(fullfile(filepath, 'slip_vs_results_*.mat'));
        case_tot = numel(flist);

        if case_tot == 0
            warning('    No slip_vs_results_*.mat in %s. Skip this (wave,part).', filepath);
            continue;
        end

        y_fit   = cell(1,case_tot);
        t_p     = cell(1,case_tot);
        v_slip  = cell(1,case_tot);
        u_slip  = cell(1,case_tot);
        v_s_all = NaN(1,case_tot);

        for j = 1:case_tot
            S = load(fullfile(filepath, flist(j).name));

            y_fit{j}   = S.y_fit(:);      % [m], PIV y coordinate, positive downward
            t_p{j}     = S.time(:);       % [s]
            v_slip{j}  = S.v_slip(:);     % [m/s]
            u_slip{j}  = S.u_slip(:);     % [m/s]
            v_s_all(j) = S.meta.v_s;      % [m/s]
        end

        % Settling velocity is nominally the same across tracks
        v_s = mean(v_s_all,'omitnan');

        % -----------------------------------------------------------------
        % 2) Define ky range and bins
        % -----------------------------------------------------------------
        k_y_min = +Inf;
        k_y_max = -Inf;

        for j = 1:case_tot
            yj = y_fit{j};
            if isempty(yj), continue; end

            k_y_min = min(k_y_min, k*min(yj));
            k_y_max = max(k_y_max, k*max(yj));
        end

        if ~isfinite(k_y_min) || ~isfinite(k_y_max)
            warning('    Invalid ky range. Skip.');
            continue;
        end

        k_y_edges   = linspace(k_y_min, k_y_max, num_bins_ky+1);
        k_y_centers = 0.5*(k_y_edges(1:end-1) + k_y_edges(2:end));

        % -----------------------------------------------------------------
        % 3) Compute ky-binned means for each track
        % -----------------------------------------------------------------
        ensemble_v_norm = NaN(case_tot, num_bins_ky);  % v_slip / v_s
        ensemble_v_unit = NaN(case_tot, num_bins_ky);  % v_slip [m/s]
        ensemble_u_unit = NaN(case_tot, num_bins_ky);  % u_slip [m/s]

        for j = 1:case_tot
            k_y_data = k * y_fit{j};  % ky
            v_data_u = v_slip{j};     % [m/s]
            u_data_u = u_slip{j};     % [m/s]

            if isempty(k_y_data) || isempty(v_data_u)
                continue;
            end

            v_data_n = v_data_u / v_s;

            [~,~,bin_idx] = histcounts(k_y_data, k_y_edges);

            binned_v_n = NaN(1, num_bins_ky);
            binned_v_u = NaN(1, num_bins_ky);
            binned_u_u = NaN(1, num_bins_ky);

            for b = 1:num_bins_ky
                mask = (bin_idx == b);

                if any(mask)
                    binned_v_n(b) = mean(v_data_n(mask),'omitnan');
                    binned_v_u(b) = mean(v_data_u(mask),'omitnan');
                    binned_u_u(b) = mean(u_data_u(mask),'omitnan');
                end
            end

            ensemble_v_norm(j,:) = binned_v_n;
            ensemble_v_unit(j,:) = binned_v_u;
            ensemble_u_unit(j,:) = binned_u_u;
        end

        % -----------------------------------------------------------------
        % 4) Ensemble mean
        % -----------------------------------------------------------------
        ensemble_binned_v_norm = mean(ensemble_v_norm, 1, 'omitnan');
        ensemble_binned_v_unit = mean(ensemble_v_unit, 1, 'omitnan');
        ensemble_binned_u_unit = mean(ensemble_u_unit, 1, 'omitnan');

        % -----------------------------------------------------------------
        % 5) Estimate per-bin SE via bootstrap
        % -----------------------------------------------------------------
        ensemble_se_v_norm = NaN(1, num_bins_ky);
        ensemble_se_v_unit = NaN(1, num_bins_ky);
        ensemble_se_u_unit = NaN(1, num_bins_ky);

        for b = 1:num_bins_ky
            vals = ensemble_v_norm(:,b);
            vals = vals(~isnan(vals));
            if numel(vals) >= 2
                res = bootstrap_std(vals, B_boot, alpha_ci);
                ensemble_se_v_norm(b) = res.sd_hat ./ sqrt(numel(vals));
            end

            vals = ensemble_v_unit(:,b);
            vals = vals(~isnan(vals));
            if numel(vals) >= 2
                res = bootstrap_std(vals, B_boot, alpha_ci);
                ensemble_se_v_unit(b) = res.sd_hat ./ sqrt(numel(vals));
            end

            vals = ensemble_u_unit(:,b);
            vals = vals(~isnan(vals));
            if numel(vals) >= 2
                res = bootstrap_std(vals, B_boot, alpha_ci);
                ensemble_se_u_unit(b) = res.sd_hat ./ sqrt(numel(vals));
            end
        end

        % -----------------------------------------------------------------
        % 6) Plot vertical profiles
        % -----------------------------------------------------------------
        ax_this = ax_profile(slotForPart(pI));
        hold(ax_this,'on');

        c_wave  = c_wI{wI};
        m_wave  = m_wI{wI};
        c_trial = lighten(c_wave, alpha_trial);

        % Track curves: light color, dashed
        for j = 1:case_tot
            h_trial = plot(ax_this, ensemble_v_norm(j,:), k_y_centers, '--', ...
                'Color', c_trial, ...
                'Marker', m_wave, ...
                'MarkerSize', 3, ...
                'LineWidth', 0.8, ...
                'HandleVisibility','off');

            set(h_trial,'UserData',struct('particle',j,'k',k));
        end

        % Ensemble + SE: dark color, solid
        h_ens = errorbar(ax_this, ensemble_binned_v_norm, k_y_centers, ...
            ensemble_se_v_norm, ensemble_se_v_norm, ...
            'horizontal', ...
            'Color', c_wave, ...
            'LineWidth', 1.3, ...
            'Marker', m_wave, ...
            'MarkerSize', 6, ...
            'CapSize', 8, ...
            'DisplayName', sprintf('%s ensemble', wave));

        set(h_ens,'UserData',struct('particle',0,'k',k));

        xline(ax_this,1,'HandleVisibility','off');

        % -----------------------------------------------------------------
        % 6-1) Draw z_p(t'=5 tau_p) as dash-dotted horizontal line
        %      tau_p is computed from particle_stag_settling.CSV.
        %      The plotted ordinate is k*y_p, consistent with the profile axis.
        % -----------------------------------------------------------------
        tau_p = tau_p_part(pI);
        t_cut = 5*tau_p;

        ky_tau_each = NaN(1,case_tot);

        for j = 1:case_tot
            tj = t_p{j};
            yj = y_fit{j};

            m_valid_tau = isfinite(tj) & isfinite(yj);

            if nnz(m_valid_tau) < 2
                continue;
            end

            tj = tj(m_valid_tau);
            yj = yj(m_valid_tau);

            t_rel = tj - min(tj);

            [t_rel, ord_t] = sort(t_rel);
            yj = yj(ord_t);

            [t_rel, ia] = unique(t_rel, 'stable');
            yj = yj(ia);

            if t_cut >= min(t_rel) && t_cut <= max(t_rel)
                y_tau = interp1(t_rel, yj, t_cut, 'linear');
                ky_tau_each(j) = k*y_tau;
            end
        end

        ky_tau5 = mean(ky_tau_each,'omitnan');
        ky_tau5_mat(wI,pI) = ky_tau5;

        if isfinite(ky_tau5)
            xl = xlim(ax_this);

            h_tau = plot(ax_this, xl, [ky_tau5 ky_tau5], '-.', ...
                'Color', c_wave, ...
                'LineWidth', 1.2, ...
                'HandleVisibility','off');

            set(h_tau,'UserData',struct( ...
                'particle',0, ...
                'k',k, ...
                'tau_p',tau_p, ...
                't_cut',t_cut, ...
                'ky_tau5',ky_tau5));
        end

        % -----------------------------------------------------------------
        % 7) ky-binning based on wave-averaged location
        % -----------------------------------------------------------------
        v_s_mat(wI,pI) = v_s;

        ky_track = NaN(1,case_tot);   % <k y_p> for each track
        v_track  = NaN(1,case_tot);   % <w_slip> for each track

        for j = 1:case_tot
            yj = y_fit{j};
            vj = v_slip{j};

            m_valid = isfinite(yj) & isfinite(vj);

            if any(m_valid)
                ky_track(j) = k * mean(yj(m_valid),'omitnan');
                v_track(j)  = mean(vj(m_valid),'omitnan');
            end
        end

        m_ok = isfinite(ky_track) & isfinite(v_track);
        ky_track = ky_track(m_ok);
        v_track  = v_track(m_ok);

        if isempty(ky_track)
            kz_bin_z_mean{wI,pI} = [];
            kz_bin_z_std{wI,pI}  = [];
            kz_bin_v_mean{wI,pI} = [];
            kz_bin_v_std{wI,pI}  = [];
        else
            [~,~,bin_idx_track] = histcounts(ky_track, k_y_edges);

            ky_mean = NaN(1,num_bins_ky);
            ky_std  = NaN(1,num_bins_ky);
            v_mean  = NaN(1,num_bins_ky);
            v_std   = NaN(1,num_bins_ky);

            for b = 1:num_bins_ky
                mb = (bin_idx_track == b);

                if any(mb)
                    ky_mean(b) = mean(ky_track(mb),'omitnan');
                    ky_std(b)  = std( ky_track(mb),0,'omitnan');
                    v_mean(b)  = mean(v_track(mb),'omitnan');
                    v_std(b)   = std( v_track(mb),0,'omitnan');
                end
            end

            kz_bin_z_mean{wI,pI} = ky_mean;
            kz_bin_z_std{wI,pI}  = ky_std;
            kz_bin_v_mean{wI,pI} = v_mean;
            kz_bin_v_std{wI,pI}  = v_std;
        end

    end % parts
end % waves

% ----------------- Legend subplot ---------------------------------------
ax_leg = ax_profile(3);
axes(ax_leg); cla(ax_leg); axis(ax_leg,'off'); hold(ax_leg,'on');

h_leg_ens   = gobjects(1,wN);
h_leg_trial = gobjects(1,wN);
h_leg_tau   = gobjects(1,wN);

for w = 1:wN
    c_wave  = c_wI{w};
    m_wave  = m_wI{w};
    c_trial = lighten(c_wave, alpha_trial);

    % Ensemble dummy: solid
    h_leg_ens(w) = plot(ax_leg, nan, nan, '-', ...
        'Color', c_wave, ...
        'Marker', m_wave, ...
        'LineWidth', 2, ...
        'MarkerSize', 6, ...
        'DisplayName', sprintf('$\\langle w_z (z_p,t) \\rangle_{[z]} \\, \\mathrm{(%s)}$',waves{w}));

    % Trial dummy: dashed
    h_leg_trial(w) = plot(ax_leg, nan, nan, '--', ...
        'Color', c_trial, ...
        'Marker', m_wave, ...
        'LineWidth', 1.2, ...
        'MarkerSize', 5, ...
        'DisplayName', sprintf('$w_z(z_p,t)\\, \\mathrm{(%s)}$',waves{w}));

    % 5 tau_p position dummy: dash-dotted
    h_leg_tau(w) = plot(ax_leg, nan, nan, '-.', ...
        'Color', c_wave, ...
        'LineWidth', 1.2, ...
        'DisplayName', sprintf('$z_p(t''=5\\tau_{p,t}) \\, \\mathrm{(%s)}$',waves{w}));
end

h_all = [h_leg_ens(:); h_leg_trial(:); h_leg_tau(:)];

legend(ax_leg, h_all, 'FontSize',16,'Interpreter','latex');

% ----------------- Panel labels and object ordering ----------------------
labels = 'ab-cdefgh';

for sp = 1:9
    if sp == 3
        continue;
    end

    ax = ax_profile(sp);
    ch = get(ax,'Children');

    if isempty(ch)
        continue;
    end

    isEns = false(size(ch));

    for i = 1:numel(ch)
        dn = get(ch(i),'DisplayName');

        if ischar(dn) && ~isempty(dn) && contains(dn,'ensemble')
            isEns(i) = true;
        end
    end

    ch_new = [ch(isEns); ch(~isEns)];
    set(ax,'Children',ch_new);

    text(ax,-0.18,1.02,sprintf('(%c)',labels(sp)), ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',18, ...
        'Clipping','off');
end

fig_profile_all.Position = [-1000 -500 1100 1200];

outpdf = 'E:\2.ParticlesInWave\Manuscript\figures\all_profiles.pdf';

set(fig_profile_all,'Renderer','painters');
% 
% exportgraphics(fig_profile_all, outpdf, ...
%     'ContentType','vector', ...             
%     'BackgroundColor','none');                

%% ----------------- Wave-averaged location profile (3x3, figure 5) ------
%   - <w_slip>_track vs <k y_p>_track (wave-averaged location)
%   - Same axes/titles as the 3x3 trial/ensemble plot (y-axis is labeled z_p but values are k y_p)

clear all; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};
parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};

wN = numel(waves); 
pN = numel(parts);

% Subplot positions for 8 particles in a 3x3 layout (3rd slot is for legend)
slotForPart = [1 2 4 5 6 7 8 9];

% Wave colors and markers (match the previous figure)
c_wI = { ...
    [0 0 0]                % W1: black
    [0.698 0.1333 0.1333]  % W2: dark reddish
    [0 0.5451 0.5451]      % W3: teal-ish
    };
m_wI = {'o','^','s'};      % W1, W2, W3

% Mean profile arrays (if needed)
v_avg = cell(wN,pN); v_std = cell(wN,pN);
u_avg = cell(wN,pN); u_std = cell(wN,pN);

% Wave-averaged (representative track) values: <k y_p>, <w_slip>
kz_bin_z_mean = cell(wN,pN);   % <k y_p> (each bin)
kz_bin_z_std  = cell(wN,pN);   % std(k y_p)
kz_bin_v_mean = cell(wN,pN);   % <w_slip>  [m/s]
kz_bin_v_std  = cell(wN,pN);   % std(w_slip) [m/s]

% Std-phase arrays (if needed)
r_avg_u = cell(wN,pN); r_avg_v = cell(wN,pN);
r_se_u  = cell(wN,pN); r_se_v  = cell(wN,pN);
thetas_std = cell(wN,pN);      % phi_std bins

% Scalar parameters
v_s = zeros(wN,pN); 
k   = v_s; 
T   = v_s; 
h   = v_s; 
Ec  = v_s; 
Es  = v_s;
v_slip = nan(wN,pN);

nu    = waterNu(24.5);
rho_w = waterRho(24.5);
g     = 9.81;

% ----------------- Particle relaxation time from still-water data --------
pInfo = readtable('particle_stag_settling.CSV');

nu    = waterNu(24.5);
rho_w = waterRho(24.5);

tau_p_part = NaN(1,pN);

for pI = 1:pN
    partinfo = parts{pI};
    case_mask = strcmp(pInfo.Case, partinfo);

    if nnz(case_mask) ~= 1
        error('particle_stag_settling.CSV must contain exactly one row for %s.', partinfo);
    end

    d_p   = pInfo.d(case_mask)     * 1e-6;   % [m]
    rho_p = pInfo.rho_p(case_mask) * 1000;   % [kg/m3]
    Re_0  = pInfo.Re(case_mask);             % [-], based on still-water settling

    SchNau = 1 + 0.15*Re_0.^0.687;

    tau_p_part(pI) = d_p^2*(rho_p/rho_w + 0.5) ...
                   / (18*nu*SchNau);         % [s]
end

wave_cond = readtable('waveInfo.CSV');
H0 = wave_cond.H;

H = zeros(wN,pN);
H(1,:) = H0(1); 
H(2,:) = H0(2); 
H(3,:) = H0(3);

for wI = 1:wN
  wave = waves{wI};
  for pI = 1:pN
    partinfo = parts{pI};
    filepath = [wave partinfo];

    S = load(fullfile(filepath, sprintf('%s%s_vh_avg_unit.mat', wave, partinfo)));

    % Phase vs slip (if needed)
    thetas_std{wI,pI} = S.thetas_std;     % phi_std [rad]
    r_avg_u{wI,pI}    = S.r_avg_u;        % [m/s]
    r_avg_v{wI,pI}    = S.r_avg_v;        % [m/s]
    r_se_u{wI,pI}     = S.r_se_u;         % [m/s]
    r_se_v{wI,pI}     = S.r_se_v;         % [m/s]

    % Parameters
    v_s(wI,pI)        = S.v_s;
    k(wI,pI)          = S.k;
    h(wI,pI)          = S.h;
    T(wI,pI)          = S.T;
    Ec(wI,pI)         = S.E_cosh_ens_out;
    Es(wI,pI)         = S.E_sinh_ens_out;

    % Wave-averaged location (already saved in file)
    kz_bin_z_mean{wI,pI} = S.kz_w_bin_z_mean;   % <k y_p>
    kz_bin_z_std{wI,pI}  = S.kz_w_bin_z_std;
    kz_bin_v_mean{wI,pI} = S.kz_w_bin_v_mean;   % <w_slip> [m/s]
    kz_bin_v_std{wI,pI}  = S.kz_w_bin_v_std;    % std(w_slip) [m/s]

  end
end

fig_waveavg = figure('Color','w');
ax_waveavg  = gobjects(1,9);
for sp = 1:9
    ax_waveavg(sp) = subplot(3,3,sp);
    if sp == 3
        axis(ax_waveavg(sp),'off');   % legend only
    else
        hold(ax_waveavg(sp),'on');
        grid(ax_waveavg(sp),'on');
        set(ax_waveavg(sp),'FontSize',16);
        xlabel(ax_waveavg(sp),'$-w_z/v_s$','Interpreter','latex');
        ylabel(ax_waveavg(sp),'$z_p$','Interpreter','latex');  % actual values are k y_p
        xlim(ax_waveavg(sp),[0.6 1.5]);
    end
end

% Titles for each particle subplot
for pI = 1:pN
    axp = ax_waveavg(slotForPart(pI));
    title(axp, parts{pI}, 'Interpreter','none');
end

% Plot wave-averaged location data (no trial/ensemble separation; only wave is distinguished)
for pI = 1:pN
    axp = ax_waveavg(slotForPart(pI));
    for wI = 1:wN
        ky_m   = kz_bin_z_mean{wI,pI};   % <k y_p>
        ky_s   = kz_bin_z_std{wI,pI};    % std(k y_p)
        v_m    = kz_bin_v_mean{wI,pI};   % <w_slip> [m/s]
        v_sdev = kz_bin_v_std{wI,pI};    % std(w_slip) [m/s]

        if isempty(ky_m) || isempty(v_m)
            continue;
        end

        m = isfinite(ky_m)   & isfinite(ky_s)   & ...
            isfinite(v_m)    & isfinite(v_sdev) & ...
            (ky_s   > 0)     & (v_sdev > 0);

        if ~any(m)
            continue;
        end

        vs_term = v_s(wI,pI);   % settling velocity for this (wave,part)
        if ~isfinite(vs_term) || vs_term == 0
            continue;
        end

        % x: <w_slip>/v_s, y: <k y_p> (labeled z_p in the paper)
        x  = v_m(m)    ./ vs_term;
        dx = v_sdev(m) ./ vs_term;
        y  = ky_m(m);
        dy = ky_s(m);

        % Center points (markers)
        plot(axp, x, y, ...
            m_wI{wI}, ...
            'MarkerSize', 6, ...
            'MarkerFaceColor','none', ...
            'Color', c_wI{wI}, ...
            'LineStyle','none');

        % Error bars: horizontal/vertical
        for ii = 1:numel(x)
            line(axp, [x(ii)-dx(ii), x(ii)+dx(ii)], [y(ii), y(ii)], ...
                 'Color', c_wI{wI}, 'LineWidth', 1.0);
            line(axp, [x(ii), x(ii)], [y(ii)-dy(ii), y(ii)+dy(ii)], ...
                 'Color', c_wI{wI}, 'LineWidth', 1.0);
        end

        % Reference line at w_slip = v_s
        xline(axp, 1, 'k-', 'HandleVisibility','off');
    end
end

% Legend in the 3rd slot (wave-specific color/marker only)
ax_leg2 = ax_waveavg(3);
axes(ax_leg2); cla(ax_leg2); axis(ax_leg2,'off'); hold(ax_leg2,'on');
H_wave2 = gobjects(1,wN);
for wI = 1:wN
    H_wave2(wI) = plot(ax_leg2, nan, nan, m_wI{wI}, ...
        'MarkerSize', 10, ...
        'MarkerFaceColor','none', ...
        'Color', c_wI{wI}, ...
        'LineStyle','none', ...
        'DisplayName', sprintf('$\\langle \\widetilde{w_z} (z_p,t) \\rangle_{[z]} \\, \\mathrm{(%s)}$',waves{wI}));
end
legend(ax_leg2, H_wave2,'FontSize',18,'Interpreter','latex');


labels = 'ab-cdefgh';

for sp = 1:9
    if sp == 3, continue; end
    ax = ax_waveavg(sp);

    text(ax,-0.18,1.02,sprintf('(%c)',labels(sp)), 'Units','normalized',...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Clipping','off')
end

fig_waveavg.Position = [-1000 -500 1100 1200];
% exportgraphics(fig_waveavg,'E:\2.ParticlesInWave\Manuscript\figures\wave_avg_profiles.pdf')

%% ----------------- Wave-averaged vertical profiles: particle velocity (fig 6)------
clear all; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};
parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};

wN = numel(waves);
pN = numel(parts);
slotForPart = [1 2 4 5 6 7 8 9];

c_wI = { ...
    [0 0 0]
    [0.698 0.1333 0.1333]
    [0 0.5451 0.5451]
    };
m_wI = {'o','^','s'};

ciFactor = 1.96;   % saved *_std are actually bootstrap SE in your processing code

% load tables
pInfo = readtable('particle_stag_settling.CSV');
wave_cond = readtable('waveInfo.CSV');

nu    = waterNu(24.5);
rho_w = waterRho(24.5);

% storage
kz_bin_z_mean  = cell(wN,pN);
kz_bin_z_std   = cell(wN,pN);
kz_bin_v_mean  = cell(wN,pN);
kz_bin_v_std   = cell(wN,pN);
kz_bin_vp_mean = cell(wN,pN);
kz_bin_vp_std  = cell(wN,pN);

v_s = zeros(wN,pN);
k   = zeros(wN,pN);
T   = zeros(wN,pN);
h   = zeros(wN,pN);
H   = zeros(wN,pN);

for wI = 1:wN
    wave = waves{wI};

    rowW = strcmp(wave_cond.Case, wave);
    H_this = wave_cond.H(rowW);
    T_this = wave_cond.T(rowW);
    h_this = wave_cond.h(rowW);

    for pI = 1:pN
        partinfo = parts{pI};
        filepath = [wave partinfo];

        S = load(fullfile(filepath, sprintf('%s%s_vh_avg_unit.mat', wave, partinfo)));

        assert(isfield(S,'kz_w_bin_vp_mean') && isfield(S,'kz_w_bin_vp_std'), ...
            'MAT file %s does not contain v_p-based wave-averaged fields. Re-run the processing script first.', ...
            sprintf('%s%s_vh_avg_unit.mat', wave, partinfo));

        v_s(wI,pI) = S.v_s;
        k(wI,pI)   = S.k;
        T(wI,pI)   = S.T;
        h(wI,pI)   = S.h;
        H(wI,pI)   = H_this;

        kz_bin_z_mean{wI,pI}  = S.kz_w_bin_z_mean;
        kz_bin_z_std{wI,pI}   = S.kz_w_bin_z_std;
        kz_bin_v_mean{wI,pI}  = S.kz_w_bin_v_mean;   % stored downward-positive = -w_z in manuscript sign
        kz_bin_v_std{wI,pI}   = S.kz_w_bin_v_std;    % bootstrap SE
        kz_bin_vp_mean{wI,pI} = S.kz_w_bin_vp_mean;  % stored downward-positive = -v_z in manuscript sign
        kz_bin_vp_std{wI,pI}  = S.kz_w_bin_vp_std;   % bootstrap SE
    end
end

inset_xl = 0.995; inset_xr = 1.005; 
inset_yl = -0.5;  inset_yr = 0; 

% ------------------------------------------------------------------------
% -<v_z>/v_s  + D22 vz-drift
% ------------------------------------------------------------------------
fig_vz = figure('Color','w');
ax_vz = gobjects(1,9);

for sp = 1:9
    ax_vz(sp) = subplot(3,3,sp);
    if sp == 3
        axis(ax_vz(sp),'off');
    else
        hold(ax_vz(sp),'on');
        grid(ax_vz(sp),'on');
        set(ax_vz(sp),'FontSize',16);
        xlabel(ax_vz(sp),'$- v_z / v_s$','Interpreter','latex');
        ylabel(ax_vz(sp),'$z_p$','Interpreter','latex');
        xlim(ax_vz(sp), [0.8 1.2]);
        ylim(ax_vz(sp), [-2 0]);
        % xlim(ax_vz(sp), [0.995 1.005]);
        % ylim(ax_vz(sp), [-0.5 0]);
    end
end

for pI = 1:pN
    axp = ax_vz(slotForPart(pI));
    title(axp, parts{pI}, 'Interpreter','none');
end

for pI = 1:pN
    partinfo = parts{pI};

    rowP = strcmp(pInfo.Case, partinfo);
    d_p   = pInfo.d(rowP) * 1e-6;
    rho_p = pInfo.rho_p(rowP) * 1000;

    for wI = 1:wN
        wave = waves{wI};

        ky_m   = kz_bin_z_mean{wI,pI};
        ky_se  = kz_bin_z_std{wI,pI};
        vp_m   = kz_bin_vp_mean{wI,pI};   % stored downward-positive = -v_z
        vp_se  = kz_bin_vp_std{wI,pI};    % bootstrap SE

        if isempty(ky_m) || isempty(vp_m)
            continue
        end

        m = isfinite(ky_m) & isfinite(ky_se) & isfinite(vp_m) & isfinite(vp_se);
        if ~any(m)
            continue
        end

        axp = ax_vz(slotForPart(pI));

        vs_term = v_s(wI,pI);
        k_case  = k(wI,pI);
        T_case  = T(wI,pI);
        h_case  = h(wI,pI);
        H_case  = H(wI,pI);

        % stored vp is downward-positive, i.e. vp = -v_z
        % therefore -<v_z>/v_s = <vp>/v_s
        x  =  vp_m(m) ./ vs_term;
        dx =  ciFactor * vp_se(m) ./ vs_term;
        y  =  ky_m(m);
        dy =  ciFactor * ky_se(m);

        % same marker style as Figure A
        plot(axp, x, y, m_wI{wI}, ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', 'none', ...
            'Color', c_wI{wI}, ...
            'LineStyle', 'none');

        for ii = 1:numel(x)
            line(axp, [x(ii)-dx(ii), x(ii)+dx(ii)], [y(ii), y(ii)], ...
                'Color', c_wI{wI}, 'LineWidth', 1.0);
            line(axp, [x(ii), x(ii)], [y(ii)-dy(ii), y(ii)+dy(ii)], ...
                'Color', c_wI{wI}, 'LineWidth', 1.0);
        end

        % dense theory grid in plotted ordinate (currently kz_p)
        y_theory = linspace(-2, 0, 400);     % plotted vertical coordinate
        z_line   = y_theory ./ k_case;      % dimensional z [m] for theory evaluation

        [~, vz_ratio_d22] = profileD22SlipVelocity( ...
            z_line, H_case, T_case, h_case, d_p, rho_p, rho_w, nu, vs_term);
        
        % D22 particle velocity: convert <v_z>/v_s -> -<v_z>/v_s
        plot(axp, -vz_ratio_d22, y_theory, '-', ...
            'Color', c_wI{wI}, 'LineWidth', 1.8);
        xline(axp, 1, 'k-', 'HandleVisibility','off');
        set(axp,'XScale','log')
    end
end

ax_leg2 = ax_vz(3);
axes(ax_leg2); cla(ax_leg2); axis(ax_leg2,'off'); hold(ax_leg2,'on');
Hleg2 = gobjects(1,2*wN);
for wI = 1:wN
    Hleg2(2*wI-1) = plot(ax_leg2, nan, nan, m_wI{wI}, ...
        'MarkerSize', 8, 'MarkerFaceColor','none', 'Color', c_wI{wI}, ...
        'LineStyle','none', ...
        'DisplayName', sprintf('$\\widetilde{v_z}(z_p,t)\\,\\mathrm{(%s)}$',waves{wI}));

    Hleg2(2*wI) = plot(ax_leg2, nan, nan, '-', ...
        'Color', c_wI{wI}, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('$v_{z-drift}$ (%s)', waves{wI}));
end
legend(ax_leg2, Hleg2, 'FontSize', 13, 'Location', 'northwest', 'Interpreter','latex');

labels = 'ab-cdefgh';
for sp = 1:9
    if sp == 3, continue; end
    ax = ax_vz(sp);
    text(ax, -0.18, 1.02, sprintf('(%c)', labels(sp)), ...
        'Units','normalized', 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',18, 'Clipping','off');
end


fig_vz.Position = [50 -500 1200 1200];


% exportgraphics(fig_vz, 'waveavg_vz_profile_with_Lagrangian_and_D22.pdf');

%% ===== Phase-ensemble + bootstrap-SE errorbar (figure 7) =====
clear; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};

% pick ONE particle
partinfo = 'G1';   % <-- change as needed

B_boot = 10000;       % bootstrap reps
nbin   = 180;         % phase bins

% wave colors (same palette)
waveColorMap.W1 = [0 0 0];
waveColorMap.W2 = [0.698 0.1333 0.1333];
waveColorMap.W3 = [0 0.5451 0.5451];

lighten = @(c,a) 1 - a*(1 - c);   % a=0 original, a=1 white
alpha_trial = 0.3;               % larger -> lighter trial curves

% ----------------- PHASE BINNING SETUP -----------------
edges_phi = linspace(0, 2*pi, nbin+1);
cent_phi  = 0.5*(edges_phi(1:end-1) + edges_phi(2:end));
phi_deg   = rad2deg(cent_phi);

wrap2pi = @(x) mod(x, 2*pi);

% ----------------- LOAD wave info (for k, T) -----------------
wave_cond = readtable('waveInfo.CSV');

% containers for plotting (per wave)
U_mean = cell(1,numel(waves));  W_mean = U_mean;
U_se   = U_mean;                W_se   = U_mean;
U_cnt  = U_mean;                W_cnt  = U_mean;

% also keep all trial curves for plotting
U_trial = cell(1,numel(waves));
W_trial = cell(1,numel(waves));

% ----------------- MAIN LOOP: W1/W2/W3 -----------------
for wI = 1:numel(waves)
    wave = waves{wI};

    % wave parameters
    T = wave_cond.T(strcmp(wave_cond.Case,wave));
    h = wave_cond.h(strcmp(wave_cond.Case,wave));
    if isempty(T) || isempty(h)
        warning('Wave %s not found in waveInfo.CSV. Skip.', wave);
        continue;
    end
    k = dispersion(h,T);  % user-provided
    omega = 2*pi/T;
    c(wI) = omega/k;

    % load trials
    filepath = fullfile([wave, partinfo]);
    flist = dir(fullfile(filepath, 'slip_vh_results_*.mat'));

    case_tot = numel(flist);
    if case_tot == 0
        warning('No slip_dn_results_*.mat for %s%s', wave, partinfo);
        continue;
    end

    % prealloc per-trial binned means + per-bin sample counts
    r_trial_u = nan(case_tot, nbin);   % trial-mean in each phase bin
    r_trial_w = nan(case_tot, nbin);
    n_samp_u  = zeros(case_tot, nbin); % raw sample counts (frames) in each bin
    n_samp_w  = zeros(case_tot, nbin);

    % determine v_s per wave/particle (mean over trials)
    v_s_all = nan(1,case_tot);

    for j = 1:case_tot
        S = load(fullfile(filepath, flist(j).name));

        % required
        u_slip = S.u_slip(:);
        v_slip = S.v_slip(:);
        v_s_all(j) = S.meta.v_s;

        % time and dt
        if isfield(S.meta,'dt') && isfinite(S.meta.dt)
            dt = S.meta.dt;
        else
            dt = 1/fps;
        end
        if isfield(S,'time')
            t = S.time(:);
        else
            t = (0:numel(u_slip)-1)'*dt;
        end

        % common orbital phase, already in theory convention
        if ~isfield(S,'orb_phi_f')
            warning('Missing orb_phi_f in %s. Re-run slip generation with phi_orb. Skip trial %d', ...
                flist(j).name, j);
            continue
        end
        
        phi_std = wrap2pi(S.orb_phi_f(:));   % [0,2π)

        % normalize by v_s (use per-wave mean later; here temporary)
        % (for trimming we only need v_slip/v_s; use trial v_s to avoid bias)
        vs_j = v_s_all(j);
        if ~isfinite(vs_j) || vs_j==0, continue; end
        vs_norm = v_slip ./ vs_j;

        % --------- TRIMMING + OUTLIER REMOVAL (use function) ----------
        [idx_keep, ~, ~] = trimTrendSigmaFilter(vs_norm(:), t(:));
        
        if isempty(idx_keep)
            warning('Filtering failed in %s (trial %d). Skip.', flist(j).name, j);
            continue
        end
        
        % trimmed + outlier-removed series
        th = phi_std(idx_keep);
        uu = u_slip(idx_keep);   % [m/s]
        ww = v_slip(idx_keep);   % [m/s]

        % --------- PHASE BINNING: per-trial mean + per-bin sample count ----------
        b = discretize(th, edges_phi);

        for kbin = 1:nbin
            mk = (b == kbin) & isfinite(uu) & isfinite(th);
            if any(mk)
                r_trial_u(j,kbin) = mean(uu(mk), 'omitnan');  % [m/s]
                n_samp_u(j,kbin)  = sum(mk);
            end

            mk = (b == kbin) & isfinite(ww) & isfinite(th);
            if any(mk)
                r_trial_w(j,kbin) = mean(ww(mk), 'omitnan');  % [m/s]
                n_samp_w(j,kbin)  = sum(mk);
            end
        end
    end

    % use wave/particle v_s as mean over trials
    v_s = mean(v_s_all,'omitnan');
    if ~isfinite(v_s) || v_s==0
        warning('Invalid v_s for %s %s', wave, partinfo);
        continue;
    end

    % convert to nondim (trial bin means)
    r_trial_u = r_trial_u ./ v_s;
    r_trial_w = r_trial_w ./ v_s;

    % store trial curves for plotting
    U_trial{wI} = r_trial_u;
    W_trial{wI} = r_trial_w;

    % ---------- ENSEMBLE MEAN ----------
    U_mean{wI} = mean(r_trial_u, 1, 'omitnan');
    W_mean{wI} = mean(r_trial_w, 1, 'omitnan');

    % ---------- SAMPLE COUNTS (sum over trials; all combined) ----------
    U_cnt{wI} = sum(n_samp_u, 1, 'omitnan');
    W_cnt{wI} = sum(n_samp_w, 1, 'omitnan');

    % ---------- BOOTSTRAP SE of the MEAN (across trials) ----------
    U_se{wI} = nan(1,nbin);
    W_se{wI} = nan(1,nbin);

    for kbin = 1:nbin
        vals = r_trial_u(:,kbin); vals = vals(isfinite(vals));
        if numel(vals) >= 2
            bootMeans = bootstrp(B_boot, @mean, vals);
            U_se{wI}(kbin) = 1.96*std(bootMeans, 0);
        end

        vals = r_trial_w(:,kbin); vals = vals(isfinite(vals));
        if numel(vals) >= 2
            bootMeans = bootstrp(B_boot, @mean, vals);
            W_se{wI}(kbin) = 1.96*std(bootMeans, 0);
        end
    end
end

% PLOT
fig = figure('Color','w','Name',sprintf('Phase ensemble (%s)',partinfo));

idxEb = 1:10:nbin;   % errorbar spacing (every 10 bins)

% ----------------- (3,2,1:2) w_x/v_s -----------------
ax1 = subplot(3,2,1:2); hold(ax1,'on'); grid(ax1,'on');

% trials (light)
for wI = 1:numel(waves)
    wave = waves{wI};
    if ~isfield(waveColorMap, wave) || isempty(U_trial{wI}), continue; end
    c_wave  = waveColorMap.(wave);
    c_trial = lighten(c_wave, alpha_trial);

    Ut = U_trial{wI};
    for j = 1:size(Ut,1)
        if j==1
        p1(wI) = plot(ax1, phi_deg, Ut(j,:), '--', 'Color', c_trial, ...
            'LineWidth', 0.8, 'DisplayName', sprintf('$\\,\\, w_x(z_p,\\phi)\\quad \\mathrm{(%s)}$', wave));
        else
        plot(ax1, phi_deg, Ut(j,:), '--', 'Color', c_trial, ...
            'LineWidth', 0.8, 'HandleVisibility','off');
        end
    end
end

% ensemble + bootstrap SE errorbar (every 10 bins)
for wI = 1:numel(waves)
    wave = waves{wI};
    if isempty(U_mean{wI}), continue; end

    em = U_mean{wI}; se = U_se{wI};
    c_wave  = waveColorMap.(wave);

    plot(ax1, phi_deg, em, '-', 'Color', c_wave, 'LineWidth', 1.5, 'HandleVisibility','off');

    p2(wI) = errorbar(ax1, phi_deg(idxEb), em(idxEb), se(idxEb), 'o', ...
        'Color', c_wave, 'MarkerFaceColor', c_wave, ...
        'LineWidth', 1.5, 'CapSize', 6, ...
        'DisplayName', sprintf('$\\langle w_x(z_p,\\phi) \\rangle_{[\\phi]} \\, \\mathrm{(%s)}$', wave));
end
xlim(ax1,[0 360]);
xlabel(ax1,'$\phi\,[\,^\circ\,]$','Interpreter','latex');
ylabel(ax1,'$w_x/v_s$','Interpreter','latex');
set(ax1,'FontSize',14);
legend(ax1,[p2(1), p1(1), p2(2), p1(2), p2(3), p1(3)],'Location','northoutside','Interpreter','latex','NumColumns',3);

% ----------------- (3,2,3:4) w_z/v_s -----------------
ax2 = subplot(3,2,3:4); hold(ax2,'on'); grid(ax2,'on');

% trials (light)
for wI = 1:numel(waves)
    wave = waves{wI};
    if ~isfield(waveColorMap, wave) || isempty(W_trial{wI}), continue; end
    c_wave  = waveColorMap.(wave);
    c_trial = lighten(c_wave, alpha_trial);

    Wt = W_trial{wI};
    for j = 1:size(Wt,1)
        if j==1
        p1(wI) = plot(ax2, phi_deg, Wt(j,:), '--', 'Color', c_trial, ...
            'LineWidth', 0.8, 'DisplayName', sprintf('$\\,\\, w_z(z_p,\\phi) \\quad \\mathrm{(%s)}$', wave));
        else
        plot(ax2, phi_deg, Wt(j,:), '--', 'Color', c_trial, ...
            'LineWidth', 0.8, 'HandleVisibility','off');
        end
    end
end

% ensemble + bootstrap SE errorbar (every 10 bins)
for wI = 1:numel(waves)
    wave = waves{wI};
    if isempty(W_mean{wI}), continue; end

    em = W_mean{wI}; se = W_se{wI};
    c_wave  = waveColorMap.(wave);

    plot(ax2, phi_deg, em, '-', 'Color', c_wave, 'LineWidth', 1.5, 'HandleVisibility','off');

    p2(wI) = errorbar(ax2, phi_deg(idxEb), em(idxEb), se(idxEb), 'o', ...
        'Color', c_wave, 'MarkerFaceColor', c_wave, ...
        'LineWidth', 1.5, 'CapSize', 6, ...
        'DisplayName', sprintf('$\\langle w_z(z_p,\\phi) \\rangle_{[\\phi]} \\, \\mathrm{(%s)}$', wave));
end
xlim(ax2,[0 360]);
xlabel(ax2,'$\phi\,[\,^\circ\,]$','Interpreter','latex');
ylabel(ax2,'$-w_z/v_s$','Interpreter','latex');
set(ax2,'FontSize',14);
legend(ax2,[p2(1), p1(1), p2(2), p1(2), p2(3), p1(3)],'Location','northoutside','Interpreter','latex','NumColumns',3);


% ----------------- (3,2,5) sample counts for w_x -----------------
ax3 = subplot(3,2,5); hold(ax3,'on'); grid(ax3,'on');
for wI = 1:numel(waves)
    wave = waves{wI};
    if ~isfield(waveColorMap, wave) || isempty(U_cnt{wI}), continue; end
    plot(ax3, phi_deg, U_cnt{wI}, '.-', 'Color', waveColorMap.(wave), 'LineWidth', 1.0);
end
xlim(ax3,[0 360]); ylim([0 100])
xlabel(ax3,'$\phi\,[\,^\circ\,]$','Interpreter','latex');
ylabel(ax3,'$N_{\phi}$','Interpreter','latex');
set(ax3,'FontSize',14);
legend({'W1','W2','W3'},'NumColumns',3)

% --------- (3,2,6) polarplot of w_z/v_s + AREA centroid (per wave) ----------
% NOTE: area centroid assumes r(phi) >= 0 for "enclosed area".
% If some bins give negative r, we clip to 0 for centroid calculation only.
m_wave = {'o','^','s'};

ax6 = subplot(3,2,6);
pos6 = get(ax6,'Position');
delete(ax6);

axp = polaraxes('Position', pos6); hold(axp,'on');
set(axp,'ThetaZeroLocation','right','ThetaDir','counterclockwise');
axp.FontSize = 14;

for wI = 1:numel(waves)
    wave = waves{wI};
    if ~isfield(waveColorMap, wave) || isempty(W_mean{wI}), continue; end

    c_wave = waveColorMap.(wave);

    th = cent_phi(:);
    r  = W_mean{wI}(:);
    m  = isfinite(th) & isfinite(r);
    th = th(m); r = r(m);

    % 1) polar curve (ensemble mean)
    polarplot(axp, th, (r-min(r))/c(wI), '-', 'Color', c_wave, 'LineWidth', 1, ...
        'DisplayName', sprintf('$r_z\\,\\mathrm{(%s)}$', wave));

    % 2) AREA centroid (using r_clipped >= 0 for area definition)
    rA = max(r - min(r), 0);  % centroid calc only
    A  = 0.5 * trapz(th, rA.^2);

    if isfinite(A) && A > 0
        xc = (1/(6*A)) * trapz(th, (rA.^3).*cos(th));
        yc = (1/(6*A)) * trapz(th, (rA.^3).*sin(th));
        
        theta_c = atan2(yc, xc);
        theta_c = mod(theta_c, 2*pi);
        r_c     = hypot(xc, yc);
        
        % centroid point
        polarplot(axp, theta_c, r_c/c(wI), m_wave{wI}, ...
            'Color', c_wave, 'MarkerSize', 8, 'MarkerFaceColor', c_wave, ...
            'DisplayName',sprintf('$\\mathbf{x_c}\\,\\mathrm{(%s)}$',wave));
    end
end
title('$r_z = \langle w_z \rangle_{[\phi]} - w_{z,\min}$','Interpreter','latex')
legend(axp,'Location','southoutside','Interpreter','latex','NumColumns',3);

% --------- panel labels (a)-(d) ----------
% ax1, ax2 were created above; subplot(3,2,5) is current axes handle from subplot call.
ax5 = subplot(3,2,5);  % get handle to panel (c)
% Put labels at upper-left of each panel (normalized coords)
text(ax1, -0.06, 1.05, '(a)', 'Units','normalized', 'FontSize', 16, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom');
text(ax2, -0.06, 1.05, '(b)', 'Units','normalized', 'FontSize', 16, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom');
text(ax5, -0.15, 1.05, '(c)', 'Units','normalized', 'FontSize', 16, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom');
text(axp, -0.15, 1.05, '(d)', 'Units','normalized', 'FontSize', 16, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom');

fig.Position = [-1000 -500 900 1500];
% exportgraphics(fig,'E:\2.ParticlesInWave\Manuscript\figures\phase_ensemble_G1.pdf')


%% ========================= Force ratio (figure 9) =======================
clear; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};
parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};

% --- Fluid properties and particle info ---
rho_w = waterRho(24.5);
g     = 9.81;
nu    = waterNu(24.5);
mu    = rho_w * nu;
C_m   = 0.5;

movwin = 125;

pInfo     = readtable('particle_stag_settling.CSV');
wave_cond = readtable('waveInfo.CSV');

% --- Base colors per wave ---
waveColorMap.W1 = [0     0       0      ];
waveColorMap.W2 = [0.698 0.1333  0.1333];
waveColorMap.W3 = [0     0.5451  0.5451];

% --- Visualization parameters ---
alpha_trial = 0.30;   % line transparency
lw_trial    = 0.8;    % trial line width
lw_ref      = 1.6;    % reference (5*omega*tau_p) line width

% === Create figure ===
figR = figure('Name','R_D & R_{wave} - all particles');
set(figR,'Position',[50 50 1200 900]);
set(figR,'Color','w');
set(figR,'Renderer','opengl');   % changed from 'painters' to enable transparency

nPart = length(parts);
nWave = length(waves);

for pI = 1:nPart
    partinfo = parts{pI};

    case_mask = strcmp(pInfo.Case, partinfo);
    d_p   = pInfo.d(case_mask) * 1e-6;
    rho_p = pInfo.rho_p(case_mask)*1000;
    gamma = rho_p / rho_w;
    V_p   = pi/6 * d_p^3;
    W_sub = (gamma - 1) * rho_w * g * V_p;
    vs    = pInfo.v_s(case_mask)/1000;
    Re_ps = vs*d_p/nu;
    tau_p = d_p^2*(gamma + C_m)/18/nu/(1+0.15*Re_ps.^0.687);

    for wI = 1:nWave
        wave = waves{wI};

        mask_w = strcmp(wave_cond.Case, wave);
        H = wave_cond.H(mask_w);
        T = wave_cond.T(mask_w);
        h = wave_cond.h(mask_w);

        omega = 2*pi/T;
        k     = dispersion(h,T);

        filepath = [wave partinfo];
        flist = dir([filepath '\slip_vh_results_*']);

        case_tot = length(flist);
        if case_tot == 0
            continue;
        end

        baseCol = waveColorMap.(wave);

        idx_RD    = 2*(pI-1) + 1;
        idx_Rwave = 2*(pI-1) + 2;

        % ----------------- R_D(t) -----------------
        ax_RD = subplot(8,2,idx_RD,'Parent',figR); hold(ax_RD,'on');

        for j = 1:case_tot
            S = load(fullfile(filepath, sprintf('slip_vh_results_%d.mat',j)));

            t    = S.time(:);
            u_sl = S.u_slip(:);
            v_sl = S.v_slip(:);

            t0   = t(1);
            t_nd = omega * (t - t0);

            Urel = hypot(u_sl,v_sl);
            Re_t  = Urel * d_p / nu;
            SchN  = 1 + 0.15 * Re_t.^0.687;
            F_D_t = 3*pi*mu*d_p .* SchN .* v_sl;

            R_D_t = abs(F_D_t) / W_sub;

            colj = trialColor(baseCol, j, case_tot);
            hp = plot(ax_RD, t_nd, movmean(R_D_t,movwin), ...
                      'Color', colj, 'LineWidth', lw_trial);
            if wI == 1
                hp.Color(4) = 0.6;
            else
                hp.Color(4) = alpha_trial;
            end
        end

        xline(ax_RD, 5*omega*tau_p, '--', 'Color', baseCol, 'LineWidth', lw_ref);

        grid(ax_RD,'on')
        if pI == nPart
            xlabel(ax_RD,'$\omega t$','Interpreter','latex')
        end
        if wI == nWave
            ylabel(ax_RD, sprintf('%s\n$\\mathcal{R}_{\\mathrm{drag}}$', partinfo), 'Interpreter','latex')
            ylim(ax_RD,[0 1.5])
        end

        % ----------------- R_wave(t) -----------------
        ax_RW = subplot(8,2,idx_Rwave,'Parent',figR); hold(ax_RW,'on');

        for j = 1:case_tot
            S = load(fullfile(filepath, sprintf('slip_vh_results_%d.mat',j)));

            t  = S.time(:);
            vf = S.v_f(:);

            t0   = t(1);
            t_nd = omega * (t - t0);

            N = numel(t);
            if N < 3
                continue;
            end

            dt = median(diff(t));

            dvfdt = nan(N,1);
            for n = 3:N-2
                dvfdt(n) = (-vf(n+2) + 8*vf(n+1) - 8*vf(n-1) + vf(n-2)) / (12*dt);
            end
            dvfdt(2)   = (vf(3) - vf(1)) / (2*dt);
            dvfdt(N-1) = (vf(N) - vf(N-2)) / (2*dt);
            dvfdt(1)   = (vf(2) - vf(1)) / dt;
            dvfdt(N)   = (vf(N) - vf(N-1)) / dt;

            R_wave_t = (1 + C_m) * abs(dvfdt) / ((gamma - 1) * g);

            colj = trialColor(baseCol, j, case_tot);
            hp = plot(ax_RW, t_nd, movmean(R_wave_t,movwin), ...
                      'Color', colj, 'LineWidth', lw_trial);
            if wI == 1
                hp.Color(4) = 0.6;
            else
                hp.Color(4) = alpha_trial;
            end
        end

        xline(ax_RW, 5*omega*tau_p, '--', 'Color', baseCol, 'LineWidth', lw_ref);

        grid(ax_RW,'on')
        if pI == nPart
            xlabel(ax_RW,'$t$','Interpreter','latex')
        end
        ylabel(ax_RW,'$\mathcal{R}_{\mathrm{wave}}$','Interpreter','latex')
        ylim(ax_RW,[0 0.15])
    end
end

set(findall(figR,'Type','axes'),'FontSize',18)

% --- Legend (opaque reference lines, full baseCol) ---
ax_legend = subplot(8,2,2); hold(ax_legend,'on');
p1 = plot(ax_legend, nan,nan, '-', 'Color', waveColorMap.W1, 'LineWidth',2, 'DisplayName','W1');
p2 = plot(ax_legend, nan,nan, '-', 'Color', waveColorMap.W2, 'LineWidth',2, 'DisplayName','W2');
p3 = plot(ax_legend, nan,nan, '-', 'Color', waveColorMap.W3, 'LineWidth',2, 'DisplayName','W3');
legend(ax_legend,[p1 p2 p3],'Location','best','NumColumns',3);

% --- Save as PDF ---
outpdf = 'E:\2.ParticlesInWave\Manuscript\figures\force_ratio.pdf';
% exportgraphics(figR, outpdf, 'ContentType','image', 'Resolution', 600);

%% ===== r_z cardioids (figure 10) =====
clear; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

% ----------------- USER SETTINGS -----------------
waves = {'W1','W2','W3'};

parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};

nbin = 180;           % phase bins

waveColorMap.W1 = [0 0 0];
waveColorMap.W2 = [0.698 0.1333 0.1333];
waveColorMap.W3 = [0 0.5451 0.5451];

% ----------------- PHASE BINNING SETUP -----------------
edges_phi = linspace(0, 2*pi, nbin+1);
cent_phi  = 0.5*(edges_phi(1:end-1) + edges_phi(2:end));

wrap2pi = @(x) mod(x, 2*pi);

% ----------------- LOAD wave info (for T, h -> trimming window) -----------------
wave_cond = readtable('waveInfo.CSV');

% ----------------- FIGURE + SUBPLOT LAYOUT -----------------
fig = figure('Color','w','Name','Polarplot (all particles)');
fig.Position = [-1000 -500 1100 1200];

% subplot index mapping: reserve (3,3,3) for legend
spOrder = [1 2 4 5 6 7 8 9];

% centroid markers per wave
m_wave = {'o','^','s'};

% ----------------- LEGEND PANEL (3,3,3) -----------------
ax_leg = subplot(3,3,3);
cla(ax_leg); axis(ax_leg,'off'); hold(ax_leg,'on');

Hleg = gobjects(1, 2*numel(waves));
legCnt = 0;
for wI = 1:numel(waves)
    wave = waves{wI};
    c_wave = waveColorMap.(wave);

    % line legend item (r_z)
    legCnt = legCnt + 1;
    Hleg(legCnt) = plot(ax_leg, nan, nan, '-', 'Color', c_wave, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('$r_z\\;\\mathrm{(%s)}$', wave));

    % % centroid legend item (x_c)
    % legCnt = legCnt + 1;
    % Hleg(legCnt) = plot(ax_leg, nan, nan, m_wave{wI}, 'Color', c_wave, ...
    %     'MarkerFaceColor', c_wave, 'MarkerSize', 8, 'LineStyle', 'none', ...
    %     'DisplayName', sprintf('$\\mathbf{x_c}\\;\\mathrm{(%s)}$', wave));

    % mean legend item (mean(r_z))
    legCnt = legCnt + 1;
    Hleg(legCnt) = plot(ax_leg, nan, nan, '--', 'Color', c_wave, ...
        'MarkerFaceColor', c_wave, 'MarkerSize', 8, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('$\\overline{r_z}\\;\\mathrm{(%s)}$', wave));
end
legend(ax_leg, Hleg, 'FontSize', 16, 'Interpreter', 'latex', 'Location', 'northwest');

% ----------------- MAIN LOOP: particles -> compute W_mean -> polarplot -----------------
for pI = 1:numel(parts)
    partinfo = parts{pI};

    % containers per particle
    W_mean  = cell(1,numel(waves));

    % --------- loop waves (W1/W2/W3) ----------
    for wI = 1:numel(waves)
        wave = waves{wI};

        % wave parameters (T,h) for trimming window
        T = wave_cond.T(strcmp(wave_cond.Case,wave));
        h = wave_cond.h(strcmp(wave_cond.Case,wave));
        if isempty(T) || isempty(h)
            warning('Wave %s not found in waveInfo.CSV. Skip.', wave);
            continue;
        end
        k = dispersion(h,T);
        omega = 2*pi/T;
        c(wI) = omega/k;

        % load trials
        filepath = [wave, partinfo];
        flist = dir(fullfile(filepath, 'slip_vh_results_*.mat'));
        case_tot = numel(flist);
        if case_tot == 0
            warning('No slip_vh_results_*.mat for %s %s', wave, partinfo);
            continue;
        end

        % prealloc per-trial binned means
        r_trial_w = nan(case_tot, nbin);
        n_samp_w  = zeros(case_tot, nbin);

        % v_s per trial -> use mean over trials
        v_s_all = nan(1,case_tot);

        for j = 1:case_tot
            S = load(fullfile(filepath, flist(j).name));

            if ~isfield(S,'u_slip') || ~isfield(S,'v_slip') || ~isfield(S,'meta') || ~isfield(S.meta,'v_s')
                warning('Missing fields in %s. Skip trial %d', flist(j).name, j);
                continue;
            end

            u_slip = S.u_slip(:);
            v_slip = S.v_slip(:);
            v_s_all(j) = S.meta.v_s;

            % dt, time
            if isfield(S.meta,'dt') && isfinite(S.meta.dt)
                dt = S.meta.dt;
            else
                dt = 1/fps;
            end
            if isfield(S,'time')
                t = S.time(:);
            else
                t = (0:numel(u_slip)-1)'*dt;
            end

            % common orbital phase, already in theory convention
            if ~isfield(S,'orb_phi_f')
                warning('Missing orb_phi_f in %s. Re-run slip generation with phi_orb. Skip trial %d', ...
                    flist(j).name, j);
                continue;
            end
            
            phi_std = wrap2pi(S.orb_phi_f(:));   % [0,2π)

            % temporary normalize for trimming (trial v_s)
            vs_j = v_s_all(j);
            if ~isfinite(vs_j) || vs_j==0, continue; end
            vs_norm = v_slip ./ vs_j;

            % --------- TRIMMING + OUTLIER REMOVAL (use function) ----------
            [idx_keep, ~, ~] = trimTrendSigmaFilter(vs_norm(:), t(:));
            
            if isempty(idx_keep)
                continue
            end
            
            % trimmed
            th = phi_std(idx_keep);
            ww = v_slip(idx_keep);

            % --------- PHASE BINNING: per-trial mean ----------
            b = discretize(th, edges_phi);
            for kbin = 1:nbin
                mk = (b == kbin) & isfinite(ww) & isfinite(th);
                if any(mk)
                    r_trial_w(j,kbin) = mean(ww(mk), 'omitnan');
                    n_samp_w(j,kbin)  = sum(mk);
                end
            end
        end

        % mean v_s over trials
        v_s = mean(v_s_all,'omitnan');
        vs(wI,pI) = v_s/c(wI);
        if ~isfinite(v_s) || v_s==0
            warning('Invalid v_s for %s %s', wave, partinfo);
            continue;
        end

        % nondim + ensemble mean (across trials)
        r_trial_w = r_trial_w ./ c(wI);
        W_mean{wI} = mean(r_trial_w, 1, 'omitnan');
    end

    % --------- PLOT: one polaraxes per particle ----------
    sp = spOrder(pI);
    ax_tmp = subplot(3,3,sp);
    pos = get(ax_tmp,'Position');
    delete(ax_tmp);

    axp = polaraxes('Position', pos); hold(axp,'on');
    ax_rz(sp) = axp;
    set(axp,'ThetaZeroLocation','right','ThetaDir','counterclockwise');
    axp.FontSize = 12;

    for wI = 1:numel(waves)
        wave = waves{wI};
        if isempty(W_mean{wI}), continue; end

        c_wave = waveColorMap.(wave);

        th = cent_phi(:);
        r  = W_mean{wI}(:);

        m = isfinite(th) & isfinite(r);
        th = th(m); r = r(m);
        if isempty(r), continue; end

        rmin = min(r); rmax = max(r);
        if ~(isfinite(rmin) && isfinite(rmax)) || (rmax - rmin) == 0
            continue;
        end

        % 1) polar curve (normalized like your original panel)
        rN = (r - rmin);
        polarplot(axp, th, rN, '-', 'Color', c_wave, 'LineWidth', 1.0);

        % % 2) AREA centroid (clip for area only)
        % rA = max(r - rmin, 0);
        % A  = 0.5 * trapz(th, rA.^2);
        % 
        % if isfinite(A) && A > 0
        %     xc = (1/(6*A)) * trapz(th, (rA.^3).*cos(th));
        %     yc = (1/(6*A)) * trapz(th, (rA.^3).*sin(th));
        % 
        %     theta_c = mod(atan2(yc, xc), 2*pi);
        %     r_c     = hypot(xc, yc);
        % 
        %     polarplot(axp, theta_c, r_c, m_wave{wI}, ...
        %         'Color', c_wave, 'MarkerSize', 7, 'MarkerFaceColor', c_wave, 'LineStyle','none');
        % end
        polarplot(axp, th, mean(r)*ones(1,length(th))-rmin,'--','Color',c_wave,'LineWidth',1.5)
    end

    title(axp, parts(pI), 'Interpreter','none');
end

labels = 'ab-cdefgh';

for sp = 1:9
    if sp == 3, continue; end              
    ax = ax_rz(sp);

    text(ax,-0.18,1.02,sprintf('(%c)',labels(sp)), 'Units','normalized',...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18,'Clipping','off')
end

% exportgraphics(fig,'E:\2.ParticlesInWave\Manuscript\figures\rz_cardioids.pdf')


%% ----------------- Spectra of slip-velocity fluctuations (fig. 14 @ App. B) ----------------

clear; clc;

scriptPath = matlab.desktop.editor.getActiveFilename;
cd(fileparts(scriptPath));

waves = {'W1','W2','W3'};
parts = {'M1','M2','G1','G2','G3','E1','E2','E3'};
pN    = numel(parts);

% 3x3 grid; slot 3 reserved for the legend.
slotMap     = [1, 2, 4, 5, 6, 7, 8, 9];
slotForPart = @(pI) slotMap(pI);

% Styling
c_wI = { ...
    [0 0 0]                % W1
    [0.698 0.1333 0.1333]  % W2
    [0 0.5451 0.5451]      % W3
    };

ls_u = '-';     % horizontal slip fluctuation
ls_v = '--';    % vertical slip fluctuation

% Signal-processing parameters (B&P-style)
minN     = 64;        % minimum samples per trial after trimming
nBoot    = 500;       % bootstrap iterations for CI
ciLevel  = 0.95;      % confidence level

% Physical parameters
g     = 9.81;
nu    = waterNu(24.5);
rho_w = waterRho(24.5);

pInfo     = readtable('particle_stag_settling.CSV');
wave_cond = readtable('waveInfo.CSV');

fig_spec = figure('Color','w');
ax_spec  = gobjects(1, 9);

% Pre-create the 3x3 subplot grid.
for sp = 1:9
    ax_spec(sp) = subplot(3, 3, sp);

    if sp == 3
        axis(ax_spec(sp), 'off');
    else
        hold(ax_spec(sp), 'on');
        grid(ax_spec(sp), 'on');
        box(ax_spec(sp), 'on');

        set(ax_spec(sp), 'FontSize', 14);

        xlabel(ax_spec(sp), '$f d_p/v_g''$', 'Interpreter','latex');
        ylabel(ax_spec(sp), '$P(f)/v_g''d_p$',  'Interpreter','latex');

        xlim(ax_spec(sp), [0 0.05]);
    end
end

for pI = 1:pN
    title(ax_spec(slotForPart(pI)), parts{pI}, 'Interpreter','none');
end

for pI = 1:pN

    partinfo = parts{pI};
    ax = ax_spec(slotForPart(pI));

    % Particle properties
    rowP = strcmp(pInfo.Case, partinfo);
    assert(any(rowP), 'Particle %s not found in particle_stag_settling.CSV.', partinfo);

    d_p   = pInfo.d(rowP)     * 1e-6;    % [m]
    rho_p = pInfo.rho_p(rowP) * 1000;    % [kg/m^3]

    v_g = sqrt((rho_p/rho_w - 1) * g * d_p);   % [m/s]

    for wI = 1:numel(waves)

        wave   = waves{wI};
        c_wave = c_wI{wI};

        % Wave properties
        rowW = strcmp(wave_cond.Case, wave);
        assert(any(rowW), 'Wave %s not found in waveInfo.CSV.', wave);

        T_wave    = wave_cond.T(rowW);     % [s]
        f_wave    = 1/T_wave;              % [Hz]
        f_wave_nd = f_wave * d_p / v_g;

        % Load trials
        filepath = [wave, partinfo];
        flist    = dir(fullfile(filepath, 'slip_vs_results_*.mat'));
        case_tot = numel(flist);

        if case_tot == 0
            warning('No slip_vs_results_*.mat in %s. Skip.', filepath);
            continue;
        end

        % Estimate v_s for tau_p from all trials
        v_s_all = NaN(1, case_tot);

        for j = 1:case_tot
            S0 = load(fullfile(filepath, flist(j).name), 'meta');
            v_s_all(j) = S0.meta.v_s;
        end

        v_s_case = mean(v_s_all, 'omitnan');

        Re_0   = v_s_case * d_p / nu;
        SchNau = 1 + 0.15 * Re_0^0.687;
        tau_p  = d_p^2 * (rho_p/rho_w + 0.5) / (18 * nu * SchNau);

        % --- Pass 1: load all trials, apply 5*tau_p exclusion, store ---
        wx_cell = {};
        wz_cell = {};
        dt_arr  = [];

        for j = 1:case_tot

            S = load(fullfile(filepath, flist(j).name));

            t  = S.time(:);
            wx = S.u_slip(:);
            wz = S.v_slip(:);

            m  = isfinite(t) & isfinite(wx) & isfinite(wz);

            t  = t(m);
            wx = wx(m);
            wz = wz(m);

            if numel(t) < minN
                continue;
            end

            % Exclude initial inertial transient.
            t0        = min(t);
            m_relaxed = (t - t0) >= 5 * tau_p;

            t  = t(m_relaxed);
            wx = wx(m_relaxed);
            wz = wz(m_relaxed);

            if numel(t) < minN
                continue;
            end

            dt_j = median(diff(t), 'omitnan');

            wx_cell{end+1} = wx;
            wz_cell{end+1} = wz;
            dt_arr(end+1)  = dt_j;
        end

        nTrial = numel(wx_cell);

        if nTrial == 0
            warning('No valid trials for %s%s.', wave, partinfo);
            continue;
        end

        % --- Match length: trim all trials to the shortest length ---
        % B&P: "we first matched the length of the slip velocity data time
        %       series for all data of a given particle ... by excluding
        %       the last few data points"
        Lcommon = min(cellfun(@numel, wx_cell));

        if Lcommon < minN
            continue;
        end

        dt = median(dt_arr);
        fs = 1 / dt;
        N  = Lcommon;

        % Native one-sided frequency grid
        f_full = (0:N-1)' * fs / N;
        kNyq   = floor(N/2) + 1;
        f      = f_full(1:kNyq);

        % --- Pass 2: per-trial FFT of zero-mean series ---
        Pux_all = zeros(nTrial, kNyq);
        Pvz_all = zeros(nTrial, kNyq);

        for j = 1:nTrial

            wx = wx_cell{j}(1:N);
            wz = wz_cell{j}(1:N);

            % Zero-mean fluctuation.
            wx_zm = wx - mean(wx, 'omitnan');
            wz_zm = wz - mean(wz, 'omitnan');

            % Direct FFT, one-sided power spectrum.
            % P(f) defined so that sum(P)*df ~ variance of the signal.
            X = fft(wx_zm);
            Z = fft(wz_zm);

            Sxx = abs(X).^2 / (fs * N);    % [m^2/s^2 / Hz]
            Szz = abs(Z).^2 / (fs * N);

            Sxx = Sxx(1:kNyq);
            Szz = Szz(1:kNyq);

            % One-sided correction.
            if kNyq >= 2
                interior = 2:(kNyq - mod(N+1, 2));   % skip Nyquist if N even
                Sxx(interior) = 2 * Sxx(interior);
                Szz(interior) = 2 * Szz(interior);
            end

            % Non-dimensionalize: P / (v_g * d_p)
            Pux_all(j,:) = Sxx(:)' / (v_g * d_p);
            Pvz_all(j,:) = Szz(:)' / (v_g * d_p);
        end

        % --- Ensemble mean + bootstrap CI ---
        Pux_mean = mean(Pux_all, 1, 'omitnan');
        Pvz_mean = mean(Pvz_all, 1, 'omitnan');

        if nTrial >= 2

            bootIdx = randi(nTrial, [nBoot, nTrial]);

            Pux_bs = zeros(nBoot, kNyq);
            Pvz_bs = zeros(nBoot, kNyq);

            for b = 1:nBoot
                Pux_bs(b,:) = mean(Pux_all(bootIdx(b,:), :), 1, 'omitnan');
                Pvz_bs(b,:) = mean(Pvz_all(bootIdx(b,:), :), 1, 'omitnan');
            end

            qlo = (1 - ciLevel) / 2;
            qhi = 1 - qlo;

            Pux_lo = quantile(Pux_bs, qlo, 1);
            Pux_hi = quantile(Pux_bs, qhi, 1);

            Pvz_lo = quantile(Pvz_bs, qlo, 1);
            Pvz_hi = quantile(Pvz_bs, qhi, 1);

        else

            Pux_lo = Pux_mean;
            Pux_hi = Pux_mean;

            Pvz_lo = Pvz_mean;
            Pvz_hi = Pvz_mean;

        end

        % Dimensionless frequency
        fhat = f * d_p / v_g;

        % patch는 x, y 좌표 크기가 정확히 같아야 하므로 row vector로 통일
        fhat_row = fhat(:)';

        % --- Plot: piecewise-linear over native frequency bins ---
        % B&P plot the spectra connecting native FFT bins directly,
        % which makes the coarse frequency resolution visible.
        plot(ax, fhat_row, Pux_mean, ls_u, ...
            'Color', c_wave, 'LineWidth', 1.5, ...
            'Marker', '.', 'MarkerSize', 8, ...
            'HandleVisibility', 'off');

        plot(ax, fhat_row, Pvz_mean, ls_v, ...
            'Color', c_wave, 'LineWidth', 1.5, ...
            'Marker', '.', 'MarkerSize', 8, ...
            'HandleVisibility', 'off');

        % CI bands
        idx_band = isfinite(Pux_lo) & isfinite(Pux_hi) & ...
                   isfinite(Pvz_lo) & isfinite(Pvz_hi) & ...
                   isfinite(fhat_row) & ...
                   (fhat_row >= 0) & (fhat_row <= 0.15);

        if any(idx_band)

            fb = fhat_row(idx_band);

            Pux_lo_b = Pux_lo(idx_band);
            Pux_hi_b = Pux_hi(idx_band);

            Pvz_lo_b = Pvz_lo(idx_band);
            Pvz_hi_b = Pvz_hi(idx_band);

            patch(ax, [fb, fliplr(fb)], ...
                      [Pux_lo_b, fliplr(Pux_hi_b)], ...
                      c_wave, ...
                      'FaceAlpha', 0.10, ...
                      'EdgeColor', 'none', ...
                      'HandleVisibility', 'off');

            patch(ax, [fb, fliplr(fb)], ...
                      [Pvz_lo_b, fliplr(Pvz_hi_b)], ...
                      c_wave, ...
                      'FaceAlpha', 0.10, ...
                      'EdgeColor', 'none', ...
                      'HandleVisibility', 'off');

        end

        % Dimensionless wave frequency
        xline(ax, f_wave_nd, ':', ...
            'Color', c_wave, ...
            'LineWidth', 1.5, ...
            'HandleVisibility', 'off');

    end
end

% --- Legend in slot 3 ---
ax_leg = ax_spec(3);
axes(ax_leg);
cla(ax_leg);
axis(ax_leg, 'off');
hold(ax_leg, 'on');

Hleg = gobjects(1, 2 * numel(waves));

for wI = 1:numel(waves)

    Hleg(2*wI-1) = plot(ax_leg, nan, nan, ls_u, ...
        'Color', c_wI{wI}, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('$(v''-u'')_x^*$ (%s)', waves{wI}));

    Hleg(2*wI) = plot(ax_leg, nan, nan, ls_v, ...
        'Color', c_wI{wI}, ...
        'LineWidth', 1.5, ...
        'Marker', '.', ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('$(v''-u'')_z^*$ (%s)', waves{wI}));

end

legend(ax_leg, Hleg, ...
    'Interpreter', 'latex', ...
    'FontSize', 13, ...
    'Location', 'northwest');

% --- Panel labels ---
labels = 'ab-cdefgh';

for sp = 1:9

    if sp == 3
        continue;
    end

    text(ax_spec(sp), -0.18, 1.02, sprintf('(%c)', labels(sp)), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'Clipping', 'off');

end

fig_spec.Position = [50 50 1200 1200];

% exportgraphics(fig_spec, ...
%     'E:\2.ParticlesInWave\Manuscript\figures\slip_fluctuation_spectra_BP_style.pdf', ...
%     'ContentType', 'vector', ...
%     'BackgroundColor', 'none');


%% ========================================================================
%  Local functions
% ========================================================================

function result = bootstrap_std(x, B, alpha)
% BOOTSTRAP_STD - Estimate standard deviation and confidence intervals via bootstrap
%
%   result = bootstrap_std(x, B, alpha)
%
% Inputs:
%   x     - data vector
%   B     - number of bootstrap resamples (e.g., 10000)
%   alpha - significance level (e.g., 0.05 -> 95% CI)
%
% Outputs:
%   result struct:
%       .sd_hat        - sample standard deviation of original data
%       .boot_mean     - mean of bootstrap estimates
%       .bias          - bias estimate
%       .boot_se       - bootstrap standard error
%       .ci_percentile - percentile CI
%       .ci_bca        - BCa CI
%       .bootStats     - bootstrap distribution of standard deviation

    if nargin < 2
        B = 10000;
    end
    if nargin < 3
        alpha = 0.05;
    end

    rng(42); % for reproducibility

    sd_hat   = std(x,0);               % sample standard deviation
    bootStats= bootstrp(B,@std,x);     % bootstrap distribution

    boot_mean = mean(bootStats);
    boot_se   = std(bootStats,0);
    bias_est  = boot_mean - sd_hat;

    ci_pct = prctile(bootStats,[100*alpha/2, 100*(1-alpha/2)]);
    ci_bca = bootci(B, {@std,x}, 'type','bca');

    result.sd_hat        = sd_hat;
    result.boot_mean     = boot_mean;
    result.bias          = bias_est;
    result.boot_se       = boot_se;
    result.ci_percentile = ci_pct;
    result.ci_bca        = ci_bca;
    result.bootStats     = bootStats;
end


function txt = ptvTip(eventObj)
% Datacursor text callback
% - Reads UserData.particle and UserData.k to identify the curve
% - From Position, returns v_slip/v_s, k*y_p, and y_p

pos = get(eventObj,'Position');   % [x, y] = [v_slip/v_s, k*y_p]
tgt = get(eventObj,'Target');
ud  = get(tgt,'UserData');

particle_num = NaN;
k = NaN;
if isstruct(ud)
    if isfield(ud,'particle'), particle_num = ud.particle; end
    if isfield(ud,'k'),        k           = ud.k;        end
end

if ~isnan(k)
    ky = pos(2);
    yp = ky / k;
else
    ky = pos(2);
    yp = NaN;
end

if isequal(particle_num,0)
    header = 'Ensemble Average';
elseif ~isnan(particle_num)
    header = sprintf('Particle #%d', particle_num);
else
    header = 'Unknown curve';
end

txt = {
    header
    sprintf('v_{slip}/v_s = %.6g', pos(1))
    sprintf('k y_p = %.6g', ky)
    sprintf('y_p = %.6g m', yp)
    };
end


function k = dispersion(h, T)

w = 2 * pi ./ T;  % angular frequency
g = 9.81;

% Ensure h and T are column/row compatible
[h_grid, T_grid] = ndgrid(h, T);
omega_grid = 2 * pi ./ T_grid;

k = zeros(size(h_grid));

% Solve dispersion for each (h, T) pair
for i = 1:numel(k)
    omega = omega_grid(i);
    h_val = h_grid(i);
    f = @(k) omega^2 - g * k * tanh(k * h_val);
    k(i) = abs(fzero(f, 0.01));
end

end

function [idx_keep, idx_trim, stats] = trimTrendSigmaFilter(x, t, opts)
% trimTrendSigmaFilter
% 1) Remove initial motion: discard everything before the first time i0 where
%    values "greater than" the lower bound (median - madMult*MAD) appear for
%    Nconsec consecutive samples. (After i0, keep samples even if they fall below the bound.)
% 2) Fit a linear trend (1st-order polyfit) on the segment after i0
% 3) Remove outliers using a 3*sigma (sigMult*sigma) threshold on residuals (2-pass refit)
%
% INPUT
%   x    : (Nx1) time series (e.g., v_slip or v_avg) [any unit]
%   t    : (Nx1) time vector (same length as x)
%   opts : struct with fields
%          .Nconsec  (default 5)
%          .madMult  (default 3)
%          .sigMult  (default 3)
%          .minN     (default 5)
%
% OUTPUT
%   idx_keep : final kept indices (w.r.t. original x,t)
%   idx_trim : indices after initial-motion trimming (= segment from i0 onward)
%   stats    : diagnostic info (threshold/fit/sigma, etc.)

    arguments
        x (:,1) double
        t (:,1) double
        opts.Nconsec (1,1) double = 5
        opts.madMult (1,1) double = 2
        opts.sigMult (1,1) double = 3
        opts.minN    (1,1) double = 5
    end

    if numel(x) ~= numel(t)
        error('x and t must have the same length.');
    end

    N = numel(x);
    finite_mask = isfinite(x) & isfinite(t);
    if nnz(finite_mask) < opts.minN
        idx_keep = [];
        idx_trim = [];
        stats = struct('reason','too_few_finite');
        return;
    end

    % Project so that the settling direction becomes positive (for detection logic)
    x0 = x; x0(~finite_mask) = NaN;
    sgn = sign(median(x0,'omitnan'));
    if sgn == 0, sgn = 1; end
    xproj = x0 * sgn;

    % --- (1) Initial-motion trim by lower bound: thr = median - madMult*MAD
    center_x = median(xproj,'omitnan');
    mad_x    = mad(xproj, 1);  % median(|x - median|)
    thr      = center_x - opts.madMult * mad_x;

    good = finite_mask & (xproj >= thr);

    d = diff([false; good(:); false]);
    rs = find(d==1);
    re = find(d==-1)-1;

    if isempty(rs)
        idx_trim = find(finite_mask);
        i0 = idx_trim(1);
    else
        runlen = re - rs + 1;
        okrun = find(runlen >= opts.Nconsec, 1, 'first');
        if isempty(okrun)
            [~, okrun] = max(runlen);
        end
        i0 = rs(okrun);
        idx_trim = find(finite_mask & ((1:N)' >= i0));
    end

    if numel(idx_trim) < opts.minN
        idx_keep = [];
        stats = struct('reason','too_few_after_trim','thr',thr,'center',center_x,'mad',mad_x,'i0',i0,'sgn',sgn);
        return;
    end

    % --- (2) Linear trend fit on trimmed segment (projected)
    tt = t(idx_trim);
    yy = xproj(idx_trim);

    p1 = polyfit(tt, yy, 1);
    yhat1 = polyval(p1, tt);
    res1 = yy - yhat1;

    % --- (3) 3-sigma residual filter (2-pass)
    sig1 = std(res1,'omitnan');
    inl1 = abs(res1) <= opts.sigMult * sig1;

    tt2 = tt(inl1);
    yy2 = yy(inl1);

    if numel(tt2) < opts.minN
        idx_keep = idx_trim;
        stats = struct('thr',thr,'center',center_x,'mad',mad_x,'i0',i0,'sgn',sgn, ...
                       'p_final',p1,'sigma_final',sig1,'pass','onepass_fallback');
        return;
    end

    p2 = polyfit(tt2, yy2, 1);
    yhat2 = polyval(p2, tt);
    res2 = yy - yhat2;

    sig2 = std(res2(inl1),'omitnan');
    inl2 = abs(res2) <= opts.sigMult * sig2;

    idx_keep = idx_trim(inl2);

    stats = struct('thr',thr,'center',center_x,'mad',mad_x,'i0',i0,'sgn',sgn, ...
                   'p1',p1,'sigma1',sig1,'p_final',p2,'sigma_final',sig2, ...
                   'n_finite',nnz(finite_mask),'n_trim',numel(idx_trim),'n_keep',numel(idx_keep));
end

function colj = trialColor(baseCol, j, n)
    if all(baseCol == 0)
        % grayscale gradient for W1 (base = black)
        gval = 0.15 + 0.5*(j-1)/max(1,n-1);
        colj = [gval gval gval];
    else
        % HSV-based gradient: vary value, keep hue, raise saturation floor
        hsv_base = rgb2hsv(baseCol);
        hsv_j = hsv_base;
        hsv_j(3) = 0.25 + 0.65*(j-1)/max(1,n-1);
        hsv_j(2) = max(0.35, hsv_base(2));
        colj = hsv2rgb(hsv_j);
    end
end

function [slip_ratio_d22, vz_ratio_d22] = profileD22SlipVelocity(z, H, T, h, d_p, rho_p, rho_w, nu, vs_dim)
    % z              : dimensional vertical coordinate [m], upward-positive, z=0 at SWL
    % vs_dim         : terminal settling speed magnitude [m/s], positive scalar
    %
    % vz_ratio_d22   = <v_z>_D22 / v_s, upward-positive convention
    % slip_ratio_d22 = -<w_z>_D22 / v_s
    %
    % We interpret the D22 vertical particle-velocity expression as
    % <v_z> = <u_z> + <w_z>. The sampled carrier-flow contribution is
    % associated with the vertical wave-velocity correction term
    % proportional to 0.5*tanh(kh)*d(u_SD)/d(kz). Therefore,
    %
    %   -<w_z>/v_s = -<v_z>/v_s + <u_z>/v_s.
    %
    % With this decomposition, the D22 slip correction retains the
    % inertia-dependent settling modification term.

    Cm = 0.5;

    omega = 2*pi/T;
    k = dispersion(h, T);
    c = omega / k;
    a = H/2;

    gamma = rho_p / rho_w;
    beta  = (1 + Cm) / (gamma + Cm);

    Re0 = vs_dim * d_p / nu;
    D0  = 1 + 0.15 * Re0^0.687;

    St_lin = omega * d_p^2 * (gamma + Cm) / (18*nu);
    St0    = St_lin / D0;

    vs_nd = vs_dim / c;

    eps = k*a / tanh(k*h);

    z_nd = k*z;
    kh   = k*h;

    A_d22 = sqrt( ...
        (1 - St0^2*(1-beta)).^2 + ...
        (St0*(1-beta)).^2 );

    uSD = eps^2 * cosh(2*(z_nd + kh)) ./ (2*cosh(kh)^2);
    duSD_dznd = eps^2 * sinh(2*(z_nd + kh)) ./ (cosh(kh)^2);

    % D22 particle vertical velocity, upward-positive:
    % <v_z>/v_s
    vz_ratio_d22 = -( ...
        1 ...
        + (A_d22.^2 ./ (1 + vs_nd^2)) .* uSD ...
        + 0.5 * tanh(kh) .* duSD_dznd );

    % Corresponding sampled carrier-flow vertical velocity contribution:
    % <u_z>/v_s
    uz_ratio_d22 = -0.5 * tanh(kh) .* duSD_dznd;

    % D22 relative vertical slip:
    % -<w_z>/v_s = -(<v_z>-<u_z>)/v_s
    slip_ratio_d22 = -vz_ratio_d22 + uz_ratio_d22;
end