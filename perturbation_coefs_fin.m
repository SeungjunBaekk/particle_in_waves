clear; clc; close all;

%-------------------------------------------------------
% 1. Wave parameters → dispersion → kh
%-------------------------------------------------------
h = [0.55 0.55 0.80];        % W1, W2, W3
T = [1.51 1.252 1.01];       % periods
for i = 1:3
k(i) = dispersion(h(i),T(i));         % user-supplied function
end
kh = k .* h;                 % dimensionless relative depth

nu = waterNu(24.5);
dp = 780*1e-6;

%-------------------------------------------------------
% 2. Parameter ranges
%-------------------------------------------------------
St_vec    = logspace(-3,0,800);   % Stokes number 10^-3 ~ 10^0
nSt       = numel(St_vec);
beta_list = [0.2 0.5 0.8];        % three density-ratio parameters
nBeta     = numel(beta_list);
nZ        = 5;                   % number of depth levels per case

% line styles for beta
lineStyles = {'-','--','-.'};     % β = 0.2, 0.5, 0.8

% grayscale colors for the D22 (small-St, no wave-shear) reference curves,
% one shade per beta (black -> light gray)
d22Colors = {[0 0 0],[0 0 0],[0 0 0]};   % β = 0.2, 0.5, 0.8
d22LineWidth = 1.5;

% for phase plots
phase_tick  = 0:0.5*pi:2*pi;
phase_label = {'0','0.5\pi','\pi','1.5\pi','2\pi'};
font_sz = 20;

%-------------------------------------------------------
% 3. Loop over 3 (h,T) cases
%-------------------------------------------------------
for icase = 1:3

    khc   = kh(icase);
    St    = St_vec;
    St2   = St.^2;

    % depth in nondimensional vertical coordinate z (kz):
    %   bottom: z = -khc, surface: z = 0
    z_vec   = linspace(-khc,0,nZ+2);       % kz
    z_vec   = z_vec(2:end-1);
    zh_vec  = z_vec / khc;               % z/h in [-1,0] for labeling
    cmap    = turbo(nZ);                 % color by z

    % preallocate: (beta, z, St)
    Ax1   = zeros(nBeta,nZ,nSt);
    Az1   = zeros(nBeta,nZ,nSt);
    Phix1 = zeros(nBeta,nZ,nSt);
    Phiz1 = zeros(nBeta,nZ,nSt);

    Ax2   = zeros(nBeta,nZ,nSt);
    Az2   = zeros(nBeta,nZ,nSt);
    Phix2 = zeros(nBeta,nZ,nSt);
    Phiz2 = zeros(nBeta,nZ,nSt);

    % --- D22 (small-St, wave-shear dropped) first-order slip coefficients ---
    % Obtained from the present slip solution by (i) removing the v_s-dependent
    % wave-shear terms and (ii) taking the small-St limit (denominator 1+St^2 -> 1),
    % consistent with Appendix A.3 (eqs. A32-A34) of the manuscript.
    Ax1_D22   = zeros(nBeta,nZ,nSt);
    Az1_D22   = zeros(nBeta,nZ,nSt);
    Phix1_D22 = zeros(nBeta,nZ,nSt);
    Phiz1_D22 = zeros(nBeta,nZ,nSt);

    % helper to pull 1D slice along St
    getSlice = @(M,ib,iz) squeeze(M(ib,iz,:)).';

    %---------------------------------------------------
    % 3-1. Coefficient calculation
    %---------------------------------------------------
    for ibeta = 1:nBeta
        beta = beta_list(ibeta);
        gamma= 0.5*(3/beta-1);
        g = 9.81;
        m = 0.687; alpha = 0.15;

        tau_p = fzero(@(x) 18*nu*(1+alpha*(dp*x*(1-beta)*g/nu)^m)*x==dp^2*(gamma + 0.5),0.1 );
        % theoretical dimensionless settling slip
        vs  = St*(1 - beta)./tanh(khc) ;    % v_s = St(1-β)/tanh(kh)
        Rep = vs*dp/nu;
        dragD0 = 1 + alpha .* Rep.^m;
        
        K_ = m.*(1-dragD0)./vs./St./dragD0;

        for iz = 1:nZ
            z = z_vec(iz);

            % linear-wave vertical structure
            A = cosh(z + khc) / cosh(khc);
            B = sinh(z + khc) / cosh(khc);

            % ----- first-order coefficients px,qx,pz,qz -----
            px = ( St.*B.*vs - St2.*A.*(1-beta) ) ./ (St2 + 1);
            qx = ( -St.*A.*(1-beta) - St2.*B.*vs ) ./ (St2 + 1);

            pz = ( St2.*A.*vs + St.*(1-beta).*B ) ./ (St2 + 1);
            qz = ( St.*A.*vs - St2.*B.*(1-beta) ) ./ (St2 + 1);

            Ax1(ibeta,iz,:)   = hypot(px,qx);
            Az1(ibeta,iz,:)   = hypot(pz,qz);
            Phix1(ibeta,iz,:) = atan2(px,qx);
            Phiz1(ibeta,iz,:) = atan2(pz,qz);

            % ----- D22 small-St slip coefficients (wave-shear removed) -----
            % w_x1 = -A(1-beta)St^2 cos - A(1-beta)St sin
            % w_z1 =  B(1-beta)St   cos - B(1-beta)St^2 sin
            pxD = -A.*(1-beta).*St2;
            qxD = -A.*(1-beta).*St;
            pzD =  B.*(1-beta).*St;
            qzD = -B.*(1-beta).*St2;

            Ax1_D22(ibeta,iz,:)   = hypot(pxD,qxD);
            Az1_D22(ibeta,iz,:)   = hypot(pzD,qzD);
            Phix1_D22(ibeta,iz,:) = atan2(pxD,qxD);
            Phiz1_D22(ibeta,iz,:) = atan2(pzD,qzD);

            % ----- second-order coefficients C,D → a1,a2,b1,b2 -----
            C0 = 0.5 * ( ...
                   A.*qx + B.*pz ...
                 + (A.*px + px.^2 + qx.^2)./St ...
                 + (1-beta).*A.*qx ...
                 - vs.*B.*(A + px) );

            C1 = 0.5 * ( ...
                  -A.*qx + B.*pz ...
                 + (A.*px + px.^2 - qx.^2)./St ...
                 - (1-beta).*A.*qx ...
                 - vs.*B.*(A + px) );

            C2 = 0.5 * ( ...
                   A.*px + B.*qz ...
                 + (A.*qx + 2*px.*qx)./St ...
                 + (1-beta).*(A.*px + B.^2) ...
                 - vs.*B.*qx );

            D0 = 0.5 * ( ...
                   B.*px + A.*qz ...
                 + (A.*pz + px.*pz + qx.*qz)./St ...
                 - (1-beta).*B.*(A + px) ...
                 - vs.*A.*qx );

            D1 = 0.5 * ( ...
                   B.*px - A.*qz ...
                 + (A.*pz + px.*pz - qx.*qz)./St ...
                 + (1-beta).*B.*(A - px) ...
                 + vs.*A.*qx );

            D2 = 0.5 * ( ...
                   B.*qx + A.*pz ...
                 + (A.*qz + px.*qz + qx.*pz)./St ...
                 + (1-beta).*B.*(2*A - qx) ...
                 - vs.*A.*(A + px) );
            
            Gamma0 = C0 + 0.5*K_.*(px.*pz + qx.*qx);
            Gamma1 = C1 + 0.5*K_.*(px.*pz - qx.*qx);
            Gamma2 = C2 + 0.5*K_.*(px.*qz + qx.*px);

            Delta0 = D0 + 0.5*K_.*(pz.^2 + qz.^2);
            Delta1 = D1 + 0.5*K_.*(pz.^2 - qz.^2);
            Delta2 = D2 +     K_.*(pz .* qz);
    
            a0 = -St .* Gamma0;
            a1 = (-Gamma1.*St - 2*Gamma2.*St2) / (1 + 4*St2);
            a2 = (-Gamma2.*St + 2*Gamma1.*St2) / (1 + 4*St2);

            b0 = -St .* Delta0;
            b1 = (-Delta1.*St - 2*Delta2.*St2) / (1 + 4*St2);
            b2 = (-Delta2.*St + 2*Delta1.*St2) / (1 + 4*St2);

            Ax2(ibeta,iz,:)   = hypot(a1,a2);
            Az2(ibeta,iz,:)   = hypot(b1,b2);
            Phix2(ibeta,iz,:) = atan2(a1,a2);
            Phiz2(ibeta,iz,:) = atan2(b1,b2);
        end
    end

    % wrap phases into [0,2π]
    Phix1 = mod(Phix1,2*pi);
    Phiz1 = mod(Phiz1,2*pi);
    Phix2 = mod(Phix2,2*pi);
    Phiz2 = mod(Phiz2,2*pi);
    Phix1_D22 = mod(Phix1_D22,2*pi);
    Phiz1_D22 = mod(Phiz1_D22,2*pi);

    %---------------------------------------------------
    % 3-2. Plot: 2 x 2 panel for this kh
    %      Panels (a,b): fixed mid-depth z'/h = -0.5, beta color-coded
    %      Present solution: solid line; D22 shear-free slip limit: dotted line
    %---------------------------------------------------
    f  = figure('Name',sprintf('Case %d: kh = %.2f',icase,khc), 'Color','w');
    f.Position = [100 100 1200 800];

    ax = gobjects(1,4);

    % Use the mid-depth curve for amplitude panels.
    % The manuscript coordinate uses z'=0 at the still-water level and z'=-h at the bed.
    % Thus the user's z/h = 0.5 corresponds to z'/h = -0.5 here.
    target_zh_amp = -0.5;
    [~, iz_amp] = min(abs(zh_vec - target_zh_amp));
    zh_amp = zh_vec(iz_amp);

    % Beta color scale for panels (a,b)
    betaMap = parula(256);
    betaMin = min(beta_list);
    betaMax = max(beta_list);
    betaColor = @(b) betaMap(max(1,min(256,round(1 + (b-betaMin)/(betaMax-betaMin)*255))),:);

    % Depth color scale retained for panels (c,d)
    depthMap = turbo(256);
    zh_norm = (zh_vec + 1);                % z'/h: -1 -> 0, 0 -> 1
    idxZ    = round( zh_norm*(size(depthMap,1)-1) ) + 1;
    idxZ    = max(min(idxZ,size(depthMap,1)),1);

    getSlice = @(M,ib,iz) squeeze(M(ib,iz,:)).';

    % ---------------- Panel (a): R_{x1}, z'/h = -0.5 ----------------
    ax(1) = subplot(2,2,1); hold(ax(1),'on');
    for ib = 1:nBeta
        thisCol = betaColor(beta_list(ib));

        loglog(ax(1), St_vec, getSlice(Ax1,ib,iz_amp), ...
            'Color',thisCol, ...
            'LineStyle','-', ...
            'LineWidth',1.2);

        loglog(ax(1), St_vec, getSlice(Ax1_D22,ib,iz_amp), ...
            'Color',thisCol, ...
            'LineStyle',':', ...
            'LineWidth',1.5);
    end
    grid(ax(1),'on'); xlim(ax(1),[0.01 1]), ylim(ax(1),[0 1]);
    xlabel(ax(1),'$St_0$','Interpreter','latex');
    ylabel(ax(1),'$R_{x1}$','Interpreter','latex');
    set(ax(1),'XScale','log','FontSize',font_sz);
    colormap(ax(1), betaMap);
    caxis(ax(1), [betaMin betaMax]);

    % ---------------- Panel (b): R_{z1}, z'/h = -0.5 ----------------
    ax(2) = subplot(2,2,2); hold(ax(2),'on');
    for ib = 1:nBeta
        thisCol = betaColor(beta_list(ib));

        loglog(ax(2), St_vec, getSlice(Az1,ib,iz_amp), ...
            'Color',thisCol, ...
            'LineStyle','-', ...
            'LineWidth',1.2);

        loglog(ax(2), St_vec, getSlice(Az1_D22,ib,iz_amp), ...
            'Color',thisCol, ...
            'LineStyle',':', ...
            'LineWidth',1.5);
    end
    grid(ax(2),'on'); xlim(ax(2),[0.01 1]), ylim(ax(2),[0 1]);
    xlabel(ax(2),'$St_0$','Interpreter','latex');
    ylabel(ax(2),'$R_{z1}$','Interpreter','latex');
    set(ax(2),'XScale','log','FontSize',font_sz);
    colormap(ax(2), betaMap);
    caxis(ax(2), [betaMin betaMax]);

    % Beta colorbar for panels (a,b)
    pos1 = get(ax(1),'Position');
    pos2 = get(ax(2),'Position');
    cb_beta = colorbar(ax(2));
    cb_beta.Label.String = '';
    cb_beta.Ticks = beta_list;
    cb_beta.TickLabels = arrayfun(@(x) sprintf('%.1f',x), beta_list, 'UniformOutput',false);
    cb_beta.FontSize = 16;
    set(ax(1),'Position',pos1);
    set(ax(2),'Position',pos2);
    cbpos = get(cb_beta,'Position');
    cbpos(1) = pos2(1) + pos2(3) + 0.02;
    cbpos(2) = pos2(2);
    cbpos(3) = 0.02;
    cbpos(4) = pos2(4);
    set(cb_beta,'Position',cbpos);
    annotation(f,'textbox', [cbpos(1)-0.01, min(cbpos(2)+cbpos(4)+0.005,0.96), cbpos(3)+0.02, 0.03], ...
        'String','$\beta$', 'Interpreter','latex', ...
        'LineStyle','none', 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',16);

    % Legend for line type in panel (b)
    h_present = plot(ax(2),NaN,NaN,'k-','Color',[0.55 0.55 0.55],'LineWidth',1.5,'DisplayName','This Study');
    h_d22     = plot(ax(2),NaN,NaN,'k:','Color',[0.55 0.55 0.55],'LineWidth',2.0,'DisplayName','DiBenedetto et al., 2022');
    legend(ax(2),[h_present h_d22], 'Location','northwest', ...
        'Interpreter','latex','FontSize',14);

    % ---------------- Panel (c): theta_{x1}; keep depth dependence --------
    ax(3) = subplot(2,2,3); hold(ax(3),'on');
    for ib = 1:nBeta
        for iz = nZ:-1:1
            thisCol = depthMap(idxZ(iz),:);
            semilogx(ax(3), St_vec, getSlice(Phix1,ib,iz), ...
                'Color',thisCol, ...
                'LineStyle',lineStyles{ib}, ...
                'LineWidth',1.2);
        end
    end
    semilogx(ax(3), St_vec, getSlice(Phix1_D22,1,1), ...
        'Color',[0 0 0],'LineStyle',':','LineWidth',d22LineWidth);
    grid(ax(3),'on'); xlim(ax(3),[0.01 1]), ylim(ax(3),[0 2*pi]);
    xlabel(ax(3),'$St_0$','Interpreter','latex');
    ylabel(ax(3),'$\theta_{x1}$','Interpreter','latex');
    set(ax(3),'XScale','log','YTick',phase_tick, ...
        'YTickLabel',phase_label,'FontSize',font_sz);

    % ---------------- Panel (d): theta_{z1}; keep depth dependence --------
    ax(4) = subplot(2,2,4); hold(ax(4),'on');
    for ib = 1:nBeta
        for iz = nZ:-1:1
            thisCol = depthMap(idxZ(iz),:);
            semilogx(ax(4), St_vec, getSlice(Phiz1,ib,iz), ...
                'Color',thisCol, ...
                'LineStyle',lineStyles{ib}, ...
                'LineWidth',1.2);
        end
    end
    semilogx(ax(4), St_vec, getSlice(Phiz1_D22,1,1), ...
        'Color',[0 0 0],'LineStyle',':','LineWidth',d22LineWidth);
    grid(ax(4),'on'); xlim(ax(4),[0.01 1]), ylim(ax(4),[0 2*pi]);
    xlabel(ax(4),'$St_0$','Interpreter','latex');
    ylabel(ax(4),'$\theta_{z1}$','Interpreter','latex');
    set(ax(4),'XScale','log','YTick',phase_tick, ...
        'YTickLabel',phase_label,'FontSize',font_sz);

    % Depth colorbar for panels (c,d)
    pos3 = get(ax(3),'Position');
    pos4 = get(ax(4),'Position');
    colormap(ax(4), depthMap);
    caxis(ax(4),[-1 0]);
    cb_depth = colorbar(ax(4));
    cb_depth.Label.String = '';
    cb_depth.Ticks = [-1, zh_vec, 0];
    cb_depth.TickLabels = arrayfun(@(x) sprintf('%.2f',x), [-1, zh_vec, 0], 'UniformOutput',false);
    cb_depth.FontSize = 16;
    set(ax(3),'Position',pos3);
    set(ax(4),'Position',pos4);
    cbpos2 = get(cb_depth,'Position');
    cbpos2(1) = pos4(1) + pos4(3) + 0.02;
    cbpos2(2) = pos4(2);
    cbpos2(3) = 0.02;
    cbpos2(4) = pos4(4);
    set(cb_depth,'Position',cbpos2);
    annotation(f,'textbox', [cbpos2(1)-0.015, min(cbpos2(2)+cbpos2(4)+0.005,0.96), cbpos2(3)+0.03, 0.03], ...
        'String','$z''/h$', 'Interpreter','latex', ...
        'LineStyle','none', 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize',16);

    % No beta legend is added to panel (c); beta is already introduced in panels (a,b).

    % Present/D22 legend for phase panels
    hold(ax(4),'on');
    Lpresent = plot(ax(4),NaN,NaN,'Color',[0.55 0.55 0.55], ...
        'LineStyle','-','LineWidth',1.4,'DisplayName','This study');
    LD22 = plot(ax(4),NaN,NaN,'Color',[0 0 0], ...
        'LineStyle',':','LineWidth',2.0,'DisplayName','DiBenedetto et al., 2022');
    legend(ax(4),[Lpresent LD22],'Location','northwest', ...
        'Interpreter','latex','FontSize',13);
    text(ax(1),0.75,0.92,'$z^\prime/h=-0.5$', ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize',14, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top');
    
    text(ax(2),0.75,0.92,'$z^\prime/h=-0.5$', ...
        'Units','normalized', ...
        'Interpreter','latex', ...
        'FontSize',14, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top');

    % Panel labels (a)--(d)
    labels = 'abcd';
    for kk = 1:4
        text(ax(kk),-0.18,1.02,sprintf('(%c)',labels(kk)), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',18,'Clipping','off');
    end

    exportgraphics(f,sprintf('pert_coefs_W%d.pdf',icase));

end

% exportgraphics(f,'pert_coefs_W3.pdf')
% exportgraphics(2,'pert_coefs_W2.pdf')
% exportgraphics(3,'pert_coefs_W3.pdf')