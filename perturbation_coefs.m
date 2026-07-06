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
d22LineWidth = 1.6;

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
    % 3-2. Plot: 4×2 panel for this kh
    %---------------------------------------------------
    f  = figure('Name',sprintf('Case %d: kh = %.2f',icase,khc));
    
    % 연속 colormap (예: 256 단계)
    cmapFull = turbo(256);
    cmapD22  = gray(256);
    colormap(f,cmapFull);
    
    ax = gobjects(1,8);
    
    % z/h ∈ [-1,0] 를 [0,1] 로 매핑해서 colormap index 계산
    %   zh_vec: 길이 nZ 의 z/h 값
    zh_norm = (zh_vec + 1);                % -1→0, 0→1
    idxZ    = round( zh_norm*(size(cmapFull,1)-1) ) + 1;
    idxZ    = max(min(idxZ,size(cmapFull,1)),1);  % safety
    
    getSlice = @(M,ib,iz) squeeze(M(ib,iz,:)).';
    
    % (a) R_{x1}
    ax(1) = subplot(2,2,1); hold on;
    for ib = 1:nBeta
        for iz = 1:nZ
            thisCol = cmapFull(idxZ(iz),:);
            loglog(St_vec, getSlice(Ax1,ib,iz), ...
                'Color',thisCol, ...
                'LineStyle',lineStyles{ib}, ...
                'LineWidth',1.0);
        end
    end
    % D22 small-St reference (grayscale, one per beta; z=surface for visibility)
    for ib = 1:nBeta
        for iz = 1:nZ
            thisCol = cmapD22(idxZ(iz),:);
    
            loglog(St_vec, getSlice(Ax1_D22,ib,iz), ...
                'Color', thisCol, ...
                'LineStyle', ':', ...
                'LineWidth', d22LineWidth);
        end
    end
    grid on; ylim([0 1]);
    xlabel('$St$','Interpreter','latex');
    ylabel('$R_{x1}$','Interpreter','latex');
    set(ax(1),'XScale','log','FontSize',font_sz);
    
    % (b) R_{z1}
    ax(2) = subplot(2,2,2); hold on;
    for ib = 1:nBeta
        for iz = 1:nZ
            thisCol = cmapFull(idxZ(iz),:);
            loglog(St_vec, getSlice(Az1,ib,iz), ...
                'Color',thisCol, ...
                'LineStyle',lineStyles{ib}, ...
                'LineWidth',1.0);
        end
    end
    % D22 small-St reference (grayscale, one per beta; z=surface for visibility)
    for ib = 1:nBeta
        for iz = 1:nZ
            thisCol = cmapD22(idxZ(iz),:);
    
            loglog(St_vec, getSlice(Az1_D22,ib,iz), ...
                'Color', thisCol, ...
                'LineStyle', ':', ...
                'LineWidth', d22LineWidth);
        end
    end
    grid on; ylim([0 1]);
    xlabel('$St$','Interpreter','latex');
    ylabel('$R_{z1}$','Interpreter','latex');
    set(ax(2),'XScale','log','FontSize',font_sz);
    
    % (c) R_{x2}
    % ax(3) = subplot(4,2,3); hold on;
    % for ib = 1:nBeta
    %     for iz = 1:nZ
    %         thisCol = cmapFull(idxZ(iz),:);
    %         loglog(St_vec, getSlice(Ax2,ib,iz), ...
    %             'Color',thisCol, ...
    %             'LineStyle',lineStyles{ib}, ...
    %             'LineWidth',1.0);
    %     end
    % end
    % grid on; ylim([0 0.3]);
    % xlabel('$St$','Interpreter','latex');
    % ylabel('$R_{x2}$','Interpreter','latex');
    % set(ax(3),'XScale','log','FontSize',font_sz);
    % 
    % % (d) R_{z2}
    % ax(4) = subplot(4,2,4); hold on;
    % for ib = 1:nBeta
    %     for iz = 1:nZ
    %         thisCol = cmapFull(idxZ(iz),:);
    %         loglog(St_vec, getSlice(Az2,ib,iz), ...
    %             'Color',thisCol, ...
    %             'LineStyle',lineStyles{ib}, ...
    %             'LineWidth',1.0);
    %     end
    % end
    % grid on; ylim([0 0.3]);
    % xlabel('$St$','Interpreter','latex');
    % ylabel('$R_{z2}$','Interpreter','latex');
    % set(ax(4),'XScale','log','FontSize',font_sz);
    
    % (e) theta_{x1}
    ax(3) = subplot(2,2,3); hold on;
    for ib = 1:nBeta
        for iz = nZ:-1:1
            thisCol = cmapFull(idxZ(iz),:);
            semilogx(St_vec, getSlice(Phix1,ib,iz), ...
                'Color',thisCol, ...
                'LineStyle',lineStyles{ib}, ...
                'LineWidth',1.0);
        end
    end
    % D22 small-St reference: theta collapses (independent of beta,z) -> single curve
    semilogx(St_vec, getSlice(Phix1_D22,1,1), ...
        'Color',[0 0 0],'LineStyle',':','LineWidth',d22LineWidth);
    grid on; ylim([0 2*pi]);
    xlabel('$St$','Interpreter','latex');
    ylabel('$\theta_{x1}$','Interpreter','latex');
    set(ax(3),'XScale','log','YTick',phase_tick,...
        'YTickLabel',phase_label,'FontSize',font_sz);
    
    % (f) theta_{z1}
    ax(4) = subplot(2,2,4); hold on;
    for ib = 1:nBeta
        for iz = nZ:-1:1
            thisCol = cmapFull(idxZ(iz),:);
            semilogx(St_vec, getSlice(Phiz1,ib,iz), ...
                'Color',thisCol, ...
                'LineStyle',lineStyles{ib}, ...
                'LineWidth',1.0);
        end
    end
    % D22 small-St reference: theta collapses (independent of beta,z) -> single curve
    hD22 = semilogx(St_vec, getSlice(Phiz1_D22,1,1), ...
        'Color',[0 0 0],'LineStyle',':','LineWidth',d22LineWidth);
    grid on; ylim([0 2*pi]);
    xlabel('$St$','Interpreter','latex');
    ylabel('$\theta_{z1}$','Interpreter','latex');
    set(ax(4),'XScale','log','YTick',phase_tick,...
        'YTickLabel',phase_label,'FontSize',font_sz);
    
    % % (g) theta_{x2}
    % ax(7) = subplot(4,2,7); hold on;
    % for ib = 1:nBeta
    %     for iz = nZ:-1:1
    %         thisCol = cmapFull(idxZ(iz),:);
    %         semilogx(St_vec, getSlice(Phix2,ib,iz), ...
    %             'Color',thisCol, ...
    %             'LineStyle',lineStyles{ib}, ...
    %             'LineWidth',1.0);
    %     end
    % end
    % grid on; ylim([0 2*pi]);
    % xlabel('$St$','Interpreter','latex');
    % ylabel('$\theta_{x2}$','Interpreter','latex');
    % set(ax(7),'XScale','log','YTick',phase_tick,...
    %     'YTickLabel',phase_label,'FontSize',font_sz);
    % 
    % % (h) theta_{z2}
    % ax(8) = subplot(4,2,8); hold on;
    % for ib = 1:nBeta
    %     for iz = nZ:-1:1
    %         thisCol = cmapFull(idxZ(iz),:);
    %         semilogx(St_vec, getSlice(Phiz2,ib,iz), ...
    %             'Color',thisCol, ...
    %             'LineStyle',lineStyles{ib}, ...
    %             'LineWidth',1.0);
    %     end
    % end
    % grid on; ylim([0 2*pi]);
    % xlabel('$St$','Interpreter','latex');
    % ylabel('$\theta_{z2}$','Interpreter','latex');
    % set(ax(8),'XScale','log','YTick',phase_tick,...
    %     'YTickLabel',phase_label,'FontSize',font_sz);
    
    %----- colorbar: depth z/h (연속, tick 은 사용한 z/h 에 맞춤) -----
    origPos = cell(1,8);
    for k = 1:4
        origPos{k} = get(ax(k),'Position');
    end
    
    cb = colorbar(ax(2));    % attach to panel (b)
    cb.Label.String      = '$z/h$';
    cb.Label.Interpreter = 'latex';
    cb.Label.Rotation    = 0;
    cb.Label.HorizontalAlignment = 'center';
    cb.Label.VerticalAlignment   = 'bottom';
    cb.Label.Position = [1, 1.05, 0];

    
    % colorbar 의 축은 [0,1] 이므로, 우리가 쓰는 z/h 들의 위치는 zh_norm
    cb.Ticks      = [0, zh_norm, 1];                   % 정확히 사용한 z/h 위치
    cb.TickLabels = arrayfun(@(x) sprintf('%.2f',x), ...% label = 실제 z/h
                             [-1, zh_vec, 0], 'uni',0);
    
    % 축 위치 복원
    for k = 1:4
        set(ax(k),'Position',origPos{k});
    end
    
    % colorbar 오른쪽으로 이동
    posB  = origPos{2};
    cbPos = get(cb,'Position');
    gap   = 0.02;
    cbPos(1) = posB(1) + posB(3) + gap;
    cbPos(2) = posB(2);
    cbPos(3) = 0.02;
    cbPos(4) = posB(4);
    set(cb,'Position',cbPos,'FontSize',16);
    
    %----- legend for beta (dummy black lines) -----
%----- (b) legend: beta line styles only -----
    hold(ax(2),'on');
    Lbeta = gobjects(1,nBeta);
    for ib = 1:nBeta
        Lbeta(ib) = plot(ax(2),NaN,NaN, ...
            'Color',[.7 .7 .7], ...
            'LineStyle',lineStyles{ib}, ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('\\beta=%.1f',beta_list(ib)));
    end
    legend(ax(2),Lbeta,'Location','northwest', ...
        'Orientation','horizontal','Interpreter','tex','FontSize',16);

    %----- (d) legend: D22 dotted entry only -----
    hold(ax(4),'on');
    Lb26 = plot(ax(4),NaN,NaN, ...
        'Color',[0.7 0.7 0.7],'LineStyle','-','LineWidth',d22LineWidth, ...
        'DisplayName','This Study');
    Ld22 = plot(ax(4),NaN,NaN, ...
        'Color',[0 0 0],'LineStyle',':','LineWidth',d22LineWidth, ...
        'DisplayName','DiBenedetto et al. (2022)');
    legend(ax(4),[Lb26 Ld22],'Location','northeast', ...
        'Interpreter','latex','FontSize',16);

    %----- (d) colorbar: depth z/h (same style as panel b) -----
    posD_before = get(ax(4),'Position');     % preserve panel size
    cb2 = colorbar(ax(4));
    cb2.Label.String      = '$z/h$';
    cb2.Label.Interpreter = 'latex';
    cb2.Label.Rotation    = 0;
    cb2.Label.HorizontalAlignment = 'center';
    cb2.Label.VerticalAlignment   = 'bottom';
    cb2.Label.Position = [1, 1.05, 0];
    cb2.Ticks      = [0, zh_norm, 1];
    cb2.TickLabels = arrayfun(@(x) sprintf('%.2f',x), ...
                              [-1, zh_vec, 0], 'uni',0);
    set(ax(4),'Position',posD_before);       % restore panel after colorbar attach

    % place colorbar to the right of panel (d), matching panel-(b) gap/width
    cbPosD = get(cb2,'Position');
    gapD   = 0.02;
    cbPosD(1) = posD_before(1) + posD_before(3) + gapD;
    cbPosD(2) = posD_before(2);
    cbPosD(3) = 0.02;
    cbPosD(4) = posD_before(4);
    set(cb2,'Position',cbPosD,'FontSize',16);

    % panel labels (a)~(h)
    labels = 'abcdefgh';
    for k = 1:4
        text(ax(k),-0.18,1.02,sprintf('(%c)',labels(k)), ...
            'Units','normalized','HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','FontSize',18,'Clipping','off');
    end
    f.Position = [100 100 1200 800];
    exportgraphics(f,sprintf('pert_coefs_W%d.pdf',icase))
end

% exportgraphics(f,'pert_coefs_W3.pdf')
% exportgraphics(2,'pert_coefs_W2.pdf')
% exportgraphics(3,'pert_coefs_W3.pdf')