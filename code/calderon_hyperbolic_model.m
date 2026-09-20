%% CALDERON_HYPERBOLIC_MODEL
% Numerical finite-section experiment for the multiplicity-normalized
% rational disk kernel
%
%   (Kf)(z) = (1/pi) int_D f(w)/|1-z*conj(w)|^2 dA(w).
%
% In angular mode m, after t=r^2, the radial block is
%
%   (K_m g)(t) = int_0^1 (st)^(|m|/2)/(1-st) g(s) ds.
%
% The exact continuous spectral multiplier is
%
%   lambda(mu) = pi*sech(pi*mu),   sigma(mu)=sqrt(lambda(mu)).
%
% The finite hyperbolic truncation x in (0,xMax), r=tanh(x/2), makes
% each block compact and produces discretization-dependent wave packets.
% The script:
%   * computes leading eigenpairs in each angular block;
%   * globally orders the distinct (|m|, radial index) modes;
%   * plots generalized singular functions in the Euclidean disk;
%   * compares the eigenvalue-inferred mu with zero-spacing estimates;
%   * compares stable numerical profiles with the hypergeometric formula.
%
% The real sine/cosine multiplicity is two for |m|>0 and one for m=0.
% The global tables and plots below list DISTINCT (|m|,j) modes unless
% explicitly stated otherwise.
%
% Requirements:
%   MATLAB R2019b or later for tiledlayout (older releases can replace it
%   by subplot).  The analytic hypergeometric comparison uses hypergeom
%   from the Symbolic Math Toolbox.  All quadrature code is included below.
%
% Working-note code, September 2026.

clearvars; close all; clc;

%% -------------------------- parameters -------------------------------
N              = 900;    % Gauss-Legendre nodes in hyperbolic radius
xMax           = 18.0;   % hyperbolic truncation
mMax           = 10;     % angular blocks |m|=0,...,mMax
nModesPerBlock = 18;     % leading finite-section modes in each block
nDisplay       = 32;     % number of disk modes displayed
displayRows    = 4;
displayCols    = 8;
zeroWindow     = 15.0;   % count radial zeros only on 0<x<zeroWindow
comparisonJ    = 4;      % radial index for analytic comparisons
comparisonM    = [0 1 2 4];
saveFigures    = false;  % true writes PNG files in the current directory

if nDisplay > displayRows*displayCols
    error('nDisplay exceeds displayRows*displayCols.');
end
if nModesPerBlock >= N
    error('nModesPerBlock must be smaller than N.');
end

%% ---------------- quadrature and coordinate transforms ---------------
% The common lgwt implementation returns nodes in descending order.
% Sorting is immaterial to quadrature, but is needed by diff, zero
% detection, interpolation and left-to-right plotting.  The weights must
% be permuted by exactly the same ordering.
[x, wx] = lgwt(N, 0, xMax);
[x, idx] = sort(x(:), 'ascend');
wx = wx(idx);

r = tanh(x/2);
t = r.^2;
omega = 1-r.^2;

% Physical radial measure is r dr.  Since dr/dx=(1-r^2)/2,
% r dr = r(1-r^2) dx/2.
wRad = wx .* r .* omega / 2;

% The physical angular block has kernel
% 2 (r rho)^m/(1-r^2 rho^2) with measure rho d rho.  Equivalently,
% after t=r^2 it is the kernel (ts)^(m/2)/(1-ts) with measure ds.
% The symmetric Nystrom matrix below uses the physical r dr weights.
sqrtW = sqrt(wRad);
TT = t*t.';

%% -------------------------- block spectra -----------------------------
blocks = cell(mMax+1,1);
modeTemplate = struct( ...
    'm', [], 'j', [], 'lambda', [], 'sigma', [], ...
    'muLambda', [], 'zeroCount', [], 'muZeros', [], ...
    'xZeros', [], 'phi', [], 'R', [], 'v', []);
allModes = repmat(modeTemplate, 0, 1);

fprintf('Computing %d blocks, N=%d, xMax=%.1f ...\n',mMax+1,N,xMax);

for m = 0:mMax
    kernel = 2 * TT.^(m/2) ./ (1-TT);
    A = (sqrtW*sqrtW.') .* kernel;
    A = (A+A.')/2;

    [V, lambda] = leading_symmetric_eigenpairs(A,nModesPerBlock);

    % Convert Euclidean Nystrom eigenvectors into radial functions in
    % L^2((0,1),r dr).  Column signs are mathematically arbitrary.
    phi = V ./ sqrtW;

    blockModes = repmat(modeTemplate,nModesPerBlock,1);
    for j = 1:nModesPerBlock
        phij = phi(:,j);

        % U f=(1-r^2)f/2 is the unitary Euclidean-to-hyperbolic
        % conjugation.  R is therefore the regular hyperbolic K-type.
        Rj = (omega/2).*phij;
        vj = sqrt(sinh(x)).*Rj;  % Liouville-stable profile

        % Fix a reproducible overall sign by making the largest stable
        % profile sample positive.  Omitting this changes no mathematics.
        [~, imax] = max(abs(vj));
        if vj(imax) < 0
            phij = -phij;
            Rj = -Rj;
            vj = -vj;
        end

        lam = lambda(j);
        if lam > 0 && lam <= pi*(1+1e-10)
            muLambda = acosh(max(1,pi/lam))/pi;
        else
            muLambda = NaN;
        end

        [xZeros, zeroCount, muZeros] = ...
            zero_spacing_estimate(x,vj,zeroWindow);

        blockModes(j).m = m;
        blockModes(j).j = j;
        blockModes(j).lambda = lam;
        blockModes(j).sigma = sqrt(max(lam,0));
        blockModes(j).muLambda = muLambda;
        blockModes(j).zeroCount = zeroCount;
        blockModes(j).muZeros = muZeros;
        blockModes(j).xZeros = xZeros;
        blockModes(j).phi = phij;
        blockModes(j).R = Rj;
        blockModes(j).v = vj;

        allModes(end+1,1) = blockModes(j); %#ok<SAGROW>
    end

    blocks{m+1} = blockModes;
    fprintf('  |m|=%2d: lambda_1=% .8e, sigma_1=% .8e\n', ...
        m,lambda(1),sqrt(lambda(1)));
end

%% ----------------------- global distinct ordering ---------------------
[~, globalOrder] = sort([allModes.sigma], 'descend');
allModes = allModes(globalOrder);

nPrint = min(20,numel(allModes));
fprintf('\nFirst %d distinct (m,j) modes in global order:\n',nPrint);
fprintf('%5s %13s %13s %4s %4s %7s %12s\n', ...
    'rank','sigma','lambda','|m|','j','zeros','muZeros');
for k = 1:nPrint
    M = allModes(k);
    if isnan(M.muZeros)
        muText = '          --';
    else
        muText = sprintf('%12.6f',M.muZeros);
    end
    fprintf('%5d %13.6e %13.6e %4d %4d %7d %s\n', ...
        k,M.sigma,M.lambda,M.m,M.j,M.zeroCount,muText);
end

%% ---------------------- figure: disk modes ----------------------------
fig1 = figure('Color','w','Position',[40 40 1800 880]);
tl = tiledlayout(displayRows,displayCols, ...
    'TileSpacing','compact','Padding','compact');

nGrid = 301;
grid1 = linspace(-1,1,nGrid);
[Xg,Yg] = meshgrid(grid1,grid1);
Rg = hypot(Xg,Yg);
Thg = atan2(Yg,Xg);
inside = Rg <= 1;

for k = 1:min(nDisplay,numel(allModes))
    nexttile;
    M = allModes(k);

    % Add endpoint values before interpolation.  At r=0 the regular
    % K-type is zero for m>0 and finite for m=0.  At r=1 it decays to zero.
    if M.m == 0
        R0 = M.R(1);
    else
        R0 = 0;
    end
    rInterp = [0; r; 1];
    RInterp = [R0; M.R; 0];

    radialField = interp1(rInterp,RInterp,Rg,'pchip',0);
    field = radialField .* cos(M.m*Thg);
    field(~inside) = NaN;
    maxAbs = max(abs(field(inside)),[],'omitnan');
    if maxAbs > 0
        field = field/maxAbs;
    end

    imagesc(grid1,grid1,field);
    axis image off xy;
    caxis([-1 1]);
    title(sprintf('#%d: |m|=%d, j=%d, z=%d, \\sigma=%.3g', ...
        k,M.m,M.j,M.zeroCount,M.sigma), ...
        'FontSize',8,'Interpreter','tex');
end
colormap(fig1,blue_white_red(257));
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Limits = [-1 1];
title(tl,{ 'Multiplicity-normalized generalized singular functions', ...
           'displayed as (1-r^2)\phi(r) cos(m\theta)/2' }, ...
           'FontWeight','normal');

if saveFigures
    export_figure(fig1,'nonradial_32_singular_functions.png');
end

%% ---------------- figure: global spectral ordering -------------------
fig2 = figure('Color','w','Position',[80 60 1150 850]);
tl2 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

rank = (1:numel(allModes)).';
sig = [allModes.sigma].';
mList = [allModes.m].';
zList = [allModes.zeroCount].';
muZ = [allModes.muZeros].';
sigPred = NaN(size(muZ));
validMu = isfinite(muZ) & muZ >= 0;
sigPred(validMu) = sqrt(pi./cosh(pi*muZ(validMu)));

nexttile;
scatter(rank,sig,34,mList,'filled'); hold on;
scatter(rank(validMu),sigPred(validMu),38,'o', ...
    'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',1.0);
grid on; box on;
ylabel('\sigma=\lambda^{1/2}','Interpreter','tex');
title('Truncated singular values, globally ordered and coloured by |m|');
legend('numerical','prediction from zero spacing','Location','southwest');
cb2 = colorbar; cb2.Label.String='|m|';
caxis([0 mMax]);
colormap(fig2,parula(max(64,mMax+1)));

nexttile;
stairs(rank,mList,'LineWidth',1.4);
grid on; box on;
ylabel('|m|'); ylim([-0.5 mMax+0.5]);
title('Angular frequency jumps in the global ordering');

nexttile;
stairs(rank,zList,'LineWidth',1.4);
grid on; box on;
xlabel('global singular-value rank');
ylabel('radial zeros');
title(sprintf('Radial zero count on 0<x<%.1f',zeroWindow));

title(tl2,sprintf('Hyperbolic truncation x_{max}=%.1f, angular truncation |m|\\leq%d', ...
    xMax,mMax),'FontWeight','normal');

if saveFigures
    export_figure(fig2,'nonradial_global_singular_values.png');
end

%% --------------- figure: analytic hypergeometric profiles ------------
fig3 = figure('Color','w','Position',[120 60 1150 850]);
tl3 = tiledlayout(numel(comparisonM),1, ...
    'TileSpacing','compact','Padding','compact');
analyticAvailable = true;

for q = 1:numel(comparisonM)
    m = comparisonM(q);
    M = blocks{m+1}(comparisonJ);
    mask = x <= zeroWindow;
    xx = x(mask);
    vNum = M.v(mask);

    nexttile;
    plot(xx,vNum/max(abs(vNum)),'k-','LineWidth',1.6); hold on;

    try
        RExact = analytic_hyperbolic_K_type(m,M.muLambda,xx);
        vExact = sqrt(sinh(xx)).*RExact;

        % Least-squares sign/amplitude alignment.  The analytic and
        % numerical generalized eigenfunctions have arbitrary scale.
        alpha = real(vExact' * vNum)/(real(vExact' * vExact));
        vExact = alpha*vExact;
        commonScale = max(abs(vNum));
        plot(xx,vExact/commonScale,'r--','LineWidth',1.4);
        legend('Nystrom','{}_2F_1 prediction','Location','northeast');
    catch ME
        analyticAvailable = false;
        warning('Analytic comparison skipped: %s',ME.message);
        legend('Nystrom','Location','northeast');
    end

    grid on; box on;
    ylabel('normalised v');
    title(sprintf('|m|=%d, radial index j=%d, \\mu_{\\lambda}=%.6g', ...
        m,comparisonJ,M.muLambda),'Interpreter','tex');
end
xlabel(tl3,'hyperbolic radius x');
ylabel(tl3,'normalised v');
title(tl3,'Stable profiles v=sqrt(sinh x)(1-r^2)\phi/2', ...
    'FontWeight','normal');

if saveFigures
    export_figure(fig3,'nonradial_analytic_comparison.png');
end

if ~analyticAvailable
    fprintf(['\nThe numerical calculation completed, but the analytic ', ...
        'hypergeometric comparison requires hypergeom.\n']);
end

%% -------------------------- explanatory note -------------------------
fprintf(['\nNote on low-|m| radial pictures:\n', ...
    'For m=0 and j=2,3, the nodal circles may be almost invisible in ', ...
    'a Euclidean disk plot.\n', ...
    'The oscillations are roughly uniform in hyperbolic x, whereas ', ...
    'r=tanh(x/2) compresses\n', ...
    'large-x oscillations into thin annuli next to r=1.  The displayed ', ...
    'quantity (1-r^2)phi/2\n', ...
    'also decays like (sinh x)^(-1/2).  The sign of every mode is arbitrary.\n']);

%% =====================================================================
%% Local functions
%% =====================================================================

function [V,lambda] = leading_symmetric_eigenpairs(A,k)
% Leading eigenpairs of a real symmetric matrix, sorted decreasingly.
    n = size(A,1);
    opts.issym = true;
    opts.isreal = true;
    opts.tol = 1e-11;
    opts.maxit = 3000;
    opts.disp = 0;

    if k >= n-1
        [V,D] = eig(A,'vector');
    else
        try
            [V,D] = eigs(A,k,'largestreal',opts);
            D = diag(D);
        catch
            try
                [V,D] = eigs(A,k,'la',opts);
                D = diag(D);
            catch
                warning('eigs failed; falling back to full eig.');
                [V,D] = eig(A,'vector');
            end
        end
    end

    [lambda,ord] = sort(real(D),'descend');
    V = real(V(:,ord));
    lambda = lambda(1:k);
    V = V(:,1:k);
end

function [xZeros,nZeros,muZeros] = zero_spacing_estimate(x,v,xWindow)
% Detect sign-changing zeros and estimate mu from mean zero spacing.
% At least two zeros are needed for one spacing, hence NaN is deliberate
% for the first radial modes.
    mask = x <= xWindow;
    xx = x(mask);
    yy = v(mask);

    amp = max(abs(yy));
    if amp == 0 || ~isfinite(amp)
        xZeros = [];
        nZeros = 0;
        muZeros = NaN;
        return;
    end

    % Ignore the numerically tiny regular tail at x=0, especially for m>0.
    first = find(abs(yy) > 1e-8*amp,1,'first');
    if isempty(first) || first >= numel(yy)
        xZeros = [];
        nZeros = 0;
        muZeros = NaN;
        return;
    end
    xx = xx(first:end);
    yy = yy(first:end);

    ind = find(yy(1:end-1).*yy(2:end) < 0);
    xZeros = xx(ind) - yy(ind).*(xx(ind+1)-xx(ind)) ./ ...
        (yy(ind+1)-yy(ind));
    nZeros = numel(xZeros);

    if nZeros >= 2
        muZeros = pi/mean(diff(xZeros));
    else
        muZeros = NaN;
    end
end

function R = analytic_hyperbolic_K_type(m,mu,x)
% Regular K-type of the hyperbolic Laplacian:
%
% R_{m,mu}(x) = tanh(x/2)^m *
%   2F1(1/2+i mu,1/2-i mu;m+1;-sinh(x/2)^2).
%
% It is real for real x and mu.  hypergeom is evaluated pointwise because
% vector behaviour varies between MATLAB releases.
    if ~isfinite(mu)
        error('mu must be finite.');
    end
    rr = tanh(x/2);
    zz = -sinh(x/2).^2;
    R = zeros(size(x));
    aa = [0.5+1i*mu, 0.5-1i*mu];
    for k = 1:numel(x)
        val = hypergeom(aa,m+1,zz(k));
        R(k) = rr(k)^m * real(double(val));
    end
end

function cmap = blue_white_red(n)
% Symmetric blue-white-red colormap without external dependencies.
    if nargin < 1, n = 257; end
    n1 = floor(n/2);
    n2 = n-n1;
    left = [linspace(0,1,n1).', linspace(0,1,n1).', ones(n1,1)];
    right = [ones(n2,1), linspace(1,0,n2).', linspace(1,0,n2).'];
    cmap = [left; right];
end

function export_figure(fig,name)
% Export with a fallback for older MATLAB releases.
    try
        exportgraphics(fig,name,'Resolution',220);
    catch
        print(fig,name,'-dpng','-r220');
    end
end

function [x,w] = lgwt(N,a,b)
%LGWT Gauss-Legendre nodes and weights on [a,b].
% Newton iteration for the roots of P_N.  This standard implementation is
% included to make the script self-contained.  It returns nodes in the
% order inherited from the cosine initial guesses, commonly descending;
% the main script sorts them and applies the same permutation to weights.

    if N < 1 || N ~= floor(N)
        error('N must be a positive integer.');
    end
    if ~(isfinite(a) && isfinite(b) && a < b)
        error('Require finite a<b.');
    end

    N0 = N-1;
    N1 = N0+1;
    N2 = N0+2;

    xu = linspace(-1,1,N1).';
    y = cos((2*(0:N0).'+1)*pi/(2*N0+2)) ...
        + (0.27/N1)*sin(pi*xu*N0/N2);

    L = zeros(N1,N2);
    Lp = zeros(N1,1);
    y0 = 2*ones(size(y));

    while max(abs(y-y0)) > 4*eps(max(1,max(abs(y))))
        L(:,1) = 1;
        L(:,2) = y;
        for k = 2:N1
            L(:,k+1) = ((2*k-1).*y.*L(:,k) ...
                -(k-1).*L(:,k-1))/k;
        end
        Lp = N2*(L(:,N1)-y.*L(:,N2))./(1-y.^2);
        y0 = y;
        y = y0-L(:,N2)./Lp;
    end

    x = (a*(1-y)+b*(1+y))/2;
    w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end
