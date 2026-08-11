function [pad, z, par] = TPJB_BinarySlip_5Pad(makePlots, allowedSlipMasks)
% =========================================================================
% 五瓦可倾瓦径向轴承二元边界滑移稳态润滑计算
% 边界条件：
%   轴颈表面无滑移；
%   五块瓦全部工作表面均具备二元滑移能力；
%   实际滑移区域由 |tau_w| > tau_c 自动判定；
%   压力负值直接截断为零。
% =========================================================================

if nargin < 1
    makePlots = true;
end

if makePlots
    clc;
    close all;
end

fprintf('============================================================\n');
fprintf(' 五瓦可倾瓦径向轴承二元边界滑移稳态计算\n');
fprintf('============================================================\n\n');

%% 1. 轴承结构参数
par.Npad = 5;
par.Db = 30.065e-3;
par.Dj = 29.9945e-3;
par.Rb = par.Db/2;
par.Rj = par.Dj/2;
par.B = 15e-3;
par.C0 = (par.Db-par.Dj)/2;
par.beta = 55*pi/180;
par.thetaPivot = [-90, -18, 54, 126, 198]*pi/180;
par.pivotRatio = 0.5;
par.preload = 0;

%% 2. 润滑油参数
par.nu = 44.1e-6;
par.rho = 871;
par.mu = par.nu*par.rho;
par.T = 40;

%% 3. 固定运行工况
par.n = 30000;
par.omega = 2*pi*par.n/60;
par.U = par.omega*par.Rj;
par.epsilon = 0.4;
par.e = par.epsilon*par.C0;
par.phi = 0;
par.thetaE = -pi/2 + par.phi;

%% 4. 二元边界滑移参数
par.b = 5e-6;
par.tauC = 50e3;
par.slipRelax = 0.30;
par.slipTolScan = 2e-4;
par.slipTolRoot = 2e-6;
par.slipTolFinal = 1e-6;
par.maxSlipIterScan = 80;
par.maxSlipIterRoot = 120;
par.maxSlipIterFinal = 200;
par.maskStableRequired = 3;

%% 5. 网格参数
par.Ntheta = 81;
par.Nz = 31;
z = linspace(-par.B/2, par.B/2, par.Nz);

% allowedSlipMasks只规定哪些网格具有滑移能力；真正发生滑移还需
% 同时满足局部剪应力超过临界剪应力。省略时五块瓦全部允许滑移。
if nargin < 2 || isempty(allowedSlipMasks)
    allowedSlipMasks = repmat( ...
        {true(par.Ntheta,par.Nz)},1,par.Npad);
end

if ~iscell(allowedSlipMasks) || numel(allowedSlipMasks) ~= par.Npad
    error('allowedSlipMasks必须是包含5个逻辑矩阵的元胞数组。');
end

for iPad = 1:par.Npad
    if ~isequal(size(allowedSlipMasks{iPad}),[par.Ntheta,par.Nz])
        error('第%d块瓦的允许滑移掩码尺寸必须为%d x %d。', ...
            iPad,par.Ntheta,par.Nz);
    end
    allowedSlipMasks{iPad} = logical(allowedSlipMasks{iPad});
end

%% 6. 压力方程数值参数
par.sorOmega = 1.60;
par.pressureTolScan = 2e-5;
par.pressureTolRoot = 2e-7;
par.pressureTolFinal = 1e-7;
par.maxPressureIter = 10000;
par.checkEvery = 10;

%% 7. 瓦块摆角求解参数
par.gammaCap = 3.0e-3;
par.gammaSafety = 0.80;
par.NgammaScan = 31;
par.maxTiltIter = 60;
par.gammaTol = 1e-11;
par.momentTol = 1e-7;
par.unloadForceTol = 1e-7;

fprintf('轴承直径              Db = %.6f mm\n', par.Db*1e3);
fprintf('轴颈直径              Dj = %.6f mm\n', par.Dj*1e3);
fprintf('轴承宽度               B = %.6f mm\n', par.B*1e3);
fprintf('半径间隙              C0 = %.6f um\n', par.C0*1e6);
fprintf('动力黏度              mu = %.8f Pa·s\n', par.mu);
fprintf('转速                   n = %.0f r/min\n', par.n);
fprintf('轴颈表面速度           U = %.6f m/s\n', par.U);
fprintf('偏心率           epsilon = %.6f\n', par.epsilon);
fprintf('滑移长度               b = %.6f um\n', par.b*1e6);
fprintf('临界剪应力         tau_c = %.6f kPa\n', par.tauC/1e3);
fprintf('周向/轴向网格           = %d x %d\n\n', par.Ntheta, par.Nz);

%% 8. 分别求解五块瓦
padCell = cell(1, par.Npad);

for iPad = 1:par.Npad
    fprintf('------------------------------------------------------------\n');
    fprintf('正在求解第 %d 块瓦，中心角 = %.3f deg\n', ...
        iPad, par.thetaPivot(iPad)*180/pi);

    padCell{iPad} = solvePadEquilibriumSlip( ...
        iPad,z,par,allowedSlipMasks{iPad});
    q = padCell{iPad};

    fprintf('状态                 = %s\n', q.status);
    fprintf('平衡摆角             = % .9f deg\n', q.gamma*180/pi);
    fprintf('最大压力             = %.6f MPa\n', q.pMax/1e6);
    fprintf('最小油膜厚度         = %.6f um\n', q.hMin*1e6);
    fprintf('瓦块承载力           = %.6f N\n', q.F);
    fprintf('实际滑移面积比例     = %.6f %%\n', 100*q.slipAreaRatio);
    fprintf('最大滑移速度         = %.6f m/s\n', q.slipSpeedMax);
    fprintf('支点无量纲力矩残差   = %.6e\n', q.momentResidual);
    fprintf('压力/滑移迭代次数    = %d / %d\n\n', ...
        q.pressureIterations, q.slipIterations);
end

pad = [padCell{:}];

%% 9. 轴承整体性能
FxTotal = sum([pad.Fx]);
FyTotal = sum([pad.Fy]);
WTotal = hypot(FxTotal, FyTotal);
loadAngle = atan2(FyTotal, FxTotal)*180/pi;

FrictionSignedTotal = sum([pad.frictionSigned]);
FrictionMagnitudeTotal = sum([pad.frictionMagnitude]);
TorqueTotal = abs(FrictionSignedTotal)*par.Rj;
PowerLoss = TorqueTotal*par.omega;

pMaxTotal = max([pad.pMax]);
hMinTotal = min([pad.hMin]);
slipAreaTotal = sum([pad.slipArea]);
padAreaTotal = par.Npad*par.Rb*par.beta*par.B;
slipAreaRatioTotal = slipAreaTotal/padAreaTotal;

fprintf('============================================================\n');
fprintf(' 五瓦轴承整体计算结果\n');
fprintf('============================================================\n');
fprintf('总油膜力Fx                 = % .8f N\n', FxTotal);
fprintf('总油膜力Fy                 = % .8f N\n', FyTotal);
fprintf('总承载力W                  = %.8f N\n', WTotal);
fprintf('承载力方向                 = %.8f deg\n', loadAngle);
fprintf('全轴承最大压力             = %.8f MPa\n', pMaxTotal/1e6);
fprintf('全轴承最小膜厚             = %.8f um\n', hMinTotal*1e6);
fprintf('全轴承实际滑移面积比例     = %.8f %%\n', 100*slipAreaRatioTotal);
fprintf('总有符号周向摩擦力         = %.8f N\n', FrictionSignedTotal);
fprintf('总局部摩擦力幅值积分       = %.8f N\n', FrictionMagnitudeTotal);
fprintf('摩擦转矩                   = %.8f N·m\n', TorqueTotal);
fprintf('摩擦功耗                   = %.8f W\n', PowerLoss);
fprintf('============================================================\n\n');

%% 10. 输出各瓦结果表
padNumber = (1:par.Npad)';
centerAngleDeg = par.thetaPivot(:)*180/pi;
status = {pad.status}';
gammaDeg = [pad.gamma]'*180/pi;
pMaxMPa = [pad.pMax]'/1e6;
hMinUm = [pad.hMin]'*1e6;
FxN = [pad.Fx]';
FyN = [pad.Fy]';
FN = [pad.F]';
slipAreaPercent = 100*[pad.slipAreaRatio]';
maxSlipSpeed = [pad.slipSpeedMax]';
frictionN = [pad.frictionMagnitude]';
momentResidual = [pad.momentResidual]';

resultTable = table( ...
    padNumber, centerAngleDeg, status, gammaDeg, pMaxMPa, hMinUm, ...
    FxN, FyN, FN, slipAreaPercent, maxSlipSpeed, frictionN, ...
    momentResidual, ...
    'VariableNames', {'Pad','CenterAngle_deg','Status','Gamma_deg', ...
    'Pmax_MPa','Hmin_um','Fx_N','Fy_N','Load_N', ...
    'SlipArea_percent','MaxSlipSpeed_mps','Friction_N', ...
    'MomentResidual'});

disp(resultTable);

%% 11. 绘制五张结果图
if makePlots
    plotFivePadSlipRegionMap(pad, z, par);
    plotFivePadCircumferentialCurvesSlip(pad, z, par, 'pressure');
    plotFivePadCircumferentialCurvesSlip(pad, z, par, 'thickness');
    plotFivePadCircumferentialCurvesSlip(pad, z, par, 'friction');
    plotFivePadCircumferentialCurvesSlip(pad, z, par, 'load');
end

end


%% ========================================================================
% 给定轴颈位置时，求解单瓦平衡摆角
% =========================================================================
function pad = solvePadEquilibriumSlip(iPad, z, par, allowedSlipMask)

thetaPivot = par.thetaPivot(iPad);
theta = linspace(thetaPivot-par.beta/2, ...
    thetaPivot+par.beta/2, par.Ntheta)';

hBase = calculateFilmThicknessSlip(theta, 0, thetaPivot, par);
sinTermMax = max(abs(sin(theta-thetaPivot)));

if sinTermMax < 1e-14
    gammaLimit = par.gammaCap;
else
    gammaLimitByFilm = par.gammaSafety*min(hBase)/(par.Rb*sinTermMax);
    gammaLimit = min(par.gammaCap, gammaLimitByFilm);
end

if gammaLimit <= 0
    error('第%d块瓦的摆角搜索范围无效。', iPad);
end

gammaScan = linspace(-gammaLimit, gammaLimit, par.NgammaScan);
Mscan = zeros(size(gammaScan));
Fscan = zeros(size(gammaScan));

for k = 1:numel(gammaScan)
    state = evaluatePadAtGammaSlip(theta, z, thetaPivot, ...
        gammaScan(k), par, 'scan', allowedSlipMask);
    Mscan(k) = state.M;
    Fscan(k) = state.F;
end

loadedLeft = Fscan(1:end-1) > par.unloadForceTol;
loadedRight = Fscan(2:end) > par.unloadForceTol;
momentChange = Mscan(1:end-1).*Mscan(2:end) <= 0;
rootIntervals = find(loadedLeft & loadedRight & momentChange);

if ~isempty(rootIntervals)
    intervalCenters = 0.5*(gammaScan(rootIntervals) + ...
        gammaScan(rootIntervals+1));
    [~, localIndex] = min(abs(intervalCenters));
    iRoot = rootIntervals(localIndex);

    gammaLow = gammaScan(iRoot);
    gammaHigh = gammaScan(iRoot+1);
    Mlow = Mscan(iRoot);

    for iterTilt = 1:par.maxTiltIter
        gammaMid = 0.5*(gammaLow+gammaHigh);
        stateMid = evaluatePadAtGammaSlip(theta, z, thetaPivot, ...
            gammaMid, par, 'root', allowedSlipMask);

        if stateMid.F > par.unloadForceTol
            residual = abs(stateMid.M)/max(stateMid.F*par.Rb, realmin);
        else
            residual = inf;
        end

        if residual < par.momentTol || ...
                abs(gammaHigh-gammaLow) < par.gammaTol
            break;
        end

        if Mlow*stateMid.M <= 0
            gammaHigh = gammaMid;
        else
            gammaLow = gammaMid;
            Mlow = stateMid.M;
        end
    end

    gammaFinal = gammaMid;
    statusFinal = '承载平衡';
else
    unloadedIndex = find(Fscan <= par.unloadForceTol);

    if ~isempty(unloadedIndex)
        [~, iNearest] = min(abs(gammaScan(unloadedIndex)));
        gammaFinal = gammaScan(unloadedIndex(iNearest));
        statusFinal = '卸载';
    else
        normalizedMoment = abs(Mscan)./max(Fscan*par.Rb, realmin);
        [~, iBest] = min(normalizedMoment);
        gammaFinal = gammaScan(iBest);
        statusFinal = '搜索边界受限';
    end
end

stateFinal = evaluatePadAtGammaSlip(theta, z, thetaPivot, ...
    gammaFinal, par, 'final', allowedSlipMask);

pad.index = iPad;
pad.thetaPivot = thetaPivot;
pad.theta = theta;
pad.gamma = gammaFinal;
pad.status = statusFinal;
pad.P = stateFinal.P;
pad.H = stateFinal.H;
pad.h = stateFinal.h;
pad.dpdx = stateFinal.dpdx;
pad.dpdz = stateFinal.dpdz;
pad.slipU = stateFinal.slipU;
pad.slipW = stateFinal.slipW;
pad.slipSpeed = stateFinal.slipSpeed;
pad.slipMask = stateFinal.slipMask;
pad.allowedSlipMask = allowedSlipMask;
pad.tauPadX = stateFinal.tauPadX;
pad.tauPadZ = stateFinal.tauPadZ;
pad.tauPadMagnitude = stateFinal.tauPadMagnitude;
pad.M = stateFinal.M;
pad.Fx = stateFinal.Fx;
pad.Fy = stateFinal.Fy;
pad.F = stateFinal.F;
pad.pMax = max(stateFinal.P(:));
pad.hMin = min(stateFinal.H(:));
pad.frictionSigned = stateFinal.frictionSigned;
pad.frictionMagnitude = stateFinal.frictionMagnitude;
pad.pressureIterations = stateFinal.pressureIterations;
pad.pressureError = stateFinal.pressureError;
pad.slipIterations = stateFinal.slipIterations;
pad.slipError = stateFinal.slipError;

dA = par.Rb*(theta(2)-theta(1))*(z(2)-z(1));
weights = ones(size(stateFinal.slipMask));
weights([1 end],:) = 0.5*weights([1 end],:);
weights(:,[1 end]) = 0.5*weights(:,[1 end]);
pad.slipArea = dA*sum(weights(:).*double(stateFinal.slipMask(:)));
pad.padArea = par.Rb*par.beta*par.B;
pad.slipAreaRatio = pad.slipArea/pad.padArea;
pad.slipSpeedMax = max(stateFinal.slipSpeed(:));

if stateFinal.F > par.unloadForceTol
    pad.momentResidual = abs(stateFinal.M)/ ...
        max(stateFinal.F*par.Rb, realmin);
else
    pad.momentResidual = 0;
end

end


%% ========================================================================
% 给定摆角时，耦合求解压力场和二元滑移速度
% =========================================================================
function state = evaluatePadAtGammaSlip( ...
    theta, z, thetaPivot, gamma, par, accuracyLevel, allowedSlipMask)

switch accuracyLevel
    case 'scan'
        pressureTol = par.pressureTolScan;
        slipTol = par.slipTolScan;
        maxSlipIter = par.maxSlipIterScan;
    case 'root'
        pressureTol = par.pressureTolRoot;
        slipTol = par.slipTolRoot;
        maxSlipIter = par.maxSlipIterRoot;
    otherwise
        pressureTol = par.pressureTolFinal;
        slipTol = par.slipTolFinal;
        maxSlipIter = par.maxSlipIterFinal;
end

h = calculateFilmThicknessSlip(theta, gamma, thetaPivot, par);

if min(h) <= 0
    error('出现非正油膜厚度：gamma=%.6e rad，hmin=%.6e m。', ...
        gamma, min(h));
end

dhdTheta = par.e*sin(theta-par.thetaE) + ...
    par.Rb*gamma*cos(theta-thetaPivot);
dhdx = dhdTheta/par.Rb;

Ntheta = numel(theta);
Nz = numel(z);
H = repmat(h, 1, Nz);

if ~isequal(size(allowedSlipMask),[Ntheta,Nz])
    error('允许滑移掩码尺寸与压力网格不一致。');
end
allowedSlipMask = logical(allowedSlipMask);

P = zeros(Ntheta, Nz);
slipU = zeros(Ntheta, Nz);
slipW = zeros(Ntheta, Nz);
oldMask = false(Ntheta, Nz);
maskStableCount = 0;
slipError = inf;

for iterSlip = 1:maxSlipIter
    PoldOuter = P;
    slipUold = slipU;
    slipWold = slipW;

    [P, ~, ~] = solvePressureSlip( ...
        theta, z, h, dhdx, slipU, slipW, P, par, pressureTol);

    gradients = calculatePressureGradients(P, theta, z, par);

    % 若尚未滑移，固定瓦面处的二维剪应力向量。
    tauTrialX = par.mu*par.U./H - 0.5*H.*gradients.dpdx;
    tauTrialZ = -0.5*H.*gradients.dpdz;
    tauTrialMagnitude = hypot(tauTrialX, tauTrialZ);

    slipMask = allowedSlipMask & (tauTrialMagnitude > par.tauC);

    % 二元滑移的矢量解析更新：
    % u_s = b*h/[mu*(h+b)] * (1-tau_c/|tau_trial|) * tau_trial
    activation = zeros(size(H));
    activation(slipMask) = ...
        1 - par.tauC./tauTrialMagnitude(slipMask);

    coefficient = par.b*H./(par.mu*(H+par.b));
    slipUcalc = coefficient.*activation.*tauTrialX;
    slipWcalc = coefficient.*activation.*tauTrialZ;

    slipU = (1-par.slipRelax)*slipUold + ...
        par.slipRelax*slipUcalc;
    slipW = (1-par.slipRelax)*slipWold + ...
        par.slipRelax*slipWcalc;

    % 未触发区域严格置零，避免欠松弛残留。
    slipU(~slipMask) = 0;
    slipW(~slipMask) = 0;

    slipDelta = hypot(slipU-slipUold, slipW-slipWold);
    slipScale = max([par.U, max(hypot(slipU,slipW),[],'all')]);
    slipError = max(slipDelta,[],'all')/max(slipScale, realmin);

    pressureOuterError = max(abs(P-PoldOuter),[],'all')/ ...
        (max(P,[],'all')+1);

    if isequal(slipMask, oldMask)
        maskStableCount = maskStableCount+1;
    else
        maskStableCount = 0;
    end
    oldMask = slipMask;

    if slipError < slipTol && pressureOuterError < 5*pressureTol && ...
            maskStableCount >= par.maskStableRequired
        break;
    end
end

% 用最终滑移速度再求一次压力，并重新计算最终场量。
[P, pressureIterations, pressureError] = solvePressureSlip( ...
    theta, z, h, dhdx, slipU, slipW, P, par, pressureTol);

metrics = calculatePadMetricsSlip(theta, z, thetaPivot, ...
    P, h, slipU, slipW, par);

state.P = P;
state.H = H;
state.h = h;
state.dpdx = metrics.dpdx;
state.dpdz = metrics.dpdz;
state.slipU = slipU;
state.slipW = slipW;
state.slipSpeed = hypot(slipU, slipW);
state.slipMask = metrics.slipMask;
state.tauPadX = metrics.tauPadX;
state.tauPadZ = metrics.tauPadZ;
state.tauPadMagnitude = metrics.tauPadMagnitude;
state.M = metrics.M;
state.Fx = metrics.Fx;
state.Fy = metrics.Fy;
state.F = metrics.F;
state.frictionSigned = metrics.frictionSigned;
state.frictionMagnitude = metrics.frictionMagnitude;
state.pressureIterations = pressureIterations;
state.pressureError = pressureError;
state.slipIterations = iterSlip;
state.slipError = slipError;

end


%% ========================================================================
% 油膜厚度
% =========================================================================
function h = calculateFilmThicknessSlip(theta, gamma, thetaPivot, par)

h = par.C0 ...
    - par.e*cos(theta-par.thetaE) ...
    + par.Rb*gamma*sin(theta-thetaPivot);

end


%% ========================================================================
% 给定滑移速度时求解扩展Reynolds方程
% =========================================================================
function [P, iter, errorP] = solvePressureSlip( ...
    theta, z, h, dhdx, slipU, slipW, Pinit, par, tolerance)

Ntheta = numel(theta);
Nz = numel(z);
dtheta = theta(2)-theta(1);
dx = par.Rb*dtheta;
dz = z(2)-z(1);

H = repmat(h, 1, Nz);
H3 = H.^3;

H3E = 0.5*(H3(2:end-1,2:end-1)+H3(3:end,2:end-1));
H3W = 0.5*(H3(2:end-1,2:end-1)+H3(1:end-2,2:end-1));
H3N = 0.5*(H3(2:end-1,2:end-1)+H3(2:end-1,3:end));
H3S = 0.5*(H3(2:end-1,2:end-1)+H3(2:end-1,1:end-2));

AE = H3E/dx^2;
AW = H3W/dx^2;
AN = H3N/dz^2;
AS = H3S/dz^2;
AP = AE+AW+AN+AS;

baseSource = repmat(6*par.mu*par.U*dhdx(2:end-1), 1, Nz-2);

HU = H.*slipU;
HW = H.*slipW;
slipDivergence = ...
    (HU(3:end,2:end-1)-HU(1:end-2,2:end-1))/(2*dx) + ...
    (HW(2:end-1,3:end)-HW(2:end-1,1:end-2))/(2*dz);

RHS = baseSource + 6*par.mu*slipDivergence;

[I,J] = ndgrid(2:Ntheta-1,2:Nz-1);
redMask = mod(I+J,2)==0;
blackMask = ~redMask;

if isequal(size(Pinit), [Ntheta,Nz]) && all(isfinite(Pinit),'all')
    P = max(Pinit,0);
else
    P = zeros(Ntheta,Nz);
end

P([1 end],:) = 0;
P(:,[1 end]) = 0;
errorP = inf;

for iter = 1:par.maxPressureIter
    if mod(iter-1,par.checkEvery)==0
        Pold = P;
    end

    Pgs = ( ...
        AE.*P(3:end,2:end-1) + AW.*P(1:end-2,2:end-1) + ...
        AN.*P(2:end-1,3:end) + AS.*P(2:end-1,1:end-2) - ...
        RHS)./AP;

    Pint = P(2:end-1,2:end-1);
    Pcand = (1-par.sorOmega)*Pint + par.sorOmega*Pgs;
    Pcand = max(Pcand,0);
    Pint(redMask) = Pcand(redMask);
    P(2:end-1,2:end-1) = Pint;

    Pgs = ( ...
        AE.*P(3:end,2:end-1) + AW.*P(1:end-2,2:end-1) + ...
        AN.*P(2:end-1,3:end) + AS.*P(2:end-1,1:end-2) - ...
        RHS)./AP;

    Pint = P(2:end-1,2:end-1);
    Pcand = (1-par.sorOmega)*Pint + par.sorOmega*Pgs;
    Pcand = max(Pcand,0);
    Pint(blackMask) = Pcand(blackMask);
    P(2:end-1,2:end-1) = Pint;

    P([1 end],:) = 0;
    P(:,[1 end]) = 0;

    if mod(iter,par.checkEvery)==0
        errorP = max(abs(P-Pold),[],'all')/(max(P,[],'all')+1);
        if errorP < tolerance
            break;
        end
    end
end

end


%% ========================================================================
% 压力梯度
% =========================================================================
function gradients = calculatePressureGradients(P, theta, z, par)

dx = par.Rb*(theta(2)-theta(1));
dz = z(2)-z(1);

dpdx = zeros(size(P));
dpdz = zeros(size(P));

dpdx(2:end-1,:) = (P(3:end,:)-P(1:end-2,:))/(2*dx);
dpdx(1,:) = (P(2,:)-P(1,:))/dx;
dpdx(end,:) = (P(end,:)-P(end-1,:))/dx;

dpdz(:,2:end-1) = (P(:,3:end)-P(:,1:end-2))/(2*dz);
dpdz(:,1) = (P(:,2)-P(:,1))/dz;
dpdz(:,end) = (P(:,end)-P(:,end-1))/dz;

gradients.dpdx = dpdx;
gradients.dpdz = dpdz;

end


%% ========================================================================
% 载荷、力矩、实际剪应力及摩擦
% =========================================================================
function metrics = calculatePadMetricsSlip( ...
    theta, z, thetaPivot, P, h, slipU, slipW, par)

Nz = numel(z);
Theta = repmat(theta,1,Nz);
H = repmat(h,1,Nz);

gradients = calculatePressureGradients(P,theta,z,par);
dpdx = gradients.dpdx;
dpdz = gradients.dpdz;

integrandFx = -P.*cos(Theta);
integrandFy = -P.*sin(Theta);

Fx = par.Rb*trapz(z,trapz(theta,integrandFx,1),2);
Fy = par.Rb*trapz(z,trapz(theta,integrandFy,1),2);
F = hypot(Fx,Fy);

momentArm = par.Rb*sin(Theta-thetaPivot);
M = par.Rb*trapz(z,trapz(theta,P.*momentArm,1),2);

% 固定瓦面实际二维剪应力。
tauPadX = par.mu*(par.U-slipU)./H - 0.5*H.*dpdx;
tauPadZ = -par.mu*slipW./H - 0.5*H.*dpdz;
tauPadMagnitude = hypot(tauPadX,tauPadZ);

slipMask = hypot(slipU,slipW) > max(1e-12,1e-10*par.U);

frictionSigned = par.Rb* ...
    trapz(z,trapz(theta,tauPadX,1),2);
frictionMagnitude = par.Rb* ...
    trapz(z,trapz(theta,tauPadMagnitude,1),2);

metrics.dpdx = dpdx;
metrics.dpdz = dpdz;
metrics.tauPadX = tauPadX;
metrics.tauPadZ = tauPadZ;
metrics.tauPadMagnitude = tauPadMagnitude;
metrics.slipMask = slipMask;
metrics.Fx = Fx;
metrics.Fy = Fy;
metrics.F = F;
metrics.M = M;
metrics.frictionSigned = frictionSigned;
metrics.frictionMagnitude = frictionMagnitude;

end


%% ========================================================================
% 图1：径向轴承周向展开的实际滑移区域
% =========================================================================
function plotFivePadSlipRegionMap(pad,z,par)

thetaZero = par.thetaPivot(1)-par.beta/2;
fig = figure('Color','w','Position',[120,80,1050,560], ...
    'Name','图1 五瓦实际滑移区域');
ax = axes('Parent',fig,'Position',[0.095,0.15,0.78,0.77]);
hold(ax,'on');

for iPad = 1:par.Npad
    thetaDeg = mod((pad(iPad).theta(:)-thetaZero)*180/pi,360);
    zStar = z/par.B;
    slipField = double(pad(iPad).slipMask).';
    imagesc(ax,thetaDeg,zStar,slipField);

    text(ax,mean(thetaDeg),0,sprintf('%d号瓦',iPad), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Microsoft YaHei','FontSize',12, ...
        'FontWeight','bold','Color','w', ...
        'BackgroundColor',[0.12 0.12 0.12],'Margin',2);
end

set(ax,'YDir','normal');
colormap(ax,[0.08 0.32 0.88; 1.00 0.82 0.06]);
clim(ax,[0 1]);
xlim(ax,[0 360]);
ylim(ax,[-0.5 0.5]);
xticks(ax,0:30:360);
yticks(ax,-0.5:0.25:0.5);
grid(ax,'on');
box(ax,'on');
set(ax,'FontName','Times New Roman','FontSize',12, ...
    'LineWidth',1.1,'TickDir','out','Layer','top');

xlabel(ax,'周向展开角度 \theta^* / (°)（1号瓦入口为0°）', ...
    'FontName','Microsoft YaHei','FontSize',14);
ylabel(ax,'轴向无量纲坐标  z/B', ...
    'FontName','Microsoft YaHei','FontSize',14);
title(ax,'五瓦可倾瓦径向轴承实际滑移区域', ...
    'FontName','Microsoft YaHei','FontSize',15,'FontWeight','bold');

cb = colorbar(ax,'Location','eastoutside');
cb.Ticks = [0 1];
cb.TickLabels = {'未滑移','实际滑移'};
cb.FontName = 'Microsoft YaHei';
cb.FontSize = 11;

annotation(fig,'textbox',[0.27,0.015,0.46,0.055], ...
    'String','(a) 五瓦实际滑移区域', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FontName','Microsoft YaHei','FontSize',14);

saveSlipFigure(fig,'TPJB_BinarySlip_Fig1_SlipRegion.png');
end


%% ========================================================================
% 图2至图5：五块径向瓦的周向二维分布曲线
% =========================================================================
function plotFivePadCircumferentialCurvesSlip(pad,z,par,fieldType)

thetaZero = par.thetaPivot(1)-par.beta/2;
thetaCurve = cell(1,par.Npad);
valueCurve = cell(1,par.Npad);
allValues = [];

for iPad = 1:par.Npad
    thetaCurve{iPad} = mod( ...
        (pad(iPad).theta(:)-thetaZero)*180/pi,360);

    switch fieldType
        case 'pressure'
            midZ = round((numel(z)+1)/2);
            valueCurve{iPad} = pad(iPad).P(:,midZ)/1e6;

        case 'thickness'
            valueCurve{iPad} = pad(iPad).h(:)*1e6;

        case 'friction'
            % 滑移模型为二维剪应力，使用切向和轴向合成幅值。
            valueCurve{iPad} = par.Rb* ...
                trapz(z,pad(iPad).tauPadMagnitude,2);

        case 'load'
            valueCurve{iPad} = par.Rb* ...
                trapz(z,pad(iPad).P,2);

        otherwise
            error('未知二维绘图变量类型：%s',fieldType);
    end

    allValues = [allValues;valueCurve{iPad}(:)]; %#ok<AGROW>
end

valueMin = min(allValues);
valueMax = max(allValues);
valueRange = valueMax-valueMin;
if valueRange <= eps(max(abs([valueMin valueMax])))
    valueRange = max(abs(valueMax),1);
end

fig = figure('Color','w','Position',[130,85,1050,620]);
ax = axes('Parent',fig,'Position',[0.105,0.16,0.76,0.76]);
hold(ax,'on');
curveColors = lines(par.Npad);
lineStyles = {'-','--','-.',':','-'};
curveHandles = gobjects(par.Npad,1);

for iPad = 1:par.Npad
    thetaDeg = thetaCurve{iPad};
    fieldValue = valueCurve{iPad};
    curveHandles(iPad) = plot(ax,thetaDeg,fieldValue, ...
        'Color',curveColors(iPad,:), ...
        'LineStyle',lineStyles{iPad}, ...
        'LineWidth',2.2, ...
        'DisplayName',sprintf('%d号瓦',iPad));

    if strcmp(fieldType,'thickness')
        labelIndex = round((numel(thetaDeg)+1)/2);
    else
        [~,labelIndex] = max(fieldValue);
    end
    labelX = thetaDeg(labelIndex);
    labelY = fieldValue(labelIndex) + ...
        (0.035+0.007*mod(iPad,2))*valueRange;
    text(ax,labelX,labelY,sprintf('%d号瓦',iPad), ...
        'Color',curveColors(iPad,:), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName','Microsoft YaHei','FontSize',10.5, ...
        'FontWeight','bold','BackgroundColor','w', ...
        'Margin',1,'Clipping','off');
end

grid(ax,'on');
box(ax,'on');
xlim(ax,[0 360]);
xticks(ax,0:30:360);
set(ax,'FontName','Times New Roman','FontSize',12, ...
    'LineWidth',1.1,'TickDir','out', ...
    'XMinorGrid','on','YMinorGrid','on');
xlabel(ax,'周向展开角度 \theta^* / (°)（1号瓦入口为0°）', ...
    'FontName','Microsoft YaHei','FontSize',14);

switch fieldType
    case 'pressure'
        ylabelText = '轴向中截面油膜压力  p / MPa';
        titleText = '二元滑移五瓦轴向中截面油膜压力分布';
        caption = '(b) 五瓦油膜压力分布';
        fileName = 'TPJB_BinarySlip_Fig2_Pressure.png';
        ylim(ax,[0 valueMax+0.14*valueRange]);

    case 'thickness'
        ylabelText = '油膜厚度  h / \mum';
        titleText = '二元滑移五瓦油膜厚度分布';
        caption = '(c) 五瓦油膜厚度分布';
        fileName = 'TPJB_BinarySlip_Fig3_FilmThickness.png';
        ylim(ax,[max(0,valueMin-0.10*valueRange), ...
            valueMax+0.16*valueRange]);

    case 'friction'
        ylabelText = '局部摩擦力幅值  dF_f/d\theta / (N·rad^{-1})';
        titleText = '二元滑移五瓦周向摩擦力分布';
        caption = '(d) 五瓦摩擦力分布';
        fileName = 'TPJB_BinarySlip_Fig4_Friction.png';
        ylim(ax,[0 valueMax+0.14*valueRange]);

    case 'load'
        ylabelText = '局部法向承载力  dW/d\theta / (N·rad^{-1})';
        titleText = '二元滑移五瓦周向承载力分布';
        caption = '(e) 五瓦承载力分布';
        fileName = 'TPJB_BinarySlip_Fig5_Load.png';
        ylim(ax,[0 valueMax+0.14*valueRange]);
end

ylabel(ax,ylabelText,'FontName','Microsoft YaHei','FontSize',14);
title(ax,titleText,'FontName','Microsoft YaHei', ...
    'FontSize',15,'FontWeight','bold');
legend(ax,curveHandles,'Location','northeastoutside', ...
    'FontName','Microsoft YaHei','FontSize',10.5,'Box','on');

annotation(fig,'textbox',[0.27,0.015,0.46,0.055], ...
    'String',caption,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','EdgeColor','none', ...
    'FontName','Microsoft YaHei','FontSize',14);

saveSlipFigure(fig,fileName);
end


%% ========================================================================
% 保存结果图到当前全滑移源文件所在目录
% =========================================================================
function saveSlipFigure(fig,fileName)

outputDirectory = fileparts(mfilename('fullpath'));
outputFigure = fullfile(outputDirectory,fileName);
exportgraphics(fig,outputFigure,'Resolution',300);
fprintf('结果图已保存至：%s\n',outputFigure);
end
