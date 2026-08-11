function [pad, z, par] = TPJB_NoSlip_5Pad(makePlots)
% =========================================================================
% 五瓦可倾瓦径向轴承无滑移稳态润滑计算
%
%   采用有限差分法和红黑投影SOR法求解Reynolds方程
%   计算压力、膜厚、承载力和瓦面摩擦剪应力
% =========================================================================

if nargin < 1
    makePlots = true;
end

if makePlots
    clc;
    close all;
end

fprintf('============================================================\n');
fprintf(' 五瓦可倾瓦径向轴承无滑移稳态润滑计算\n');
fprintf('============================================================\n\n');

%% 1. 轴承结构参数

par.Npad = 5;                    % 瓦块数量

par.Db = 30.065e-3;              % 轴承直径，m
par.Dj = 29.9945e-3;             % 轴颈直径，m

par.Rb = par.Db / 2;             % 轴承半径，m
par.Rj = par.Dj / 2;             % 轴颈半径，m

par.B  = 15e-3;                  % 轴承宽度，m
par.C0 = (par.Db - par.Dj) / 2;  % 半径间隙，m

par.beta = 55 * pi / 180;        % 单瓦张角，rad

% 五块瓦中心角
par.thetaPivot = [-90, -18, 54, 126, 198] * pi / 180;

% 中心支点、无预载
par.pivotRatio = 0.5;
par.preload    = 0;

%% 2. 润滑油参数

par.nu  = 44.1e-6;               % 运动黏度，m^2/s
par.rho = 871;                    % 密度，kg/m^3
par.mu  = par.nu * par.rho;       % 动力黏度，Pa·s
par.T   = 40;                     % 等温参考温度，摄氏度

%% 3. 固定运行工况

par.n = 30000;                    % 转速，r/min
par.omega = 2*pi*par.n/60;        % 角速度，rad/s
par.U = par.omega * par.Rj;       % 轴颈表面线速度，m/s

par.epsilon = 0.4;                % 固定偏心率
par.e = par.epsilon * par.C0;     % 偏心距，m

% 本程序规定：
% phi=0表示轴颈中心沿着1号瓦中心方向偏移，即沿全局-y方向偏移
par.phi = 0;
par.thetaE = -pi/2 + par.phi;     % 偏心方向角，rad

%% 4. 网格参数

par.Ntheta = 81;                  % 每块瓦周向网格数
par.Nz     = 31;                  % 每块瓦轴向网格数

z = linspace(-par.B/2, par.B/2, par.Nz);

%% 5. Reynolds方程数值参数

par.sorOmega = 1.60;              % 红黑SOR松弛因子

par.pressureTolScan  = 1e-5;      % 摆角粗扫描时压力精度
par.pressureTolRoot  = 1e-6;      % 摆角二分求根时压力精度
par.pressureTolFinal = 1e-7;      % 最终压力场精度

par.maxPressureIter = 10000;      % 压力最大迭代次数
par.checkEvery      = 10;         % 每隔多少步检查压力收敛

%% 6. 瓦块摆角求解参数

par.gammaCap = 3.0e-3;            % 最大搜索摆角绝对值，rad

% 保证摆角搜索时不会出现负膜厚
par.gammaSafety = 0.80;

par.NgammaScan = 31;              % 摆角粗扫描点数
par.maxTiltIter = 60;             % 摆角二分求根最大次数

par.gammaTol  = 1e-10;            % 摆角区间收敛精度，rad
par.momentTol = 1e-7;             % 无量纲力矩残差

% 小于该值认为瓦块不承载
par.unloadForceTol = 1e-7;        % N

%% 7. 输出参数信息

fprintf('轴承直径              Db = %.6f mm\n', par.Db*1e3);
fprintf('轴颈直径              Dj = %.6f mm\n', par.Dj*1e3);
fprintf('轴承宽度               B = %.6f mm\n', par.B*1e3);
fprintf('半径间隙              C0 = %.6f um\n', par.C0*1e6);
fprintf('单瓦张角            beta = %.6f deg\n', par.beta*180/pi);
fprintf('动力黏度              mu = %.8f Pa·s\n', par.mu);
fprintf('转速                   n = %.0f r/min\n', par.n);
fprintf('轴颈表面速度           U = %.6f m/s\n', par.U);
fprintf('偏心率           epsilon = %.6f\n', par.epsilon);
fprintf('偏心距                 e = %.6f um\n', par.e*1e6);
fprintf('周向网格数        Ntheta = %d\n', par.Ntheta);
fprintf('轴向网格数            Nz = %d\n\n', par.Nz);

%% 8. 分别求解五块瓦

% 每块瓦的求解结果先存入元胞数组。直接把带字段的返回结构体
% 赋给 struct([]) 在部分 MATLAB 版本中会触发
% “在不同结构体之间进行下标赋值”。
padCell = cell(1, par.Npad);

for iPad = 1:par.Npad

    fprintf('------------------------------------------------------------\n');
    fprintf('正在求解第 %d 块瓦，中心角 = %.3f deg\n', ...
        iPad, par.thetaPivot(iPad)*180/pi);

    padCell{iPad} = solvePadEquilibrium(iPad, z, par);

    fprintf('状态                 = %s\n', padCell{iPad}.status);
    fprintf('平衡摆角             = % .9f deg\n', ...
        padCell{iPad}.gamma*180/pi);
    fprintf('最大压力             = %.6f MPa\n', ...
        padCell{iPad}.pMax/1e6);
    fprintf('最小油膜厚度         = %.6f um\n', ...
        padCell{iPad}.hMin*1e6);
    fprintf('瓦块承载力           = %.6f N\n', padCell{iPad}.F);
    fprintf('支点无量纲力矩残差   = %.6e\n', ...
        padCell{iPad}.momentResidual);
    fprintf('最终压力迭代次数     = %d\n\n', ...
        padCell{iPad}.pressureIterations);
end

% 求解完成后再一次性合并。此时五个结构体具有完全相同的字段。
pad = [padCell{:}];

%% 9. 轴承整体性能

FxTotal = sum([pad.Fx]);
FyTotal = sum([pad.Fy]);

WTotal = sqrt(FxTotal^2 + FyTotal^2);
loadAngle = atan2(FyTotal, FxTotal) * 180/pi;

FrictionSignedTotal = sum([pad.frictionSigned]);
FrictionAbsTotal    = sum([pad.frictionMagnitude]);

TorqueTotal = abs(FrictionSignedTotal) * par.Rj;
PowerLoss   = TorqueTotal * par.omega;

pMaxTotal = max([pad.pMax]);
hMinTotal = min([pad.hMin]);

fprintf('============================================================\n');
fprintf(' 五瓦轴承整体计算结果\n');
fprintf('============================================================\n');
fprintf('总油膜力Fx                 = % .8f N\n', FxTotal);
fprintf('总油膜力Fy                 = % .8f N\n', FyTotal);
fprintf('总承载力W                  = %.8f N\n', WTotal);
fprintf('承载力方向                 = %.8f deg\n', loadAngle);
fprintf('全轴承最大压力             = %.8f MPa\n', pMaxTotal/1e6);
fprintf('全轴承最小膜厚             = %.8f um\n', hMinTotal*1e6);
fprintf('总有符号周向摩擦力         = %.8f N\n', FrictionSignedTotal);
fprintf('总局部摩擦力幅值积分       = %.8f N\n', FrictionAbsTotal);
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
frictionN = [pad.frictionMagnitude]';
momentResidual = [pad.momentResidual]';

resultTable = table( ...
    padNumber, ...
    centerAngleDeg, ...
    status, ...
    gammaDeg, ...
    pMaxMPa, ...
    hMinUm, ...
    FxN, ...
    FyN, ...
    FN, ...
    frictionN, ...
    momentResidual, ...
    'VariableNames', { ...
    'Pad', ...
    'CenterAngle_deg', ...
    'Status', ...
    'Gamma_deg', ...
    'Pmax_MPa', ...
    'Hmin_um', ...
    'Fx_N', ...
    'Fy_N', ...
    'Load_N', ...
    'Friction_N', ...
    'MomentResidual'});

disp(resultTable);

%% 11. 绘制五张结果图

if makePlots
    plotFivePadSlipRegionMapNoSlip(pad, z, par);
    plotFivePadCircumferentialCurvesNoSlip(pad, z, par, 'pressure');
    plotFivePadCircumferentialCurvesNoSlip(pad, z, par, 'thickness');
    plotFivePadCircumferentialCurvesNoSlip(pad, z, par, 'friction');
    plotFivePadCircumferentialCurvesNoSlip(pad, z, par, 'load');
end

end


%% ========================================================================
%  求解单块瓦的平衡摆角
% =========================================================================
function pad = solvePadEquilibrium(iPad, z, par)

thetaPivot = par.thetaPivot(iPad);

theta = linspace( ...
    thetaPivot - par.beta/2, ...
    thetaPivot + par.beta/2, ...
    par.Ntheta)';

%% 计算安全摆角搜索范围

hBase = calculateFilmThickness(theta, 0, thetaPivot, par);

sinTermMax = max(abs(sin(theta - thetaPivot)));

if sinTermMax < 1e-14
    gammaLimit = par.gammaCap;
else
    gammaLimitByFilm = ...
        par.gammaSafety * min(hBase) / (par.Rb*sinTermMax);

    gammaLimit = min(par.gammaCap, gammaLimitByFilm);
end

if gammaLimit <= 0
    error('第%d块瓦的摆角搜索范围无效。', iPad);
end

%% 摆角粗扫描

gammaScan = linspace(-gammaLimit, gammaLimit, par.NgammaScan);

Mscan = zeros(size(gammaScan));
Fscan = zeros(size(gammaScan));

for k = 1:length(gammaScan)

    state = evaluatePadAtGamma( ...
        theta, ...
        z, ...
        thetaPivot, ...
        gammaScan(k), ...
        par, ...
        par.pressureTolScan);

    Mscan(k) = state.M;
    Fscan(k) = state.F;
end

%% 查找承载状态下的力矩变号区间

loadedLeft  = Fscan(1:end-1) > par.unloadForceTol;
loadedRight = Fscan(2:end)   > par.unloadForceTol;

momentChange = ...
    Mscan(1:end-1).*Mscan(2:end) <= 0;

rootIntervals = find(loadedLeft & loadedRight & momentChange);

if ~isempty(rootIntervals)

    % 如果存在多个根，选择距离零摆角最近的根
    intervalCenters = ...
        0.5*(gammaScan(rootIntervals) + ...
        gammaScan(rootIntervals+1));

    [~, localIndex] = min(abs(intervalCenters));
    iRoot = rootIntervals(localIndex);

    gammaLow  = gammaScan(iRoot);
    gammaHigh = gammaScan(iRoot+1);

    Mlow = Mscan(iRoot);

    for iterTilt = 1:par.maxTiltIter

        gammaMid = 0.5*(gammaLow + gammaHigh);

        stateMid = evaluatePadAtGamma( ...
            theta, ...
            z, ...
            thetaPivot, ...
            gammaMid, ...
            par, ...
            par.pressureTolRoot);

        if stateMid.F > par.unloadForceTol

            residual = abs(stateMid.M) / ...
                max(stateMid.F*par.Rb, realmin);

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

    %% 未找到承载力矩平衡根，判断是否为卸载瓦

    unloadedIndex = find(Fscan <= par.unloadForceTol);

    if ~isempty(unloadedIndex)

        % 自由卸载瓦的摆角在纯流体模型中并不唯一。
        % 本程序选择距离零摆角最近、且油膜压力消失的状态。
        [~, iNearest] = min(abs(gammaScan(unloadedIndex)));

        gammaFinal = gammaScan(unloadedIndex(iNearest));
        statusFinal = '卸载';

    else

        % 极少数情况下，在搜索范围内既无力矩根也未完全卸载。
        % 此时选择无量纲力矩残差最小的位置。
        normalizedMoment = abs(Mscan) ./ ...
            max(Fscan*par.Rb, realmin);

        [~, iBest] = min(normalizedMoment);

        gammaFinal = gammaScan(iBest);
        statusFinal = '搜索边界受限';
    end
end

%% 使用最终摆角重新高精度求解

stateFinal = evaluatePadAtGamma( ...
    theta, ...
    z, ...
    thetaPivot, ...
    gammaFinal, ...
    par, ...
    par.pressureTolFinal);

%% 保存结果

pad.index = iPad;
pad.thetaPivot = thetaPivot;
pad.theta = theta;
pad.gamma = gammaFinal;
pad.status = statusFinal;

pad.P = stateFinal.P;
pad.H = stateFinal.H;
pad.h = stateFinal.h;

pad.dpdx = stateFinal.dpdx;

pad.tauPad = stateFinal.tauPad;
pad.tauJournal = stateFinal.tauJournal;

pad.M = stateFinal.M;
pad.Fx = stateFinal.Fx;
pad.Fy = stateFinal.Fy;
pad.F  = stateFinal.F;

pad.pMax = max(stateFinal.P(:));
pad.hMin = min(stateFinal.H(:));

pad.frictionSigned = stateFinal.frictionSigned;
pad.frictionMagnitude = stateFinal.frictionMagnitude;

pad.pressureIterations = stateFinal.pressureIterations;
pad.pressureError = stateFinal.pressureError;

if stateFinal.F > par.unloadForceTol

    pad.momentResidual = ...
        abs(stateFinal.M) / ...
        max(stateFinal.F*par.Rb, realmin);

else

    pad.momentResidual = 0;
end

% 无滑移工况下所有网格点的滑移指示量均为0
pad.slipIndicator = zeros(size(stateFinal.P));

end


%% ========================================================================
%  在给定摆角下计算单块瓦
% =========================================================================
function state = evaluatePadAtGamma( ...
    theta, z, thetaPivot, gamma, par, pressureTolerance)

%% 油膜厚度

h = calculateFilmThickness(theta, gamma, thetaPivot, par);

if min(h) <= 0

    error(['出现非正油膜厚度：gamma = %.6e rad，' ...
        'hmin = %.6e m。'], gamma, min(h));
end

%% 油膜厚度沿周向弧长的导数

% h = C0 - e*cos(theta-thetaE)
%           + Rb*gamma*sin(theta-thetaPivot)
%
% dh/dtheta =
%       e*sin(theta-thetaE)
%       + Rb*gamma*cos(theta-thetaPivot)

dhdTheta = ...
    par.e*sin(theta-par.thetaE) + ...
    par.Rb*gamma*cos(theta-thetaPivot);

dhdx = dhdTheta / par.Rb;

%% 求解Reynolds方程

[P, pressureIterations, pressureError] = ...
    solvePressureNoSlip(theta, z, h, dhdx, par, pressureTolerance);

%% 计算载荷、力矩和摩擦

metrics = calculatePadMetrics( ...
    theta, z, thetaPivot, P, h, par);

%% 保存状态

state.P = P;
state.H = metrics.H;
state.h = h;

state.dpdx = metrics.dpdx;

state.tauPad = metrics.tauPad;
state.tauJournal = metrics.tauJournal;

state.M = metrics.M;
state.Fx = metrics.Fx;
state.Fy = metrics.Fy;
state.F  = metrics.F;

state.frictionSigned = metrics.frictionSigned;
state.frictionMagnitude = metrics.frictionMagnitude;

state.pressureIterations = pressureIterations;
state.pressureError = pressureError;

end


%% ========================================================================
%  油膜厚度
% =========================================================================
function h = calculateFilmThickness(theta, gamma, thetaPivot, par)

% 基础圆柱间隙加轴颈偏心，再叠加瓦块绕支点的小角度倾转
%
% phi=0时，偏心方向为thetaE=-90°，即朝向1号瓦中心
%
% gamma>0时：
%   theta>thetaPivot一侧的膜厚增加；
%   theta<thetaPivot一侧的膜厚减小。

h = ...
    par.C0 ...
    - par.e*cos(theta-par.thetaE) ...
    + par.Rb*gamma*sin(theta-thetaPivot);

end


%% ========================================================================
%  无滑移Reynolds方程：红黑投影SOR法
% =========================================================================
function [P, iter, errorP] = solvePressureNoSlip( ...
    theta, z, h, dhdx, par, tolerance)

Ntheta = length(theta);
Nz = length(z);

dtheta = theta(2)-theta(1);
dx = par.Rb*dtheta;
dz = z(2)-z(1);

H = repmat(h, 1, Nz);
H3 = H.^3;

%% 内部节点系数

H3E = 0.5*( ...
    H3(2:end-1, 2:end-1) + ...
    H3(3:end,   2:end-1));

H3W = 0.5*( ...
    H3(2:end-1, 2:end-1) + ...
    H3(1:end-2, 2:end-1));

H3N = 0.5*( ...
    H3(2:end-1, 2:end-1) + ...
    H3(2:end-1, 3:end));

H3S = 0.5*( ...
    H3(2:end-1, 2:end-1) + ...
    H3(2:end-1, 1:end-2));

AE = H3E / dx^2;
AW = H3W / dx^2;
AN = H3N / dz^2;
AS = H3S / dz^2;

AP = AE + AW + AN + AS;

%% Reynolds方程右端项

% d/dx(h^3 dp/dx) + d/dz(h^3 dp/dz)
%                       = 6*mu*U*dh/dx

source1D = 6*par.mu*par.U*dhdx(2:end-1);

RHS = repmat(source1D, 1, Nz-2);

%% 红黑网格

[I, J] = ndgrid(2:Ntheta-1, 2:Nz-1);

redMask = mod(I+J, 2) == 0;
blackMask = ~redMask;

%% 压力初始化

P = zeros(Ntheta, Nz);

errorP = inf;

for iter = 1:par.maxPressureIter

    if mod(iter-1, par.checkEvery) == 0
        Pold = P;
    end

    %% 红色节点更新

    Pgs = ( ...
        AE.*P(3:end,   2:end-1) + ...
        AW.*P(1:end-2, 2:end-1) + ...
        AN.*P(2:end-1, 3:end)   + ...
        AS.*P(2:end-1, 1:end-2) - ...
        RHS) ./ AP;

    Pinterior = P(2:end-1, 2:end-1);

    Pcandidate = ...
        (1-par.sorOmega).*Pinterior + ...
        par.sorOmega.*Pgs;

    % 负压截断
    Pcandidate = max(Pcandidate, 0);

    Pinterior(redMask) = Pcandidate(redMask);

    P(2:end-1, 2:end-1) = Pinterior;

    %% 黑色节点更新

    Pgs = ( ...
        AE.*P(3:end,   2:end-1) + ...
        AW.*P(1:end-2, 2:end-1) + ...
        AN.*P(2:end-1, 3:end)   + ...
        AS.*P(2:end-1, 1:end-2) - ...
        RHS) ./ AP;

    Pinterior = P(2:end-1, 2:end-1);

    Pcandidate = ...
        (1-par.sorOmega).*Pinterior + ...
        par.sorOmega.*Pgs;

    % 负压截断
    Pcandidate = max(Pcandidate, 0);

    Pinterior(blackMask) = Pcandidate(blackMask);

    P(2:end-1, 2:end-1) = Pinterior;

    %% 明确设置全部压力边界为0

    P(1, :)   = 0;
    P(end, :) = 0;
    P(:, 1)   = 0;
    P(:, end) = 0;

    %% 收敛判断

    if mod(iter, par.checkEvery) == 0

        numerator = max(abs(P(:)-Pold(:)));
        denominator = max(P(:)) + 1;

        errorP = numerator / denominator;

        if errorP < tolerance
            break;
        end
    end
end

end


%% ========================================================================
%  单块瓦载荷、力矩和摩擦计算
% =========================================================================
function metrics = calculatePadMetrics( ...
    theta, z, thetaPivot, P, h, par)

Nz = length(z);

Theta = repmat(theta, 1, Nz);
H = repmat(h, 1, Nz);

dtheta = theta(2)-theta(1);
dx = par.Rb*dtheta;

%% 压力沿周向弧长的梯度

dpdx = zeros(size(P));

dpdx(2:end-1, :) = ...
    (P(3:end, :) - P(1:end-2, :)) / (2*dx);

dpdx(1, :) = ...
    (P(2, :) - P(1, :)) / dx;

dpdx(end, :) = ...
    (P(end, :) - P(end-1, :)) / dx;

%% 油膜力

% 压力对轴颈的作用方向沿径向向内，因此有负号

integrandFx = -P.*cos(Theta);
integrandFy = -P.*sin(Theta);

Fx = par.Rb * ...
    trapz(z, trapz(theta, integrandFx, 1), 2);

Fy = par.Rb * ...
    trapz(z, trapz(theta, integrandFy, 1), 2);

F = sqrt(Fx^2 + Fy^2);

%% 瓦块支点力矩

% 压力作用点相对支点的周向力臂：
%
% l = Rb*sin(theta-thetaPivot)
%
% 力矩为0表示压力中心通过支点。

momentArm = par.Rb*sin(Theta-thetaPivot);

M = par.Rb * ...
    trapz(z, trapz(theta, P.*momentArm, 1), 2);

%% 无滑移条件下的壁面剪应力

% 速度分布：
%
% u(y) = [1/(2mu)]*dp/dx*(y^2-hy) + U*y/h
%
% 固定瓦面y=0处：
%
% tau_pad = mu*U/h - h/2*dp/dx
%
% 运动轴颈表面y=h处：
%
% tau_journal = mu*U/h + h/2*dp/dx

tauPad = ...
    par.mu*par.U./H ...
    - 0.5*H.*dpdx;

tauJournal = ...
    par.mu*par.U./H ...
    + 0.5*H.*dpdx;

%% 摩擦力

% 有符号摩擦力反映整体周向合力
frictionSigned = par.Rb * ...
    trapz(z, trapz(theta, tauPad, 1), 2);

% 局部幅值积分便于评价总摩擦强度
frictionMagnitude = par.Rb * ...
    trapz(z, trapz(theta, abs(tauPad), 1), 2);

%% 输出

metrics.H = H;

metrics.dpdx = dpdx;

metrics.tauPad = tauPad;
metrics.tauJournal = tauJournal;

metrics.Fx = Fx;
metrics.Fy = Fy;
metrics.F  = F;

metrics.M = M;

metrics.frictionSigned = frictionSigned;
metrics.frictionMagnitude = frictionMagnitude;

end


%% ========================================================================
%  图1：径向轴承周向展开滑移区域图
% =========================================================================
function plotFivePadSlipRegionMapNoSlip(pad, z, par)

% 横坐标从1号瓦入口开始，五块径向瓦按旋转方向依次展开。
% 纵坐标为轴向无量纲坐标z/B。无滑移=蓝色，滑移=黄色。
thetaZero = par.thetaPivot(1)-par.beta/2;

fig = figure('Color','w','Position',[120,80,1050,560], ...
    'Name','图1 五瓦滑移区域');
ax = axes('Parent',fig,'Position',[0.095,0.15,0.78,0.77]);
hold(ax,'on');

for iPad = 1:par.Npad
    thetaDeg = mod((pad(iPad).theta(:)-thetaZero)*180/pi,360);
    zStar = z/par.B;
    slipField = double(pad(iPad).slipIndicator).';

    imagesc(ax,thetaDeg,zStar,slipField);

    thetaCenter = mean(thetaDeg);
    text(ax,thetaCenter,0,sprintf('%d号瓦',iPad), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontName','Microsoft YaHei', ...
        'FontSize',12, ...
        'FontWeight','bold', ...
        'Color','w', ...
        'BackgroundColor',[0.12 0.12 0.12], ...
        'Margin',2);
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
title(ax,'五瓦可倾瓦径向轴承滑移区域', ...
    'FontName','Microsoft YaHei','FontSize',15,'FontWeight','bold');

cb = colorbar(ax,'Location','eastoutside');
cb.Ticks = [0 1];
cb.TickLabels = {'无滑移','滑移'};
cb.FontName = 'Microsoft YaHei';
cb.FontSize = 11;

annotation(fig,'textbox',[0.27,0.015,0.46,0.055], ...
    'String','(a) 五瓦滑移区域（无滑移工况）', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FontName','Microsoft YaHei','FontSize',14);

saveNoSlipFigure(fig,'TPJB_NoSlip_Fig1_SlipRegion.png');
end


%% ========================================================================
%  图2至图5：五块径向瓦的周向二维分布曲线
% =========================================================================
function plotFivePadCircumferentialCurvesNoSlip(pad, z, par, fieldType)

thetaZero = par.thetaPivot(1)-par.beta/2;
thetaCurve = cell(1,par.Npad);
valueCurve = cell(1,par.Npad);
allValues = [];

for iPad = 1:par.Npad
    thetaCurve{iPad} = mod( ...
        (pad(iPad).theta(:)-thetaZero)*180/pi,360);

    switch fieldType
        case 'pressure'
            % 径向轴承常用轴向中截面压力曲线。
            midZ = round((numel(z)+1)/2);
            valueCurve{iPad} = pad(iPad).P(:,midZ)/1e6;

        case 'thickness'
            % 当前刚性、无不对中模型的膜厚沿轴向不变。
            valueCurve{iPad} = pad(iPad).h(:)*1e6;

        case 'friction'
            % dF_f/dtheta = Rb*integral(|tau_pad| dz)，单位N/rad。
            valueCurve{iPad} = par.Rb* ...
                trapz(z,abs(pad(iPad).tauPad),2);

        case 'load'
            % dW/dtheta = Rb*integral(p dz)，单位N/rad。
            % 这是局部法向承载力分布；积分后得到单瓦压力合力尺度。
            valueCurve{iPad} = par.Rb* ...
                trapz(z,pad(iPad).P,2);

        otherwise
            error('未知二维绘图变量类型：%s',fieldType);
    end

    allValues = [allValues; valueCurve{iPad}(:)]; %#ok<AGROW>
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
    labelY = fieldValue(labelIndex) + (0.035+0.007*mod(iPad,2))*valueRange;
    text(ax,labelX,labelY,sprintf('%d号瓦',iPad), ...
        'Color',curveColors(iPad,:), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName','Microsoft YaHei', ...
        'FontSize',10.5, ...
        'FontWeight','bold', ...
        'BackgroundColor','w', ...
        'Margin',1, ...
        'Clipping','off');
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
        titleText = '无滑移五瓦轴向中截面油膜压力分布';
        caption = '(b) 五瓦油膜压力分布';
        fileName = 'TPJB_NoSlip_Fig2_Pressure.png';
        ylim(ax,[0 valueMax+0.14*valueRange]);

    case 'thickness'
        ylabelText = '油膜厚度  h / \mum';
        titleText = '无滑移五瓦油膜厚度分布';
        caption = '(c) 五瓦油膜厚度分布';
        fileName = 'TPJB_NoSlip_Fig3_FilmThickness.png';
        ylim(ax,[max(0,valueMin-0.10*valueRange), ...
            valueMax+0.16*valueRange]);

    case 'friction'
        ylabelText = '局部摩擦力幅值  dF_f/d\theta / (N·rad^{-1})';
        titleText = '无滑移五瓦周向摩擦力分布';
        caption = '(d) 五瓦摩擦力分布';
        fileName = 'TPJB_NoSlip_Fig4_Friction.png';
        ylim(ax,[0 valueMax+0.14*valueRange]);

    case 'load'
        ylabelText = '局部法向承载力  dW/d\theta / (N·rad^{-1})';
        titleText = '无滑移五瓦周向承载力分布';
        caption = '(e) 五瓦承载力分布';
        fileName = 'TPJB_NoSlip_Fig5_Load.png';
        ylim(ax,[0 valueMax+0.14*valueRange]);
end

ylabel(ax,ylabelText,'FontName','Microsoft YaHei','FontSize',14);
title(ax,titleText,'FontName','Microsoft YaHei', ...
    'FontSize',15,'FontWeight','bold');

legend(ax,curveHandles,'Location','northeastoutside', ...
    'FontName','Microsoft YaHei','FontSize',10.5,'Box','on');

annotation(fig,'textbox',[0.27,0.015,0.46,0.055], ...
    'String',caption, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FontName','Microsoft YaHei','FontSize',14);

saveNoSlipFigure(fig,fileName);
end


%% ========================================================================
%  将结果图保存到当前m文件所在的outputs目录
% =========================================================================
function saveNoSlipFigure(fig,fileName)

outputDirectory = fileparts(mfilename('fullpath'));
outputFigure = fullfile(outputDirectory,fileName);
exportgraphics(fig,outputFigure,'Resolution',300);
fprintf('结果图已保存至：%s\n',outputFigure);
end
