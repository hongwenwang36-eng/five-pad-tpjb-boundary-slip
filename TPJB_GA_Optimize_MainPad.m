function result = TPJB_GA_Optimize_MainPad(replotOnly)
% =========================================================================
% 五瓦可倾瓦径向轴承主承载瓦边界滑移区域遗传算法优化
%   1. 其他四块瓦保持全表面具有滑移能力；
%   2. 自动识别主承载瓦，并对其允许滑移区域进行0/1编码；
%   3. GA阶段固定全滑移基准平衡摆角，采用粗网格快速评价；
%   4. 最优区域映射回81x31细网格后，重新求解五瓦摆角平衡；
%   5. 优化目标为降低全轴承摩擦因数，同时限制承载力损失并抑制
%      过度零散的滑移区域。
% =========================================================================

if nargin < 1
    replotOnly = false;
end

clc;
close all;
rng(20260806,'twister');

outputDir = fullfile(fileparts(mfilename('fullpath')), ...
    'GA_Optimization_Results');

% Recreate publication figures from the saved verified result without
% rerunning the expensive genetic algorithm.
if replotOnly
    loaded = load(fullfile(outputDir,'GA_Optimization_Result.mat'),'result');
    result = loaded.result;
    parFine.Npad = numel(result.basePad);
    parFine.Ntheta = size(result.basePad(1).P,1);
    parFine.Nz = size(result.basePad(1).P,2);
    parFine.Rb = 15.0325e-3;
    parFine.B = 15e-3;
    parFine.beta = 55*pi/180;
    parFine.thetaPivot = [result.basePad.thetaPivot];
    parOpt = parFine;
    zFine = linspace(-parFine.B/2,parFine.B/2,parFine.Nz);
    zFineOpt = zFine;
    plotOptimizationResults(result,zFine,zFineOpt, ...
        parFine,parOpt,outputDir);
    fprintf('Figures regenerated in: %s\n',outputDir);
    return;
end

fprintf('============================================================\n');
fprintf(' 五瓦可倾瓦径向轴承主承载瓦滑移区域GA优化\n');
fprintf('============================================================\n\n');

%% 1. 全滑移基准工况
fprintf('步骤1/4：求解五瓦全滑移基准工况...\n');
[basePad,zFine,parFine] = TPJB_BinarySlip_5Pad(false);
baseMetrics = systemMetrics(basePad,parFine);

[~,mainPad] = max([basePad.F]);
fprintf('\n自动识别的主承载瓦：%d号瓦，承载力 %.6f N\n', ...
    mainPad,basePad(mainPad).F);

%% 2. 遗传算法参数
ga.designNtheta = 8;
ga.designNz = 6;
ga.nGenes = ga.designNtheta*ga.designNz;
ga.coverage = 0.50;
ga.nOnes = round(ga.coverage*ga.nGenes);
ga.populationSize = 20;
ga.maxGenerations = 25;
ga.eliteCount = 2;
ga.tournamentSize = 3;
ga.crossoverProbability = 0.85;
ga.mutationProbability = 0.05;
ga.maxStallGenerations = 10;
ga.capacityFloor = 0.98;
ga.capacityPenaltyWeight = 8.0;
ga.fragmentPenaltyWeight = 0.015;
ga.balanceMomentTol = 2e-5;
ga.balanceMaxIter = 8;
ga.balanceStep = 2.5e-5;

% GA粗网格与收敛参数
coarse = parFine;
coarse.Ntheta = 41;
coarse.Nz = 17;
coarse.sorOmega = 1.55;
coarse.pressureTol = 8e-5;
coarse.slipTol = 8e-5;
coarse.maxPressureIter = 1800;
coarse.maxSlipIter = 35;
coarse.checkEvery = 5;
coarse.maskStableRequired = 2;

thetaPivot = parFine.thetaPivot(mainPad);
thetaCoarse = linspace(thetaPivot-parFine.beta/2, ...
    thetaPivot+parFine.beta/2,coarse.Ntheta)';
zCoarse = linspace(-parFine.B/2,parFine.B/2,coarse.Nz);

ctx.mainPad = mainPad;
ctx.theta = thetaCoarse;
ctx.z = zCoarse;
ctx.thetaPivot = thetaPivot;
ctx.gammaStart = basePad(mainPad).gamma;
ctx.par = coarse;
ctx.otherFx = sum([basePad(setdiff(1:parFine.Npad,mainPad)).Fx]);
ctx.otherFy = sum([basePad(setdiff(1:parFine.Npad,mainPad)).Fy]);
ctx.otherFriction = sum( ...
    [basePad(setdiff(1:parFine.Npad,mainPad)).frictionMagnitude]);
ctx.ga = ga;

% The GA reference must use exactly the same coarse grid and the same
% pad-moment equilibrium process as every chromosome.  Mixing a fine-grid
% reference with a coarse, fixed-tilt candidate creates a false improvement.
coarseFullSlipMask = true(coarse.Ntheta,coarse.Nz);
coarseReference = solveBalancedPadGA(thetaCoarse,zCoarse,thetaPivot, ...
    ctx.gammaStart,coarse,coarseFullSlipMask,ga);
refFx = ctx.otherFx+coarseReference.Fx;
refFy = ctx.otherFy+coarseReference.Fy;
ctx.refW = hypot(refFx,refFy);
ctx.refFriction = ctx.otherFriction+coarseReference.frictionMagnitude;
ctx.refMu = ctx.refFriction/ctx.refW;
ctx.coarseReference = coarseReference;

fprintf(['GA coarse reference: gamma = %.7f deg, W = %.6f N, ' ...
    'mu = %.9f\n'],coarseReference.gamma*180/pi,ctx.refW,ctx.refMu);

fprintf('\n步骤2/4：开始GA搜索（%d个变量，固定%d个滑移单元）...\n', ...
    ga.nGenes,ga.nOnes);

[bestChromosome,history,gaSummary] = runBinaryGA(ctx);
bestDesignMask = reshape(bestChromosome, ...
    ga.designNtheta,ga.designNz);

%% 3. 细网格五瓦重新平衡验证
fprintf('\n步骤3/4：将最优区域映射到细网格并重新求解五瓦平衡...\n');
fineMainMask = mapDesignMask(bestDesignMask, ...
    parFine.Ntheta,parFine.Nz);

allowedMasks = repmat({true(parFine.Ntheta,parFine.Nz)}, ...
    1,parFine.Npad);
allowedMasks{mainPad} = fineMainMask;

[optPad,zFineOpt,parOpt] = TPJB_BinarySlip_5Pad(false,allowedMasks);
optMetrics = systemMetrics(optPad,parOpt);

%% 4. 输出、保存和绘图
fprintf('\n步骤4/4：整理优化结果并输出图形...\n');

result.mainPad = mainPad;
result.ga = ga;
result.gaSummary = gaSummary;
result.bestChromosome = bestChromosome;
result.bestDesignMask = bestDesignMask;
result.fineAllowedMask = fineMainMask;
result.basePad = basePad;
result.optimizedPad = optPad;
result.baseMetrics = baseMetrics;
result.optimizedMetrics = optMetrics;
result.history = history;

result.coveragePercent = 100*nnz(fineMainMask)/numel(fineMainMask);
result.actualSlipPercentMain = 100*optPad(mainPad).slipAreaRatio;
result.loadChangePercent = 100*(optMetrics.W/baseMetrics.W-1);
result.frictionChangePercent = 100*(optMetrics.friction/ ...
    baseMetrics.friction-1);
result.frictionCoefficientChangePercent = 100*( ...
    optMetrics.frictionCoefficient/baseMetrics.frictionCoefficient-1);

fprintf('\n============================================================\n');
fprintf(' GA优化最终结果\n');
fprintf('============================================================\n');
fprintf('主承载瓦                         = %d号瓦\n',mainPad);
fprintf('GA实际完成代数                   = %d\n',gaSummary.generations);
fprintf('GA唯一方案评价次数               = %d\n',gaSummary.uniqueEvaluations);
fprintf('主承载瓦允许滑移面积比例         = %.6f %%\n', ...
    result.coveragePercent);
fprintf('主承载瓦最终实际滑移面积比例     = %.6f %%\n', ...
    result.actualSlipPercentMain);
fprintf('全滑移基准总承载力               = %.8f N\n',baseMetrics.W);
fprintf('优化后总承载力                   = %.8f N\n',optMetrics.W);
fprintf('承载力变化                       = %+.6f %%\n', ...
    result.loadChangePercent);
fprintf('全滑移基准摩擦力                 = %.8f N\n', ...
    baseMetrics.friction);
fprintf('优化后摩擦力                     = %.8f N\n', ...
    optMetrics.friction);
fprintf('摩擦力变化                       = %+.6f %%\n', ...
    result.frictionChangePercent);
fprintf('全滑移基准摩擦因数               = %.10f\n', ...
    baseMetrics.frictionCoefficient);
fprintf('优化后摩擦因数                   = %.10f\n', ...
    optMetrics.frictionCoefficient);
fprintf('摩擦因数变化                     = %+.6f %%\n', ...
    result.frictionCoefficientChangePercent);
fprintf('全滑移基准最大压力               = %.8f MPa\n', ...
    baseMetrics.pMax/1e6);
fprintf('优化后最大压力                   = %.8f MPa\n', ...
    optMetrics.pMax/1e6);
fprintf('全滑移基准最小膜厚               = %.8f um\n', ...
    baseMetrics.hMin*1e6);
fprintf('优化后最小膜厚                   = %.8f um\n', ...
    optMetrics.hMin*1e6);
fprintf('============================================================\n');

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

save(fullfile(outputDir,'GA_Optimization_Result.mat'),'result','-v7.3');
writeResultTable(result,outputDir);
plotOptimizationResults(result,zFine,zFineOpt,parFine,parOpt,outputDir);

fprintf('\n结果文件夹：%s\n',outputDir);
fprintf('运行完成。\n');

end


%% ========================================================================
%  自编二进制遗传算法
% =========================================================================
function [bestChromosome,history,summary] = runBinaryGA(ctx)

ga = ctx.ga;
population = initializePopulation(ga);
cache = containers.Map('KeyType','char','ValueType','any');

history.bestObjective = nan(ga.maxGenerations,1);
history.meanObjective = nan(ga.maxGenerations,1);
history.bestLoadRatio = nan(ga.maxGenerations,1);
history.bestFrictionRatio = nan(ga.maxGenerations,1);
history.bestMuRatio = nan(ga.maxGenerations,1);
history.uniqueEvaluations = nan(ga.maxGenerations,1);

globalBestObjective = inf;
bestChromosome = population(1,:);
stallCount = 0;

for generation = 1:ga.maxGenerations
    [objective,stats,cache] = evaluatePopulation(population,ctx,cache);
    [objective,order] = sort(objective,'ascend');
    population = population(order,:);
    stats = stats(order);

    generationBest = objective(1);
    if generationBest < globalBestObjective-1e-9
        globalBestObjective = generationBest;
        bestChromosome = population(1,:);
        globalBestStats = stats(1); %#ok<NASGU>
        stallCount = 0;
    else
        stallCount = stallCount+1;
    end

    history.bestObjective(generation) = generationBest;
    history.meanObjective(generation) = mean(objective);
    history.bestLoadRatio(generation) = stats(1).W/ctx.refW;
    history.bestFrictionRatio(generation) = ...
        stats(1).friction/ctx.refFriction;
    history.bestMuRatio(generation) = stats(1).mu/ctx.refMu;
    history.uniqueEvaluations(generation) = cache.Count;

    fprintf(['第%2d代：J=%.6f，W/W0=%.5f，Ff/Ff0=%.5f，' ...
        'mu/mu0=%.5f，唯一评价=%d\n'],generation,generationBest, ...
        history.bestLoadRatio(generation), ...
        history.bestFrictionRatio(generation), ...
        history.bestMuRatio(generation),cache.Count);

    if stallCount >= ga.maxStallGenerations
        fprintf('最佳目标连续%d代未改善，提前停止。\n',stallCount);
        break;
    end

    newPopulation = false(size(population));
    newPopulation(1:ga.eliteCount,:) = ...
        population(1:ga.eliteCount,:);

    iChild = ga.eliteCount+1;
    while iChild <= ga.populationSize
        parent1 = tournamentSelect(population,objective, ...
            ga.tournamentSize);
        parent2 = tournamentSelect(population,objective, ...
            ga.tournamentSize);

        child1 = parent1;
        child2 = parent2;
        if rand < ga.crossoverProbability
            crossMask = rand(1,ga.nGenes) < 0.5;
            child1(crossMask) = parent2(crossMask);
            child2(crossMask) = parent1(crossMask);
        end

        child1 = mutateAndRepair(child1,ga);
        child2 = mutateAndRepair(child2,ga);

        newPopulation(iChild,:) = child1;
        if iChild+1 <= ga.populationSize
            newPopulation(iChild+1,:) = child2;
        end
        iChild = iChild+2;
    end
    population = newPopulation;
end

lastGeneration = generation;
fields = fieldnames(history);
for k = 1:numel(fields)
    history.(fields{k}) = history.(fields{k})(1:lastGeneration,:);
end

summary.generations = lastGeneration;
summary.uniqueEvaluations = cache.Count;
summary.bestObjective = globalBestObjective;

end


%% ========================================================================
%  初始种群：含入口、出口、中心、轴向半区和随机设计
% =========================================================================
function population = initializePopulation(ga)

population = false(ga.populationSize,ga.nGenes);
seedMasks = cell(1,4);

mask = false(ga.designNtheta,ga.designNz);
mask(1:round(ga.designNtheta/2),:) = true;
seedMasks{1} = mask;

mask = false(ga.designNtheta,ga.designNz);
mask(end-round(ga.designNtheta/2)+1:end,:) = true;
seedMasks{2} = mask;

mask = false(ga.designNtheta,ga.designNz);
startRow = floor((ga.designNtheta-round(ga.designNtheta/2))/2)+1;
mask(startRow:startRow+round(ga.designNtheta/2)-1,:) = true;
seedMasks{3} = mask;

mask = false(ga.designNtheta,ga.designNz);
mask(:,1:round(ga.designNz/2)) = true;
seedMasks{4} = mask;

for i = 1:min(numel(seedMasks),ga.populationSize)
    population(i,:) = repairCoverage(seedMasks{i}(:)',ga.nOnes);
end

for i = numel(seedMasks)+1:ga.populationSize
    chromosome = false(1,ga.nGenes);
    chromosome(randperm(ga.nGenes,ga.nOnes)) = true;
    population(i,:) = chromosome;
end

end


%% ========================================================================
%  评价种群，缓存重复染色体
% =========================================================================
function [objective,stats,cache] = evaluatePopulation(population,ctx,cache)

n = size(population,1);
objective = zeros(n,1);
stats = repmat(struct('W',0,'friction',0,'mu',0,'pMax',0, ...
    'actualSlipRatio',0,'gamma',0,'momentResidual',inf),n,1);

for i = 1:n
    key = char(double(population(i,:))+'0');
    if isKey(cache,key)
        value = cache(key);
    else
        [obj,st] = evaluateChromosome(population(i,:),ctx);
        value.objective = obj;
        value.stats = st;
        cache(key) = value;
    end
    objective(i) = value.objective;
    stats(i) = value.stats;
end

end


%% ========================================================================
%  单个染色体粗网格评价
% =========================================================================
function [objective,stats] = evaluateChromosome(chromosome,ctx)

designMask = reshape(logical(chromosome), ...
    ctx.ga.designNtheta,ctx.ga.designNz);
allowedMask = mapDesignMask(designMask, ...
    ctx.par.Ntheta,ctx.par.Nz);

state = solveBalancedPadGA(ctx.theta,ctx.z,ctx.thetaPivot, ...
    ctx.gammaStart,ctx.par,allowedMask,ctx.ga);

FxTotal = ctx.otherFx+state.Fx;
FyTotal = ctx.otherFy+state.Fy;
stats.W = hypot(FxTotal,FyTotal);
stats.friction = ctx.otherFriction+state.frictionMagnitude;
stats.mu = stats.friction/max(stats.W,realmin);
stats.pMax = max(state.P,[],'all');
stats.actualSlipRatio = nnz(state.slipMask)/numel(state.slipMask);
stats.gamma = state.gamma;
stats.momentResidual = state.momentResidual;

muRatio = stats.mu/ctx.refMu;
loadRatio = stats.W/ctx.refW;
capacityPenalty = ctx.ga.capacityPenaltyWeight* ...
    max(0,ctx.ga.capacityFloor-loadRatio)^2;
fragmentPenalty = ctx.ga.fragmentPenaltyWeight* ...
    normalizedFragmentation(designMask);

objective = muRatio+capacityPenalty+fragmentPenalty;

if ~isfinite(objective) || stats.W <= 0 || ...
        stats.momentResidual > 10*ctx.ga.balanceMomentTol
    objective = 1e6;
end

end


%% ========================================================================
%  固定摆角下压力-滑移耦合求解（GA粗网格）
% =========================================================================
function state = evaluateFixedPad(theta,z,thetaPivot,gamma,par,allowedMask)

h = par.C0-par.e*cos(theta-par.thetaE)+ ...
    par.Rb*gamma*sin(theta-thetaPivot);
if min(h) <= 0
    error('GA评价出现非正油膜厚度。');
end

dhdTheta = par.e*sin(theta-par.thetaE)+ ...
    par.Rb*gamma*cos(theta-thetaPivot);
dhdx = dhdTheta/par.Rb;

H = repmat(h,1,numel(z));
P = zeros(size(H));
slipU = zeros(size(H));
slipW = zeros(size(H));
oldMask = false(size(H));
stableCount = 0;

for iterSlip = 1:par.maxSlipIter
    Pold = P;
    slipUold = slipU;
    slipWold = slipW;

    [P,~,~] = solvePressureGA(theta,z,h,dhdx, ...
        slipU,slipW,P,par);
    [dpdx,dpdz] = pressureGradientsGA(P,theta,z,par);

    tauTrialX = par.mu*par.U./H-0.5*H.*dpdx;
    tauTrialZ = -0.5*H.*dpdz;
    tauTrialMagnitude = hypot(tauTrialX,tauTrialZ);
    activeMask = logical(allowedMask) & ...
        (tauTrialMagnitude > par.tauC);

    activation = zeros(size(H));
    activation(activeMask) = 1-par.tauC./ ...
        tauTrialMagnitude(activeMask);
    coefficient = par.b*H./(par.mu*(H+par.b));

    slipUcalc = coefficient.*activation.*tauTrialX;
    slipWcalc = coefficient.*activation.*tauTrialZ;
    slipU = (1-par.slipRelax)*slipUold+par.slipRelax*slipUcalc;
    slipW = (1-par.slipRelax)*slipWold+par.slipRelax*slipWcalc;
    slipU(~activeMask) = 0;
    slipW(~activeMask) = 0;

    slipError = max(hypot(slipU-slipUold,slipW-slipWold), ...
        [],'all')/max(par.U,realmin);
    pressureError = max(abs(P-Pold),[],'all')/(max(P,[],'all')+1);

    if isequal(activeMask,oldMask)
        stableCount = stableCount+1;
    else
        stableCount = 0;
    end
    oldMask = activeMask;

    if slipError < par.slipTol && ...
            pressureError < 5*par.pressureTol && ...
            stableCount >= par.maskStableRequired
        break;
    end
end

[P,~,~] = solvePressureGA(theta,z,h,dhdx, ...
    slipU,slipW,P,par);
[dpdx,dpdz] = pressureGradientsGA(P,theta,z,par);

Theta = repmat(theta,1,numel(z));
tauPadX = par.mu*(par.U-slipU)./H-0.5*H.*dpdx;
tauPadZ = -par.mu*slipW./H-0.5*H.*dpdz;
tauMagnitude = hypot(tauPadX,tauPadZ);

state.P = P;
state.Fx = par.Rb*trapz(z,trapz(theta,-P.*cos(Theta),1),2);
state.Fy = par.Rb*trapz(z,trapz(theta,-P.*sin(Theta),1),2);
state.F = hypot(state.Fx,state.Fy);
momentArm = par.Rb*sin(Theta-thetaPivot);
state.M = par.Rb*trapz(z,trapz(theta,P.*momentArm,1),2);
state.frictionMagnitude = par.Rb* ...
    trapz(z,trapz(theta,tauMagnitude,1),2);
state.slipMask = hypot(slipU,slipW) > max(1e-12,1e-10*par.U);
state.hMin = min(H,[],'all');

end


%% ========================================================================
%  Rebalance the main pad for every GA chromosome (coarse grid)
% =========================================================================
function state = solveBalancedPadGA(theta,z,thetaPivot,gammaStart,par, ...
    allowedMask,ga)

hBase = par.C0-par.e*cos(theta-par.thetaE);
sinMax = max(abs(sin(theta-thetaPivot)));
gammaLimit = min(par.gammaCap, ...
    par.gammaSafety*min(hBase)/(par.Rb*max(sinMax,realmin)));

gamma0 = min(max(gammaStart,-gammaLimit),gammaLimit);
state0 = evaluateFixedPad(theta,z,thetaPivot,gamma0,par,allowedMask);
r0 = normalizedMomentGA(state0,par);

if abs(r0) <= ga.balanceMomentTol
    state = state0;
    state.gamma = gamma0;
    state.momentResidual = abs(r0);
    return;
end

% A safeguarded secant iteration is much cheaper than a full tilt scan for
% every chromosome.  The step direction is arbitrary; the second point and
% subsequent secant updates determine the local moment slope.
gamma1 = min(max(gamma0+ga.balanceStep,-gammaLimit),gammaLimit);
if gamma1 == gamma0
    gamma1 = min(max(gamma0-ga.balanceStep,-gammaLimit),gammaLimit);
end
state1 = evaluateFixedPad(theta,z,thetaPivot,gamma1,par,allowedMask);
r1 = normalizedMomentGA(state1,par);

bestState = state0;
bestGamma = gamma0;
bestResidual = abs(r0);
if abs(r1) < bestResidual
    bestState = state1;
    bestGamma = gamma1;
    bestResidual = abs(r1);
end

for k = 1:ga.balanceMaxIter
    if abs(r1) <= ga.balanceMomentTol
        bestState = state1;
        bestGamma = gamma1;
        bestResidual = abs(r1);
        break;
    end

    denominator = r1-r0;
    if abs(denominator) < 1e-12
        gamma2 = 0.5*(gamma0+gamma1);
    else
        gamma2 = gamma1-r1*(gamma1-gamma0)/denominator;
    end

    % Keep the secant step physical and prevent an over-aggressive jump.
    maxJump = 4*ga.balanceStep;
    gamma2 = min(max(gamma2,gamma1-maxJump),gamma1+maxJump);
    gamma2 = min(max(gamma2,-gammaLimit),gammaLimit);
    if abs(gamma2-gamma1) < 1e-11
        gamma2 = min(max(gamma1-sign(r1)*ga.balanceStep, ...
            -gammaLimit),gammaLimit);
    end

    state2 = evaluateFixedPad(theta,z,thetaPivot,gamma2,par,allowedMask);
    r2 = normalizedMomentGA(state2,par);
    if abs(r2) < bestResidual
        bestState = state2;
        bestGamma = gamma2;
        bestResidual = abs(r2);
    end

    gamma0 = gamma1;
    r0 = r1;
    gamma1 = gamma2;
    r1 = r2;
    state1 = state2;
end

state = bestState;
state.gamma = bestGamma;
state.momentResidual = bestResidual;

end


function residual = normalizedMomentGA(state,par)

residual = state.M/max(state.F*par.Rb,realmin);

end


%% ========================================================================
%  GA粗网格扩展Reynolds方程
% =========================================================================
function [P,iter,errorP] = solvePressureGA( ...
    theta,z,h,dhdx,slipU,slipW,Pinit,par)

Ntheta = numel(theta);
Nz = numel(z);
dx = par.Rb*(theta(2)-theta(1));
dz = z(2)-z(1);
H = repmat(h,1,Nz);
H3 = H.^3;

AE = 0.5*(H3(2:end-1,2:end-1)+H3(3:end,2:end-1))/dx^2;
AW = 0.5*(H3(2:end-1,2:end-1)+H3(1:end-2,2:end-1))/dx^2;
AN = 0.5*(H3(2:end-1,2:end-1)+H3(2:end-1,3:end))/dz^2;
AS = 0.5*(H3(2:end-1,2:end-1)+H3(2:end-1,1:end-2))/dz^2;
AP = AE+AW+AN+AS;

baseSource = repmat(6*par.mu*par.U*dhdx(2:end-1),1,Nz-2);
HU = H.*slipU;
HW = H.*slipW;
slipDivergence = ...
    (HU(3:end,2:end-1)-HU(1:end-2,2:end-1))/(2*dx)+ ...
    (HW(2:end-1,3:end)-HW(2:end-1,1:end-2))/(2*dz);
RHS = baseSource+6*par.mu*slipDivergence;

[I,J] = ndgrid(2:Ntheta-1,2:Nz-1);
redMask = mod(I+J,2)==0;
blackMask = ~redMask;

P = max(Pinit,0);
P([1 end],:) = 0;
P(:,[1 end]) = 0;
errorP = inf;

for iter = 1:par.maxPressureIter
    if mod(iter-1,par.checkEvery)==0
        Pold = P;
    end

    Pgs = (AE.*P(3:end,2:end-1)+AW.*P(1:end-2,2:end-1)+ ...
        AN.*P(2:end-1,3:end)+AS.*P(2:end-1,1:end-2)-RHS)./AP;
    Pint = P(2:end-1,2:end-1);
    Pcand = max((1-par.sorOmega)*Pint+par.sorOmega*Pgs,0);
    Pint(redMask) = Pcand(redMask);
    P(2:end-1,2:end-1) = Pint;

    Pgs = (AE.*P(3:end,2:end-1)+AW.*P(1:end-2,2:end-1)+ ...
        AN.*P(2:end-1,3:end)+AS.*P(2:end-1,1:end-2)-RHS)./AP;
    Pint = P(2:end-1,2:end-1);
    Pcand = max((1-par.sorOmega)*Pint+par.sorOmega*Pgs,0);
    Pint(blackMask) = Pcand(blackMask);
    P(2:end-1,2:end-1) = Pint;

    P([1 end],:) = 0;
    P(:,[1 end]) = 0;

    if mod(iter,par.checkEvery)==0
        errorP = max(abs(P-Pold),[],'all')/(max(P,[],'all')+1);
        if errorP < par.pressureTol
            break;
        end
    end
end

end


%% ========================================================================
%  压力梯度
% =========================================================================
function [dpdx,dpdz] = pressureGradientsGA(P,theta,z,par)

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

end


%% ========================================================================
%  遗传操作
% =========================================================================
function parent = tournamentSelect(population,objective,tournamentSize)

indices = randi(size(population,1),tournamentSize,1);
[~,bestLocal] = min(objective(indices));
parent = population(indices(bestLocal),:);

end


function chromosome = mutateAndRepair(chromosome,ga)

mutationMask = rand(size(chromosome)) < ga.mutationProbability;
chromosome(mutationMask) = ~chromosome(mutationMask);
chromosome = repairCoverage(chromosome,ga.nOnes);

end


function chromosome = repairCoverage(chromosome,nOnes)

chromosome = logical(chromosome);
nCurrent = nnz(chromosome);
if nCurrent > nOnes
    oneIndex = find(chromosome);
    chromosome(oneIndex(randperm(numel(oneIndex),nCurrent-nOnes))) = false;
elseif nCurrent < nOnes
    zeroIndex = find(~chromosome);
    chromosome(zeroIndex(randperm(numel(zeroIndex),nOnes-nCurrent))) = true;
end

end


%% ========================================================================
%  设计网格映射到压力网格
% =========================================================================
function pressureMask = mapDesignMask(designMask,Ntheta,Nz)

[Mtheta,Mz] = size(designMask);
thetaIndex = min(floor(linspace(0,1,Ntheta)*Mtheta)+1,Mtheta);
zIndex = min(floor(linspace(0,1,Nz)*Mz)+1,Mz);
pressureMask = logical(designMask(thetaIndex,zIndex));

end


%% ========================================================================
%  区域零散程度
% =========================================================================
function value = normalizedFragmentation(mask)

transitions = sum(abs(diff(double(mask),1,1)),'all')+ ...
    sum(abs(diff(double(mask),1,2)),'all');
maximum = (size(mask,1)-1)*size(mask,2)+ ...
    size(mask,1)*(size(mask,2)-1);
value = transitions/max(maximum,1);

end


%% ========================================================================
%  系统性能汇总
% =========================================================================
function metrics = systemMetrics(pad,par)

metrics.Fx = sum([pad.Fx]);
metrics.Fy = sum([pad.Fy]);
metrics.W = hypot(metrics.Fx,metrics.Fy);
metrics.friction = sum([pad.frictionMagnitude]);
metrics.frictionCoefficient = metrics.friction/max(metrics.W,realmin);
metrics.pMax = max([pad.pMax]);
metrics.hMin = min([pad.hMin]);
metrics.powerLoss = abs(sum([pad.frictionSigned]))*par.Rj*par.omega;

end


%% ========================================================================
%  输出CSV结果表
% =========================================================================
function writeResultTable(result,outputDir)

name = {'总承载力_N';'总摩擦力_N';'摩擦因数';'最大压力_MPa'; ...
    '最小膜厚_um';'摩擦功耗_W'};
baseline = [result.baseMetrics.W;result.baseMetrics.friction; ...
    result.baseMetrics.frictionCoefficient;result.baseMetrics.pMax/1e6; ...
    result.baseMetrics.hMin*1e6;result.baseMetrics.powerLoss];
optimized = [result.optimizedMetrics.W;result.optimizedMetrics.friction; ...
    result.optimizedMetrics.frictionCoefficient; ...
    result.optimizedMetrics.pMax/1e6;result.optimizedMetrics.hMin*1e6; ...
    result.optimizedMetrics.powerLoss];
changePercent = 100*(optimized./baseline-1);
T = table(name,baseline,optimized,changePercent);
writetable(T,fullfile(outputDir,'GA_Performance_Comparison.csv'));

end


%% ========================================================================
%  绘制并保存优化结果图
% =========================================================================
%{
旧版绘图实现仅保留为注释，当前程序统一调用下方的新绘图函数。

mainPad = result.mainPad;
basePad = result.basePad;
optPad = result.optimizedPad;
history = result.history;

%% 图1：GA收敛
fig = figure('Color','w','Position',[100,80,980,680]);
tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
ax1 = nexttile;
plot(ax1,history.bestObjective,'b-o','LineWidth',1.8, ...
    'MarkerSize',4,'DisplayName','最优目标函数');
grid(ax1,'on');
xlabel(ax1,'遗传代数');
ylabel(ax1,'目标函数 J');
legend(ax1,'Location','best');
bestRange = [min(history.bestObjective),max(history.bestObjective)];
bestMargin = max(0.02*(bestRange(2)-bestRange(1)),0.002);
ylim(ax1,[bestRange(1)-bestMargin,bestRange(2)+bestMargin]);
title(ax1,'遗传算法收敛过程');

ax2 = nexttile;
plot(ax2,100*(history.bestLoadRatio-1),'r-o','LineWidth',1.8, ...
    'MarkerSize',4,'DisplayName','承载力变化');
hold(ax2,'on');
plot(ax2,100*(history.bestFrictionRatio-1),'g-s','LineWidth',1.8, ...
    'MarkerSize',4,'DisplayName','摩擦力变化');
plot(ax2,100*(history.bestMuRatio-1),'m-^','LineWidth',1.8, ...
    'MarkerSize',4,'DisplayName','摩擦因数变化');
yline(ax2,0,'k:','HandleVisibility','off');
grid(ax2,'on');
xlabel(ax2,'遗传代数');
ylabel(ax2,'相对全滑移基准变化 / %');
legend(ax2,'Location','best');
set([ax1 ax2],'FontName','Microsoft YaHei','FontSize',11);
exportgraphics(fig,fullfile(outputDir,'GA_01_Convergence.png'), ...
    'Resolution',300);

%% 图2：主承载瓦允许滑移区和实际滑移区
fig = figure('Color','w','Position',[100,80,1100,460]);
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
thetaLocal = linspace(0,100,parFine.Ntheta);
zmm = zFine*1e3;

ax1 = nexttile;
imagesc(ax1,thetaLocal,zmm,double(result.fineAllowedMask)');
axis(ax1,'xy');
clim(ax1,[0 1]);
title(ax1,sprintf('%d号瓦GA优化允许滑移区域',mainPad));
xlabel(ax1,'入口 → 出口的周向位置 / %');
ylabel(ax1,'轴向位置 z / mm');

ax2 = nexttile;
imagesc(ax2,thetaLocal,zFineOpt*1e3, ...
    double(optPad(mainPad).slipMask)');
axis(ax2,'xy');
clim(ax2,[0 1]);
title(ax2,sprintf('%d号瓦最终实际滑移区域',mainPad));
xlabel(ax2,'入口 → 出口的周向位置 / %');
ylabel(ax2,'轴向位置 z / mm');

colormap(fig,[1 1 1;0.00 0.45 0.74]);
cb = colorbar(ax2);
cb.Ticks = [0 1];
cb.TickLabels = {'无滑移','滑移'};
set([ax1 ax2],'FontName','Microsoft YaHei','FontSize',11);
exportgraphics(fig,fullfile(outputDir,'GA_02_MainPad_SlipRegion.png'), ...
    'Resolution',300);

%% 图3：主承载瓦中截面压力比较
fig = figure('Color','w','Position',[120,100,900,580]);
ax = axes(fig);
midBase = round((size(basePad(mainPad).P,2)+1)/2);
midOpt = round((size(optPad(mainPad).P,2)+1)/2);
plot(ax,thetaLocal,basePad(mainPad).P(:,midBase)/1e6, ...
    'k--','LineWidth',2.2,'DisplayName','五瓦全滑移基准');
hold(ax,'on');
plot(ax,thetaLocal,optPad(mainPad).P(:,midOpt)/1e6, ...
    'b-','LineWidth',2.2,'DisplayName','主承载瓦GA优化');
grid(ax,'on');
xlabel(ax,'入口 → 出口的周向位置 / %');
ylabel(ax,'油膜压力 p / MPa');
title(ax,sprintf('%d号主承载瓦轴向中截面压力对比',mainPad));
legend(ax,'Location','best');
set(ax,'FontName','Microsoft YaHei','FontSize',12);
exportgraphics(fig,fullfile(outputDir,'GA_03_Pressure_Comparison.png'), ...
    'Resolution',300);

%% 图4：五瓦周向摩擦力分布比较
fig = figure('Color','w','Position',[80,70,1120,680]);
ax = axes(fig,'Position',[0.09 0.14 0.68 0.79]);
hold(ax,'on');
colors = lines(parFine.Npad);
legendText = cell(2*parFine.Npad,1);
handles = gobjects(2*parFine.Npad,1);
thetaDisplayZero = parFine.thetaPivot(1)-parFine.beta/2;

for iPad = 1:parFine.Npad
    thetaDegree = mod((basePad(iPad).theta-thetaDisplayZero)*180/pi,360);
    baseCurve = parFine.Rb*trapz(zFine, ...
        abs(basePad(iPad).tauPadX),2);
    optCurve = parOpt.Rb*trapz(zFineOpt, ...
        abs(optPad(iPad).tauPadX),2);

    handles(2*iPad-1) = plot(ax,thetaDegree,baseCurve, ...
        'Color',colors(iPad,:),'LineStyle','--','LineWidth',2);
    handles(2*iPad) = plot(ax,thetaDegree,optCurve, ...
        'Color',colors(iPad,:),'LineStyle','-','LineWidth',2.2);
    legendText{2*iPad-1} = sprintf('%d号瓦-全滑移',iPad);
    legendText{2*iPad} = sprintf('%d号瓦-GA优化后',iPad);
end

%}


%% ========================================================================
% 统一的遗传算法结果绘图
% 图1为五瓦实际滑移区域；图2~图5均按1号瓦入口为0度展开。
% 全滑移基准采用虚线，GA优化后采用实线，同一瓦块保持同一颜色。
% =========================================================================
function plotOptimizationResults(result,~,zFineOpt,~,parOpt,outputDir)

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

plotGAConvergenceUnified(result.history,outputDir);
plotGAFivePadSlipRegion(result.optimizedPad,zFineOpt,parOpt, ...
    result.mainPad,outputDir);
plotGAMainPadDesignRegion(result,zFineOpt,parOpt,outputDir);
plotGAFivePadOptimized(result.optimizedPad,zFineOpt,parOpt, ...
    'pressure',outputDir);
plotGAFivePadOptimized(result.optimizedPad,zFineOpt,parOpt, ...
    'thickness',outputDir);
plotGAFivePadOptimized(result.optimizedPad,zFineOpt,parOpt, ...
    'friction',outputDir);
plotGAFivePadOptimized(result.optimizedPad,zFineOpt,parOpt, ...
    'load',outputDir);

end


%% ========================================================================
% 附图：遗传算法收敛过程
% =========================================================================
function plotGAConvergenceUnified(history,outputDir)

fig = figure('Color','w','Position',[110,80,980,680], ...
    'Name','GA收敛过程');
tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
plot(ax1,1:numel(history.bestObjective),history.bestObjective, ...
    'b-o','LineWidth',1.8,'MarkerSize',4, ...
    'DisplayName','最优目标函数');
grid(ax1,'on');
box(ax1,'on');
xlabel(ax1,'遗传代数');
ylabel(ax1,'目标函数 J');
title(ax1,'遗传算法收敛过程');
legend(ax1,'Location','best');

ax2 = nexttile;
generation = 1:numel(history.bestLoadRatio);
plot(ax2,generation,100*(history.bestLoadRatio-1), ...
    'r-o','LineWidth',1.8,'MarkerSize',4, ...
    'DisplayName','承载力变化');
hold(ax2,'on');
plot(ax2,generation,100*(history.bestFrictionRatio-1), ...
    'g-s','LineWidth',1.8,'MarkerSize',4, ...
    'DisplayName','摩擦力变化');
plot(ax2,generation,100*(history.bestMuRatio-1), ...
    'm-^','LineWidth',1.8,'MarkerSize',4, ...
    'DisplayName','摩擦因数变化');
yline(ax2,0,'k:','HandleVisibility','off');
grid(ax2,'on');
box(ax2,'on');
xlabel(ax2,'遗传代数');
ylabel(ax2,'相对全滑移基准变化 / %');
legend(ax2,'Location','best');
set([ax1 ax2],'FontName','Microsoft YaHei','FontSize',11, ...
    'LineWidth',1.05,'TickDir','out');

saveGAFigure(fig,outputDir,'GA_00_Convergence.png');

end


%% ========================================================================
% 图1：GA优化后五瓦实际滑移区域
% 蓝色=未实际滑移；黄色=实际滑移。
% =========================================================================
function plotGAFivePadSlipRegion(pad,z,par,mainPad,outputDir)

thetaZero = pad(1).theta(1);
fig = figure('Color','w','Position',[100,80,1120,590], ...
    'Name','图1 GA优化后五瓦实际滑移区域');
ax = axes('Parent',fig,'Position',[0.09,0.15,0.79,0.76]);
hold(ax,'on');

for iPad = 1:par.Npad
    thetaDeg = (pad(iPad).theta(:)-thetaZero)*180/pi;
    zStar = z/par.B;
    imagesc(ax,thetaDeg,zStar,double(pad(iPad).slipMask).');

    if iPad == mainPad
        labelText = sprintf('%d号瓦（主承载瓦）',iPad);
        labelColor = [0.72 0.05 0.05];
    else
        labelText = sprintf('%d号瓦',iPad);
        labelColor = [0.12 0.12 0.12];
    end
    text(ax,mean(thetaDeg),0,labelText, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Microsoft YaHei','FontSize',11, ...
        'FontWeight','bold','Color','w', ...
        'BackgroundColor',labelColor,'Margin',2);
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
xlabel(ax,'周向展开角度 \theta^* / (°)，1号瓦入口为0°', ...
    'FontName','Microsoft YaHei','FontSize',14);
ylabel(ax,'轴向无量纲坐标 z/B', ...
    'FontName','Microsoft YaHei','FontSize',14);
title(ax,'主承载瓦滑移区域GA优化后的五瓦实际滑移分布', ...
    'FontName','Microsoft YaHei','FontSize',15,'FontWeight','bold');

cb = colorbar(ax,'Location','eastoutside');
cb.Ticks = [0 1];
cb.TickLabels = {'未滑移','实际滑移'};
cb.FontName = 'Microsoft YaHei';
cb.FontSize = 11;

saveGAFigure(fig,outputDir,'GA_01_SlipRegion.png');

end


%% ========================================================================
% 附图：主承载瓦GA设计区域与最终实际滑移区域
% =========================================================================
function plotGAMainPadDesignRegion(result,z,par,outputDir)

mainPad = result.mainPad;
pad = result.optimizedPad(mainPad);
thetaLocal = linspace(0,100,size(result.fineAllowedMask,1));
zStar = z/par.B;

fig = figure('Color','w','Position',[100,90,1120,470], ...
    'Name','主承载瓦GA设计区域与实际滑移区域');
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
imagesc(ax1,thetaLocal,zStar,double(result.fineAllowedMask).');
axis(ax1,'xy');
clim(ax1,[0 1]);
title(ax1,sprintf('%d号主承载瓦：GA设计允许滑移区域',mainPad));
xlabel(ax1,'入口 → 出口的周向位置 / %');
ylabel(ax1,'轴向坐标 z/B');

ax2 = nexttile;
imagesc(ax2,thetaLocal,zStar,double(pad.slipMask).');
axis(ax2,'xy');
clim(ax2,[0 1]);
title(ax2,sprintf('%d号主承载瓦：最终实际滑移区域',mainPad));
xlabel(ax2,'入口 → 出口的周向位置 / %');
ylabel(ax2,'轴向坐标 z/B');

colormap(fig,[0.08 0.32 0.88; 1.00 0.82 0.06]);
cb = colorbar(ax2,'Location','eastoutside');
cb.Ticks = [0 1];
cb.TickLabels = {'无滑移','滑移'};
set([ax1 ax2],'FontName','Microsoft YaHei','FontSize',11, ...
    'LineWidth',1.05,'TickDir','out');

saveGAFigure(fig,outputDir,'GA_01B_MainPad_Design_vs_Actual.png');

end


%% ========================================================================
% 图2~图5：只绘制GA优化后的五瓦二维周向分布
% =========================================================================
function plotGAFivePadOptimized(pad,z,par,fieldType,outputDir)

thetaZero = pad(1).theta(1);
colors = lines(par.Npad);
curves = cell(1,par.Npad);
thetaCurves = cell(1,par.Npad);
allValues = [];

for iPad = 1:par.Npad
    thetaCurves{iPad} = (pad(iPad).theta(:)-thetaZero)*180/pi;

    switch fieldType
        case 'pressure'
            midZ = round((size(pad(iPad).P,2)+1)/2);
            curves{iPad} = pad(iPad).P(:,midZ)/1e6;
            yLabelText = '轴向中截面油膜压力 p / MPa';
            titleText = 'GA优化后的五瓦油膜压力分布';
            fileName = 'GA_02_Pressure_Optimized.png';

        case 'thickness'
            curves{iPad} = pad(iPad).h(:)*1e6;
            yLabelText = '油膜厚度 h / μm';
            titleText = 'GA优化后的五瓦油膜厚度分布';
            fileName = 'GA_03_FilmThickness_Optimized.png';

        case 'friction'
            curves{iPad} = par.Rb* ...
                trapz(z,pad(iPad).tauPadMagnitude,2);
            yLabelText = '局部摩擦力 dF_f/dθ / (N·rad^{-1})';
            titleText = 'GA优化后的五瓦摩擦力分布';
            fileName = 'GA_04_Friction_Optimized.png';

        case 'load'
            curves{iPad} = par.Rb*trapz(z,pad(iPad).P,2);
            yLabelText = '局部承载力 dW/dθ / (N·rad^{-1})';
            titleText = 'GA优化后的五瓦承载力分布';
            fileName = 'GA_05_Load_Optimized.png';

        otherwise
            error('未知GA绘图变量类型：%s',fieldType);
    end

    allValues = [allValues;curves{iPad}(:)]; %#ok<AGROW>
end

valueMin = min(allValues);
valueMax = max(allValues);
valueRange = valueMax-valueMin;
if valueRange <= eps(max(abs([valueMin valueMax])))
    valueRange = max(abs(valueMax),1);
end

fig = figure('Color','w','Position',[95,70,1080,660], ...
    'Name',titleText);
ax = axes('Parent',fig,'Position',[0.10,0.14,0.74,0.79]);
hold(ax,'on');
handles = gobjects(par.Npad,1);

for iPad = 1:par.Npad
    thetaDeg = thetaCurves{iPad};
    handles(iPad) = plot(ax,thetaDeg,curves{iPad},'-', ...
        'Color',colors(iPad,:),'LineWidth',2.3, ...
        'DisplayName',sprintf('%d号瓦',iPad));

    iLabel = round((numel(thetaDeg)+1)/2);
    if strcmp(fieldType,'pressure') || strcmp(fieldType,'load')
        [~,iLabel] = max(curves{iPad});
    end
    text(ax,thetaDeg(iLabel),curves{iPad}(iLabel)+0.025*valueRange, ...
        sprintf('%d号瓦',iPad), ...
        'Color',colors(iPad,:),'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontName','Microsoft YaHei', ...
        'FontSize',10.5,'FontWeight','bold','BackgroundColor','w', ...
        'Margin',1,'Clipping','off');
end

xlim(ax,[0 360]);
xticks(ax,0:30:360);
if valueMin >= 0
    ylim(ax,[0,valueMax+0.12*valueRange]);
else
    ylim(ax,[valueMin-0.08*valueRange,valueMax+0.12*valueRange]);
end
grid(ax,'on');
box(ax,'on');
set(ax,'FontName','Microsoft YaHei','FontSize',11, ...
    'LineWidth',1.1,'TickDir','out');
xlabel(ax,'周向展开角度 \theta^* / (°)，1号瓦入口为0°', ...
    'FontSize',13);
ylabel(ax,yLabelText,'FontSize',13);
title(ax,titleText,'FontSize',14,'FontWeight','bold');
legend(ax,handles,'Location','eastoutside','FontSize',10);

saveGAFigure(fig,outputDir,fileName);

end


%% ========================================================================
% 统一保存GA图形
% =========================================================================
function saveGAFigure(fig,outputDir,fileName)

% 批处理导出前关闭坐标区交互工具栏，避免工具栏图标进入论文图片。
axesList = findall(fig,'Type','axes');
for iAx = 1:numel(axesList)
    try
        axesList(iAx).Toolbar.Visible = 'off';
        disableDefaultInteractivity(axesList(iAx));
    catch
        % 某些特殊坐标对象不支持工具栏属性，不影响图形导出。
    end
end
drawnow;
exportgraphics(fig,fullfile(outputDir,fileName),'Resolution',300);

end

%{

xlim(ax,[0 360]);
xticks(ax,0:30:360);
grid(ax,'on');
box(ax,'on');
xlabel(ax,'周向展开角度 \theta^* / (°)（1号瓦入口为0°）');
ylabel(ax,'局部周向摩擦力幅值 dF_f/d\theta / (N·rad^{-1})');
title(ax,'五瓦全滑移与主承载瓦GA优化后摩擦力分布对比');
legend(ax,handles,legendText,'Location','eastoutside', ...
    'FontSize',9);
set(ax,'FontName','Microsoft YaHei','FontSize',11);
exportgraphics(fig,fullfile(outputDir,'GA_04_Friction_Comparison.png'), ...
    'Resolution',300);

%% 图5：综合性能归一化比较
fig = figure('Color','w','Position',[160,100,900,560]);
labels = {'总承载力','总摩擦力','摩擦因数','最大压力','最小膜厚'};
normalized = [1, ...
    result.optimizedMetrics.W/result.baseMetrics.W; ...
    1,result.optimizedMetrics.friction/result.baseMetrics.friction; ...
    1,result.optimizedMetrics.frictionCoefficient/ ...
    result.baseMetrics.frictionCoefficient; ...
    1,result.optimizedMetrics.pMax/result.baseMetrics.pMax; ...
    1,result.optimizedMetrics.hMin/result.baseMetrics.hMin]*100;
bar(normalized,'grouped');
yline(100,'k:','基准','HandleVisibility','off');
grid on;
xticklabels(labels);
ylabel('相对全滑移基准 / %');
title('主承载瓦滑移区域GA优化前后综合性能');
legend({'五瓦全滑移','主承载瓦GA优化'},'Location','best');
set(gca,'FontName','Microsoft YaHei','FontSize',11);
exportgraphics(fig,fullfile(outputDir,'GA_05_Performance_Comparison.png'), ...
    'Resolution',300);

%% 图6：主承载瓦中截面油膜厚度比较
fig = figure('Color','w','Position',[120,100,900,580]);
ax = axes(fig);
plot(ax,thetaLocal,basePad(mainPad).H(:,midBase)*1e6, ...
    'k--','LineWidth',2.2,'DisplayName','五瓦全滑移基准');
hold(ax,'on');
plot(ax,thetaLocal,optPad(mainPad).H(:,midOpt)*1e6, ...
    'b-','LineWidth',2.2,'DisplayName','主承载瓦GA优化');
grid(ax,'on');
xlabel(ax,'入口 → 出口的周向位置 / %');
ylabel(ax,'油膜厚度 h / \mum');
title(ax,sprintf('%d号主承载瓦轴向中截面油膜厚度对比',mainPad));
legend(ax,'Location','best');
set(ax,'FontName','Microsoft YaHei','FontSize',12);
exportgraphics(fig,fullfile(outputDir,'GA_06_FilmThickness_Comparison.png'), ...
    'Resolution',300);

end
%}
