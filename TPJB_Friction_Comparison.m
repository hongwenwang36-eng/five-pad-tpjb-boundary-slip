function TPJB_Friction_Comparison
% =========================================================================
% 五瓦可倾瓦径向轴承三种边界条件二维对比
% =========================================================================

clc;
close all;

fprintf('\n============================================================\n');
fprintf(' 正在计算无滑移、全滑移和GA优化滑移三种工况\n');
fprintf('============================================================\n');

%% 1. 无滑移和全滑移基准
[padNoSlip,zNoSlip,parNoSlip] = TPJB_NoSlip_5Pad(false);
[padFullSlip,zFullSlip,parFullSlip] = TPJB_BinarySlip_5Pad(false);

%% 2. 读取已经完成的GA优化结果
gaFile = fullfile(fileparts(mfilename('fullpath')), ...
    'GA_Optimization_Results','GA_Optimization_Result.mat');

if ~isfile(gaFile)
    error(['没有找到GA优化结果：%s\n' ...
        '请先运行 result = TPJB_GA_Optimize_MainPad;'],gaFile);
end

loaded = load(gaFile,'result');
if ~isfield(loaded,'result') || ...
        ~isfield(loaded.result,'optimizedPad')
    error('GA结果文件中缺少result.optimizedPad。');
end

padGA = loaded.result.optimizedPad;
parGA = parFullSlip;
parGA.Npad = numel(padGA);
parGA.Ntheta = size(padGA(1).P,1);
parGA.Nz = size(padGA(1).P,2);
zGA = linspace(-parGA.B/2,parGA.B/2,parGA.Nz);

%% 3. 检查三种工况是否能够直接比较
validateThreeCases(padNoSlip,zNoSlip,parNoSlip, ...
    padFullSlip,zFullSlip,parFullSlip,padGA,zGA,parGA);

caseData(1) = makeCaseData('无滑移',':',2.25, ...
    padNoSlip,zNoSlip,parNoSlip,'noSlip');
caseData(2) = makeCaseData('全滑移','--',2.25, ...
    padFullSlip,zFullSlip,parFullSlip,'slip');
caseData(3) = makeCaseData('GA优化滑移','-',2.60, ...
    padGA,zGA,parGA,'slip');

%% 4. 绘制两张二维对比图
frictionIntegral = plotThreeCaseField(caseData,'friction', ...
    'TPJB_ThreeCases_Friction_Comparison.png');
loadIntegral = plotThreeCaseField(caseData,'load', ...
    'TPJB_ThreeCases_Load_Comparison.png');

%% 5. 输出数值汇总
padNumber = (1:parNoSlip.Npad)';
frictionTable = table(padNumber,frictionIntegral(:,1), ...
    frictionIntegral(:,2),frictionIntegral(:,3), ...
    'VariableNames',{'Pad','NoSlip_N','FullSlip_N','GAOptimized_N'});

loadTable = table(padNumber,loadIntegral(:,1), ...
    loadIntegral(:,2),loadIntegral(:,3), ...
    'VariableNames',{'Pad','NoSlipNormalLoad_N', ...
    'FullSlipNormalLoad_N','GAOptimizedNormalLoad_N'});

fprintf('\n各瓦周向摩擦力（由二维曲线积分）：\n');
disp(frictionTable);
fprintf('三种工况摩擦力合计：无滑移=%.8f N，全滑移=%.8f N，GA=%.8f N\n', ...
    sum(frictionIntegral,1));

fprintf('\n各瓦局部法向承载力（由二维曲线积分）：\n');
disp(loadTable);

WNoSlip = hypot(sum([padNoSlip.Fx]),sum([padNoSlip.Fy]));
WFullSlip = hypot(sum([padFullSlip.Fx]),sum([padFullSlip.Fy]));
WGA = hypot(sum([padGA.Fx]),sum([padGA.Fy]));
fprintf('三种工况总油膜力矢量幅值：无滑移=%.8f N，全滑移=%.8f N，GA=%.8f N\n', ...
    WNoSlip,WFullSlip,WGA);

end


%% ========================================================================
% 构造统一工况结构
% =========================================================================
function data = makeCaseData(name,lineStyle,lineWidth,pad,z,par,type)

data.name = name;
data.lineStyle = lineStyle;
data.lineWidth = lineWidth;
data.pad = pad;
data.z = z;
data.par = par;
data.type = type;

end


%% ========================================================================
% 检查几何、工况和网格
% =========================================================================
function validateThreeCases(padNoSlip,zNoSlip,parNoSlip, ...
    padFullSlip,zFullSlip,parFullSlip,padGA,zGA,parGA)

referenceFields = {'Npad','Rb','Rj','B','C0','n','epsilon','phi'};
parList = {parFullSlip,parGA};

for iCase = 1:numel(parList)
    parOther = parList{iCase};
    for k = 1:numel(referenceFields)
        name = referenceFields{k};
        a = double(parNoSlip.(name));
        b = double(parOther.(name));
        tolerance = 1e-12*max(1,max(abs(a(:))));
        if ~isequal(size(a),size(b)) || max(abs(a(:)-b(:))) > tolerance
            error('三种工况的参数%s不一致，不能直接比较。',name);
        end
    end
end

padList = {padNoSlip,padFullSlip,padGA};
zList = {zNoSlip,zFullSlip,zGA};

for iCase = 1:numel(padList)
    pads = padList{iCase};
    zNow = zList{iCase};
    if numel(pads) ~= parNoSlip.Npad
        error('第%d种工况的瓦块数量不一致。',iCase);
    end
    if numel(zNow) ~= numel(zNoSlip) || ...
            max(abs(zNow(:)-zNoSlip(:))) > 1e-14
        error('第%d种工况的轴向网格不一致。',iCase);
    end
    for iPad = 1:parNoSlip.Npad
        if numel(pads(iPad).theta) ~= numel(padNoSlip(iPad).theta) || ...
                max(abs(pads(iPad).theta(:)- ...
                padNoSlip(iPad).theta(:))) > 1e-14
            error('第%d种工况第%d块瓦的周向网格不一致。', ...
                iCase,iPad);
        end
    end
end

end


%% ========================================================================
% 绘制三种工况的五瓦二维分布
% =========================================================================
function integralValues = plotThreeCaseField(caseData,fieldType,fileName)

nCase = numel(caseData);
nPad = caseData(1).par.Npad;
colors = lines(nPad);
thetaZero = caseData(1).pad(1).theta(1);

thetaCurve = cell(nPad,1);
fieldCurve = cell(nPad,nCase);
integralValues = zeros(nPad,nCase);
allValues = [];

for iPad = 1:nPad
    thetaCurve{iPad} = ...
        (caseData(1).pad(iPad).theta(:)-thetaZero)*180/pi;

    for iCase = 1:nCase
        pad = caseData(iCase).pad(iPad);
        z = caseData(iCase).z;
        par = caseData(iCase).par;

        switch fieldType
            case 'friction'
                if strcmp(caseData(iCase).type,'noSlip')
                    tauX = abs(pad.tauPad);
                else
                    tauX = abs(pad.tauPadX);
                end
                fieldCurve{iPad,iCase} = par.Rb*trapz(z,tauX,2);

            case 'load'
                fieldCurve{iPad,iCase} = par.Rb*trapz(z,pad.P,2);

            otherwise
                error('未知绘图变量类型：%s',fieldType);
        end

        integralValues(iPad,iCase) = trapz( ...
            pad.theta,fieldCurve{iPad,iCase});
        allValues = [allValues;fieldCurve{iPad,iCase}(:)]; %#ok<AGROW>
    end
end

valueMin = min(allValues);
valueMax = max(allValues);
valueRange = valueMax-valueMin;
if valueRange <= eps(max(abs([valueMin valueMax])))
    valueRange = max(abs(valueMax),1);
end

switch fieldType
    case 'friction'
        figureName = '三种滑移状态五瓦摩擦力分布对比';
        yLabelText = '局部摩擦力 dF_f/d\theta / (N·rad^{-1})';
        caption = '(a) 无滑移、全滑移与GA优化滑移摩擦力对比';
    case 'load'
        figureName = '三种滑移状态五瓦承载力分布对比';
        yLabelText = '局部承载力 dW/d\theta / (N·rad^{-1})';
        caption = '(b) 无滑移、全滑移与GA优化滑移承载力对比';
end

fig = figure('Color','w','Position',[75,60,1220,720], ...
    'Name',figureName,'NumberTitle','on');
ax = axes('Parent',fig,'Position',[0.085,0.15,0.72,0.78]);
hold(ax,'on');

for iPad = 1:nPad
    for iCase = 1:nCase
        plot(ax,thetaCurve{iPad},fieldCurve{iPad,iCase}, ...
            'Color',colors(iPad,:), ...
            'LineStyle',caseData(iCase).lineStyle, ...
            'LineWidth',caseData(iCase).lineWidth, ...
            'HandleVisibility','off');
    end

    % 在GA优化曲线上标注瓦号。
    if strcmp(fieldType,'load')
        [~,labelIndex] = max(fieldCurve{iPad,3});
    else
        labelIndex = round((numel(thetaCurve{iPad})+1)/2);
    end
    labelY = fieldCurve{iPad,3}(labelIndex)+ ...
        (0.025+0.006*mod(iPad,2))*valueRange;
    text(ax,thetaCurve{iPad}(labelIndex),labelY, ...
        sprintf('%d号瓦',iPad),'Color',colors(iPad,:), ...
        'FontName','Microsoft YaHei','FontSize',10.5, ...
        'FontWeight','bold','HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','BackgroundColor','w', ...
        'Margin',1,'Clipping','off');
end

% 紧凑图例：5种颜色表示瓦号，3种线型表示工况。
padLegend = gobjects(nPad,1);
for iPad = 1:nPad
    padLegend(iPad) = plot(ax,nan,nan,'-', ...
        'Color',colors(iPad,:),'LineWidth',2.5, ...
        'DisplayName',sprintf('%d号瓦',iPad));
end

caseLegend = gobjects(nCase,1);
for iCase = 1:nCase
    caseLegend(iCase) = plot(ax,nan,nan, ...
        'Color','k','LineStyle',caseData(iCase).lineStyle, ...
        'LineWidth',caseData(iCase).lineWidth, ...
        'DisplayName',caseData(iCase).name);
end

xlim(ax,[0 360]);
xticks(ax,0:30:360);
if valueMin >= 0
    ylim(ax,[0,valueMax+0.13*valueRange]);
else
    ylim(ax,[valueMin-0.08*valueRange,valueMax+0.13*valueRange]);
end
grid(ax,'on');
box(ax,'on');
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
set(ax,'FontName','Times New Roman','FontSize',12, ...
    'LineWidth',1.1,'TickDir','out');

xlabel(ax,'周向展开角度 \theta^* / (°)，1号瓦入口为0°', ...
    'FontName','Microsoft YaHei','FontSize',14);
ylabel(ax,yLabelText,'FontName','Microsoft YaHei','FontSize',14);
title(ax,figureName,'FontName','Microsoft YaHei', ...
    'FontSize',15,'FontWeight','bold');
legend(ax,[padLegend;caseLegend],'Location','eastoutside', ...
    'NumColumns',1,'FontName','Microsoft YaHei','FontSize',10);

annotation(fig,'textbox',[0.20,0.012,0.55,0.055], ...
    'String',caption,'HorizontalAlignment','center', ...
    'VerticalAlignment','middle','EdgeColor','none', ...
    'FontName','Microsoft YaHei','FontSize',14);

axesList = findall(fig,'Type','axes');
for iAx = 1:numel(axesList)
    try
        axesList(iAx).Toolbar.Visible = 'off';
        disableDefaultInteractivity(axesList(iAx));
    catch
    end
end
drawnow;

outputFile = fullfile(fileparts(mfilename('fullpath')),fileName);
exportgraphics(fig,outputFile,'Resolution',300);
fprintf('%s已保存至：%s\n',figureName,outputFile);

end
