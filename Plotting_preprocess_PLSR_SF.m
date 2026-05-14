clear; clc; close all;

% Script for plotting single-response PLSR model results for shake flask
% Celine + Copilot, updated 12.05.26
% User input, change files for wanted shake flask for plotting

resultFiles = {
    '/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/SF_VB_70_30_SinglePLSR_SF7_xylose_asls_noSNV.mat'
    '/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/SF_VB_70_30_SinglePLSR_SF7_xylose_asls_withSNV.mat'
    '/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/SF_VB_70_30_SinglePLSR_SF7_xylose_emsc.mat'
    '/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/SF_VB_70_30_SinglePLSR_SF7_xylose_MS_noSNV.mat'
    '/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/SF_VB_70_30_SinglePLSR_SF7_xylose_MS_withSNV.mat'
    '/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/SF_VB_70_30_SinglePLSR_SF7_xylose_MS_withSNV_maxLV20.mat'
};

methodNames = {
    'AsLS + SG'
    'AsLS + SG + SNV'
    'EMSC'
    'msbackadj + SG'
    'msbackadj + SG + SNV'
    'msbackadj + SG + SNV, maxLV = 20'
};

bestMethodIdx = 5; % mark the best method in the plot
targetSugar = "xylose"; % change to wanted sugar

relErrThreshold = 0.1; % threshold for error analysis to avoid large error when HPLC values are zero
showErrorAsLine = true;

% Load the results from the provided files

nMethods = numel(resultFiles);

allPred = [];
timeNumRef = [];
timeLabelRef = [];
measuredRef = [];

for k = 1:nMethods

    S = load(resultFiles{k});
    out = S.out;

    T = out.predTablePredict;

    [~, ord] = sort(T.timeNum);
    T = T(ord,:);

    if isempty(timeNumRef)
        timeNumRef = T.timeNum;
        timeLabelRef = string(T.timeLabel);
        measuredRef = T.measuredConc;
        allPred = nan(numel(timeNumRef), nMethods);
    end

    allPred(:,k) = T.predictedConc;
end

% Calculating the error

measuredRef = measuredRef(:);
allPred = double(allPred);

signedError = allPred - measuredRef;
absError = abs(signedError);
sqError  = signedError.^2;

% Overall RMSEP across all timepoints for each preprocessing method
RMSEP_method = sqrt(mean(sqError, 1, 'omitnan'));

% Error per timepoint
RMSEP_timepoint = absError;

% Relative error only where true HPLC value is sufficiently high
relError = nan(size(allPred));
validRel = measuredRef > relErrThreshold;

for k = 1:nMethods
    relError(validRel,k) = signedError(validRel,k) ./ measuredRef(validRel);
end

% Plotting

figure('Name', sprintf('Preprocessing comparison for Exp 3 SF7 - %s', targetSugar), ...
    'Units','normalized','Position',[0.05 0.08 0.90 0.78]);

tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

C = lines(nMethods);

% 1) Time-series: all methods + true HPLC

nexttile;
hold on;

plot(timeNumRef, measuredRef, 'k-o', ...
    'LineWidth', 2.5, ...
    'MarkerSize', 6, ...
    'DisplayName', 'HPLC measured');

for k = 1:nMethods

    if k == bestMethodIdx
        lw = 2.8;
        ms = 7;
        marker = 'o';
    else
        lw = 1.2;
        ms = 5;
        marker = '.';
    end

    plot(timeNumRef, allPred(:,k), ...
        '-', ...
        'Color', C(k,:), ...
        'LineWidth', lw, ...
        'Marker', marker, ...
        'MarkerSize', ms, ...
        'DisplayName', methodNames{k});
end

grid on;
xlabel('Time (h)');
ylabel(sprintf('%s concentration (mM)', targetSugar));
title('Predicted and measured concentration over time');
xticks(timeNumRef);
xticklabels(timeLabelRef);
xtickangle(45);
legend('Location','eastoutside','Interpreter','none');

% 2) Overall RMSEP per preprocessing method

nexttile;
b = bar(RMSEP_method);
grid on;
ylabel('RMSEP (mM)');
title('Overall RMSEP per preprocessing method on external dataset');

xticks(1:nMethods);
xticklabels(methodNames);
xtickangle(35);

% Mark best method
hold on;
plot(bestMethodIdx, RMSEP_method(bestMethodIdx), 'kp', ...
    'MarkerSize', 14, ...
    'MarkerFaceColor','k');

% 3) Prediction error per timepoint

nexttile;
hold on;

for k = 1:nMethods
    if k == bestMethodIdx
        lw = 2.8;
        ms = 7;
        marker = 'o';
    else
        lw = 1.2;
        ms = 5;
        marker = '.';
    end

    plot(timeNumRef, RMSEP_timepoint(:,k), ...
        '-', ...
        'Color', C(k,:), ...
        'LineWidth', lw, ...
        'Marker', marker, ...
        'MarkerSize', ms, ...
        'DisplayName', methodNames{k});
end

grid on;
ylabel('|Predicted - measured| (mM)');
xlabel('Time (h)');
title('Absolute prediction error per timepoint');

xticks(timeNumRef);
xticklabels(timeLabelRef);
xtickangle(45);
legend(methodNames, 'Location','eastoutside','Interpreter','none');

% 4) Relative error normalized by true value

nexttile;
hold on;

for k = 1:nMethods
    if k == bestMethodIdx
        lw = 2.8;
        ms = 7;
        marker = 'o';
    else
        lw = 1.2;
        ms = 5;
        marker = '.';
    end

    plot(timeNumRef, 100*relError(:,k), ...
        '-', ...
        'Color', C(k,:), ...
        'LineWidth', lw, ...
        'Marker', marker, ...
        'MarkerSize', ms, ...
        'DisplayName', methodNames{k});
end

grid on;
yline(0,'k--');
ylabel('Relative error (%)');
xlabel('Time (h)');
title(sprintf('Relative error, shown only when measured xylose > %.1f mM', relErrThreshold));

xticks(timeNumRef);
xticklabels(timeLabelRef);
xtickangle(45);
legend(methodNames, 'Location','eastoutside','Interpreter','none');

title(tl, sprintf('Comparison of preprocessing methods for Exp 3 SF7 %s', targetSugar));