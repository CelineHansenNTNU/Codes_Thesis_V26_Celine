% Calculate RPD for intermediate multi-response PLSR model

clear; clc;

% Select .mat file

matFile = 'Multi_VB_res_msbackadj_SNVon.mat';

S = load(matFile);

% The struct is assumed to be called "results"
results = S.results;

% Sugar names

if isfield(results, 'sugarNames') && ~isempty(results.sugarNames)
    sugarNames = string(results.sugarNames);
else
    sugarNames = "Sugar_" + string(1:size(results.Y_trainval, 2));
end

% CV RPD: training/validation set

Yref_cv = results.Y_trainval;
Ypred_cv = results.Yhat_trainCV;

RMSEP_CV_calc = sqrt(mean((Yref_cv - Ypred_cv).^2, 1, 'omitnan'));
SD_CV = std(Yref_cv, 0, 1, 'omitnan');
RPD_CV = SD_CV ./ RMSEP_CV_calc;

% Held-out test RPD: test groups / test samples

Yref_test = results.Y_test;
Ypred_test = results.Yhat_test;

RMSEP_test_calc = sqrt(mean((Yref_test - Ypred_test).^2, 1, 'omitnan'));
SD_test = std(Yref_test, 0, 1, 'omitnan');
RPD_test = SD_test ./ RMSEP_test_calc;

%% 5) Collect sample-level results in table

sampleLevelRPD = table( ...
    sugarNames(:), ...
    RMSEP_CV_calc(:), SD_CV(:), RPD_CV(:), ...
    RMSEP_test_calc(:), SD_test(:), RPD_test(:), ...
    'VariableNames', { ...
        'Sugar', ...
        'RMSEP_CV_calc', 'SD_CV', 'RPD_CV', ...
        'RMSEP_test_calc', 'SD_test', 'RPD_test'} ...
);

disp('Sample-level RPD results:')
disp(sampleLevelRPD)

