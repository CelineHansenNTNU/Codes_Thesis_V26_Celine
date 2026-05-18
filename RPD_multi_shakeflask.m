% Calculated RPD from saved multi-response PLSR .mat file

clear; clc;

% Path to .mat file
matFile = 'SF_VB_70_30_MultiPLSR_SF8_asls_noSNV';

% Load file
S = load(matFile);
out = S.out;

% Define sugar names in same order as Y columns
sugarNames = ["Glucose", "Xylose", "Mannitol"];

fprintf('\nRPD results for: %s\n\n', matFile);

% CV RPD: Venetian blinds CV on 70% training/calibration set

Yref_cv = out.Y_trainval;
Ypred_cv = out.Yhat_trainCV;

RMSE_CV = sqrt(mean((Yref_cv - Ypred_cv).^2, 1, 'omitnan'));
SD_CV = std(Yref_cv, 0, 1, 'omitnan');
RPD_CV = SD_CV ./ RMSE_CV;

% Internal test RPD: 30% held-out test set from experiment 2

Yref_test = out.Y_test;
Ypred_test = out.Yhat_test;

RMSEP_test = sqrt(mean((Yref_test - Ypred_test).^2, 1, 'omitnan'));
SD_test = std(Yref_test, 0, 1, 'omitnan');
RPD_test = SD_test ./ RMSEP_test;

% External test RPD: experiment 3

if isfield(out, 'Y_external') && isfield(out, 'Yhat_external') ...
        && ~isempty(out.Y_external) && ~isempty(out.Yhat_external)

    Yref_ext = out.Y_external;
    Ypred_ext = out.Yhat_external;

    RMSEP_external = sqrt(mean((Yref_ext - Ypred_ext).^2, 1, 'omitnan'));
    SD_external = std(Yref_ext, 0, 1, 'omitnan');
    RPD_external = SD_external ./ RMSEP_external;

else
    RMSEP_external = nan(1, numel(sugarNames));
    SD_external = nan(1, numel(sugarNames));
    RPD_external = nan(1, numel(sugarNames));
end


% Collect sample-level RPD results in table

sampleLevelTable = table( ...
    sugarNames', ...
    RMSE_CV', SD_CV', RPD_CV', ...
    RMSEP_test', SD_test', RPD_test', ...
    RMSEP_external', SD_external', RPD_external', ...
    'VariableNames', { ...
        'Sugar', ...
        'RMSE_CV', 'SD_CV', 'RPD_CV', ...
        'RMSEP_test', 'SD_test', 'RPD_test', ...
        'RMSEP_external', 'SD_external', 'RPD_external'} ...
);

disp('Sample-level RPD results:')
disp(sampleLevelTable)
