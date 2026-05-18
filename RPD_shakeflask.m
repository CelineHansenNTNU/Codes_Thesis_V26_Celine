% Calculate RPD from saved single-response PLSR .mat file

clear; clc;

% Path to .mat file
matFile = 'SF_VB_70_30_MultiPLSR_SF9_asls_noSNV';

% Load file
S = load(matFile);
out = S.out;

fprintf('\nRPD results for: %s\n', matFile);
fprintf('Target sugar: %s\n\n', out.targetSugar);

% CV RPD: Venetian Blinds CV on 70% training/calibration set
Yref_cv = out.Y_trainval;
Ypred_cv = out.Yhat_trainCV;
valid_cv = ~isnan(Yref_cv) & ~isnan(Ypred_cv);
RMSE_CV = sqrt(mean((Yref_cv(valid_cv) - Ypred_cv(valid_cv)).^2));
SD_CV = std(Yref_cv(valid_cv));

RPD_CV = SD_CV/RMSE_CV;

% Internal test RPD: 30% held-out test set from experiment 2
Yref_test = out.Y_test;
Ypred_test = out.Yhat_test;
valid_test = ~isnan(Yref_test) & ~isnan(Ypred_test);
RMSEP_test = sqrt(mean((Yref_test(valid_test) - Ypred_test(valid_test)).^2));
SD_test = std(Yref_test(valid_test));

RPD_test = SD_test / RMSEP_test;

% External test RPD: experiment 3
if isfield(out, 'Y_external') && isfield(out, 'Yhat_external') && ~isempty(out.Y_external) && ~isempty(out.Yhat_external)
    Yref_external = out.Y_external;
    Ypred_external = out.Yhat_external;
    valid_external = ~isnan(Yref_external) & ~isnan(Ypred_external);
    RMSEP_external = sqrt(mean((Yref_external(valid_external) - Ypred_external(valid_external)).^2));
    SD_external = std(Yref_external(valid_external));

    RPD_external = SD_external / RMSEP_external;
else
    RMSEP_external = NaN;
    SD_external = NaN;
    RPD_external = NaN;
end

% Collect results in a table
RPD_table = table(...
    ["CV train/valiation"; "Internal test"; "External test"], ...
    [RMSE_CV; RMSEP_test; RMSEP_external], ...
    [SD_CV; SD_test; SD_external], ...
    [RPD_CV; RPD_test; RPD_external], ...
    'VariableNames',{'Evaluation', 'RMSEP', 'SD_reference', 'RPD'});

disp(RPD_table);