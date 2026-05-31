function out = SF_AllExp_VB_70_30_UnivariatePLSR_MeanSpectra(rootDirs, hplcExcelPath, opts)
% SF_AllExp_VB_70_30_UnivariatePLSR_MeanSpectra
%
% Single-response PLSR for shake flask mean spectra using ALL shake flask experiments together:
%   - one mean RAW spectrum per shakeflask per timepoint (SF/T)
%   - build datasets from several Raman root folders and HPLC sheets
%   - mix all experiments into one dataset
%   - random 70/30 split on SF/T mean samples
%   - Venetian blinds CV on the 70% train/validation set
%   - fit final model on all 70%
%   - predict held-out internal 30%
%
% Example:
%   rootDirs = {fullfile(pwd,'Shakeflasks_v1'), ...
%               fullfile(pwd,'Shakeflasks_v2'), ...
%               fullfile(pwd,'Shakeflasks_v3')};
%   opts.hplcSheets = {'EXP1','EXP2','EXP3'};
%   opts.targetSugar = "glucose";
%   out = SF_AllExp_VB_70_30_UnivariatePLSR_MeanSpectra(rootDirs, hplcExcelPath, opts);
%
% Celine Hansen + adapted for pooled internal testing across all SF experiments + ChatGPT
% Date: 21.05.26

    if nargin < 3
        opts = struct();
    end

    if ischar(rootDirs) || isstring(rootDirs)
        rootDirs = cellstr(rootDirs);
    end

    if ~iscell(rootDirs)
        error('rootDirs must be a cell array, string array, or char/string path.');
    end

    nExp = numel(rootDirs);

    if ~isfield(opts,'shiftRange'),         opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'power_mw'),           opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),               opts.t_ms = 1000; end
    if ~isfield(opts,'snvOn'),              opts.snvOn = false; end
    if ~isfield(opts,'maxLV_wish'),         opts.maxLV_wish = 10; end
    if ~isfield(opts,'lvRule'),             opts.lvRule = "min"; end
    if ~isfield(opts,'tolFactor'),          opts.tolFactor = 1.05; end
    if ~isfield(opts,'hplcSheets'),         opts.hplcSheets = arrayfun(@(i) sprintf('EXP%d',i), 1:nExp, 'UniformOutput', false); end
    if ~isfield(opts,'expNames'),           opts.expNames = opts.hplcSheets; end
    if ~isfield(opts,'naMeansZero'),        opts.naMeansZero = true; end
    if ~isfield(opts,'saveMat'),            opts.saveMat = true; end
    if ~isfield(opts,'doPlot'),             opts.doPlot = true; end
    if ~isfield(opts,'trainFraction'),      opts.trainFraction = 0.70; end
    if ~isfield(opts,'randomSeed'),         opts.randomSeed = 1; end
    if ~isfield(opts,'nBlinds'),            opts.nBlinds = 10; end
    if ~isfield(opts,'internalCV_folds'),   opts.internalCV_folds = 10; end
    if ~isfield(opts,'preprocMethod'),      opts.preprocMethod = "asls"; end
    if ~isfield(opts,'emscPolyOrder'),      opts.emscPolyOrder = 1; end
    if ~isfield(opts,'targetSugar'),        opts.targetSugar = "glucose"; end
    if ~isfield(opts,'tagDuplicateKeys'),   opts.tagDuplicateKeys = true; end

    hplcSheets = cellstr(opts.hplcSheets);
    expNames   = cellstr(opts.expNames);

    if numel(hplcSheets) ~= nExp
        error('opts.hplcSheets must have the same length as rootDirs.');
    end
    if numel(expNames) ~= nExp
        error('opts.expNames must have the same length as rootDirs.');
    end

    targetSugar = lower(string(opts.targetSugar));
    sugarNames = ["glucose","xylose","mannitol"];

    if ~ismember(targetSugar, sugarNames)
        error('opts.targetSugar must be "glucose", "xylose", or "mannitol".');
    end

    sugarIdx = find(sugarNames == targetSugar, 1);

    msgs = {};
    msgs{end+1} = sprintf('\n=== SF_AllExp_VB_70_30_UnivariatePLSR_MeanSpectra ===');
    msgs{end+1} = sprintf('HPLC Excel:         %s', hplcExcelPath);
    msgs{end+1} = sprintf('Target sugar:       %s', targetSugar);
    msgs{end+1} = sprintf('Preprocessing:      %s | SNV = %d', string(opts.preprocMethod), opts.snvOn);
    msgs{end+1} = sprintf('Train fraction:     %.2f | Test fraction: %.2f', opts.trainFraction, 1-opts.trainFraction);
    msgs{end+1} = sprintf('Number of pooled experiments: %d', nExp);

    %% PART A: BUILD AND MIX ALL EXPERIMENT DATASETS
    X_all = [];
    Xraw_all = [];
    Y_all = [];
    G_all = [];
    keys_all = {};
    timeNum_all = [];
    timeLabel_all = strings(0,1);
    nRep_all = [];
    expID_all = [];
    expName_all = strings(0,1);
    commonShift = [];

    for e = 1:nExp
        rootDir_e = rootDirs{e};
        sheet_e = hplcSheets{e};
        expName_e = string(expNames{e});

        msgs{end+1} = sprintf('--- Building %s from root: %s | HPLC sheet: %s', expName_e, rootDir_e, sheet_e);

        hplcMap_e = buildHPLCMap_Shakeflask_Exact(hplcExcelPath, sheet_e, opts.naMeansZero);

        [X_e, Y_e_all, G_e, keys_e, timeNum_e, timeLabel_e, shift_e, nRep_e, Xraw_e] = ...
            build_mean_dataset_from_root_local(rootDir_e, hplcMap_e, opts);

        assert(~isempty(X_e), 'No mean spectra could be built for %s.', rootDir_e);

        shift_e = shift_e(:).';
        if isempty(commonShift)
            commonShift = shift_e;
        else
            if numel(shift_e) ~= numel(commonShift) || any(abs(shift_e - commonShift) > 1e-10)
                error('Raman shift grid mismatch between experiments.');
            end
        end

        if opts.tagDuplicateKeys
            keys_e = string(expName_e) + "_" + string(keys_e);
            keys_e = cellstr(keys_e(:));
        end

        X_all = [X_all; X_e]; %#ok<AGROW>
        Xraw_all = [Xraw_all; Xraw_e]; %#ok<AGROW>
        Y_all = [Y_all; Y_e_all]; %#ok<AGROW>
        G_all = [G_all; G_e(:)]; %#ok<AGROW>
        keys_all = [keys_all; keys_e(:)]; %#ok<AGROW>
        timeNum_all = [timeNum_all; timeNum_e(:)]; %#ok<AGROW>
        timeLabel_all = [timeLabel_all; timeLabel_e(:)]; %#ok<AGROW>
        nRep_all = [nRep_all; nRep_e(:)]; %#ok<AGROW>
        expID_all = [expID_all; e*ones(size(X_e,1),1)]; %#ok<AGROW>
        expName_all = [expName_all; repmat(expName_e, size(X_e,1), 1)]; %#ok<AGROW>

        msgs{end+1} = sprintf('  Added %d SF/T mean spectra from %s.', size(X_e,1), expName_e);
    end

    assert(~isempty(X_all), 'No pooled mean spectra could be built.');

    Y_target = Y_all(:, sugarIdx);

    msgs{end+1} = sprintf('Pooled dataset: %d SF/T samples from %d experiments.', size(X_all,1), nExp);
    for e = 1:nExp
        msgs{end+1} = sprintf('  %s: %d samples', string(expNames{e}), sum(expID_all == e));
    end

    %% PART B: RANDOM 70/30 SPLIT ACROSS ALL POOLED SF/T SAMPLES
    nSamples = size(X_all,1);
    rng(opts.randomSeed);

    perm = randperm(nSamples);
    nTrain = max(2, round(opts.trainFraction * nSamples));
    nTrain = min(nTrain, nSamples-1);

    idxTrain = false(nSamples,1);
    idxTrain(perm(1:nTrain)) = true;
    idxTest = ~idxTrain;

    X_trainval = X_all(idxTrain,:);
    Y_trainval = Y_target(idxTrain,:);
    G_trainval = G_all(idxTrain);
    keys_trainval = keys_all(idxTrain);
    timeNum_trainval = timeNum_all(idxTrain);
    timeLabel_trainval = timeLabel_all(idxTrain);
    nRep_trainval = nRep_all(idxTrain);
    expID_trainval = expID_all(idxTrain);
    expName_trainval = expName_all(idxTrain);

    X_test = X_all(idxTest,:);
    Y_test = Y_target(idxTest,:);
    G_test = G_all(idxTest);
    keys_test = keys_all(idxTest);
    timeNum_test = timeNum_all(idxTest);
    timeLabel_test = timeLabel_all(idxTest);
    nRep_test = nRep_all(idxTest);
    expID_test = expID_all(idxTest);
    expName_test = expName_all(idxTest);

    msgs{end+1} = sprintf('Train/validation SF/T samples: %d | Internal test SF/T samples: %d', size(X_trainval,1), size(X_test,1));
    for e = 1:nExp
        msgs{end+1} = sprintf('  %s split: train/val=%d | test=%d', string(expNames{e}), sum(expID_trainval == e), sum(expID_test == e));
    end

    %% PART C: OPTIONAL EMSC AFTER SPLIT, FITTED ONLY ON TRAIN/VALIDATION
    EMSC_model = [];
    if lower(string(opts.preprocMethod)) == "emsc"
        msgs{end+1} = sprintf('Applying EMSC preprocessing fitted on train/validation set only.');
        [X_trainval, EMSC_model, EMSC_train_residuals, EMSC_train_params] = fit_emsc_local(X_trainval, commonShift, keys_trainval, opts.emscPolyOrder);
        [X_test, EMSC_test_residuals, EMSC_test_params] = apply_emsc_model_local(X_test, commonShift, keys_test, EMSC_model);
    else
        EMSC_train_residuals = [];
        EMSC_train_params = [];
        EMSC_test_residuals = [];
        EMSC_test_params = [];
    end

    %% PART D: VENETIAN BLINDS ON THE 70% TRAIN/VALIDATION SET
    nBlinds = min(opts.nBlinds, size(X_trainval,1));
    foldID_train = make_venetian_blinds_folds_simple(size(X_trainval,1), nBlinds);

    Yhat_trainCV = nan(size(Y_trainval));
    bestLV_perFold = nan(nBlinds,1);
    RMSECV_all = nan(opts.maxLV_wish, nBlinds);
    RMSECV_perFold = nan(nBlinds,1);

    msgs{end+1} = sprintf('Venetian blinds CV on train/validation set: %d blinds | n=%d mean spectra', nBlinds, size(X_trainval,1));

    for b = 1:nBlinds
        valIdx = (foldID_train == b);
        trainIdx = ~valIdx;

        if ~any(valIdx)
            warning('Blind %d has no validation samples; skipping.', b);
            continue;
        end

        X_train = X_trainval(trainIdx,:);
        Y_train = Y_trainval(trainIdx);
        X_val = X_trainval(valIdx,:);
        Y_val = Y_trainval(valIdx);

        maxLV_try = min([opts.maxLV_wish, size(X_train,1)-2, size(X_train,2)]);
        maxLV_try = max(maxLV_try, 1);

        RMSECV = nan(maxLV_try,1);
        innerFolds = min(opts.internalCV_folds, size(X_train,1)-1);
        innerFolds = max(innerFolds, 2);

        for a = 1:maxLV_try
            try
                [~,~,~,~,~,~,MSE_CV] = plsregress(X_train, Y_train, a, 'CV', innerFolds);
                RMSECV(a) = sqrt(MSE_CV(2, a+1));
            catch ME
                warning('plsregress CV failed for LV=%d (blind %d): %s', a, b, ME.message);
                RMSECV(a) = Inf;
            end
        end

        RMSECV_all(1:maxLV_try,b) = RMSECV;

        [minCV, idxMin] = min(RMSECV);
        bestLV = idxMin;

        if opts.lvRule == "tol"
            cand = find(RMSECV <= opts.tolFactor * minCV, 1, 'first');
            if ~isempty(cand)
                bestLV = cand;
            end
        end

        bestLV_perFold(b) = bestLV;

        [~,~,~,~,BETA_best] = plsregress(X_train, Y_train, bestLV);
        Ypred_val = [ones(size(X_val,1),1) X_val] * BETA_best;
        Ypred_val = max(Ypred_val, 0);

        Yhat_trainCV(valIdx) = Ypred_val;

        rmse_fold = sqrt(mean((Y_val - Ypred_val).^2, 'omitnan'));
        RMSECV_perFold(b) = rmse_fold;

        msgs{end+1} = sprintf('Blind %d/%d | nTrain = %d nVal = %d | bestLV = %d | RMSECV = %.3f mM', b, nBlinds, size(X_train,1), size(X_val,1), bestLV, rmse_fold);
    end

    %% PART E: INTERNAL CV METRICS ON TRAIN/VALIDATION SET
    rmse_cv = sqrt(mean((Y_trainval - Yhat_trainCV).^2, 'omitnan'));
    R2_cv = calcR2_local(Y_trainval, Yhat_trainCV);

    msgs{end+1} = sprintf('Internal Venetian blinds CV on 70%% train/validation set:');
    msgs{end+1} = sprintf('  RMSECV = %.3f mM', rmse_cv);
    msgs{end+1} = sprintf('  R^2    = %.3f', R2_cv);

    %% PART F: SELECT FINAL LV FROM MEAN CV CURVE
    meanRMSECV = mean(RMSECV_all, 2, 'omitnan');
    validLV = find(~isnan(meanRMSECV) & isfinite(meanRMSECV));

    if isempty(validLV)
        finalLV = 1;
    else
        [minCV, idxMin] = min(meanRMSECV(validLV));
        finalLV = validLV(idxMin);

        if opts.lvRule == "tol"
            cand = validLV(meanRMSECV(validLV) <= opts.tolFactor * minCV);
            if ~isempty(cand)
                finalLV = cand(1);
            end
        end
    end

    finalLV = max(1, min(finalLV, min([opts.maxLV_wish, size(X_trainval,1)-2, size(X_trainval,2)])));
    msgs{end+1} = sprintf('Selected final LV = %d', finalLV);

    %% PART G: FIT FINAL MODEL ON 70%, PREDICT HELD-OUT INTERNAL 30%
    [XL, YL, XS_train, YS_train, BETA_locked, PCTVAR_locked, MSE_locked, stats_locked] = ...
        plsregress(X_trainval, Y_trainval, finalLV); %#ok<ASGLU>

    Yhat_train_fit = [ones(size(X_trainval,1),1) X_trainval] * BETA_locked;
    Yhat_test = [ones(size(X_test,1),1) X_test] * BETA_locked;
    Yhat_train_fit = max(Yhat_train_fit, 0);
    Yhat_test = max(Yhat_test, 0);

    rmse_fit = sqrt(mean((Y_trainval - Yhat_train_fit).^2, 'omitnan'));
    R2_fit = calcR2_local(Y_trainval, Yhat_train_fit);
    rmse_test = sqrt(mean((Y_test - Yhat_test).^2, 'omitnan'));
    R2_test = calcR2_local(Y_test, Yhat_test);

    msgs{end+1} = sprintf('Prediction on held-out internal 30%% set:');
    msgs{end+1} = sprintf('  RMSEP = %.3f mM', rmse_test);
    msgs{end+1} = sprintf('  R^2   = %.3f', R2_test);

    %% PART H: PER-EXPERIMENT TEST METRICS
    testMetrics = table();
    for e = 1:nExp
        idxE = (expID_test == e);
        if any(idxE)
            rmseE = sqrt(mean((Y_test(idxE) - Yhat_test(idxE)).^2, 'omitnan'));
            R2E = calcR2_local(Y_test(idxE), Yhat_test(idxE));
            testMetrics = [testMetrics; table(string(expNames{e}), sum(idxE), rmseE, R2E, ...
                'VariableNames', {'Experiment','nTest','RMSEP','R2'})]; %#ok<AGROW>
            msgs{end+1} = sprintf('  %s internal test: n = %d | RMSEP = %.3f mM | R^2 = %.3f', string(expNames{e}), sum(idxE), rmseE, R2E);
        end
    end

    %% PART I: VIP
    VIP = computeVIP_univariate(X_trainval, Y_trainval, XS_train, stats_locked.W, finalLV);

    %% PART J: PLOTS
    if opts.doPlot
        plot_all_exp_70_30_univariate_results( ...
            Y_trainval, Yhat_trainCV, expID_trainval, ...
            Y_test, Yhat_test, expID_test, ...
            rmse_cv, R2_cv, rmse_test, R2_test, ...
            meanRMSECV, finalLV, RMSECV_perFold, bestLV_perFold, ...
            commonShift, XL, VIP, XS_train, targetSugar, expNames);
    end

    %% PART K: OUTPUT
    out = struct();
    out.targetSugar = targetSugar;
    out.sugarIdx = sugarIdx;
    out.sugarNames = sugarNames;
    out.rootDirs = rootDirs;
    out.hplcExcelPath = hplcExcelPath;
    out.hplcSheets = hplcSheets;
    out.expNames = expNames;
    out.commonShift = commonShift;

    out.finalLV = finalLV;
    out.meanRMSECV = meanRMSECV;
    out.bestLV_perFold = bestLV_perFold;
    out.RMSECV_all = RMSECV_all;
    out.RMSECV_perFold = RMSECV_perFold;

    out.rmse_cv = rmse_cv;
    out.R2_cv = R2_cv;
    out.rmse_fit = rmse_fit;
    out.R2_fit = R2_fit;
    out.rmse_test = rmse_test;
    out.R2_test = R2_test;
    out.testMetrics = testMetrics;

    out.BETA_locked = BETA_locked;
    out.PCTVAR_locked = PCTVAR_locked;
    out.MSE_locked = MSE_locked;
    out.stats_locked = stats_locked;
    out.VIP = VIP;
    out.XL = XL;
    out.XS_train = XS_train;
    out.YL = YL;
    out.YS_train = YS_train;

    out.trainFraction = opts.trainFraction;
    out.randomSeed = opts.randomSeed;
    out.nBlinds = nBlinds;
    out.internalCV_folds = opts.internalCV_folds;

    out.X_all = X_all;
    out.Xraw_all = Xraw_all;
    out.Y_all_sugars = Y_all;
    out.Y_all = Y_target;
    out.flaskID_all = G_all;
    out.keys_all = keys_all;
    out.timeNum_all = timeNum_all;
    out.timeLabel_all = timeLabel_all;
    out.nRep_all = nRep_all;
    out.expID_all = expID_all;
    out.expName_all = expName_all;

    out.idxTrain = idxTrain;
    out.idxTest = idxTest;

    out.X_trainval = X_trainval;
    out.Y_trainval = Y_trainval;
    out.G_trainval = G_trainval;
    out.keys_trainval = keys_trainval;
    out.timeNum_trainval = timeNum_trainval;
    out.timeLabel_trainval = timeLabel_trainval;
    out.nRep_trainval = nRep_trainval;
    out.expID_trainval = expID_trainval;
    out.expName_trainval = expName_trainval;

    out.X_test = X_test;
    out.Y_test = Y_test;
    out.G_test = G_test;
    out.keys_test = keys_test;
    out.timeNum_test = timeNum_test;
    out.timeLabel_test = timeLabel_test;
    out.nRep_test = nRep_test;
    out.expID_test = expID_test;
    out.expName_test = expName_test;

    out.Yhat_trainCV = Yhat_trainCV;
    out.Yhat_train_fit = Yhat_train_fit;
    out.Yhat_test = Yhat_test;

    out.EMSC_model = EMSC_model;
    out.EMSC_train_residuals = EMSC_train_residuals;
    out.EMSC_train_params = EMSC_train_params;
    out.EMSC_test_residuals = EMSC_test_residuals;
    out.EMSC_test_params = EMSC_test_params;
    out.emscPolyOrder = opts.emscPolyOrder;

    out.messages = msgs;

    for k = 1:numel(msgs)
        fprintf('%s\n', msgs{k});
    end
    drawnow;

    if opts.saveMat
        save(sprintf('SF_AllExp_VB_70_30_SinglePLSR_%s_%s_noSNV.mat', targetSugar, string(opts.preprocMethod)), 'out');
    end
end

function R2 = calcR2_local(y, yhat)
    ok = isfinite(y) & isfinite(yhat);
    if sum(ok) < 2
        R2 = NaN;
        return;
    end
    denom = sum((y(ok) - mean(y(ok))).^2);
    if denom > eps
        R2 = 1 - sum((y(ok) - yhat(ok)).^2) / denom;
    else
        R2 = NaN;
    end
end

function VIP = computeVIP_univariate(X, Y, XS, W, A)
    [~, p] = size(X);
    A = min(A, size(XS,2));
    W = W(:,1:A);
    T = XS(:,1:A);

    Wn = W ./ sqrt(sum(W.^2,1));
    SSY = zeros(A,1);

    for a = 1:A
        t = T(:,a);
        denom = (t' * t);
        if denom == 0
            SSY(a) = 0;
        else
            q = (t' * Y) / denom;
            SSY(a) = (q.^2) * denom;
        end
    end

    SSYtot = sum(SSY);
    if SSYtot == 0
        VIP = zeros(p,1);
        return;
    end

    VIP = sqrt(p * ((Wn.^2) * SSY) / SSYtot);
    VIP = VIP(:);
end

function plot_all_exp_70_30_univariate_results( ...
    Y_trainval, Yhat_trainCV, expID_trainval, ...
    Y_test, Yhat_test, expID_test, ...
    rmse_cv, R2_cv, rmse_test, R2_test, ...
    meanRMSECV, finalLV, RMSECV_perFold, bestLV_perFold, ...
    ramanShift, XL, VIP, XS_train, targetSugar, expNames)

    figure('Name', sprintf('All experiments pooled 70/30 - %s', targetSugar), ...
           'Units','normalized','Position',[0.05 0.07 0.88 0.80]);
    tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    scatter(Y_trainval, Yhat_trainCV, 32, expID_trainval, 'filled'); hold on;
    mn = min([Y_trainval; Yhat_trainCV]);
    mx = max([Y_trainval; Yhat_trainCV]);
    plot([mn mx],[mn mx],'k--','LineWidth',1);
    axis equal; grid on;
    xlabel('Measured (mM)');
    ylabel('CV predicted (mM)');
    title(sprintf('%s CV | RMSECV = %.3f | R^2 = %.3f', targetSugar, rmse_cv, R2_cv));
    cb = colorbar; cb.Label.String = 'Experiment ID';

    nexttile;
    scatter(Y_test, Yhat_test, 32, expID_test, 'filled'); hold on;
    mn = min([Y_test; Yhat_test]);
    mx = max([Y_test; Yhat_test]);
    plot([mn mx],[mn mx],'k--','LineWidth',1);
    axis equal; grid on;
    xlabel('Measured (mM)');
    ylabel('Predicted (mM)');
    title(sprintf('%s internal test | RMSEP = %.3f mM | R^2 = %.3f', targetSugar, rmse_test, R2_test));
    cb = colorbar; cb.Label.String = 'Experiment ID';

    nexttile;
    yyaxis left;
    bar(RMSECV_perFold);
    ylabel('RMSECV per blind');
    yyaxis right;
    plot(bestLV_perFold, '-ok', 'LineWidth', 1.4, 'MarkerFaceColor','w');
    ylabel('Best LV');
    xlabel('Blind');
    title('RMSECV and selected LV per blind');
    grid on;

    nexttile;
    plot(1:numel(meanRMSECV), meanRMSECV, '-o', 'LineWidth', 1.4); hold on;
    xline(finalLV, 'k--', sprintf('finalLV=%d', finalLV));
    xlabel('LV');
    ylabel('Mean RMSECV (mM)');
    title('Mean RMSECV vs LV');
    grid on;

    nexttile;
    plot(ramanShift, XL(:,1), 'k-', 'LineWidth', 1.2); hold on;
    yyaxis right;
    plot(ramanShift, VIP, 'r--', 'LineWidth', 1.0);
    xlabel('Raman shift');
    ylabel('VIP');
    yyaxis left;
    ylabel('Loading 1');
    title('Loading 1 and VIP');
    grid on;

    nexttile;
    if size(XS_train,2) >= 2
        scatter(XS_train(:,1), XS_train(:,2), 35, expID_trainval, 'filled');
        xlabel('Score 1');
        ylabel('Score 2');
        title('Train/validation scores, colored by experiment');
        cb = colorbar; cb.Label.String = 'Experiment ID';
        grid on;
    else
        text(0.2,0.5,'Not enough score components'); axis off;
    end

    nexttile;
    scatter(Yhat_test, Yhat_test - Y_test, 24, expID_test, 'filled'); hold on;
    yline(0,'k--');
    xlabel('Predicted');
    ylabel('Residual');
    title('Test residual vs predicted');
    cb = colorbar; cb.Label.String = 'Experiment ID';
    grid on;

    nexttile;
    scatter(Y_test, Yhat_test - Y_test, 24, expID_test, 'filled'); hold on;
    yline(0,'k--');
    xlabel('Measured');
    ylabel('Residual');
    title('Test residual vs measured');
    cb = colorbar; cb.Label.String = 'Experiment ID';
    grid on;

    nexttile;
    nExp = numel(expNames);
    countsTrain = zeros(nExp,1);
    countsTest = zeros(nExp,1);
    for e = 1:nExp
        countsTrain(e) = sum(expID_trainval == e);
        countsTest(e) = sum(expID_test == e);
    end
    bar([countsTrain countsTest]);
    xticks(1:nExp);
    xticklabels(string(expNames));
    xtickangle(30);
    ylabel('Number of SF/T samples');
    legend({'Train/validation','Test'}, 'Location','best');
    title('Split composition');
    grid on;
end

function [XmeanProc, Y, G, keysOut, timeNum, timeLabel, commonShift, nRep, XmeanRaw] = build_mean_dataset_from_root_local(rootDir, hplcMap, opts)

    [files, keys, flaskID, timeNum_all, timeLabel_all] = collect_shakeflask_files_all_local(rootDir, opts.power_mw, opts.t_ms);

    assert(~isempty(files), 'No Raman files found in %s', rootDir);

    keysS = string(keys);
    [keysUniqueS, ~, ic] = unique(keysS, 'stable');

    XmeanRaw  = [];
    XmeanProc = [];
    Y         = [];
    G         = [];
    keysOut   = {};
    timeNum   = [];
    timeLabel = strings(0,1);
    nRep      = [];
    commonShift = [];

    kOut = 0;

    for i = 1:numel(keysUniqueS)
        idx = (ic == i);
        key_i = char(keysUniqueS(i));

        if ~isempty(hplcMap)
            if ~isKey(hplcMap, key_i)
                continue;
            end
        end

        files_i = files(idx);
        flask_i = mode(flaskID(idx));
        tnum_i  = mode(timeNum_all(idx));
        tlab_i  = unique(timeLabel_all(idx));
        tlab_i  = tlab_i(1);

        switch lower(string(opts.preprocMethod))
            case "asls"
                [~, shift_i, Xraw_i] = PreprocessingSpectras_AsLS(files_i, opts.shiftRange, false);
            case "msbackadj"
                [~, shift_i, Xraw_i] = PreprocessingSpectras2(files_i, opts.shiftRange, false);
            case "peaknorm"
                [~, shift_i, Xraw_i] = PreprocessingSpectrasNew(files_i, opts.shiftRange, false);
            case "emsc"
                [~, shift_i, Xraw_i] = PreprocessingSpectras_AsLS(files_i, opts.shiftRange, false);
            otherwise
                error('Unknown preprocessing method: %s', string(opts.preprocMethod));
        end

        shift_i = shift_i(:).';

        if isempty(commonShift)
            commonShift = shift_i;
        else
            if numel(shift_i) ~= numel(commonShift) || any(abs(shift_i - commonShift) > 1e-10)
                error('Raman shift grid mismatch between keys while building mean spectra.');
            end
        end

        xmean_raw_i = mean(Xraw_i, 1);

        switch lower(string(opts.preprocMethod))
            case "asls"
                xmean_proc_i = PreprocessingSpectras_AsLS(xmean_raw_i, commonShift, opts.snvOn);
            case "msbackadj"
                xmean_proc_i = PreprocessingSpectras2(xmean_raw_i, commonShift, opts.snvOn);
            case "peaknorm"
                xmean_proc_i = PreprocessingSpectrasNew(xmean_raw_i, commonShift, opts.snvOn);
            case "emsc"
                xmean_proc_i = xmean_raw_i;
            otherwise
                error('Unknown preprocessing method: %s', string(opts.preprocMethod));
        end

        kOut = kOut + 1;
        XmeanRaw(kOut,:)  = xmean_raw_i; %#ok<AGROW>
        XmeanProc(kOut,:) = xmean_proc_i; %#ok<AGROW>
        G(kOut,1)         = flask_i; %#ok<AGROW>
        keysOut{kOut,1}   = key_i; %#ok<AGROW>
        timeNum(kOut,1)   = tnum_i; %#ok<AGROW>
        timeLabel(kOut,1) = tlab_i; %#ok<AGROW>
        nRep(kOut,1)      = numel(files_i); %#ok<AGROW>

        if ~isempty(hplcMap)
            Y(kOut,:) = hplcMap(key_i); %#ok<AGROW>
        end
    end
end

function [files, keys, flaskID, timeNum, timeLabel] = collect_shakeflask_files_all_local(rootDir, power_mw, t_ms)

    d = dir(fullfile(rootDir, '**', '*.spc'));

    files = {};
    keys = {};
    flaskID = [];
    timeNum = [];
    timeLabel = strings(0,1);

    for i = 1:numel(d)
        fp = fullfile(d(i).folder, d(i).name);

        hasPower = true;
        hasTime  = true;

        if ~isempty(power_mw)
            hasPower = contains(d(i).name, sprintf('%dmw', power_mw), 'IgnoreCase', true);
        end
        if ~isempty(t_ms)
            hasTime = contains(d(i).name, sprintf('%dms', t_ms), 'IgnoreCase', true);
        end

        if ~(hasPower && hasTime)
            continue;
        end

        parts = strsplit(d(i).folder, filesep);
        if numel(parts) < 2
            continue;
        end

        flaskPart = string(parts{end});
        timePart  = string(parts{end-1});

        tokF = regexp(char(flaskPart), 'shakeflask(\d+)', 'tokens', 'once', 'ignorecase');
        tokT = regexp(char(timePart),  'T(\d+)',         'tokens', 'once', 'ignorecase');

        if isempty(tokF) || isempty(tokT)
            continue;
        end

        fnum = str2double(tokF{1});
        tnum = str2double(tokT{1});

        files{end+1,1} = fp; %#ok<AGROW>
        keys{end+1,1} = sprintf('SF%d_T%d', fnum, tnum); %#ok<AGROW>
        flaskID(end+1,1) = fnum; %#ok<AGROW>
        timeNum(end+1,1) = tnum; %#ok<AGROW>
        timeLabel(end+1,1) = "T" + string(tnum); %#ok<AGROW>
    end
end

function hplcMap = buildHPLCMap_Shakeflask_Exact(excelPath, sheetName, naMeansZero)

    raw = readcell(excelPath, 'Sheet', sheetName);
    hplcMap = containers.Map('KeyType','char','ValueType','any');

    nFoundLabels = 0;

    for i = 1:size(raw,1)

        labelCol = [];
        tok = [];

        for j = 1:size(raw,2)
            v = raw{i,j};
            if ischar(v) || isstring(v)
                s = strtrim(string(v));
                s = regexprep(s, '\s+', ' ');
                tok_here = regexp(char(s), '^T\s*(\d+)\s+Flask\s*(\d+)$', 'tokens', 'once', 'ignorecase');
                if ~isempty(tok_here)
                    labelCol = j;
                    tok = tok_here;
                    break;
                end
            end
        end

        if isempty(tok)
            continue;
        end

        nFoundLabels = nFoundLabels + 1;

        tnum = str2double(tok{1});
        fnum = str2double(tok{2});

        glu = parseHPLCValue_local(raw{i,labelCol+1}, naMeansZero);
        xyl = parseHPLCValue_local(raw{i,labelCol+2}, naMeansZero);
        mtl = parseHPLCValue_local(raw{i,labelCol+3}, naMeansZero);

        key = sprintf('SF%d_T%d', fnum, tnum);
        hplcMap(key) = [glu, xyl, mtl];
    end

    fprintf('Found %d HPLC label rows.\n', nFoundLabels);
    fprintf('Built HPLC map with %d entries from sheet %s.\n', hplcMap.Count, sheetName);
end

function v = parseHPLCValue_local(x, naMeansZero)
    if isnumeric(x)
        if isnan(x)
            v = NaN;
        else
            v = x;
        end
        return;
    end

    if isstring(x) || ischar(x)
        s = strtrim(string(x));
        if s == ""
            v = NaN;
        elseif strcmpi(s,'n.a.') || strcmpi(s,'n.a') || strcmpi(s,'na')
            if naMeansZero
                v = 0;
            else
                v = NaN;
            end
        else
            v = str2double(s);
        end
        return;
    end

    v = NaN;
end

function foldID = make_venetian_blinds_folds_simple(nSamples, nBlinds)
    nBlinds = min(nBlinds, nSamples);
    foldID = mod((1:nSamples)' - 1, nBlinds) + 1;
end

function [Xcorr, EMSC_model, residuals, params] = fit_emsc_local(X, shift, sampleNames, polyOrder)
    if nargin < 4 || isempty(polyOrder)
        polyOrder = 1;
    end
    Z = make_saisir_local(X, shift, sampleNames);
    EMSC_model = make_emsc_modfunc(Z, polyOrder);
    [Zcorr, Zres, params] = cal_emsc(Z, EMSC_model);
    Xcorr = Zcorr.d;
    residuals = Zres.d;
end

function [Xcorr, residuals, params] = apply_emsc_model_local(X, shift, sampleNames, EMSC_model)
    Z = make_saisir_local(X, shift, sampleNames);
    [Zcorr, Zres, params] = cal_emsc(Z, EMSC_model);
    Xcorr = Zcorr.d;
    residuals = Zres.d;
end

function Z = make_saisir_local(X, shift, sampleNames)
    Z = [];
    Z.d = X;
    Z.v = num2str(shift(:));

    if nargin < 3 || isempty(sampleNames)
        sampleNames = "sample" + string((1:size(X,1))');
    end

    sampleNames = string(sampleNames(:));
    Z.i = char(pad(sampleNames));
end
