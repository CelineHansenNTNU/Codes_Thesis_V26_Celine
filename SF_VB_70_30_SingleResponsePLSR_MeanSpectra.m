function out = SF_VB_70_30_SingleResponsePLSR_MeanSpectra(trainRootDir, hplcExcelPath, predictRootDir, opts)
% SF_VB_70_30_SingleResponsePLSR_MeanSpectra
%
% Single-response PLSR for shake flask mean spectra:
%   - one mean RAW spectrum per shakeflask per timepoint (SF/T)
%   - preprocess mean spectrum
%   - 70/30 split on SF/T samples
%   - Venetian blinds CV on 70% train/validation
%   - fit final model on all 70%
%   - predict held-out 30%
%   - optional application to another experiment/root
%
% Choose one sugar with:
%   opts.targetSugar = "glucose" or "xylose" or "mannitol"
%
% Optional:
%   opts.wantedPredictFlask = 1   % only plot/predict one flask in predictRootDir
%
% Celine Hansen + adapted from multivariate version

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts,'shiftRange'),         opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'power_mw'),           opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),               opts.t_ms = 1000; end
    if ~isfield(opts,'snvOn'),              opts.snvOn = false; end
    if ~isfield(opts,'maxLV_wish'),         opts.maxLV_wish = 10; end
    if ~isfield(opts,'lvRule'),             opts.lvRule = "min"; end
    if ~isfield(opts,'tolFactor'),          opts.tolFactor = 1.05; end
    if ~isfield(opts,'hplcSheet'),          opts.hplcSheet = 'EXP2'; end
    if ~isfield(opts,'naMeansZero'),        opts.naMeansZero = true; end
    if ~isfield(opts,'saveMat'),            opts.saveMat = true; end
    if ~isfield(opts,'doPlot'),             opts.doPlot = true; end
    if ~isfield(opts,'trainFraction'),      opts.trainFraction = 0.70; end
    if ~isfield(opts,'randomSeed'),         opts.randomSeed = 1; end
    if ~isfield(opts,'nBlinds'),            opts.nBlinds = 10; end
    if ~isfield(opts,'internalCV_folds'),   opts.internalCV_folds = 10; end
    if ~isfield(opts,'preprocMethod'),      opts.preprocMethod = "asls"; end
    if ~isfield(opts,'emscPolyOrder'),      opts.emscPolyOrder = 1; end % 1 = basic EMSC
    if ~isfield(opts,'wantedPredictFlask'), opts.wantedPredictFlask = []; end
    if ~isfield(opts,'targetSugar'),        opts.targetSugar = "glucose"; end
    if ~isfield(opts, 'predictHPLCSheet'),  opts.predictHPLCSheet = ""; end
    if ~isfield(opts, 'evaluateExternalRMSEP'), opts.evaluateExternalRMSEP = false; end

    targetSugar = lower(string(opts.targetSugar));
    sugarNames = ["glucose","xylose","mannitol"];

    if ~ismember(targetSugar, sugarNames)
        error('opts.targetSugar must be "glucose", "xylose", or "mannitol".');
    end

    sugarIdx = find(sugarNames == targetSugar, 1);

    msgs = {};
    msgs{end+1} = sprintf('\n=== SF_VB_70_30_UnivariatePLSR_MeanSpectra ===');
    msgs{end+1} = sprintf('Train Raman root:   %s', trainRootDir);
    msgs{end+1} = sprintf('HPLC Excel:         %s', hplcExcelPath);
    msgs{end+1} = sprintf('Predict Raman root: %s', string(predictRootDir));
    msgs{end+1} = sprintf('Target sugar:       %s', targetSugar);
    msgs{end+1} = sprintf('Preprocessing:      %s | SNV = %d', string(opts.preprocMethod), opts.snvOn);
    msgs{end+1} = sprintf('Train fraction:     %.2f | Test fraction: %.2f', opts.trainFraction, 1-opts.trainFraction);

    %% PART A: BUILD TRAIN DATASET
    hplcMap = buildHPLCMap_Shakeflask_Exact(hplcExcelPath, opts.hplcSheet, opts.naMeansZero);

    [X_exp, Y_all, G_exp, keys_exp, timeNum_exp, timeLabel_exp, commonShift, nRep_exp, XrawMean_exp] = ...
        build_mean_dataset_from_root_local(trainRootDir, hplcMap, opts);

    assert(~isempty(X_exp), 'No train mean spectra could be built.');

    % Pick one sugar only
    Y_exp = Y_all(:, sugarIdx);

    msgs{end+1} = sprintf('Built train mean-spectra dataset with %d SF/T samples across %d flasks.', ...
        size(X_exp,1), numel(unique(G_exp)));

    %% PART B: 70/30 SPLIT
    nSamples = size(X_exp,1);
    rng(opts.randomSeed);

    perm = randperm(nSamples);
    nTrain = max(2, round(opts.trainFraction * nSamples));
    nTrain = min(nTrain, nSamples-1);

    idxTrain = false(nSamples,1);
    idxTrain(perm(1:nTrain)) = true;
    idxTest = ~idxTrain;

    X_trainval = X_exp(idxTrain,:);
    Y_trainval = Y_exp(idxTrain,:);
    G_trainval = G_exp(idxTrain);
    keys_trainval = keys_exp(idxTrain);
    timeNum_trainval = timeNum_exp(idxTrain);
    timeLabel_trainval = timeLabel_exp(idxTrain);
    nRep_trainval = nRep_exp(idxTrain);

    X_test = X_exp(idxTest,:);
    Y_test = Y_exp(idxTest,:);
    G_test = G_exp(idxTest);
    keys_test = keys_exp(idxTest);
    timeNum_test = timeNum_exp(idxTest);
    timeLabel_test = timeLabel_exp(idxTest);
    nRep_test = nRep_exp(idxTest);

    msgs{end+1} = sprintf('Train/val SF/T samples: %d | Test SF/T samples: %d', ...
        size(X_trainval,1), size(X_test,1));

    % Optional EMSC preprocessing after split
    EMSC_model = [];
    if lower(string(opts.preprocMethod)) == "emsc"
        msgs{end+1} = sprintf('Applying EMSC preprocessing fitted on train/val set only.');
        % Fit EMSC on train/val spectra
        [X_trainval, EMSC_model, EMSC_train_residuals, EMSC_train_params] = fit_emsc_local(X_trainval, commonShift, keys_trainval, opts.emscPolyOrder);
        % apply same EMSC model to held-out spectra
        [X_test, EMSC_test_residuals, EMSC_test_params] = apply_emsc_model_local(X_test, commonShift, keys_test, EMSC_model);
    else
        EMSC_train_residuals = [];
        EMSC_train_params = [];
        EMSC_test_residuals = [];
        EMSC_test_params = [];
    end

    %% PART C: VENETIAN BLINDS ON TRAIN/VAL
    nBlinds = min(opts.nBlinds, size(X_trainval,1));
    foldID_train = make_venetian_blinds_folds_simple(size(X_trainval,1), nBlinds);

    Yhat_trainCV = nan(size(Y_trainval));
    bestLV_perFold = nan(nBlinds,1);
    RMSECV_all = nan(opts.maxLV_wish, nBlinds);
    RMSECV_perFold = nan(nBlinds,1);

    msgs{end+1} = sprintf('Venetian blinds on train/validation set: %d blinds | n=%d mean spectra', ...
        nBlinds, size(X_trainval,1));

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

        Yhat_trainCV(valIdx) = Ypred_val;

        rmse_fold = sqrt(mean((Y_val - Ypred_val).^2, 'omitnan'));
        RMSECV_perFold(b) = rmse_fold;

        msgs{end+1} = sprintf('Blind %d/%d | nTrain=%d nVal=%d | bestLV=%d | RMSECV=%.3f mM', ...
            b, nBlinds, size(X_train,1), size(X_val,1), bestLV, rmse_fold);
    end

    %% PART D: OVERALL INTERNAL CV METRICS
    rmse_cv = sqrt(mean((Y_trainval - Yhat_trainCV).^2, 'omitnan'));

    ok = isfinite(Y_trainval) & isfinite(Yhat_trainCV);
    denom = sum((Y_trainval(ok) - mean(Y_trainval(ok))).^2);
    if denom > eps
        R2_cv = 1 - sum((Y_trainval(ok) - Yhat_trainCV(ok)).^2) / denom;
    else
        R2_cv = NaN;
    end

    msgs{end+1} = sprintf('Internal Venetian blinds CV on 70%% set:');
    msgs{end+1} = sprintf('  RMSECV = %.3f mM', rmse_cv);
    msgs{end+1} = sprintf('  R^2    = %.3f', R2_cv);

    %% PART E: CHOOSE FINAL LV
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

    %% PART F: FIT FINAL MODEL ON 70%, PREDICT HELD-OUT 30%
    [XL, YL, XS_train, YS_train, BETA_locked, PCTVAR_locked, MSE_locked, stats_locked] = ...
        plsregress(X_trainval, Y_trainval, finalLV); %#ok<ASGLU>

    Yhat_train_fit = [ones(size(X_trainval,1),1) X_trainval] * BETA_locked;
    Yhat_test = [ones(size(X_test,1),1) X_test] * BETA_locked;
    Yhat_train_fit = max(Yhat_train_fit, 0);
    Yhat_test = max(Yhat_test, 0);

    rmse_fit = sqrt(mean((Y_trainval - Yhat_train_fit).^2, 'omitnan'));
    rmse_test = sqrt(mean((Y_test - Yhat_test).^2, 'omitnan'));

    ok = isfinite(Y_test) & isfinite(Yhat_test);
    denom = sum((Y_test(ok) - mean(Y_test(ok))).^2);
    if denom > eps
        R2_test = 1 - sum((Y_test(ok) - Yhat_test(ok)).^2) / denom;
    else
        R2_test = NaN;
    end

    msgs{end+1} = sprintf('Prediction on held-out 30%% set:');
    msgs{end+1} = sprintf('  RMSEP = %.3f mM', rmse_test);
    msgs{end+1} = sprintf('  R^2   = %.3f', R2_test);

    %% PART G: TEST SET GROUP-MEAN VIEW
    % Already one row per SF/T
    mu_true_test = Y_test;
    mu_pred_test = Yhat_test;
    sd_pred_test = nan(size(Yhat_test));

    RMSEP_test_group = sqrt(mean((mu_true_test - mu_pred_test).^2, 'omitnan'));
    %msgs{end+1} = sprintf('  RMSEP on prediction SF/T means = %.3f mM', RMSEP_test_group);

    %% PART H: OPTIONAL APPLY TO ANOTHER EXPERIMENT, ONE CHOSEN FLASK + Optinal external RMSEP
    X_exp_rep = [];
    XrawMean_exp_rep = [];
    Yhat_exp_rep = [];
    Y_exp_rep = [];
    keys_exp_rep = {};
    G_exp_rep = [];
    timeNum_exp_rep = [];
    timeLabel_exp_rep = strings(0,1);
    nRep_exp_rep = [];
    predTableExp_rep = table();

    if nargin >= 3 && ~isempty(predictRootDir)
        % If external HPLC should be used, build HPLC map for predict set too  
        hplcMap_predict = [];
        if opts.evaluateExternalRMSEP
            if strlength(string(opts.predictHPLCSheet)) == 0
                error('opts.predictHPLCSheet must be provided when opts.evaluateExternalRMSEP = true.');
            end
            hplcMap_predict = buildHPLCMap_Shakeflask_Exact(hplcExcelPath, opts.predictHPLCSheet, opts.naMeansZero);
        end

        [X_exp_rep, Yrep_all, G_exp_rep, keys_exp_rep, timeNum_exp_rep, timeLabel_exp_rep, shift_exp_rep, nRep_exp_rep, XrawMean_exp_rep] = build_mean_dataset_from_root_local(predictRootDir, hplcMap_predict, opts);

        assert(~isempty(X_exp_rep), 'No replicate mean spectra could be built.');

        shift_exp_rep = shift_exp_rep(:).';
        if numel(shift_exp_rep) ~= numel(commonShift) || any(abs(shift_exp_rep - commonShift) > 1e-10)
            error('Predict Raman shift grid does not match training grid.');
        end

        % If HPLC was included, keep only chosen sugar
        if opts.evaluateExternalRMSEP
            Y_exp_rep = Yrep_all(:, sugarIdx);
        else
            Y_exp_rep = [];
        end

        if lower(string(opts.preprocMethod)) == "emsc"
            [X_exp_rep, EMSC_pred_residuals, EMSC_pred_params] = apply_emsc_model_local(X_exp_rep, commonShift, keys_exp_rep, EMSC_model);
        else
            EMSC_pred_residuals = [];
            EMSC_pred_params = [];
        end

        Yhat_exp_rep = [ones(size(X_exp_rep, 1), 1) X_exp_rep] * BETA_locked;
        Yhat_exp_rep = max(Yhat_exp_rep, 0);

        % External RMSEP / R2
        if opts.evaluateExternalRMSEP
            rmse_external = sqrt(mean((Y_exp_rep -  Yhat_exp_rep).^2, 'omitnan'));

            ok = isfinite(Y_exp_rep) & isfinite(Yhat_exp_rep);
            denom = sum((Y_exp_rep(ok) - mean(Y_exp_rep(ok))).^2);
            if denom > eps
                R2_external = 1 - sum((Y_exp_rep(ok) - Yhat_exp_rep(ok)).^2) / denom;
            else
                R2_external = NaN;
            end

            RMSEP_external_group = rmse_external;

            msgs{end+1} = sprintf('External prediction on %s:', string(opts.predictHPLCSheet));
            msgs{end+1} = sprintf('  RMSEP_external = %.3f mM', rmse_external);
            msgs{end+1} = sprintf('  R^2_external   = %.3f', R2_external);
        end

        if isempty(opts.wantedPredictFlask)
            flasksToPlot = unique(G_exp_rep);
        else
            flasksToPlot = opts.wantedPredictFlask(:);
        end

        if ~isempty(flasksToPlot)
            thisFlask = flasksToPlot(1);
            idxF = (G_exp_rep == thisFlask);

            if any(idxF)
                [tSorted, ord] = sort(timeNum_exp_rep(idxF));
                tLabel    = timeLabel_exp_rep(idxF); tLabel = tLabel(ord);
                nRepThis  = nRep_exp_rep(idxF);      nRepThis = nRepThis(ord);
                YhatThis  = Yhat_exp_rep(idxF);      YhatThis = YhatThis(ord);

                predTableExp_rep = table(tSorted(:), tLabel(:), nRepThis(:), YhatThis(:), ...
                    'VariableNames', {'timeNum','timeLabel','nRep','predictedConc'});

                if opts.evaluateExternalRMSEP
                    YtrueThis = Y_exp_rep(idxF);
                    YtrueThis = YtrueThis(ord);
                    predTableExp_rep.measuredConc = YtrueThis(:);
                    predTableExp_rep.residual = YhatThis(:) - YtrueThis(:);
                end

                if opts.doPlot
                    figure('Name', sprintf('%s predictions - shakeflask%d', targetSugar, thisFlask), ...
                           'Units','normalized','Position',[0.18 0.18 0.62 0.42]);
                
                    if opts.evaluateExternalRMSEP
                        plot(tSorted, YhatThis, '-o', ...
                            'LineWidth', 1.5, 'MarkerSize', 6, ...
                            'Color', [0 0.4470 0.7410], ...
                            'DisplayName','Predicted'); hold on;
                
                        plot(tSorted, YtrueThis, '-s', ...
                            'LineWidth', 1.2, 'MarkerSize', 6, ...
                            'Color', [0.8500 0.3250 0.0980], ...
                            'DisplayName','Measured');
                
                        % Blue predicted values
                        for ii = 1:numel(tSorted)
                            text(tSorted(ii), YhatThis(ii), sprintf(' %.2f', YhatThis(ii)), ...
                                'Color', [0 0.4470 0.7410], ...
                                'FontSize', 8, ...
                                'VerticalAlignment','bottom', ...
                                'HorizontalAlignment','left');
                        end
                
                        % Red measured values
                        for ii = 1:numel(tSorted)
                            text(tSorted(ii), YtrueThis(ii), sprintf(' %.2f', YtrueThis(ii)), ...
                                'Color', [0.8500 0.3250 0.0980], ...
                                'FontSize', 8, ...
                                'VerticalAlignment','top', ...
                                'HorizontalAlignment','left');
                        end
                
                        legend('Location','best');
                    else
                        plot(tSorted, YhatThis, '-o', ...
                            'LineWidth', 1.5, 'MarkerSize', 6, ...
                            'Color', [0 0.4470 0.7410]); hold on;
                
                        for ii = 1:numel(tSorted)
                            text(tSorted(ii), YhatThis(ii), sprintf(' %.2f', YhatThis(ii)), ...
                                'Color', [0 0.4470 0.7410], ...
                                'FontSize', 8, ...
                                'VerticalAlignment','bottom', ...
                                'HorizontalAlignment','left');
                        end
                    end
                
                    grid on;
                    xlabel('Time (h)');
                    ylabel(sprintf('%s concentration (mM)', targetSugar));
                    title(sprintf('Predicted %s vs time - shakeflask%d', targetSugar, thisFlask));
                    xticks(tSorted);
                    xticklabels(tLabel);
                    xtickangle(45);
                end
            else
                warning('No mean spectra found for shakeflask%d in predictRootDir.', thisFlask);
            end
        end
    end

    %% PART I: VIP
    VIP = computeVIP_univariate(X_trainval, Y_trainval, XS_train, stats_locked.W, finalLV);

    %% PART J: PLOTS
    if opts.doPlot
        plot_shakeflask_70_30_univariate_results( ...
            Y_trainval, Yhat_trainCV, ...
            Y_test, Yhat_test, ...
            rmse_cv, R2_cv, ...
            rmse_test, R2_test, ...
            meanRMSECV, finalLV, ...
            RMSECV_perFold, bestLV_perFold, ...
            commonShift, XL, VIP, XS_train, ...
            mu_true_test, mu_pred_test, RMSEP_test_group, targetSugar);

        if opts.doPlot && ~isempty(Y_exp_rep) && opts.evaluateExternalRMSEP
            figure('Name', sprintf('External test - %s', targetSugar), ...
                   'Units','normalized','Position',[0.15 0.15 0.70 0.45]);
            tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

            nexttile;
            plot(Y_exp_rep, Yhat_exp_rep, '.', 'MarkerSize', 12); hold on;
            mn = min([Y_exp_rep; Yhat_exp_rep]);
            mx = max([Y_exp_rep; Yhat_exp_rep]);
            plot([mn mx],[mn mx],'k--','LineWidth',1);
            axis equal; grid on;
            xlabel('Measured (mM)');
            ylabel('Predicted (mM)');
            title(sprintf('External | RMSEP=%.3f mM | R^2=%.3f', rmse_external, R2_external));

            nexttile;
            scatter(Yhat_exp_rep, Yhat_exp_rep - Y_exp_rep, 20, 'filled'); hold on;
            yline(0,'k--');
            xlabel('Predicted');
            ylabel('Residual');
            title('Residual vs predicted');
            grid on;

            nexttile;
            scatter(Y_exp_rep, Yhat_exp_rep - Y_exp_rep, 20, 'filled'); hold on;
            yline(0,'k--');
            xlabel('Measured');
            ylabel('Residual');
            title('Residual vs measured');
            grid on;
        end
    end

    %% PART K: OUTPUT
    out = struct();

    out.targetSugar = targetSugar;
    out.sugarIdx = sugarIdx;

    out.commonShift = commonShift;
    out.finalLV = finalLV;
    out.meanRMSECV = meanRMSECV;
    out.bestLV_perFold = bestLV_perFold;
    out.RMSECV_all = RMSECV_all;
    out.RMSECV_perFold = RMSECV_perFold;

    out.rmse_cv = rmse_cv;
    out.R2_cv = R2_cv;
    out.rmse_fit = rmse_fit;
    out.rmse_test = rmse_test;
    out.R2_test = R2_test;
    out.RMSEP_test_group = RMSEP_test_group;

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

    out.X_exp = X_exp;
    out.XrawMean_exp = XrawMean_exp;
    out.Y_exp = Y_exp;
    out.keys_exp = keys_exp;
    out.flaskID_exp = G_exp;
    out.timeNum_exp = timeNum_exp;
    out.timeLabel_exp = timeLabel_exp;
    out.nRep_exp = nRep_exp;

    out.idxTrain = idxTrain;
    out.idxTest = idxTest;

    out.X_trainval = X_trainval;
    out.Y_trainval = Y_trainval;
    out.G_trainval = G_trainval;
    out.keys_trainval = keys_trainval;
    out.timeNum_trainval = timeNum_trainval;
    out.timeLabel_trainval = timeLabel_trainval;
    out.nRep_trainval = nRep_trainval;

    out.X_test = X_test;
    out.Y_test = Y_test;
    out.G_test = G_test;
    out.keys_test = keys_test;
    out.timeNum_test = timeNum_test;
    out.timeLabel_test = timeLabel_test;
    out.nRep_test = nRep_test;

    out.Yhat_trainCV = Yhat_trainCV;
    out.Yhat_train_fit = Yhat_train_fit;
    out.Yhat_test = Yhat_test;

    out.mu_true_test = mu_true_test;
    out.mu_pred_test = mu_pred_test;
    out.sd_pred_test = sd_pred_test;

    out.X_predict = X_exp_rep;
    out.XrawMean_predict = XrawMean_exp_rep;
    out.Yhat_predict = Yhat_exp_rep;
    out.keys_predict = keys_exp_rep;
    out.flaskID_predict = G_exp_rep;
    out.timeNum_predict = timeNum_exp_rep;
    out.timeLabel_predict = timeLabel_exp_rep;
    out.nRep_predict = nRep_exp_rep;
    out.predTablePredict = predTableExp_rep;

    out.Y_external = Y_exp_rep;
    out.Yhat_external = Yhat_exp_rep;
    out.rmse_external = rmse_external;
    out.R2_external = R2_external;
    out.RMSEP_external_group = RMSEP_external_group;
    out.predictHPLCSheet = opts.predictHPLCSheet;
    out.evaluateExternalRMSEP = opts.evaluateExternalRMSEP;

    out.EMSC_model = EMSC_model;
    out.EMSC_train_residuals = EMSC_train_residuals;
    out.EMSC_train_params = EMSC_train_params;
    out.EMSC_test_residuals = EMSC_test_residuals;
    out.EMSC_test_params = EMSC_test_params;
    out.EMSC_pred_residuals = EMSC_pred_residuals;
    out.EMSC_pred_params = EMSC_pred_params;
    out.emscPolyOrder = opts.emscPolyOrder;

    out.messages = msgs;

    for k = 1:numel(msgs)
        fprintf('%s\n', msgs{k});
    end
    drawnow;

    if opts.saveMat
        save(sprintf('SF_VB_70_30_SinglePLSR_SF7_%s_emsc.mat', targetSugar), 'out');
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

function plot_shakeflask_70_30_univariate_results( ...
    Y_trainval, Yhat_trainCV, ...
    Y_test, Yhat_test, ...
    rmse_cv, R2_cv, ...
    rmse_test, R2_test, ...
    meanRMSECV, finalLV, ...
    RMSECV_perFold, bestLV_perFold, ...
    ramanShift, XL, VIP, XS_train, ...
    mu_true_test, mu_pred_test, RMSEP_test_group, targetSugar)

    figure('Name', sprintf('Univariate 70/30 - %s', targetSugar), ...
           'Units','normalized','Position',[0.06 0.08 0.86 0.78]);
    tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(Y_trainval, Yhat_trainCV, '.', 'MarkerSize', 12); hold on;
    mn = min([Y_trainval; Yhat_trainCV]);
    mx = max([Y_trainval; Yhat_trainCV]);
    plot([mn mx],[mn mx],'k--','LineWidth',1);
    axis equal; grid on;
    xlabel('Measured (mM)');
    ylabel('CV predicted (mM)');
    title(sprintf('%s CV | RMSECV=%.3f | R^2=%.3f', targetSugar, rmse_cv, R2_cv));

    nexttile;
    plot(Y_test, Yhat_test, '.', 'MarkerSize', 12); hold on;
    mn = min([Y_test; Yhat_test]);
    mx = max([Y_test; Yhat_test]);
    plot([mn mx],[mn mx],'k--','LineWidth',1);
    axis equal; grid on;
    xlabel('Measured (mM)');
    ylabel('Predicted (mM)');
    title(sprintf('%s test | RMSEP=%.3f mM | R^2=%.3f', targetSugar, rmse_test, R2_test));

    nexttile;
    yyaxis left;
    bar(RMSECV_perFold);
    ylabel('RMSECV per fold');
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
        scatter(XS_train(:,1), XS_train(:,2), 35, Y_trainval, 'filled');
        xlabel('Score 1');
        ylabel('Score 2');
        title('Scores');
        colorbar;
        grid on;
    else
        text(0.2,0.5,'Not enough score components'); axis off;
    end

    nexttile;
    scatter(Yhat_test, Yhat_test - Y_test, 20, 'filled'); hold on;
    yline(0,'k--');
    xlabel('Predicted');
    ylabel('Residual');
    title('Residual vs predicted');
    grid on;

    nexttile;
    scatter(Y_test, Yhat_test - Y_test, 20, 'filled'); hold on;
    yline(0,'k--');
    xlabel('Measured');
    ylabel('Residual');
    title('Residual vs measured');
    grid on;

    % nexttile;
    % plot(mu_true_test, mu_pred_test, 'o', 'LineWidth', 1.2); hold on;
    % mn = min([mu_true_test; mu_pred_test]);
    % mx = max([mu_true_test; mu_pred_test]);
    % plot([mn mx],[mn mx],'k--','LineWidth',1);
    % axis equal; grid on;
    % xlabel('Measured SF/T mean');
    % ylabel('Predicted SF/T mean');
    % title(sprintf('Test means | RMSEP=%.3f', RMSEP_test_group));
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

        % Read raw replicate spectra, only raw part needed here
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

        % Mean over RAW replicates
        xmean_raw_i = mean(Xraw_i, 1);

        % Preprocess mean raw spectrum
        switch lower(string(opts.preprocMethod))
            case "asls"
                xmean_proc_i = PreprocessingSpectras_AsLS(xmean_raw_i, commonShift, opts.snvOn);
            case "msbackadj"
                xmean_proc_i = PreprocessingSpectras2(xmean_raw_i, commonShift, opts.snvOn);
            case "peaknorm"
                xmean_proc_i = PreprocessingSpectrasNew(xmean_raw_i, commonShift, opts.snvOn); % For this, snvOn is really doMeanAfter
            case "emsc"
                xmean_proc_i = xmean_raw_i;
            otherwise
                error('Unknown preprocessing method: %s', string(opts.preprocMethod));
        end

        kOut = kOut + 1;
        XmeanRaw(kOut,:)  = xmean_raw_i;    %#ok<AGROW>
        XmeanProc(kOut,:) = xmean_proc_i;   %#ok<AGROW>
        G(kOut,1)         = flask_i;        %#ok<AGROW>
        keysOut{kOut,1}   = key_i;          %#ok<AGROW>
        timeNum(kOut,1)   = tnum_i;         %#ok<AGROW>
        timeLabel(kOut,1) = tlab_i;         %#ok<AGROW>
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