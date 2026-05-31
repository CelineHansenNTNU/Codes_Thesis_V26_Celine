function out = SF_VB_70_30_MultiSugarPLSR_MeanSpec(trainRootDir, hplcExcelPath, predictRootDir, opts)
% Shake flask multivariate PLSR using:
%   - one mean RAW spectrum per shakeflask per timepoint (SF/T)
%   - preprocessing of that mean spectrum
%   - 70/30 split on SF/T groups
%   - Venetian blinds CV on the 70% train/validation part only
%   - final model fitted on all 70%
%   - prediction on held-out 30%
%
% Celine Hansen + merged from LOLO + Venetian blinds setup

    if nargin < 4
        opts = struct();
    end

    if ~isfield(opts,'shiftRange'),      opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'power_mw'),        opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),            opts.t_ms = 1000; end
    if ~isfield(opts,'snvOn'),           opts.snvOn = false; end
    if ~isfield(opts,'maxLV_wish'),      opts.maxLV_wish = 10; end
    if ~isfield(opts,'lvRule'),          opts.lvRule = "min"; end
    if ~isfield(opts,'tolFactor'),       opts.tolFactor = 1.05; end
    if ~isfield(opts,'hplcSheet'),       opts.hplcSheet = 'EXP2'; end
    if ~isfield(opts,'naMeansZero'),     opts.naMeansZero = true; end
    if ~isfield(opts,'saveMat'),         opts.saveMat = true; end
    if ~isfield(opts,'doPlot'),          opts.doPlot = true; end
    if ~isfield(opts,'trainFraction'),   opts.trainFraction = 0.70; end
    if ~isfield(opts,'randomSeed'),      opts.randomSeed = 1; end
    if ~isfield(opts,'nBlinds'),         opts.nBlinds = 10; end
    if ~isfield(opts,'internalCV_folds'),opts.internalCV_folds = 10; end
    if ~isfield(opts,'preprocMethod'),   opts.preprocMethod = "asls"; end
    if ~isfield(opts,'wantedExp1Flask'), opts.wantedExp1Flask = []; end
    if ~isfield(opts,'predictHPLCSheet'), opts.predictHPLCSheet = ""; end
    if ~isfield(opts,'evaluateExternalRMSEP'), opts.evaluateExternalRMSEP = false; end
    if ~isfield(opts,'plotSugarsPredict'), opts.plotSugarsPredict = [1 2 3]; end % 1 = glucose, 2 = xylose, 3 = mannitol

    msgs = {};

    msgs{end+1} = sprintf('\n=== SF_HPLC_VenetianBlinds_70_30_MultiSugarPLSR_MeanSpectra ===\n');
    if iscell(trainRootDir)
        msgs{end+1} = sprintf('Training Raman roots:\n %s\n', strjoin(string(trainRootDir), newline + "  "));
    else
        msgs{end+1} = sprintf('Training Raman root:   %s\n', string(trainRootDir));
    end
    msgs{end+1} = sprintf('HPLC Excel:        %s\n', hplcExcelPath);
    msgs{end+1} = sprintf('Exp1 Raman root:   %s\n', string(predictRootDir));
    msgs{end+1} = sprintf('Preprocessing:     %s | SNV = %d\n', string(opts.preprocMethod), opts.snvOn);
    msgs{end+1} = sprintf('Train fraction:    %.2f | Test fraction: %.2f\n', opts.trainFraction, 1-opts.trainFraction);

    % PART A: BUILD EXP MEAN-SPECTRUM DATASET
    trainRootDirs = trainRootDir;

    if ischar(trainRootDirs) || isstring(trainRootDirs)
        trainRootDirs = cellstr(trainRootDirs);
    end
    
    hplcSheets = string(opts.hplcSheet);
    
    if numel(hplcSheets) ~= numel(trainRootDirs)
        error('Number of training folders must match number of HPLC sheets.');
    end
    
    X_exp = [];
    Y_all = [];
    G_exp = [];
    keys_exp = {};
    timeNum_exp = [];
    timeLabel_exp = strings(0,1);
    nRep_exp = [];
    XrawMean_exp = [];
    commonShift = [];
    
    for e = 1:numel(trainRootDirs)
    
        hplcMap_e = buildHPLCMap_Shakeflask_Exact( ...
            hplcExcelPath, hplcSheets(e), opts.naMeansZero);
    
        [X_e, Y_e, G_e, keys_e, timeNum_e, timeLabel_e, shift_e, nRep_e, XrawMean_e] = ...
            build_mean_dataset_from_root_local(trainRootDirs{e}, hplcMap_e, opts);
    
        if isempty(commonShift)
            commonShift = shift_e;
        else
            if numel(shift_e) ~= numel(commonShift) || any(abs(shift_e - commonShift) > 1e-10)
                error('Raman shift grid mismatch between training experiments.');
            end
        end
    
        keys_e = strcat("EXP", string(e), "_", string(keys_e));
    
        X_exp = [X_exp; X_e];
        Y_all = [Y_all; Y_e];
        G_exp = [G_exp; G_e + 100*e];
        keys_exp = [keys_exp; cellstr(keys_e(:))];
        timeNum_exp = [timeNum_exp; timeNum_e];
        timeLabel_exp = [timeLabel_exp; timeLabel_e];
        nRep_exp = [nRep_exp; nRep_e];
        XrawMean_exp = [XrawMean_exp; XrawMean_e];
    end

    assert(~isempty(X_exp), 'No Exp mean spectra could be built.');

    msgs{end+1} = sprintf('Built Exp mean-spectra dataset with %d SF/T samples across %d flasks.\n', ...
        size(X_exp,1), numel(unique(G_exp)));

    % PART B: 70/30 SPLIT ON SF/T GROUPS
    nSamples = size(X_exp,1);
    rng(opts.randomSeed);

    perm = randperm(nSamples);
    nTrain = max(2, round(opts.trainFraction * nSamples));
    nTrain = min(nTrain, nSamples-1);

    idxTrain = false(nSamples,1);
    idxTrain(perm(1:nTrain)) = true;
    idxTest = ~idxTrain;

    X_trainval = X_exp(idxTrain,:);
    Y_trainval = Y_all(idxTrain,:);
    G_trainval = G_exp(idxTrain);
    keys_trainval = keys_exp(idxTrain);
    timeNum_trainval = timeNum_exp(idxTrain);
    timeLabel_trainval = timeLabel_exp(idxTrain);
    nRep_trainval = nRep_exp(idxTrain);

    X_test = X_exp(idxTest,:);
    Y_test = Y_all(idxTest,:);
    G_test = G_exp(idxTest);
    keys_test = keys_exp(idxTest);
    timeNum_test = timeNum_exp(idxTest);
    timeLabel_test = timeLabel_exp(idxTest);
    nRep_test = nRep_exp(idxTest);

    msgs{end+1} = sprintf('Train/val SF/T samples: %d | Test SF/T samples: %d\n', size(X_trainval,1), size(X_test,1));

    % PART C: VENETIAN BLINDS ON TRAIN/VAL SET ONLY
    nBlinds = min(opts.nBlinds, size(X_trainval,1));
    foldID_train = make_venetian_blinds_folds_simple(size(X_trainval,1), nBlinds);

    Yhat_trainCV = nan(size(Y_trainval));
    bestLV_perFold = nan(nBlinds,1);
    RMSECV_all = nan(opts.maxLV_wish, nBlinds);
    RMSECV_perFold_perSugar = nan(nBlinds,3);
    RMSECV_perFold_mean = nan(nBlinds,1);

    msgs{end+1} = sprintf('\nVenetian blinds on train/validation set: %d blinds | n=%d mean spectra\n', ...
        nBlinds, size(X_trainval,1));

    for b = 1:nBlinds
        valIdx = (foldID_train == b);
        trainIdx = ~valIdx;

        if ~any(valIdx)
            warning('Blind %d has no validation samples; skipping.', b);
            continue;
        end

        X_train = X_trainval(trainIdx,:);
        Y_train = Y_trainval(trainIdx,:);
        X_val = X_trainval(valIdx,:);
        Y_val = Y_trainval(valIdx,:);

        maxLV_try = min([opts.maxLV_wish, size(X_train,1)-2, size(X_train,2)]);
        maxLV_try = max(maxLV_try, 1);

        RMSECV = nan(maxLV_try,1);
        innerFolds = min(opts.internalCV_folds, size(X_train,1)-1);
        innerFolds = max(innerFolds, 2);

        for a = 1:maxLV_try
            try
                [~,~,~,~,~,~,MSE_CV] = plsregress(X_train, Y_train, a, 'CV', innerFolds);
                mseY = MSE_CV(2:end, a+1);
                rmseY = sqrt(mseY(:));
                RMSECV(a) = mean(rmseY, 'omitnan');
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

        Yhat_trainCV(valIdx,:) = Ypred_val;

        [rmse_fold_sug, rmse_fold_mean] = rmse_multi_local(Y_val, Ypred_val);
        RMSECV_perFold_perSugar(b,:) = rmse_fold_sug;
        RMSECV_perFold_mean(b) = rmse_fold_mean;

        msgs{end+1} = sprintf('Blind %d/%d | nTrain = %d nVal = %d | bestLV = %d | RMSECVmean = %.3f mM | glu = %.3f mM xyl = %.3f mM mtl = %.3f mM\n', ...
            b, nBlinds, size(X_train,1), size(X_val,1), bestLV, rmse_fold_mean, ...
            rmse_fold_sug(1), rmse_fold_sug(2), rmse_fold_sug(3));
    end

    % PART D: OVERALL INTERNAL CV METRICS
    [rmse_cv_perSugar, rmse_cv_mean] = rmse_multi_local(Y_trainval, Yhat_trainCV);

    R2_cv_perSugar = nan(1,3);
    for j = 1:3
        Yj = Y_trainval(:,j);
        Yhatj = Yhat_trainCV(:,j);
        ok = isfinite(Yj) & isfinite(Yhatj);
        denom = sum((Yj(ok) - mean(Yj(ok))).^2);
        if denom > eps
            R2_cv_perSugar(j) = 1 - sum((Yj(ok) - Yhatj(ok)).^2) / denom;
        end
    end
    R2_cv_mean = mean(R2_cv_perSugar, 'omitnan');

    msgs{end+1} = sprintf('\nInternal Venetian blinds CV on 70%% set:\n');
    msgs{end+1} = sprintf('  RMSECV mean = %.3f mM\n', rmse_cv_mean);
    msgs{end+1} = sprintf('  RMSECV per sugar: glu = %.3f mM | xyl = %.3f mM | mtl = %.3f mM\n', rmse_cv_perSugar(1), rmse_cv_perSugar(2), rmse_cv_perSugar(3));
    msgs{end+1} = sprintf('  R^2 per sugar:    glu = %.3f | xyl = %.3f | mtl = %.3f\n', R2_cv_perSugar(1), R2_cv_perSugar(2), R2_cv_perSugar(3));

    % PART E: CHOOSE FINAL LV
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

    msgs{end+1} = sprintf('Selected final LV from 70%% Venetian blinds = %d\n', finalLV);

    % PART F: FIT FINAL MODEL ON ALL TRAIN/VAL, PREDICT TEST
    [XL, YL, XS_train, YS_train, BETA_locked, PCTVAR_locked, MSE_locked, stats_locked] = plsregress(X_trainval, Y_trainval, finalLV); %#ok<ASGLU>

    Yhat_train_fit = [ones(size(X_trainval,1),1) X_trainval] * BETA_locked;
    Yhat_test = [ones(size(X_test,1),1) X_test] * BETA_locked;

    % optional practical clipping: concentrations cannot be negative
    Yhat_train_fit = max(Yhat_train_fit, 0);
    Yhat_test = max(Yhat_test, 0);

    [rmse_fit_perSugar, rmse_fit_mean] = rmse_multi_local(Y_trainval, Yhat_train_fit);
    [rmse_test_perSugar, rmse_test_mean] = rmse_multi_local(Y_test, Yhat_test);

    R2_test_perSugar = nan(1,3);
    for j = 1:3
        Yj = Y_test(:,j);
        Yhatj = Yhat_test(:,j);
        ok = isfinite(Yj) & isfinite(Yhatj);
        denom = sum((Yj(ok) - mean(Yj(ok))).^2);
        if denom > eps
            R2_test_perSugar(j) = 1 - sum((Yj(ok) - Yhatj(ok)).^2) / denom;
        end
    end
    R2_test_mean = mean(R2_test_perSugar, 'omitnan');

    msgs{end+1} = sprintf('\nPrediction on held-out 30%% set:\n');
    msgs{end+1} = sprintf('  RMSEP mean = %.3f mM\n', rmse_test_mean);
    msgs{end+1} = sprintf('  RMSEP per sugar: glu = %.3f mM | xyl = %.3f mM | mtl = %.3f mM\n', rmse_test_perSugar(1), rmse_test_perSugar(2), rmse_test_perSugar(3));
    msgs{end+1} = sprintf('  R^2 per sugar:   glu = %.3f | xyl = %.3f | mtl = %.3f\n', R2_test_perSugar(1), R2_test_perSugar(2), R2_test_perSugar(3));

    % PART G: TEST SET GROUP-MEAN VIEW
    % Here each row already is one mean spectrum per SF/T, so this is already the group-mean level.
    mu_true_test = Y_test;
    mu_pred_test = Yhat_test;
    sd_pred_test = nan(size(Yhat_test));

    [RMSEP_test_group_perSugar, RMSEP_test_group_mean] = rmse_multi_local(mu_true_test, mu_pred_test);

    %msgs{end+1} = sprintf('  RMSEP on prediction SF/T means: mean = %.3f mM | glu = %.3f mM xyl = %.3f mM mtl = %.3f mM\n', ...
        %RMSEP_test_group_mean, RMSEP_test_group_perSugar(1), RMSEP_test_group_perSugar(2), RMSEP_test_group_perSugar(3));


    % PART H: OPTIONAL APPLY TO ANOTHER TECHNICAL REP. EXP + OPTIONAL EXTERNAL RMSEP
    X_exp_rep = [];
    XrawMean_exp_rep = [];
    Yhat_exp_rep = [];
    Yrep_all = [];
    Y_exp_rep = [];
    keys_exp_rep = {};
    G_exp_rep = [];
    timeNum_exp_rep = [];
    timeLabel_exp_rep = strings(0,1);
    nRep_exp_rep = [];
    predTablesExp_rep = {};

    rmse_external_perSugar = nan(1,3);
    rmse_external_mean = NaN;
    R2_external_perSugar = nan(1,3);
    R2_external_mean = NaN;
    RMSEP_external_group_perSugar = nan(1,3);
    RMSEP_external_group_mean = NaN;

    if nargin >= 3 && ~isempty(predictRootDir)

        % If external HPLC should be used, build HPLC map for predict set too
        hplcMap_predict = [];
        if opts.evaluateExternalRMSEP
            if strlength(string(opts.predictHPLCSheet)) == 0
                error('opts.predictHPLCSheet must be provided when opts.evaluateExternalRMSEP = true.');
            end
            hplcMap_predict = buildHPLCMap_Shakeflask_Exact(hplcExcelPath, opts.predictHPLCSheet, opts.naMeansZero);
        end

        [X_exp_rep, Yrep_all, G_exp_rep, keys_exp_rep, timeNum_exp_rep, timeLabel_exp_rep, shift_exp_rep, nRep_exp_rep, XrawMean_exp_rep] = ...
            build_mean_dataset_from_root_local(predictRootDir, hplcMap_predict, opts);

        assert(~isempty(X_exp_rep), 'No Exp tech. rep. mean spectra could be built.');

        shift_exp_rep = shift_exp_rep(:).';
        if numel(shift_exp_rep) ~= numel(commonShift) || any(abs(shift_exp_rep - commonShift) > 1e-10)
            error('Exp tech. rep. Raman shift grid does not match first Exp training grid.');
        end

        if opts.evaluateExternalRMSEP
            Y_exp_rep = Yrep_all;
        else
            Y_exp_rep = [];
        end

        Yhat_exp_rep = [ones(size(X_exp_rep,1),1) X_exp_rep] * BETA_locked;

        % optional practical clipping: concentrations cannot be negative
        Yhat_exp_rep = max(Yhat_exp_rep, 0);

        % ===== External RMSEP / R² =====
        if opts.evaluateExternalRMSEP
            [rmse_external_perSugar, rmse_external_mean] = rmse_multi_local(Y_exp_rep, Yhat_exp_rep);

            for j = 1:3
                Yj = Y_exp_rep(:,j);
                Yhat_j = Yhat_exp_rep(:,j);
                ok = isfinite(Yj) & isfinite(Yhat_j);
                denom = sum((Yj(ok) - mean(Yj(ok))).^2);
                if denom > eps
                    R2_external_perSugar(j) = 1 - sum((Yj(ok) - Yhat_j(ok)).^2) / denom;
                else
                    R2_external_perSugar(j) = NaN;
                end
            end
            R2_external_mean = mean(R2_external_perSugar, 'omitnan');

            % Already one mean spectrum per SF/T, so external group means = same rows
            RMSEP_external_group_perSugar = rmse_external_perSugar;
            RMSEP_external_group_mean = rmse_external_mean;

            msgs{end+1} = sprintf('\nExternal prediction on %s:', string(opts.predictHPLCSheet));
            msgs{end+1} = sprintf('  RMSEP_external mean = %.3f mM', rmse_external_mean);
            msgs{end+1} = sprintf('  RMSEP_external per sugar: glu = %.3f mM | xyl = %.3f mM | mtl = %.3f mM', ...
                rmse_external_perSugar(1), rmse_external_perSugar(2), rmse_external_perSugar(3));
            msgs{end+1} = sprintf('  R^2_external per sugar:   glu = %.3f | xyl = %.3f | mtl = %.3f', ...
                R2_external_perSugar(1), R2_external_perSugar(2), R2_external_perSugar(3));
        end

        if isempty(opts.wantedExp1Flask)
            flasksToPlot = unique(G_exp_rep);
        else
            flasksToPlot = opts.wantedExp1Flask(:);
        end

        predTablesExp_rep = cell(numel(flasksToPlot),1);
        sugarNames = {'Glucose','Xylose','Mannitol'};
        sugarsToPlot = opts.plotSugarsPredict(:)';

        for f = 1:numel(flasksToPlot)
            thisFlask = flasksToPlot(f);
            idxF = (G_exp_rep == thisFlask);

            if ~any(idxF)
                warning('No Exp tech. rep. mean spectra found for shakeflask%d.', thisFlask);
                continue;
            end

            [tSorted, ord] = sort(timeNum_exp_rep(idxF));
            tLabel    = timeLabel_exp_rep(idxF); tLabel = tLabel(ord);
            nRepThis  = nRep_exp_rep(idxF);      nRepThis = nRepThis(ord);
            YhatThis  = Yhat_exp_rep(idxF,:);    YhatThis = YhatThis(ord,:);

            predTable = table(tSorted(:), tLabel(:), nRepThis(:), ...
                              YhatThis(:,1), YhatThis(:,2), YhatThis(:,3), ...
                'VariableNames', {'timeNum','timeLabel','nRep','glu_pred','xyl_pred','mtl_pred'});

            if opts.evaluateExternalRMSEP
                YtrueThis = Y_exp_rep(idxF,:); YtrueThis = YtrueThis(ord,:);
                predTable.glu_meas = YtrueThis(:,1);
                predTable.xyl_meas = YtrueThis(:,2);
                predTable.mtl_meas = YtrueThis(:,3);
                predTable.glu_res  = YhatThis(:,1) - YtrueThis(:,1);
                predTable.xyl_res  = YhatThis(:,2) - YtrueThis(:,2);
                predTable.mtl_res  = YhatThis(:,3) - YtrueThis(:,3);
            end

            predTablesExp_rep{f} = predTable;

            figure('Name', sprintf('Exp predictions - shakeflask%d', thisFlask), ...
                   'Units','normalized','Position',[0.12 0.12 0.78 0.72]);
            tiledlayout(numel(sugarsToPlot),1,'TileSpacing','compact','Padding','compact');

            for jj = 1:numel(sugarsToPlot)
                j = sugarsToPlot(jj);
                nexttile;

                if opts.evaluateExternalRMSEP
                    YtrueThis = Y_exp_rep(idxF,:); 
                    YtrueThis = YtrueThis(ord,:);

                    plot(tSorted, YhatThis(:,j), '-o', 'LineWidth', 1.4, 'MarkerSize', 6, ...
                        'Color', [0 0.4470 0.7410], 'DisplayName','Predicted'); hold on;
                    plot(tSorted, YtrueThis(:,j), '-s', 'LineWidth', 1.2, 'MarkerSize', 6, ...
                        'Color', [0.8500 0.3250 0.0980], 'DisplayName','Measured');

                    % Add blue/red numeric labels next to points
                    for ii = 1:numel(tSorted)
                        text(tSorted(ii), YhatThis(ii,j), sprintf(' %.2f', YhatThis(ii,j)), ...
                            'Color', [0 0.4470 0.7410], 'FontSize', 8, ...
                            'VerticalAlignment','bottom', 'HorizontalAlignment','left');

                        text(tSorted(ii), YtrueThis(ii,j), sprintf(' %.2f', YtrueThis(ii,j)), ...
                            'Color', [0.8500 0.3250 0.0980], 'FontSize', 8, ...
                            'VerticalAlignment','top', 'HorizontalAlignment','left');
                    end

                    legend('Location','best');
                else
                    plot(tSorted, YhatThis(:,j), '-o', 'LineWidth', 1.4, 'MarkerSize', 6, ...
                        'Color', [0 0.4470 0.7410]); hold on;

                    for ii = 1:numel(tSorted)
                        text(tSorted(ii), YhatThis(ii,j), sprintf(' %.2f', YhatThis(ii,j)), ...
                            'Color', [0 0.4470 0.7410], 'FontSize', 8, ...
                            'VerticalAlignment','bottom', 'HorizontalAlignment','left');
                    end
                end

                grid on;
                xlabel('Time (h)');
                ylabel('Concentration (mM)');
                title(sprintf('%s - Exp shakeflask%d', sugarNames{j}, thisFlask));
                xticks(tSorted);
                xticklabels(tLabel);
                xtickangle(45);
            end
        end
    end

    % PART I: VIP
    VIP = computeVIP_multi(X_trainval, Y_trainval, XS_train, stats_locked.W, finalLV);

    % PART J: PLOTS
    if opts.doPlot
        plot_shakeflask_70_30_results( ...
            Y_trainval, Yhat_trainCV, ...
            Y_test, Yhat_test, ...
            rmse_cv_perSugar, R2_cv_perSugar, ...
            rmse_test_perSugar, R2_test_perSugar, ...
            meanRMSECV, finalLV, ...
            RMSECV_perFold_mean, bestLV_perFold, ...
            commonShift, XL, VIP, XS_train, ...
            mu_true_test, mu_pred_test, RMSEP_test_group_perSugar);

        if opts.doPlot && ~isempty(Y_exp_rep) && opts.evaluateExternalRMSEP
            sugarNames = {'Glucose','Xylose','Mannitol'};
    
            figure('Name','External test - multivariate','Units','normalized','Position',[0.10 0.12 0.82 0.72]);
            tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
    
            for j = 1:3
                nexttile;
                plot(Y_exp_rep(:,j), Yhat_exp_rep(:,j), '.', 'MarkerSize', 12); hold on;
                mn = min([Y_exp_rep(:,j); Yhat_exp_rep(:,j)]);
                mx = max([Y_exp_rep(:,j); Yhat_exp_rep(:,j)]);
                plot([mn mx],[mn mx],'k--','LineWidth',1);
                axis equal; grid on;
                xlabel('Measured (mM)');
                ylabel('Predicted (mM)');
                title(sprintf('%s external | RMSEP=%.3f mM | R^2=%.3f', ...
                    sugarNames{j}, rmse_external_perSugar(j), R2_external_perSugar(j)));
            end
    
            for j = 1:3
                nexttile;
                res = Yhat_exp_rep(:,j) - Y_exp_rep(:,j);
                scatter(Y_exp_rep(:,j), res, 20, 'filled'); hold on;
                yline(0,'k--');
                xlabel('Measured');
                ylabel('Residual');
                title(sprintf('%s external residuals', sugarNames{j}));
                grid on;
            end
        end
    end

    % PART K: OUTPUT
    out = struct();

    out.commonShift = commonShift;
    out.finalLV = finalLV;
    out.meanRMSECV = meanRMSECV;
    out.bestLV_perFold = bestLV_perFold;
    out.RMSECV_all = RMSECV_all;
    out.RMSECV_perFold_perSugar = RMSECV_perFold_perSugar;
    out.RMSECV_perFold_mean = RMSECV_perFold_mean;

    out.rmse_cv_perSugar = rmse_cv_perSugar;
    out.rmse_cv_mean = rmse_cv_mean;
    out.R2_cv_perSugar = R2_cv_perSugar;
    out.R2_cv_mean = R2_cv_mean;

    out.rmse_fit_perSugar = rmse_fit_perSugar;
    out.rmse_fit_mean = rmse_fit_mean;

    out.rmse_test_perSugar = rmse_test_perSugar;
    out.rmse_test_mean = rmse_test_mean;
    out.R2_test_perSugar = R2_test_perSugar;
    out.R2_test_mean = R2_test_mean;

    out.RMSEP_test_group_perSugar = RMSEP_test_group_perSugar;
    out.RMSEP_test_group_mean = RMSEP_test_group_mean;

    out.BETA_locked = BETA_locked;
    out.PCTVAR_locked = PCTVAR_locked;
    out.MSE_locked = MSE_locked;
    out.stats_locked = stats_locked;
    out.VIP = VIP;
    out.XL = XL;
    out.XS_train = XS_train;

    out.trainFraction = opts.trainFraction;
    out.randomSeed = opts.randomSeed;
    out.nBlinds = nBlinds;
    out.internalCV_folds = opts.internalCV_folds;

    out.exp2_X = X_exp;
    out.exp2_XrawMean = XrawMean_exp;
    out.exp2_Y = Y_all;
    out.exp2_keys = keys_exp;
    out.exp2_flaskID = G_exp;
    out.exp2_timeNum = timeNum_exp;
    out.exp2_timeLabel = timeLabel_exp;
    out.exp2_nRep = nRep_exp;

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

    out.exp1_X = X_exp_rep;
    out.exp1_XrawMean = XrawMean_exp_rep;
    out.exp1_Yhat = Yhat_exp_rep;
    out.exp1_keys = keys_exp_rep;
    out.exp1_flaskID = G_exp_rep;
    out.exp1_timeNum = timeNum_exp_rep;
    out.exp1_timeLabel = timeLabel_exp_rep;
    out.exp1_nRep = nRep_exp_rep;
    out.predTablesExp1 = predTablesExp_rep;

    out.Y_external = Y_exp_rep;
    out.Yhat_external = Yhat_exp_rep;

    out.rmse_external_perSugar = rmse_external_perSugar;
    out.rmse_external_mean = rmse_external_mean;
    out.R2_external_perSugar = R2_external_perSugar;
    out.R2_external_mean = R2_external_mean;

    out.RMSEP_external_group_perSugar = RMSEP_external_group_perSugar;
    out.RMSEP_external_group_mean = RMSEP_external_group_mean;

    out.predictHPLCSheet = opts.predictHPLCSheet;
    out.evaluateExternalRMSEP = opts.evaluateExternalRMSEP;

    for k = 1:numel(msgs)
        fprintf('%s\n', msgs{k});
    end
    drawnow;

    if opts.saveMat
        save('SF_HPLC_VB_70_30_MultiPLSR_meanSpec_msback_SNVon_maxLV20.mat', 'out');
    end
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

function foldID = make_venetian_blinds_folds_simple(nSamples, nBlinds)
    nBlinds = min(nBlinds, nSamples);
    foldID = mod((1:nSamples)' - 1, nBlinds) + 1;
end

function [rmse_perSugar, rmse_mean] = rmse_multi_local(Y, Yhat)
    e = Y - Yhat;
    rmse_perSugar = sqrt(mean(e.^2, 1, 'omitnan'));
    rmse_mean = mean(rmse_perSugar, 'omitnan');
end

function VIP = computeVIP_multi(X, Y, XS, W, A)
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
            SSY(a) = sum((q.^2) * denom);
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

function plot_shakeflask_70_30_results( ...
    Y_trainval, Yhat_trainCV, ...
    Y_test, Yhat_test, ...
    rmse_cv_perSugar, R2_cv_perSugar, ...
    rmse_test_perSugar, R2_test_perSugar, ...
    meanRMSECV, finalLV, ...
    RMSECV_perFold_mean, bestLV_perFold, ...
    ramanShift, XL, VIP, XS_train, ...
    mu_true_test, mu_pred_test, RMSEP_test_group_perSugar)

    sugarNames = {'Glucose','Xylose','Mannitol'};

    figure('Name','Shake flask 70/30 - internal CV','Units','normalized','Position',[0.05 0.08 0.88 0.78]);
    tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        nexttile;
        plot(Y_trainval(:,j), Yhat_trainCV(:,j), '.', 'MarkerSize', 12); hold on;
        mn = min([Y_trainval(:,j); Yhat_trainCV(:,j)]);
        mx = max([Y_trainval(:,j); Yhat_trainCV(:,j)]);
        plot([mn mx],[mn mx],'k--','LineWidth',1);
        axis equal; grid on;
        xlabel('Measured (mM)');
        ylabel('CV predicted (mM)');
        title(sprintf('%s CV | RMSECV = %.3f mM | R^2 = %.3f', sugarNames{j}, rmse_cv_perSugar(j), R2_cv_perSugar(j)));
    end

    nexttile;
    yyaxis left;
    bar(RMSECV_perFold_mean);
    ylabel('RMSECV mean');
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
        scatter(XS_train(:,1), XS_train(:,2), 35, Y_trainval(:,1), 'filled');
        xlabel('Score 1');
        ylabel('Score 2');
        title('Scores');
        colorbar;
        grid on;
    else
        text(0.2,0.5,'Not enough score components'); axis off;
    end

    nexttile;
    histogram(bestLV_perFold);
    xlabel('Selected LV');
    ylabel('Count');
    title('Distribution of selected LV');
    grid on;

    figure('Name','Shake flask 70/30 - held-out test','Units','normalized','Position',[0.08 0.10 0.82 0.72]);
    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        nexttile;
        plot(Y_test(:,j), Yhat_test(:,j), '.', 'MarkerSize', 12); hold on;
        mn = min([Y_test(:,j); Yhat_test(:,j)]);
        mx = max([Y_test(:,j); Yhat_test(:,j)]);
        plot([mn mx],[mn mx],'k--','LineWidth',1);
        axis equal; grid on;
        xlabel('Measured (mM)');
        ylabel('Predicted (mM)');
        title(sprintf('%s test | RMSEP = %.3f mM | R^2 = %.3f', sugarNames{j}, rmse_test_perSugar(j), R2_test_perSugar(j)));
    end
end

