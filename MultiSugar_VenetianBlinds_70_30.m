function results = MultiSugar_VenetianBlinds_70_30(baseDir, tablePath, nBlinds, internalCV_folds, opts)
% Multivariat PLSR for intermediate mixed-sugar dataset.
%
% Workflow:
%   1) Read same filesystem and same Excel/CSV mixture table as before
%   2) Split mixture groups into 70% train/validation and 30% prediction
%   3) Use Venetian blinds CV ONLY on the 70% train/validation part
%   4) Choose final LV from train/validation CV
%   5) Fit final model on all 70% train/validation data
%   6) Predict the held-out 30% prediction set
%
% Same preprocessing options as before:
%   opts.preprocMethod = "msbackadj" or "asls" or "peaknorm"
%   opts.useSNV = true/false
%   for peaknorm, useSNV if to use mean spectrum or not
%
% Needs access to:
%   - PreprocessingSpectras2.m
%   - PreprocessingSpectras_AsLS.m
%
% Celine Hansen + Copilot Chat 15.04.26

    if nargin < 3 || isempty(nBlinds), nBlinds = 10; end
    if nargin < 4 || isempty(internalCV_folds), internalCV_folds = 10; end
    if nargin < 5, opts = struct(); end

    if ~isfield(opts,'shiftRange'),      opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'power_mw'),        opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),            opts.t_ms = 1000; end
    if ~isfield(opts,'maxLV_wish'),      opts.maxLV_wish = 10; end
    if ~isfield(opts,'lvRule'),          opts.lvRule = "min"; end
    if ~isfield(opts,'tolFactor'),       opts.tolFactor = 1.05; end
    if ~isfield(opts,'useAchieved'),     opts.useAchieved = true; end
    if ~isfield(opts,'firstID'),         opts.firstID = 2; end
    if ~isfield(opts,'doPlot'),          opts.doPlot = false; end
    if ~isfield(opts,'preprocMethod'),   opts.preprocMethod = "msbackadj"; end
    if ~isfield(opts,'useSNV'),          opts.useSNV = false; end
    if ~isfield(opts,'trainFraction'),   opts.trainFraction = 0.70; end
    if ~isfield(opts,'randomSeed'),      opts.randomSeed = 1; end
    if ~isfield(opts,'useMeanSpec'),      opts.useMeanSpec = false; end
    if ~isfield(opts,'saveMat'),         opts.saveMat = true; end

    fprintf('\n=== MultiSugar_VenetianBlinds_70_30 ===\n');
    fprintf('BaseDir: %s\n', baseDir);
    fprintf('Table:   %s\n', tablePath);
    fprintf('Preprocessing: %s | useSNV=%d\n', string(opts.preprocMethod), opts.useSNV);
    fprintf('Train/val fraction: %.2f | Prediction fraction: %.2f\n', ...
        opts.trainFraction, 1-opts.trainFraction);

    %% 1) Build groupID -> [glu xyl mtl]
    mixMap = buildMixtureMapFromTable(tablePath, opts.firstID, opts.useAchieved);

    %% 2) Collect spc files
    [files, G] = collect_spc_files_by_group(baseDir, opts.power_mw, opts.t_ms);

    inMap = isKey(mixMap, num2cell(G));
    files = files(inMap);
    G = G(inMap);

    assert(~isempty(files), 'No .spc files found after filtering by map IDs.');
    fprintf('Found %d spectra across %d groups.\n', numel(files), numel(unique(G)));

    %% 3) Preprocess spectra using selected method
    [X, ramanShift, ~] = run_selected_preprocessing(files, opts.shiftRange, opts.preprocMethod, opts.useSNV);

    %% 4) Build Y = [glu xyl mtl]
    n = numel(files);
    Y = nan(n,3);
    for i = 1:n
        Y(i,:) = mixMap(G(i));
    end

    %% 4b) Optionally average preprocessed spectra per group
    if opts.useMeanSpec
        [X, Y, G, groupcounts] = average_preprocessed(X, Y, G);
        fprintf('Using mean preprocessed spectrum.\n');
        fprintf('Reduced dataset to %d groups.\n', numel(G));
    else
        groupcounts = [];
    end

    sugarNames = {'Glucose','Xylose','Mannitol'};

    %% 5) Group-based 70/30 split
    groups = unique(G, 'stable');
    nGroups = numel(groups);

    rng(opts.randomSeed);
    perm = randperm(nGroups);

    nTrainGroups = max(1, round(opts.trainFraction * nGroups));
    nTrainGroups = min(nTrainGroups, nGroups-1); % ensure at least 1 test group

    trainGroups = groups(perm(1:nTrainGroups));
    testGroups  = groups(perm(nTrainGroups+1:end));

    trainMask = ismember(G, trainGroups);
    testMask  = ismember(G, testGroups);

    X_trainval = X(trainMask,:);
    Y_trainval = Y(trainMask,:);
    G_trainval = G(trainMask);

    X_test = X(testMask,:);
    Y_test = Y(testMask,:);
    G_test = G(testMask);

    assert(~isempty(X_trainval), 'Training/validation set is empty.');
    assert(~isempty(X_test), 'Prediction/test set is empty.');

    fprintf('Train/val groups: %d | Prediction groups: %d\n', numel(trainGroups), numel(testGroups));
    fprintf('Train/val spectra: %d | Prediction spectra: %d\n', size(X_trainval,1), size(X_test,1));

    %% 6) Venetian blinds assignment ONLY within training groups
    if opts.useMeanSpec
        nBlinds = min(nBlinds, numel(unique(G_trainval)));
        foldID_train = make_venetian_blinds_folds_simple(numel(G_trainval), nBlinds);
    else
        foldID_train = make_venetian_blinds_folds_by_group(G_trainval, nBlinds);
    end

    Yhat_trainCV = nan(size(Y_trainval));
    bestLV_perFold = nan(nBlinds,1);
    RMSEP_perFold_perSugar = nan(nBlinds,3);
    RMSEP_perFold_mean = nan(nBlinds,1);
    RMSECV_all = nan(opts.maxLV_wish, nBlinds);

    fprintf('\nVenetian blinds on train/validation set: %d blinds | %d groups | n=%d spectra\n', ...
        nBlinds, numel(unique(G_trainval)), size(X_trainval,1));

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
        %innerFolds = max(2, min(internalCV_folds, size(X_train,1)));
        innerFolds = min(internalCV_folds, size(X_train, 1) -1);
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

        [rmse_fold_sug, rmse_fold_mean] = rmse_multi(Y_val, Ypred_val);
        RMSEP_perFold_perSugar(b,:) = rmse_fold_sug;
        RMSEP_perFold_mean(b) = rmse_fold_mean;

        fprintf(['Blind %d/%d | nTrain=%d nVal=%d | bestLV=%d | RMSECVmean=%.3f | ' ...
                 'glu=%.3f xyl=%.3f mtl=%.3f\n'], ...
            b, nBlinds, size(X_train,1), size(X_val,1), bestLV, rmse_fold_mean, ...
            rmse_fold_sug(1), rmse_fold_sug(2), rmse_fold_sug(3));
    end

    %% 7) Internal train/validation CV metrics
    [RMSEP_CV_perSugar, RMSEP_CV_mean] = rmse_multi(Y_trainval, Yhat_trainCV);
    R2_CV_perSugar = r2_multi(Y_trainval, Yhat_trainCV);
    R2_CV_mean = mean(R2_CV_perSugar, 'omitnan');

    fprintf('\nInternal Venetian blinds CV on 70%% set:\n');
    fprintf('  RMSECV mean = %.3f mM\n', RMSEP_CV_mean);
    fprintf('  RMSECV per sugar: glu=%.3f | xyl=%.3f | mtl=%.3f\n', ...
        RMSEP_CV_perSugar(1), RMSEP_CV_perSugar(2), RMSEP_CV_perSugar(3));
    fprintf('  R^2 per sugar:  glu=%.3f | xyl=%.3f | mtl=%.3f\n', ...
        R2_CV_perSugar(1), R2_CV_perSugar(2), R2_CV_perSugar(3));

    %% 8) Choose final LV from mean RMSECV across folds
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

    fprintf('Selected final LV from 70%% Venetian blinds = %d\n', finalLV);

    %% 9) Final model on all 70% train/validation data -> predict 30% test
    [XL, YL, XS_train, YS_train, BETA, PCTVAR, MSE, stats] = plsregress(X_trainval, Y_trainval, finalLV); %#ok<ASGLU>

    Yhat_test = [ones(size(X_test,1),1) X_test] * BETA;
    Yhat_train_fit = [ones(size(X_trainval,1),1) X_trainval] * BETA;

    [RMSEP_test_perSugar, RMSEP_test_mean] = rmse_multi(Y_test, Yhat_test);
    R2_test_perSugar = r2_multi(Y_test, Yhat_test);
    R2_test_mean = mean(R2_test_perSugar, 'omitnan');

    fprintf('\nPrediction on held-out 30%% set:\n');
    fprintf('  RMSEP mean = %.3f mM\n', RMSEP_test_mean);
    fprintf('  RMSEP per sugar: glu=%.3f | xyl=%.3f | mtl=%.3f\n', ...
        RMSEP_test_perSugar(1), RMSEP_test_perSugar(2), RMSEP_test_perSugar(3));
    fprintf('  R^2 per sugar:  glu=%.3f | xyl=%.3f | mtl=%.3f\n', ...
        R2_test_perSugar(1), R2_test_perSugar(2), R2_test_perSugar(3))

    %% 10) Group means for test/prediction set
    if opts.useMeanSpec
        % already one row per group  
        mu_true_test = Y_test;
        mu_pred_test = Yhat_test;
        sd_pred_test = nan(size(Yhat_test));

        [RMSEP_test_group_perSugar, RMSEP_test_group_mean] = rmse_multi(mu_true_test, mu_pred_test);
        fprintf(' RMSEP on pred-group means: mean = %.3f | glu = %.3f xyl = %.3f mtl = %.3f\n', RMSEP_test_group_mean, RMSEP_test_group_perSugar(1), RMSEP_test_group_perSugar(2), RMSEP_test_group_perSugar(3));
    else
        testGroupsUnique = unique(G_test);
        mu_true_test = nan(numel(testGroupsUnique),3);
        mu_pred_test = nan(numel(testGroupsUnique),3);
        sd_pred_test = nan(numel(testGroupsUnique),3);

        for i = 1:numel(testGroupsUnique)
            idx = (G_test == testGroupsUnique(i));
            mu_true_test(i,:) = mean(Y_test(idx,:), 1, 'omitnan');
            mu_pred_test(i,:) = mean(Yhat_test(idx,:), 1, 'omitnan');
            sd_pred_test(i,:) = std(Yhat_test(idx,:), 0, 1, 'omitnan');
        end

        [RMSEP_test_group_perSugar, RMSEP_test_group_mean] = rmse_multi(mu_true_test, mu_pred_test);

        fprintf('  RMSEP on prediction-group means: mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', ...
            RMSEP_test_group_mean, RMSEP_test_group_perSugar(1), RMSEP_test_group_perSugar(2), RMSEP_test_group_perSugar(3));
    end

    %% 11) VIP from final model
    VIP = computeVIP_multi(X_trainval, Y_trainval, XS_train, stats.W, finalLV);

    %% 12) Pack results
    results.files = files;
    results.G = G;
    results.Y = Y;

    results.trainGroups = trainGroups;
    results.testGroups = testGroups;
    results.trainMask = trainMask;
    results.testMask = testMask;

    results.X_trainval = X_trainval;
    results.Y_trainval = Y_trainval;
    results.G_trainval = G_trainval;

    results.X_test = X_test;
    results.Y_test = Y_test;
    results.G_test = G_test;

    results.Yhat_trainCV = Yhat_trainCV;
    results.Yhat_train_fit = Yhat_train_fit;
    results.Yhat_test = Yhat_test;

    results.foldID_train = foldID_train;
    results.ramanShift = ramanShift;
    results.sugarNames = sugarNames;

    results.RMSEP_CV_perSugar = RMSEP_CV_perSugar;
    results.RMSEP_CV_mean = RMSEP_CV_mean;
    results.R2_CV_perSugar = R2_CV_perSugar;
    results.R2_CV_mean = R2_CV_mean;

    results.RMSEP_test_perSugar = RMSEP_test_perSugar;
    results.RMSEP_test_mean = RMSEP_test_mean;
    results.R2_test_perSugar = R2_test_perSugar;
    results.R2_test_mean = R2_test_mean;

    results.mu_true_test = mu_true_test;
    results.mu_pred_test = mu_pred_test;
    results.sd_pred_test = sd_pred_test;
    results.RMSEP_test_group_perSugar = RMSEP_test_group_perSugar;
    results.RMSEP_test_group_mean = RMSEP_test_group_mean;

    results.bestLV_perFold = bestLV_perFold;
    results.RMSEP_perFold_perSugar = RMSEP_perFold_perSugar;
    results.RMSEP_perFold_mean = RMSEP_perFold_mean;
    results.RMSECV_all = RMSECV_all;
    results.meanRMSECV = meanRMSECV;
    results.finalLV = finalLV;

    results.XL = XL;
    results.YL = YL;
    results.XS_train = XS_train;
    results.YS_train = YS_train;
    results.BETA = BETA;
    results.PCTVAR = PCTVAR;
    results.MSE = MSE;
    results.stats = stats;
    results.VIP = VIP;
    results.preprocessing = opts.preprocMethod;
    results.useSNV = opts.useSNV;
    results.trainFraction = opts.trainFraction;
    results.randomSeed = opts.randomSeed;
    results.useMeanSpec = opts.useMeanSpec;
    results.groupcounts = groupcounts;

    %% 13) Plotting
    if opts.doPlot
        plot_multisugar_70_30_results(results);

        [meas_sorted, idx] = sort(results.Y_test(:,1));
        pred_sorted = results.Yhat_test(idx,1);
        figure('Name', 'Sorted measured vs predicted (glucose)');
        plot(meas_sorted, '-o', 'DisplayName', 'Measured'); hold on;
        plot(pred_sorted, '-s', 'DisplayName', 'Predicted');
        legend('Location', 'best'); xlabel('Sorted sample index'); ylabel('Glucose (mM)');
        title('Glucose: measured vs predicted (test set)');
        grid on;

        [meas_sorted, idx] = sort(results.Y_test(:,2));
        pred_sorted = results.Yhat_test(idx,2);
        figure('Name', 'Sorted measured vs predicted (xylose)');
        plot(meas_sorted, '-o', 'DisplayName', 'Measured'); hold on;
        plot(pred_sorted, '-s', 'DisplayName', 'Predicted');
        legend('Location', 'best'); xlabel('Sorted sample index'); ylabel('Xylose (mM)');
        title('Xylose: measured vs predicted (test set)');
        grid on;

        [meas_sorted, idx] = sort(results.Y_test(:,3));
        pred_sorted = results.Yhat_test(idx,3);
        figure('Name', 'Sorted measured vs predicted (mannitol)');
        plot(meas_sorted, '-o', 'DisplayName', 'Measured'); hold on;
        plot(pred_sorted, '-s', 'DisplayName', 'Predicted');
        legend('Location', 'best'); xlabel('Sorted sample index'); ylabel('Mannitol (mM)');
        title('Mannitol: measured vs predicted (test set)');
        grid on;
    end

    %% Print compact result table
    fprintf('\nPredictions on test set:\n');
    tbl_print = table(results.files(results.testMask), results.G_test, results.Y_test, results.Yhat_test, 'VariableNames', {'File', 'Group', 'Measured_mM', 'Predicted_mM'});
    disp(tbl_print);

    if opts.saveMat
        save('Multi_VB_res_msbackadj_noSNV_2', 'results');
    end
end

%% ========================= Helpers =========================

function [X_processed, ramanShift, X_raw] = run_selected_preprocessing(inputData, shiftOrAxis, preprocMethod, useSNV)

    preprocMethod = lower(string(preprocMethod));

    switch preprocMethod
        case "msbackadj"
            [X_processed, ramanShift, X_raw] = PreprocessingSpectras2(inputData, shiftOrAxis, useSNV);

        case "asls"
            [X_processed, ramanShift, X_raw] = PreprocessingSpectras_AsLS(inputData, shiftOrAxis, useSNV);

        case "peaknorm"
            [X_processed, ramanShift, X_raw] = PreprocessingSpectrasNew(inputData, shiftOrAxis, useSNV);

        otherwise
            error('Unknown preprocessing method: %s. Use "msbackadj" or "asls".', preprocMethod);
    end
end

function mixMap = buildMixtureMapFromTable(tablePath, firstID, useAchieved)
    T = readtable(tablePath);

    if any(strcmpi(T.Properties.VariableNames, 'Origin'))
        newIdx = strcmpi(string(T.Origin), "new");
        Tn = T(newIdx,:);
    else
        error('Table must contain an Origin column with value "new".');
    end

    ids = (firstID : firstID + height(Tn) - 1).';
    vn = string(Tn.Properties.VariableNames);

    if useAchieved
        gluCol = pickColFlexible(vn, ["cglu","glu"], ["ach"]);
        xylCol = pickColFlexible(vn, ["cxyl","xyl"], ["ach"]);
        mtlCol = pickColFlexible(vn, ["cmtl","mtl"], ["ach"]);
    else
        gluCol = pickColFlexible(vn, ["cglu","glu"], []);
        xylCol = pickColFlexible(vn, ["cxyl","xyl"], []);
        mtlCol = pickColFlexible(vn, ["cmtl","mtl"], []);
    end

    C = [Tn.(gluCol), Tn.(xylCol), Tn.(mtlCol)];

    mixMap = containers.Map('KeyType','double','ValueType','any');
    for i = 1:numel(ids)
        mixMap(ids(i)) = C(i,:);
    end
end

function colName = pickColFlexible(varNames, mustContainAny, alsoContainAll)
    namesLower = lower(varNames);
    idx = false(size(varNames));

    for k = 1:numel(mustContainAny)
        idx = idx | contains(namesLower, lower(mustContainAny(k)));
    end

    if ~isempty(alsoContainAll)
        for k = 1:numel(alsoContainAll)
            idx = idx & contains(namesLower, lower(alsoContainAll(k)));
        end
    end

    hits = varNames(idx);

    if isempty(hits) && ~isempty(alsoContainAll)
        idx = false(size(varNames));
        for k = 1:numel(mustContainAny)
            idx = idx | contains(namesLower, lower(mustContainAny(k)));
        end
        hits = varNames(idx);
    end

    if isempty(hits)
        error('Could not find matching column.');
    end

    colName = hits(1);
end

function [files, G] = collect_spc_files_by_group(baseDir, power_mw, t_ms)
    d = dir(fullfile(baseDir, '**', '*.spc'));
    files = strings(numel(d),1);
    G = nan(numel(d),1);

    for i = 1:numel(d)
        fp = fullfile(d(i).folder, d(i).name);
        files(i) = fp;

        [parentFolder, ~] = fileparts(fp);
        [~, groupName] = fileparts(parentFolder);
        gid = str2double(groupName);

        if ~isfinite(gid)
            tok = regexp(d(i).name, '^(?<id>\d+)_', 'names', 'once');
            if ~isempty(tok)
                gid = str2double(tok.id);
            end
        end
        G(i) = gid;
    end

    keep = true(size(files));
    if ~isempty(power_mw)
        keep = keep & contains(files, sprintf('%dmw', power_mw), 'IgnoreCase', true);
    end
    if ~isempty(t_ms)
        keep = keep & contains(files, sprintf('%dms', t_ms), 'IgnoreCase', true);
    end

    files = cellstr(files(keep));
    G = G(keep);
end

function foldID = make_venetian_blinds_folds_by_group(G, nBlinds)
    G = G(:);
    foldID = nan(size(G));
    groups = unique(G, 'stable');

    for i = 1:numel(groups)
        idx = find(G == groups(i));
        f = mod((1:numel(idx))' - 1, nBlinds) + 1;
        foldID(idx) = f;
    end
end

function [rmse_perSugar, rmse_mean] = rmse_multi(Y, Yhat)
    E = Y - Yhat;
    rmse_perSugar = sqrt(mean(E.^2, 1, 'omitnan'));
    rmse_mean = mean(rmse_perSugar, 'omitnan');
end

function R2 = r2_multi(Y, Yhat)
    nResp = size(Y,2);
    R2 = nan(1,nResp);
    for j = 1:nResp
        ok = isfinite(Y(:,j)) & isfinite(Yhat(:,j));
        if sum(ok) > 1 && numel(unique(Y(ok,j))) > 1
            R2(j) = corr(Y(ok,j), Yhat(ok,j)).^2;
        end
    end
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

function plot_multisugar_70_30_results(results)

    sugarNames = results.sugarNames;

    %% Figure 1: Internal CV on 70% set
    figure('Name','Internal Venetian blinds on 70% train/validation set', ...
           'Units','normalized','Position',[0.05 0.08 0.88 0.78]);

    tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        nexttile;
        plot(results.Y_trainval(:,j), results.Yhat_trainCV(:,j), '.', 'MarkerSize', 10); hold on;
        mn = min([results.Y_trainval(:,j); results.Yhat_trainCV(:,j)]);
        mx = max([results.Y_trainval(:,j); results.Yhat_trainCV(:,j)]);
        plot([mn mx],[mn mx],'k--','LineWidth',1);
        grid on; axis equal;
        xlabel('Measured (mM)');
        ylabel('CV predicted (mM)');
        title(sprintf('%s CV | RMSECV = %.3f mM | R^2 = %.3f', ...
            sugarNames{j}, results.RMSEP_CV_perSugar(j), results.R2_CV_perSugar(j)));
    end

    nexttile;
    yyaxis left;
    bar(results.RMSEP_perFold_mean);
    ylabel('RMSECV mean');
    yyaxis right;
    plot(results.bestLV_perFold, '-ok', 'LineWidth', 1.5, 'MarkerFaceColor','w');
    ylabel('Best LV');
    xlabel('Blind');
    grid on;
    title('Train/val: RMSECV and best LV per fold');

    nexttile;
    plot(1:numel(results.meanRMSECV), results.meanRMSECV, '-o', 'LineWidth', 1.5); hold on;
    xline(results.finalLV, 'k--', sprintf('finalLV=%d', results.finalLV), ...
        'LabelVerticalAlignment','bottom');
    grid on;
    xlabel('LV');
    ylabel('Mean RMSECV');
    title('Train/val: RMSECV vs LV');

    nexttile;
    plot(results.ramanShift, results.XL(:,1), 'k-', 'LineWidth', 1.2); hold on;
    yyaxis right;
    plot(results.ramanShift, results.VIP, 'r--', 'LineWidth', 1.0);
    xlabel('Raman shift');
    ylabel('VIP');
    yyaxis left;
    ylabel('Loading 1');
    grid on;
    title('Final model: Loading 1 and VIP');

    nexttile;
    if size(results.XS_train,2) >= 2
        gscatter(results.XS_train(:,1), results.XS_train(:,2), results.Y_trainval(:,1));
        xlabel('Score 1'); ylabel('Score 2');
        title('Scores (colored by glucose)');
        grid on;
    else
        text(0.1,0.5,'Not enough score components','FontAngle','italic'); axis off;
    end

    nexttile;
    histogram(results.bestLV_perFold);
    xlabel('Selected LV');
    ylabel('Count');
    title('Distribution of selected LV');
    grid on;

    %% Figure 2: External prediction on 30% set
    figure('Name','Prediction on held-out 30% set', ...
           'Units','normalized','Position',[0.08 0.10 0.82 0.72]);

    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        nexttile;
        plot(results.Y_test(:,j), results.Yhat_test(:,j), '.', 'MarkerSize', 10); hold on;
        mn = min([results.Y_test(:,j); results.Yhat_test(:,j)]);
        mx = max([results.Y_test(:,j); results.Yhat_test(:,j)]);
        plot([mn mx],[mn mx],'k--','LineWidth',1);
        grid on; axis equal;
        xlabel('Measured (mM)');
        ylabel('Predicted (mM)');
        title(sprintf('%s test | RMSEP = %.3f mM | R^2 = %.3f', ...
            sugarNames{j}, results.RMSEP_test_perSugar(j), results.R2_test_perSugar(j)));
    end

    for j = 1:3
        nexttile;
        errorbar(results.mu_true_test(:,j), results.mu_pred_test(:,j), results.sd_pred_test(:,j), ...
            'o', 'LineWidth', 1.2); hold on;
        mn = min([results.mu_true_test(:,j); results.mu_pred_test(:,j)]);
        mx = max([results.mu_true_test(:,j); results.mu_pred_test(:,j)]);
        plot([mn mx],[mn mx],'k--','LineWidth',1);
        grid on; axis equal;
        xlabel('Measured group mean (mM)');
        ylabel('Predicted group mean (mM)');
        title(sprintf('%s test group means | RMSEP = %.3f mM', ...
            sugarNames{j}, results.RMSEP_test_group_perSugar(j)));
    end

    %% Figure 3: Residuals on test set
    figure('Name','Residual diagnostics on held-out 30% set', ...
           'Units','normalized','Position',[0.12 0.12 0.78 0.65]);

    tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        res = results.Yhat_test(:,j) - results.Y_test(:,j);

        nexttile;
        scatter(results.Yhat_test(:,j), res, 18, 'filled'); hold on;
        yline(0,'k--');
        xlabel('Predicted');
        ylabel('Residual');
        title(sprintf('%s: Residual vs predicted', sugarNames{j}));
        grid on;

        nexttile;
        scatter(results.Y_test(:,j), res, 18, 'filled'); hold on;
        yline(0,'k--');
        xlabel('Measured');
        ylabel('Residual');
        title(sprintf('%s: Residual vs measured', sugarNames{j}));
        grid on;
    end
end

function [Xg, Yg, Gg, counts] = average_preprocessed(X, Y, G)
    G = G(:);
    groups = unique(G, 'stable');
    nGroups = numel(groups);
    Xg = zeros(nGroups, size(X, 2));
    Yg = zeros(nGroups, size(Y, 2));
    Gg = zeros(nGroups,1);
    counts = zeros(nGroups,1);

    for i = 1:nGroups
        idx = (G == groups(i));
        Xg(i,:) = mean(X(idx,:), 1, 'omitnan');
        Yg(i,:) = mean(Y(idx,:), 1, 'omitnan');
        Gg(i) = groups(i);
        counts(i) = sum(idx);
    end
end

function foldID = make_venetian_blinds_folds_simple(nSamples, nBlinds)
    nBlinds = min(nBlinds, nSamples);
    foldID = mod((1:nSamples)' - 1, nBlinds) + 1;
end