function results = MultiSugar_LOGO_Holdout(baseDir, csvPath, heldGroups, opts)
% MultiSugar_LOGO_Holdout
% - Multi-response PLSR (3 output: glucose, xylose, mannitol) from Raman spectra
% - Group-aware: every sugar mixture has a groupID (42..93)
% - External test where 1 or 2 mixtures are completely held out
% - Intern choice of model with LOGO (leave-one-group-out) on training set to choose bestLV
%
% Folder structure expected:
%   baseDir/
%     88/
%       88_450mw_1000ms_00091.spc
%       ...
%     89/
%       ...
% GroupID is taken from map-name (88,89,...).
%
% CSV expected:
%   DOE_3_sugars_output...csv with Origin column and concentration columns.
%   For Origin=="new", rows correspond to groupIDs 42..93 in order.
%
% Needs access to:
%   PreprocessingSpectras2.m
%
% Celine Hansen + ChatGPT, 2026-02-25

    if nargin < 3 || isempty(heldGroups)
        heldGroups = [52 72]; % example
    end
    heldGroups = heldGroups(:);

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'shiftRange'),    opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'power_mw'),      opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),          opts.t_ms = 1000; end
    if ~isfield(opts,'snvOn'),         opts.snvOn = false; end
    if ~isfield(opts,'maxLV_wish'),    opts.maxLV_wish = 10; end
    if ~isfield(opts,'lvRule'),        opts.lvRule = "min"; end            % "min" or "tol"
    if ~isfield(opts,'tolFactor'),     opts.tolFactor = 1.05; end          % used if lvRule="tol"
    if ~isfield(opts,'useAchieved'),   opts.useAchieved = true; end        % achieved vs nominal from CSV
    if ~isfield(opts,'firstID'),       opts.firstID = 2; end              % start ID for "new" rows

    fprintf('\n=== MultiSugar_LOGO_Holdout ===\n');
    fprintf('BaseDir: %s\n', baseDir);
    fprintf('CSV: %s\n', csvPath);
    fprintf('Held-out groups: %s\n', mat2str(heldGroups.'));

    % 1) Build mapping groupID -> [glu xyl mtl]
    mixMap = buildMixtureMapFromCSV(csvPath, opts.firstID, opts.useAchieved);

    % 2) Collect all SPC files + their groupIDs
    [files, G] = collect_spc_files_by_group(baseDir, opts.power_mw, opts.t_ms);

    % Keep only files whose group exists in map (42..93)
    inMap = isKey(mixMap, num2cell(G));
    files = files(inMap);
    G = G(inMap);

    assert(~isempty(files), 'No .spc files found after filtering by map IDs.');

    fprintf('Found %d spectra across %d groups.\n', numel(files), numel(unique(G)));

    % 3) Read + Hampel + preprocess directly on native grid
    [X, ramanShift, Xraw] = PreprocessingSpectras2(files, opts.shiftRange, opts.snvOn);
    ramanShift = ramanShift(:).';

    % 5) Build Y (n x 3) from map
    n = numel(files);
    Y = nan(n,3);
    for i = 1:n
        Y(i,:) = mixMap(G(i));
    end

    % 6) External split: hold out groups
    testIdx = ismember(G, heldGroups);
    trainIdx = ~testIdx;

    X_train = X(trainIdx,:);
    Y_train = Y(trainIdx,:);
    G_train = G(trainIdx);

    X_test = X(testIdx,:);
    Y_test = Y(testIdx,:);
    G_test = G(testIdx);

    assert(~isempty(X_test), 'No test samples matched heldGroups. Check heldGroups.');
    fprintf('Train samples: %d | Test samples: %d\n', size(X_train,1), size(X_test,1));

    % 7) Choose bestLV by INTERNAL LOGO on TRAINING groups
    maxLV_try = min([opts.maxLV_wish, size(X_train,1)-2, size(X_train,2)]);
    maxLV_try = max(maxLV_try, 1);

    [bestLV, rmse_logo_curve] = chooseLV_LOGO_multiY(X_train, Y_train, G_train, maxLV_try, opts.lvRule, opts.tolFactor);

    fprintf('Internal LOGO (train only): bestLV=%d | RMSE_LOGO(mean sugars)=%.4f mM\n', ...
        bestLV, rmse_logo_curve(bestLV));

    % 8) Fit final multivariate PLSR on full training with bestLV and predict test
    [~,~,~,~,BETA_best,PCTVAR_best,MSE_best,stats_best] = plsregress(X_train, Y_train, bestLV);

    Yhat_train = [ones(size(X_train,1),1) X_train] * BETA_best;   % nTrain x 3
    Yhat_test  = [ones(size(X_test,1),1)  X_test ] * BETA_best;   % nTest x 3

    % 9) Metrics
    [rmse_train_sug, rmse_train_mean] = rmse_multi(Y_train, Yhat_train);
    [rmse_test_sug,  rmse_test_mean ] = rmse_multi(Y_test,  Yhat_test );

    fprintf('Final model RMSE (train): mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', ...
        rmse_train_mean, rmse_train_sug(1), rmse_train_sug(2), rmse_train_sug(3));
    fprintf('External test RMSE:       mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', ...
        rmse_test_mean, rmse_test_sug(1), rmse_test_sug(2), rmse_test_sug(3));

    % 10) Plots
    %plot_multisugar_results(Y_train, Yhat_train, Y_test, Yhat_test, rmse_logo_curve, bestLV, heldGroups);

    % 11) Pack results
    results.heldGroups = heldGroups;
    results.ramanShift = ramanShift;
    results.bestLV = bestLV;
    results.rmse_logo_curve = rmse_logo_curve;

    results.rmse_train_perSugar = rmse_train_sug;
    results.rmse_train_mean = rmse_train_mean;
    results.rmse_test_perSugar = rmse_test_sug;
    results.rmse_test_mean = rmse_test_mean;

    results.BETA = BETA_best;
    results.stats = stats_best;
    results.PCTVAR = PCTVAR_best;
    results.MSE = MSE_best;

    results.G_train = G_train;
    results.G_test  = G_test;
    results.Y_train = Y_train;
    results.Yhat_train = Yhat_train;
    results.Y_test = Y_test;
    results.Yhat_test = Yhat_test;

    % --- plotting: per-held-out-group means ± std (test) ---
    %G = results.G_test;
    %Y = results.Y_test;
    %Yhat = results.Yhat_test;

    %groups = unique(G);
    %mu_true = zeros(numel(groups),3);
    %mu_pred = zeros(numel(groups),3);
    %sd_pred = zeros(numel(groups),3);

    %for i = 1:numel(groups)
    %    idx = (G==groups(i));
    %    mu_true(i,:) = mean(Y(idx,:),1);
    %    mu_pred(i,:) = mean(Yhat(idx,:),1);
    %    sd_pred(i,:) = std(Yhat(idx,:),0,1);
    %end

    % Recompute per-unique measured means (train/test) robustly
    [uniqY_train_rows, ~, ic_tr] = unique(Y_train, 'rows', 'stable');
    meanPred_train_rows = zeros(size(uniqY_train_rows));
    for k = 1:size(uniqY_train_rows,1)
        meanPred_train_rows(k,:) = mean(Yhat_train(ic_tr==k,:), 1);
    end

    [uniqY_test_rows, ~, ic_te] = unique(Y_test, 'rows', 'stable');
    meanPred_test_rows = zeros(size(uniqY_test_rows));
    for k = 1:size(uniqY_test_rows,1)
        meanPred_test_rows(k,:) = mean(Yhat_test(ic_te==k,:), 1);
    end

    % Per-held-out-group means ± std (test)
    groups = unique(G_test);
    mu_true = zeros(numel(groups),3);
    mu_pred = zeros(numel(groups),3);
    sd_pred = zeros(numel(groups),3);
    for i = 1:numel(groups)
        idx = (G_test==groups(i));
        mu_true(i,:) = mean(Y_test(idx,:),1);
        mu_pred(i,:) = mean(Yhat_test(idx,:),1);
        sd_pred(i,:) = std(Yhat_test(idx,:),0,1);
    end

    sug = {'Glu','Xyl','Mtl'};
    figure('Name','Test results per held-out mixture (means ± std)','Units','normalized','Position',[0.2 0.2 0.75 0.35]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        nexttile;
        errorbar(mu_true(:,j), mu_pred(:,j), sd_pred(:,j), 'o', 'LineWidth', 1.5); hold on;
        mn = min([mu_true(:,j); mu_pred(:,j)]); mx = max([mu_true(:,j); mu_pred(:,j)]);
        plot([mn mx],[mn mx],'k--','LineWidth',1);
        grid on; axis equal;
        lim = [mn-0.1*(mx-mn), mx+0.1*(mx-mn)];
        xlim(lim); ylim(lim);
        % plot unique-sample mean predictions too
        plot(uniqY_test_rows(:,j), meanPred_test_rows(:,j), 's', 'MarkerSize', 8, 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'DisplayName','mean per mixture');
        text(mu_true(:,j), mu_pred(:,j), "  ID " + string(groups), 'VerticalAlignment','bottom');
        xlabel('Measured (mM)'); ylabel('Predicted (mM)');
        title(sprintf('%s (held-out means)', sug{j}));
        legend('Location','best');
    end

    % Summary figure with LV curve and train scatter
    figure('Name','Multi-sugar PLSR summary','Units','normalized','Position',[0.1 0.1 0.85 0.75]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % LV curve (use rmse_logo_curve)
    nexttile;
    plot(1:numel(rmse_logo_curve), rmse_logo_curve, '-o'); grid on;
    xline(bestLV,'k--',sprintf('bestLV=%d',bestLV),'LabelVerticalAlignment','bottom');
    xlabel('LV'); ylabel('RMSE (mean over sugars)');
    title('Internal LOGO on training');

    sugarNames = {'Glucose','Xylose','Mannitol'};
    for j = 1:3
        nexttile;
        plot(Y_train(:,j), Yhat_train(:,j), 'o'); hold on;
        mn = min([Y_train(:,j); Yhat_train(:,j)]); mx = max([Y_train(:,j); Yhat_train(:,j)]);
        plot([mn mx],[mn mx],'k--');
        grid on; axis equal;
        lim = [mn-0.1*(mx-mn), mx+0.1*(mx-mn)];
        xlim(lim); ylim(lim);
        % plot per-unique-sample mean predictions for train
        plot(uniqY_train_rows(:,j), meanPred_train_rows(:,j), 's', 'MarkerSize', 8, 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'DisplayName','mean per mixture');
        xlabel('Measured'); ylabel('Predicted');
        title(sprintf('Train: %s', sugarNames{j}));
        legend('Location','best');
    end

    % RMSE på mixture-mean (ikke per spektrum)
    groups = unique(results.G_test);
    Ytrue_m = zeros(numel(groups),3);
    Ypred_m = zeros(numel(groups),3);

    for i = 1:numel(groups)
        idx = results.G_test == groups(i);
        Ytrue_m(i,:) = mean(results.Y_test(idx,:),1);
        Ypred_m(i,:) = mean(results.Yhat_test(idx,:),1);
    end

    rmse_mix = sqrt(mean((Ytrue_m - Ypred_m).^2,1));   % per sugar
    rmse_mix_mean = mean(rmse_mix);

    fprintf('External test RMSE on mixture-means: mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', rmse_mix_mean, rmse_mix(1), rmse_mix(2), rmse_mix(3));

    save('MultiPred_results.mat', 'results');
end

%% ========================= Helpers =========================

function mixMap = buildMixtureMapFromCSV(csvPath, firstID, useAchieved)
% Returns containers.Map where key=groupID (double) value=[glu xyl mtl]

    T = readtable(csvPath);

    % Filter to "new"
    if any(strcmpi(T.Properties.VariableNames, 'Origin'))
        newIdx = strcmpi(string(T.Origin), "new");
        Tn = T(newIdx,:);
    else
        error('CSV must contain an Origin column (with "new").');
    end

    ids = (firstID : firstID + height(Tn) - 1).';

    % Robust column lookup
    vn = string(Tn.Properties.VariableNames);

    if useAchieved
        gluCol = pickCol(vn, ["Cglu","glu"], "achieved");
        xylCol = pickCol(vn, ["Cxyl","xyl"], "achieved");
        mtlCol = pickCol(vn, ["Cmtl","mtl"], "achieved");
    else
        gluCol = pickCol(vn, ["Cglu","glu"], []);
        xylCol = pickCol(vn, ["Cxyl","xyl"], []);
        mtlCol = pickCol(vn, ["Cmtl","mtl"], []);
    end

    C = [Tn.(gluCol), Tn.(xylCol), Tn.(mtlCol)];

    mixMap = containers.Map('KeyType','double','ValueType','any');
    for i = 1:numel(ids)
        mixMap(ids(i)) = C(i,:);
    end

    fprintf('Built mixture map for %d groups (%d..%d)\n', numel(ids), ids(1), ids(end));
end

function colName = pickCol(varNames, mustContainAny, alsoContainAll)
% varNames: string array of column names
% mustContainAny: string array; at least one must match
% alsoContainAll: string array; all must be present (optional)

    idx = false(size(varNames));
    for k = 1:numel(mustContainAny)
        idx = idx | contains(lower(varNames), lower(mustContainAny(k)));
    end
    if ~isempty(alsoContainAll)
        for k = 1:numel(alsoContainAll)
            idx = idx & contains(lower(varNames), lower(alsoContainAll(k)));
        end
    end

    hits = varNames(idx);
    if isempty(hits)
        error('Could not find column matching %s (+%s).', mat2str(mustContainAny), mat2str(alsoContainAll));
    end
    colName = hits(1); % first match
end

function [files, G] = collect_spc_files_by_group(baseDir, power_mw, t_ms)
% Recursively find .spc files and extract groupID from parent folder name (42..93)
% Also optionally filter by power/t_ms in filename if present.

    d = dir(fullfile(baseDir, '**', '*.spc'));
    files = strings(numel(d),1);
    G = nan(numel(d),1);

    for i = 1:numel(d)
        fp = fullfile(d(i).folder, d(i).name);
        files(i) = fp;

        % groupID from parent folder name
        [parentFolder, ~] = fileparts(fp);
        [~, groupName] = fileparts(parentFolder);
        gid = str2double(groupName);
        if ~isfinite(gid)
            % fallback: parse from filename (e.g., 88_450mw_1000ms_00091.spc)
            tok = regexp(d(i).name, '^(?<id>\d+)_', 'names', 'once');
            if ~isempty(tok), gid = str2double(tok.id); end
        end
        G(i) = gid;
    end

    % Filter by power/time tags if present in filename
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

function [bestLV, rmse_curve] = chooseLV_LOGO_multiY(Xtr, Ytr, Gtr, maxLV_try, lvRule, tolFactor)
% Internal LOGO: hold out one mixture group at a time on training.
% Choose LV that minimizes mean RMSE across sugars AND held-out groups.

    groups = unique(Gtr);
    ng = numel(groups);
    rmse_perLV = nan(maxLV_try, ng);

    for a = 1:maxLV_try
        rmse_g = nan(ng,1);
        for gi = 1:ng
            gHold = groups(gi);
            te = (Gtr == gHold);
            tr = ~te;

            [~,~,~,~,BETA] = plsregress(Xtr(tr,:), Ytr(tr,:), a);
            Yhat = [ones(sum(te),1) Xtr(te,:)] * BETA;

            [~, rmseMean] = rmse_multi(Ytr(te,:), Yhat);
            rmse_g(gi) = rmseMean;
        end
        rmse_perLV(a,:) = rmse_g;
    end

    rmse_curve = mean(rmse_perLV, 2, 'omitnan');

    [minRMSE, idxMin] = min(rmse_curve);
    bestLV = idxMin;

    if lvRule == "tol"
        cand = find(rmse_curve <= tolFactor*minRMSE, 1, 'first');
        if ~isempty(cand), bestLV = cand; end
    end
end

function [rmse_perSugar, rmse_mean] = rmse_multi(Y, Yhat)
% RMSE per column + mean RMSE across outputs
    e = Y - Yhat;
    rmse_perSugar = sqrt(mean(e.^2, 1));
    rmse_mean = mean(rmse_perSugar);
end

%function plot_multisugar_results(Ytr, Yhtr, Yte, Yhte, rmse_curve, bestLV, heldGroups)
% Simple summary plots for 3 sugars

    %sugarNames = {'Glucose','Xylose','Mannitol'};

    %figure('Name', 'Choice of optimal number of LVs', 'Units', 'normalized', 'Position', [0.1 0.1 0.85 0.75]);

    

    % Test scatter for each sugar
    %for j = 1:3
     %   nexttile;
      % mn = min(Yte(:,j)); mx = max(Yte(:,j));
       % plot([mn mx],[mn mx],'k--'); grid on; axis equal;
       % xlabel('Measured'); ylabel('Predicted');
       % title(sprintf('Test (held %s): %s', mat2str(heldGroups.'), sugarNames{j}));
    %end
%end

