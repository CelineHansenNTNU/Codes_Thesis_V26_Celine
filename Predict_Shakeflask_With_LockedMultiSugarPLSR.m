function out = Predict_Shakeflask_With_LockedMultiSugarPLSR(mixBaseDir, csvPath, heldGroups, shakeRootDir, wantedFlask, opts)
% What this code deos:
% 1) Trains a multi-response PLSR model (glucose, xylose, mannitol) on mixture standards
%    using a external holdout of standard groups
% 2) Estimates RMSEP on held-out standard mixtures
% 3) Locks the model trained on the training standards
% 4) Applies that locked model to one chosen shakeflask across all timepoints
%
% Folder structure expected for standard mixtures:
%   mixBaseDir/
%     42/
%       42_450mw_1000ms_00001.spc
%       ...
%     43/
%       ...
%
% Folder structure expected for shakeflasks:
%   shakeRootDir/
%     T0/
%       shakeflask1/
%         SF1_450mw_1000ms_00000.spc
%         ...
%       shakeflask2/
%         ...
%     T3/
%       shakeflask1/
%         ...
%
% Needs access to:
%   PreprocessingSpectras2.m
%
% Notes:
% - RMSEP is computed ONLY on held-out standard mixtures
% - Predictions on shakeflasks are external model application; they do NOT
%   have RMSEP unless you have reference concentrations for those samples
%
% This script is based closely on MultiSugar_LOGO_Holdout.m
%
% Celine Hansen + Copilot, 09.03.26

    if nargin < 3 || isempty(heldGroups)
        heldGroups = [52 72];
    end
    heldGroups = heldGroups(:);

    if nargin < 4 || isempty(shakeRootDir)
        shakeRootDir = fullfile(pwd, 'Shakeflask');
    end

    if nargin < 5 || isempty(wantedFlask)
        wantedFlask = 1;
    end

    if nargin < 6
        opts = struct();
    end

    % Keep options as similar as possible to your existing script
    if ~isfield(opts,'shiftRange'),    opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'power_mw'),      opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),          opts.t_ms = 1000; end
    if ~isfield(opts,'snvOn'),         opts.snvOn = false; end
    if ~isfield(opts,'maxLV_wish'),    opts.maxLV_wish = 10; end
    if ~isfield(opts,'lvRule'),        opts.lvRule = "tol"; end
    if ~isfield(opts,'tolFactor'),     opts.tolFactor = 1.05; end
    if ~isfield(opts,'useAchieved'),   opts.useAchieved = true; end
    if ~isfield(opts,'firstID'),       opts.firstID = 42; end

    fprintf('\n=== Predict_Shakeflask_With_LockedMultiSugarPLSR ===\n');
    fprintf('Mixture baseDir:   %s\n', mixBaseDir);
    fprintf('CSV:               %s\n', csvPath);
    fprintf('Held-out groups:   %s\n', mat2str(heldGroups.'));
    fprintf('Shake root dir:    %s\n', shakeRootDir);
    fprintf('Wanted shakeflask: %d\n', wantedFlask);

    %% PART A: BUILD / EVALUATE LOCKED MODEL FROM STANDARD MIXTURES

    % 1) Mixture map: groupID -> [glu xyl mtl]
    mixMap = buildMixtureMapFromCSV_local(csvPath, opts.firstID, opts.useAchieved);

    % 2) Collect all mixture SPC files + group IDs
    [mixFiles, Gmix] = collect_spc_files_by_group_local(mixBaseDir, opts.power_mw, opts.t_ms);

    % Keep only groups present in map
    inMap = isKey(mixMap, num2cell(Gmix));
    mixFiles = mixFiles(inMap);
    Gmix = Gmix(inMap);

    assert(~isempty(mixFiles), 'No mixture .spc files found after filtering by map IDs.');

    msgs = {};
    msgs{end+1} = sprintf('Found %d mixture spectra across %d groups.\n', numel(mixFiles), numel(unique(Gmix)));

    % 3) Read + Hampel + preprocess mixture data directly on native grid
    [Xmix, commonShift, Xmix_raw] = PreprocessingSpectras2(mixFiles, opts.shiftRange, opts.snvOn);
    commonShift = commonShift(:).'; 

    % 4) Plot raw mixture spectra 
    nPerFig = 8;                      % max subplots per figure
    groups = unique(Gmix);
    nG = numel(groups);
    nPages = ceil(nG / nPerFig);
    
    % One distinct color per group (consistent across pages)
    groupCmap = lines(nG);   
    
    for p = 1:nPages
        startIdx = (p-1)*nPerFig + 1;
        endIdx = min(p*nPerFig, nG);
        pageGroups = groups(startIdx:endIdx);
        nThis = numel(pageGroups);
    
        % Layout: sensible grid
        nCol = min(4, nThis);
        nRow = ceil(nThis / nCol);
    
        figure('Name', sprintf('Raw mixture spectra (individual) - groups %d..%d', pageGroups(1), pageGroups(end)), ...
               'Units','normalized','Position',[0.12 0.12 0.76 0.72]);
    
        for k = 1:nThis
            gi = pageGroups(k);
            idx = find(Gmix == gi);
            nSpec = numel(idx);
            baseColor = groupCmap(startIdx + k - 1, :);  % color for this group
    
            subplot(nRow, nCol, k);
            hold on;
    
            % Plot individual spectra with shaded variations of the base color
            for sIdx = 1:nSpec
                t = (sIdx-1) / max(1, nSpec-1);            % 0..1 along spectra
                % fade towards white for earlier spectra and darker for later ones
                fadeFactor = 0.35 + 0.65 * t;             % in [0.35, 1.0]
                lineColor = 1 - (1 - baseColor) * fadeFactor;
                plot(commonShift, Xmix_raw(idx(sIdx),:), 'Color', lineColor, 'LineWidth', 0.8);
            end
    
            % Plot group mean on top for contrast
            mu = mean(Xmix_raw(idx,:), 1);
            plot(commonShift, mu, 'Color', baseColor * 0.15 + 0.85 * [0 0 0], 'LineWidth', 1.6);
    
            hold off;
            xlabel('Raman shift (cm^{-1})');
            ylabel('Intensity');
            title(sprintf('Group %d (n=%d)', gi, nSpec));
            xlim([min(commonShift) max(commonShift)]);
            grid on;
        end
    
        sgtitle(sprintf('Raw sugar mixture spectra (page %d/%d)', p, nPages));
    end


    % 5) Build Y for mixtures
    nMix = numel(mixFiles);
    Ymix = nan(nMix, 3);
    for i = 1:nMix
        Ymix(i,:) = mixMap(Gmix(i));
    end

    % 6) External split on standards
    testIdx_std  = ismember(Gmix, heldGroups);
    trainIdx_std = ~testIdx_std;

    X_train = Xmix(trainIdx_std,:);
    Y_train = Ymix(trainIdx_std,:);
    G_train = Gmix(trainIdx_std);

    X_test = Xmix(testIdx_std,:);
    Y_test = Ymix(testIdx_std,:);
    G_test = Gmix(testIdx_std);

    assert(~isempty(X_test), 'No test samples matched heldGroups. Check heldGroups.');
    msgs{end+1} = sprintf('Standard train samples: %d | Standard holdout test samples: %d\n', size(X_train,1), size(X_test,1));

    % 7) Choose best LV using INTERNAL LOGO only on training standards
    maxLV_try = min([opts.maxLV_wish, size(X_train,1)-2, size(X_train,2)]);
    maxLV_try = max(maxLV_try, 1);

    [bestLV, rmse_logo_curve] = chooseLV_LOGO_multiY_local(X_train, Y_train, G_train, maxLV_try, opts.lvRule, opts.tolFactor);

    msgs{end+1} = sprintf('Internal LOGO on training standards: bestLV = %d | RMSE_LOGO(mean sugars) = %.4f\n', bestLV, rmse_logo_curve(bestLV));

    % 8) Fit LOCKED model on training standards only
    [~,~,~,~,BETA_locked,PCTVAR_locked,MSE_locked,stats_locked] = plsregress(X_train, Y_train, bestLV);

    % 9) Evaluate on standard train and standard holdout
    Yhat_train = [ones(size(X_train,1),1) X_train] * BETA_locked;
    Yhat_test  = [ones(size(X_test,1),1)  X_test ] * BETA_locked;

    [rmse_train_sug, rmse_train_mean] = rmse_multi_local(Y_train, Yhat_train);
    [rmse_test_sug,  rmse_test_mean ] = rmse_multi_local(Y_test,  Yhat_test);

    msgs{end+1} = sprintf('LOCKED model RMSE (train standards): mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', rmse_train_mean, rmse_train_sug(1), rmse_train_sug(2), rmse_train_sug(3));
    msgs{end+1} = sprintf('LOCKED model RMSEP (holdout stds):   mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', rmse_test_mean, rmse_test_sug(1), rmse_test_sug(2), rmse_test_sug(3));

    % RMSEP on mixture-means as well
    heldUnique = unique(G_test);
    Ytrue_m = zeros(numel(heldUnique),3);
    Ypred_m = zeros(numel(heldUnique),3);

    for i = 1:numel(heldUnique)
        idx = (G_test == heldUnique(i));
        Ytrue_m(i,:) = mean(Y_test(idx,:), 1);
        Ypred_m(i,:) = mean(Yhat_test(idx,:), 1);
    end

    rmse_mixmeans = sqrt(mean((Ytrue_m - Ypred_m).^2, 1));
    rmse_mixmeans_mean = mean(rmse_mixmeans);

    msgs{end+1} = sprintf('LOCKED model RMSEP on held-out mixture means: mean=%.3f | glu=%.3f xyl=%.3f mtl=%.3f\n', rmse_mixmeans_mean, rmse_mixmeans(1), rmse_mixmeans(2), rmse_mixmeans(3));

    %% PART B: APPLY LOCKED MODEL TO ONE SHAKEFLASK

    % 10) Collect shakeflask files for the chosen flask only
    shakeInfo = collect_shakeflask_files_oneflask_local(shakeRootDir, wantedFlask, opts.power_mw, opts.t_ms);

    assert(~isempty(shakeInfo.files), 'No shakeflask files found for shakeflask%d.', wantedFlask);

    msgs{end+1} = sprintf('Found %d shakeflask spectra across %d timepoints for shakeflask%d.\n', ...
        numel(shakeInfo.files), numel(unique(shakeInfo.timeNum)), wantedFlask);

    % % 11) Read + Hampel + preprocess shakeflask data directly on native grid
    [Xshake, shakeShift, Xshake_raw] = PreprocessingSpectras2(shakeInfo.files, opts.shiftRange, opts.snvOn);
    shakeShift = shakeShift(:).';

    % Safety check: shakeflask grid must match training grid
    if numel(shakeShift) ~= numel(commonShift) || any(abs(shakeShift - commonShift) > 1e-6)
        error('Shakeflask spectra do not match the training Raman shift grid.');
    end

    % 13) Predict with LOCKED model
    Yhat_shake = [ones(size(Xshake,1),1) Xshake] * BETA_locked;   % [nSpectra x 3]

    % 14) Aggregate by timepoint
    tUnique = unique(shakeInfo.timeNum);
    nT = numel(tUnique);

    mu_pred = zeros(nT, 3);
    sd_pred = zeros(nT, 3);
    n_rep   = zeros(nT, 1);
    timeLabel_sorted = strings(nT,1);

    for i = 1:nT
        idx = (shakeInfo.timeNum == tUnique(i));
        mu_pred(i,:) = mean(Yhat_shake(idx,:), 1);
        sd_pred(i,:) = std(Yhat_shake(idx,:), 0, 1);
        n_rep(i) = sum(idx);
        lab = unique(string(shakeInfo.timeLabel(idx)));
        timeLabel_sorted(i) = lab(1);
    end

    % 15) Plot predicted concentrations vs time
    sugarNames = {'Glucose','Xylose','Mannitol'};

    figure('Name', sprintf('Predicted sugars vs time - shakeflask%d', wantedFlask), 'Units','normalized','Position',[0.12 0.12 0.78 0.72]);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    for j = 1:3
        nexttile;
        errorbar(tUnique, mu_pred(:,j), sd_pred(:,j), '-o', 'LineWidth', 1.4, 'MarkerSize', 6);
        grid on;
        xlabel('Time (h)');
        ylabel('Predicted concentration (mM)');
        title(sprintf('%s - shakeflask%d', sugarNames{j}, wantedFlask));

        xticks(tUnique);
        xticklabels(timeLabel_sorted);

        % uncomment if you prefer numeric time on x-axis:
        % xticklabels(string(tUnique));

        if j == 3
            xtickangle(45);
        end
    end

    sgtitle(sprintf('Locked multi-sugar PLSR predictions for shakeflask%d', wantedFlask));

    % 16) Optional combined plot
    figure('Name', sprintf('Predicted sugars combined - shakeflask%d', wantedFlask), 'Units','normalized','Position',[0.18 0.18 0.72 0.42]);
    hold on;
    for j = 1:3
        errorbar(tUnique, mu_pred(:,j), sd_pred(:,j), '-o', 'LineWidth', 1.4, 'DisplayName', sugarNames{j});
    end
    grid on;
    xlabel('Time (h)');
    ylabel('Predicted concentration (mM)');
    title(sprintf('Predicted sugars vs time - shakeflask%d', wantedFlask));
    xticks(tUnique);
    xticklabels(timeLabel_sorted);
    xtickangle(45);
    legend('Location','best');
    hold off;

    % 17) Print compact table
    predTable = table(tUnique(:), timeLabel_sorted(:), n_rep(:), ...
                      mu_pred(:,1), sd_pred(:,1), ...
                      mu_pred(:,2), sd_pred(:,2), ...
                      mu_pred(:,3), sd_pred(:,3), ...
        'VariableNames', {'timeNum','timeLabel','nRep', ...
                          'glu_mean','glu_std', ...
                          'xyl_mean','xyl_std', ...
                          'mtl_mean','mtl_std'});

    fprintf('\nPredicted concentrations for shakeflask%d:\n', wantedFlask);
    disp(predTable);

    %% PART C: PACK OUTPUT

    out = struct();

    % Locked model / standard performance
    out.commonShift = commonShift;
    out.bestLV = bestLV;
    out.rmse_logo_curve = rmse_logo_curve;

    out.rmse_train_perSugar = rmse_train_sug;
    out.rmse_train_mean = rmse_train_mean;

    out.rmse_holdout_perSugar = rmse_test_sug;
    out.rmse_holdout_mean = rmse_test_mean;

    out.rmse_holdout_mixmeans_perSugar = rmse_mixmeans;
    out.rmse_holdout_mixmeans_mean = rmse_mixmeans_mean;

    out.BETA_locked = BETA_locked;
    out.PCTVAR_locked = PCTVAR_locked;
    out.MSE_locked = MSE_locked;
    out.stats_locked = stats_locked;

    out.G_train = G_train;
    out.Y_train = Y_train;
    out.Yhat_train = Yhat_train;

    out.G_test = G_test;
    out.Y_test = Y_test;
    out.Yhat_test = Yhat_test;

    % Shakeflask predictions
    out.wantedFlask = wantedFlask;
    out.shake_files = shakeInfo.files;
    out.shake_timeNum = shakeInfo.timeNum;
    out.shake_timeLabel = shakeInfo.timeLabel;
    out.Yhat_shake = Yhat_shake;
    out.predTable = predTable;

    out.messages = msgs;

    for k = 1:numel(msgs)
        fprintf('%s\n', msgs{k});
    end
    drawnow;

    save('ShakeflaskPred_results.mat', 'out');
end

%% LOCAL HELPERS

function mixMap = buildMixtureMapFromCSV_local(csvPath, firstID, useAchieved)
% same logic as your existing function

    T = readtable(csvPath);

    if any(strcmpi(T.Properties.VariableNames, 'Origin'))
        newIdx = strcmpi(string(T.Origin), "new");
        Tn = T(newIdx,:);
    else
        error('CSV must contain an Origin column (with "new").');
    end

    ids = (firstID : firstID + height(Tn) - 1).';

    vn = string(Tn.Properties.VariableNames);

    if useAchieved
        gluCol = pickCol_local(vn, ["Cglu","glu"], "achieved");
        xylCol = pickCol_local(vn, ["Cxyl","xyl"], "achieved");
        mtlCol = pickCol_local(vn, ["Cmtl","mtl"], "achieved");
    else
        gluCol = pickCol_local(vn, ["Cglu","glu"], []);
        xylCol = pickCol_local(vn, ["Cxyl","xyl"], []);
        mtlCol = pickCol_local(vn, ["Cmtl","mtl"], []);
    end

    C = [Tn.(gluCol), Tn.(xylCol), Tn.(mtlCol)];

    mixMap = containers.Map('KeyType','double','ValueType','any');
    for i = 1:numel(ids)
        mixMap(ids(i)) = C(i,:);
    end

    fprintf('Built mixture map for %d groups (%d..%d)\n', numel(ids), ids(1), ids(end));
end

function colName = pickCol_local(varNames, mustContainAny, alsoContainAll)

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
        error('Could not find column matching requested pattern.');
    end
    colName = hits(1);
end

function [files, G] = collect_spc_files_by_group_local(baseDir, power_mw, t_ms)
% same logic as your existing function

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

function [bestLV, rmse_curve] = chooseLV_LOGO_multiY_local(Xtr, Ytr, Gtr, maxLV_try, lvRule, tolFactor)

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

            [~, rmseMean] = rmse_multi_local(Ytr(te,:), Yhat);
            rmse_g(gi) = rmseMean;
        end

        rmse_perLV(a,:) = rmse_g;
    end

    rmse_curve = mean(rmse_perLV, 2, 'omitnan');

    [minRMSE, idxMin] = min(rmse_curve);
    bestLV = idxMin;

    if lvRule == "tol"
        cand = find(rmse_curve <= tolFactor * minRMSE, 1, 'first');
        if ~isempty(cand)
            bestLV = cand;
        end
    end
end

function [rmse_perSugar, rmse_mean] = rmse_multi_local(Y, Yhat)

    e = Y - Yhat;
    rmse_perSugar = sqrt(mean(e.^2, 1));
    rmse_mean = mean(rmse_perSugar);
end

function shakeInfo = collect_shakeflask_files_oneflask_local(shakeRootDir, wantedFlask, power_mw, t_ms)
% Collect all .spc files for one chosen shakeflask across all T folders.
%
% Expected:
%   shakeRootDir/T0/shakeflask1/*.spc
%   shakeRootDir/T3/shakeflask1/*.spc
%   ...

    d = dir(fullfile(shakeRootDir, '**', '*.spc'));

    files = {};
    timeNum = [];
    timeLabel = strings(0,1);
    flaskNum = [];

    for i = 1:numel(d)
        fp = fullfile(d(i).folder, d(i).name);

        % filename filters
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

        % expecting ... / Txx / shakeflaskN
        if numel(parts) < 2
            continue;
        end

        lastPart = string(parts{end});       % shakeflaskN
        prevPart = string(parts{end-1});     % Txx

        tokF = regexp(char(lastPart), 'shakeflask(\d+)', 'tokens', 'once', 'ignorecase');
        tokT = regexp(char(prevPart), 'T(\d+)', 'tokens', 'once', 'ignorecase');

        if isempty(tokF) || isempty(tokT)
            continue;
        end

        fnum = str2double(tokF{1});
        tnum = str2double(tokT{1});

        if fnum ~= wantedFlask
            continue;
        end

        files{end+1,1} = fp; %#ok<AGROW>
        flaskNum(end+1,1) = fnum; %#ok<AGROW>
        timeNum(end+1,1) = tnum; %#ok<AGROW>
        timeLabel(end+1,1) = "T" + string(tnum); %#ok<AGROW>
    end

    if isempty(files)
        shakeInfo = struct('files', {{}}, 'flaskNum', [], 'timeNum', [], 'timeLabel', strings(0,1));
        return;
    end

    % Sort by time
    [~, ord] = sort(timeNum);
    files = files(ord);
    flaskNum = flaskNum(ord);
    timeNum = timeNum(ord);
    timeLabel = timeLabel(ord);

    shakeInfo = struct();
    shakeInfo.files = files;
    shakeInfo.flaskNum = flaskNum;
    shakeInfo.timeNum = timeNum;
    shakeInfo.timeLabel = timeLabel;
end

