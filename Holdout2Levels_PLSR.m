function results = Holdout2Levels_PLSR(HPLCmap, heldLevels, opts)
% Hold out two concentration levels completely (external test),
% choose #LV using internal CV on training only, refit, and compute RMSEP on external test set.
%
% Requires:
%   find_all_files_by_folder.m
%   PreprocessingSpectras2.m
%
% Celine Hansen + Copilot, updated 09.03.26

    % Defaults
    if nargin < 1 || isempty(HPLCmap)
        %% Glucose + MQ
        %keys = {'glu100mm', 'glu50mm', 'glu25mm', 'glu12.5mm', 'glu6.25mm', 'glu3.125mm', 'glu1.5625mm'}; %'glu250mm', 'glu100mm', 
        %vals = [101.31 50.014 24.952 12.487 6.3219 3.2119 1.5088]; % Actual measured HPLC values for specific sugar, 246.91 101.31 
        
        %% Glucose + CGXII
        keys = {'CGXII50mM', 'CGXII25mM', 'CGXII12.5mM', 'CGXII6.25mM', 'CGXII3.125mM', 'CGXII1.5625mM'}; % 'CGXII40mM', 'CGXII30mM',, 'CGXII20mM'
        vals = [49.2094 27.6768 10.0836 6.5985 3.5054 1.2093]; %161.3888,  110.7403, 25.2327, 18.9182, 12.6749 , 10.0836
        
        %% Mannitol + MQ
        %keys = {'mtl100mM2', 'mtl50mM2', 'mtl25mM2', 'mtl12.5mM', 'mtl6.25mM', 'mtl3.125mM', 'mtl1.5625mM'};
        %vals = [48.6179 25.3189 15.5550 9.0447 4.5253 2.2671];

        %% Xylose + MQ
        %keys = {'Xyl100mM', 'Xyl50mM', 'Xyl25mM', 'Xyl12.5mM', 'Xyl6.25mM', 'Xyl3.125mM', 'Xyl1.5625mM'};
        %vals = [100.2628 49.7676 25.6101 12.385 5.9981 2.9972 1.4225];

        %% Arabinose + MQ
        %keys = {'ara100mM', 'ara50mM', 'ara25mM', 'ara12.5mM', 'ara6.25mM', 'ara3.125mM', 'ara1.5625mM'};
        %vals = [99.990 50.9224 22.9899 12.6587 6.4765 3.2980 1.6872];

        %% First T0 of exp 2
        %keys = {'shakeflask1', 'shakeflask2', 'shakeflask3', 'shakeflask4', 'shakeflask5', 'shakeflask6', 'shakeflask7', 'shakeflask8', 'shakeflask9'}; 
        %vals = [83.26 83.26 33.3 0 83.26 33.3 0.0 0 33.3]; 

        %vals = [100 50 25 12.5 6.25 3.125 1.5625];
        HPLCmap = containers.Map(keys, num2cell(vals));
    end

    if nargin < 2 || isempty(heldLevels)
        % Default: hold out mid-levels as example
        heldLevels = [12.5 3.125];
    end
    heldLevels = heldLevels(:).'; % row

    if nargin < 3, opts = struct(); end

    % Data/source settings
    if ~isfield(opts,'baseDir'),           opts.baseDir = fullfile(pwd,'glu+CGXII_CH/'); end
    if ~isfield(opts,'power_mw'),          opts.power_mw = 450; end
    if ~isfield(opts,'t_ms'),              opts.t_ms = 1000; end
    if ~isfield(opts,'shiftRange'),        opts.shiftRange = [400 1800]; end
    if ~isfield(opts,'maxLV_wish'),        opts.maxLV_wish = 10; end
    if ~isfield(opts,'internalCV_folds'),  opts.internalCV_folds = 7; end
    if ~isfield(opts,'matchTol'),          opts.matchTol = 1e-2; end %1e-2

    % LV selection rule
    % "min"  -> pick LV with minimum CV RMSE
    % "1se"  -> pick simplest LV within 1 standard error of min (more conservative)
    if ~isfield(opts,'lvRule'),            opts.lvRule = "min"; end

    fprintf('\nHoldout2Levels_PLSR: holding out levels: %s mM\n', mat2str(heldLevels,4));

    % Load spectra and Y
    [allFiles, Y] = find_all_files_by_folder(opts.baseDir, opts.power_mw, opts.t_ms, HPLCmap);
    assert(~isempty(allFiles), 'Found no .spc-files.');

    % Read files + preprocessing pipeline
    [X, ramanShift, Xraw] = PreprocessingSpectras2(allFiles, opts.shiftRange, true);
    ramanShift = ramanShift(:).';

    % Split into train/test by held-out levels
    testIdx = false(size(Y));
    for k = 1:numel(heldLevels)
        testIdx = testIdx | (abs(Y - heldLevels(k)) < opts.matchTol);
    end
    trainIdx = ~testIdx;

    X_train = X(trainIdx,:);
    Y_train = Y(trainIdx);
    X_test  = X(testIdx,:);
    Y_test  = Y(testIdx);

    assert(~isempty(X_test), 'No samples matched heldLevels. Check heldLevels or matchTol.');

    fprintf('Train samples: %d | Test samples (held-out): %d\n', size(X_train,1), size(X_test,1));

    % Safe max LV for this training set
    maxLV_try = min([opts.maxLV_wish, size(X_train,1)-2, size(X_train,2)]);
    maxLV_try = max(maxLV_try, 1);

    % Internal CV on training ONLY: choose bestLV
    % plsregress with CV returns MSE with size [2 x (maxLV+1)].
    % For univariate Y, row 2 corresponds to response. Col j corresponds to model with (j-1) LVs.
    [~,~,~,~,~,~,MSE_CV] = plsregress(X_train, Y_train, maxLV_try, 'CV', opts.internalCV_folds);

    rmse_cv = sqrt(MSE_CV(2, 2:end));  % 1..maxLV_try
    [minRMSE, idxMin] = min(rmse_cv);

    bestLV = idxMin;

    % Optional 1-SE rule (more conservative)
    if opts.lvRule == "1se"
        % Estimate SE of CV RMSE using fold-averaged MSE is not directly accessible from plsregress.
        % Practical alternative: use a simple "within 10% of min" rule if you want conservative LV.
        tol = 1.05; % 5% within min
        cand = find(rmse_cv <= tol*minRMSE, 1, 'first');
        if ~isempty(cand), bestLV = cand; end
    end

    fprintf('Internal CV (training only): bestLV = %d | RMSECV = %.4f mM\n', bestLV, rmse_cv(bestLV));

    % Fit final model on full training with bestLV
    [~,~,~,~,BETA_best,~,~,~] = plsregress(X_train, Y_train, bestLV);

    plot_beta_stability(X_train, Y_train, ramanShift, maxLV_try, bestLV);

    % Predictions
    Yhat_train = [ones(size(X_train,1),1) X_train] * BETA_best;
    Yhat_test  = [ones(size(X_test,1),1)  X_test]  * BETA_best;

    % Calculating the mean prediction
    [uniqY_train, ~, ic_train] = unique(Y_train);
    meanPred_train = accumarray(ic_train, Yhat_train, [], @mean);
    [uniqY_test, ~, ic_test] = unique(Y_test);
    meanPred_test = accumarray(ic_test, Yhat_test, [], @mean);

    % Metrics
    RMSEC  = sqrt(mean((Y_train - Yhat_train).^2));
    RMSEP  = sqrt(mean((Y_test  - Yhat_test ).^2));

    R2_train = corr(Y_train, Yhat_train).^2;
    R2_test  = corr(Y_test,  Yhat_test ).^2;

    bias_test = mean(Yhat_test - Y_test);

    fprintf('Final model: RMSEC = %.4f | R2_train = %.4f\n', RMSEC, R2_train);
    fprintf('External test (held-out levels): RMSEP = %.4f | R2_test = %.4f | bias = %.4f\n', RMSEP, R2_test, bias_test);

    % Collect results
    results.heldLevels = heldLevels;
    results.nTrain = numel(Y_train);
    results.nTest  = numel(Y_test);

    results.bestLV = bestLV;
    results.rmse_cv_curve = rmse_cv(:);
    results.RMSECV_best = rmse_cv(bestLV);

    results.RMSEC = RMSEC;
    results.RMSEP = RMSEP;

    results.R2_train = R2_train;
    results.R2_test  = R2_test;

    results.bias_test = bias_test;

    results.Y_train = Y_train;
    results.Yhat_train = Yhat_train;
    results.Y_test = Y_test;
    results.Yhat_test = Yhat_test;
    results.BETA = BETA_best; % BETA from plsregress

    % Plots
    figure('Name','Hold-out 2 levels: LV selection + test performance', ...
           'Units','normalized','Position',[0.15 0.15 0.75 0.65]);
    t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    % 1) RMSECV curve on training
    nexttile;
    plot(1:maxLV_try, rmse_cv, '-o'); grid on;
    xline(bestLV,'k--',sprintf('bestLV=%d',bestLV),'LabelVerticalAlignment','bottom');
    xlabel('Number of LVs'); ylabel('RMSECV (mM)');
    title(sprintf('Internal CV on training (K=%d)', opts.internalCV_folds));
    lims = xlim; xlim([1 lims(2)]);

    % 2) Train measured vs predicted
    nexttile;
    plot(Y_train, Yhat_train, 'o'); hold on;
    plot([min(Y_train) max(Y_train)], [min(Y_train) max(Y_train)], 'k--');
    plot(uniqY_train, meanPred_train, 's', 'MarkerSize', 8, 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
    for k = 1:numel(uniqY_train)
        txt = sprintf('%.2f', meanPred_train(k));       % one number only
        text(uniqY_train(k), meanPred_train(k), txt, 'VerticalAlignment','bottom', 'FontSize',10, 'Color','r');
    end
    grid on; axis equal;
    xlabel('Measured (mM)'); ylabel('Predicted (mM)');
    title(sprintf('Training fit: R^2=%.3f, RMSEC=%.3f', R2_train, RMSEC));

    % 3) Test measured vs predicted
    nexttile;
    plot(Y_test, Yhat_test, 'o'); hold on;
    plot([min(Y_test) max(Y_test)], [min(Y_test) max(Y_test)], 'k--');
    plot(uniqY_test, meanPred_test, 's', 'MarkerSize', 8, 'MarkerFaceColor','r', 'MarkerEdgeColor','k');
    for k = 1:numel(uniqY_test)
        txt = sprintf('%.2f', meanPred_test(k));       % one number only
        text(uniqY_test(k), meanPred_test(k), txt, 'VerticalAlignment','bottom', 'FontSize',10, 'Color','r');
    end
    grid on; axis equal;
    xlabel('Measured (mM)'); ylabel('Predicted (mM)');
    title(sprintf('External test (held levels): R^2=%.3f, RMSEP=%.3f', R2_test, RMSEP));

    % 4) Residuals on test (by level)
    nexttile;
    res = Yhat_test - Y_test;
    scatter(Y_test, res, 40, 'filled'); grid on;
    yline(0,'k--');
    xlabel('Test level (mM)'); ylabel('Residual (Pred - True)');
    title(sprintf('Test residuals (bias=%.3f)', bias_test));

    save('PLSR_results_xylMQ.mat', 'results');
end

function plot_beta_stability(X_train, Y_train, ramanShift, maxLV_try, bestLV)
% Plot beta curves across LV choices + smoothing for bestLV.
% X_train: [n x p], Y_train: [n x 1]
% ramanShift: [1 x p] or [p x 1]
% maxLV_try: max LV to evaluate
% bestLV: selected LV

    ramanShift = ramanShift(:);
    p = size(X_train,2);

    % Fit betas for LV=1..maxLV_try
    B = nan(p, maxLV_try);   % store beta(2:end) columns per LV

    for a = 1:maxLV_try
        [~,~,~,~,BETA] = plsregress(X_train, Y_train, a);
        B(:,a) = BETA(2:end);
    end

    % 1) All beta curves (raw) 
    figure('Name','PLSR beta stability across #LV','Units','normalized','Position',[0.1 0.1 0.85 0.7]);
    tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(ramanShift, B, 'LineWidth', 0.8); grid on;
    hold on;
    xline(ramanShift(1), 'k:'); % dummy to keep style consistent
    title('\beta(\nu) for LV = 1..maxLV');
    xlabel('Raman shift (cm^{-1})');
    ylabel('\beta (mM / a.u.)');
    legend(arrayfun(@(k)sprintf('LV%d',k),1:maxLV_try,'UniformOutput',false), ...
        'Location','eastoutside');

    % Highlight chosen bestLV
    plot(ramanShift, B(:,bestLV), 'k', 'LineWidth', 2.0);
    hold off;

    %  2) Zoomed view (optional): fingerprint region 
    nexttile;
    mask = ramanShift >= 400 & ramanShift <= 1800;
    plot(ramanShift(mask), B(mask,bestLV), 'k', 'LineWidth', 1.6); grid on;
    title(sprintf('Chosen model: LV%d (raw \\beta)', bestLV));
    xlabel('Raman shift (cm^{-1})'); ylabel('\beta');

    %  3) Smoothed beta for bestLV 
    nexttile;
    beta_best = B(:,bestLV);

    % Savitzky-Golay smoothing:
    % window length must be odd and <= p. Pick something sensible.
    win = min(101, p - mod(p+1,2)); % nearest odd <= p
    if mod(win,2)==0, win = win-1; end
    if win < 11, win = 11; end % safety

    polyOrder = 3;
    beta_smooth = sgolayfilt(beta_best, polyOrder, win);

    plot(ramanShift(mask), beta_best(mask), 'Color',[0.6 0.6 0.6], 'LineWidth', 0.8); hold on;
    plot(ramanShift(mask), beta_smooth(mask), 'r', 'LineWidth', 1.8);
    grid on; hold off;
    title(sprintf('Chosen LV%d: raw vs smoothed (sgolay win=%d)', bestLV, win));
    xlabel('Raman shift (cm^{-1})'); ylabel('\beta');
    legend('raw','smoothed','Location','best');

    % 4) Absolute beta (importance proxy) 
    nexttile;
    plot(ramanShift(mask), abs(beta_best(mask)), 'k', 'LineWidth', 1.3); grid on;
    title(sprintf('|\\beta| for chosen LV%d (importance proxy)', bestLV));
    xlabel('Raman shift (cm^{-1})'); ylabel('| \beta |');

end
