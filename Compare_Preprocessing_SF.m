%% Compare_Preprocessing_Shakeflask_IndividualSpectra
% Compare preprocessing pipelines on INDIVIDUAL Raman spectra
% for one selected shakeflask across all timepoints.
%
% Folder structure expected:
%   Shakeflask/dateFolder/shakeflaskX/*.spc
%
% Requires:
%   tgspcread
%
% Celine Hansen + Copilot Chat, spring 2026

clear; clc; close all;

%% USER INPUT

baseDir      = fullfile(pwd, 'Shakeflasks_v2/');
wantedFlask  = 1;              % choose shakeflask here
power_mw     = 450;
t_ms         = 1000;
shiftWin     = [400 1800];

% IMPORTANT:
% Old noise window [2400 2600] does NOT fit when shiftWin = [400 1800].
% Choose a quieter region inside the cropped range.
noiseWin     = [1180 1220];

% Analysis windows where real Raman peaks are expected
analWins     = [ ...
    1050 1400;   % sugar / fingerprint region
    900 1050;  % strong peak region
    1550 1700  % broad high-wavenumber feature in your data
    ];

% Which pipeline to use for detailed baseline example plot
examplePipelineIdx = 1;

% Plot labels for mean spectra
labelXpos = 1180;

%% DEFINE PREPROCESSING PIPELINES

pipelines = {};

% 1) Current-style pipeline
pipelines{end+1} = struct( ...
    'name', 'Current_msbackadj_SNV_SG', ...
    'baseline', 'msbackadj', ...
    'useSNV', true, ...
    'useVectorNorm', false, ...
    'sgolayOrder', 2, ...
    'sgolayFrame', 11, ...
    'sgolayDeriv', 0, ...
    'aslsLambda', [], ...
    'aslsP', [], ...
    'aslsNIter', [] );

% 2) Current-style without SNV
pipelines{end+1} = struct( ...
    'name', 'Current_msbackadj_noSNV_SG', ...
    'baseline', 'msbackadj', ...
    'useSNV', false, ...
    'useVectorNorm', true, ...
    'sgolayOrder', 2, ...
    'sgolayFrame', 11, ...
    'sgolayDeriv', 0, ...
    'aslsLambda', [], ...
    'aslsP', [], ...
    'aslsNIter', [] );

% 3) AsLS + SNV + SG
pipelines{end+1} = struct( ...
    'name', 'AsLS_SNV_SG', ...
    'baseline', 'asls', ...
    'useSNV', true, ...
    'useVectorNorm', false, ...
    'sgolayOrder', 2, ...
    'sgolayFrame', 11, ...
    'sgolayDeriv', 0, ...
    'aslsLambda', 1e6, ...
    'aslsP', 0.001, ...
    'aslsNIter', 15 );

% 4) AsLS + vector norm + SG
pipelines{end+1} = struct( ...
    'name', 'AsLS_vecnorm_SG', ...
    'baseline', 'asls', ...
    'useSNV', false, ...
    'useVectorNorm', true, ...
    'sgolayOrder', 2, ...
    'sgolayFrame', 11, ...
    'sgolayDeriv', 0, ...
    'aslsLambda', 1e6, ...
    'aslsP', 0.001, ...
    'aslsNIter', 15 );

% 5) AsLS + SG first derivative + SNV
pipelines{end+1} = struct( ...
    'name', 'AsLS_SNV_SG1stDer', ...
    'baseline', 'asls', ...
    'useSNV', true, ...
    'useVectorNorm', false, ...
    'sgolayOrder', 2, ...
    'sgolayFrame', 15, ...
    'sgolayDeriv', 1, ...
    'aslsLambda', 1e6, ...
    'aslsP', 0.001, ...
    'aslsNIter', 15 );

%% LOAD INDIVIDUAL SPECTRA

D = collectShakeflaskSpectra(baseDir, wantedFlask, power_mw, t_ms, shiftWin);

if isempty(D)
    error('No spectra found for shakeflask%d.', wantedFlask);
end

fprintf('Found %d timepoints for shakeflask%d:\n', numel(D), wantedFlask);
disp(string({D.timeLabel})');

x = D(1).ramanShift(:);
nTime = numel(D);

%% PLOT RAW MEAN SPECTRA

figure('Name', sprintf('RAW mean spectra - shakeflask%d', wantedFlask));
ax = axes; hold(ax,'on');
C = lines(nTime);

for t = 1:nTime
    mu = mean(D(t).Xraw, 1);
    plot(x, mu, 'Color', C(t,:), 'LineWidth', 1.2);

    [~, ix] = min(abs(x - labelXpos));
    text(x(ix)+5, mu(ix), D(t).timeLabel, ...
        'Color', C(t,:), 'FontSize', 9, 'Interpreter', 'none');
end

grid on;
xlabel('Raman shift (cm^{-1})');
ylabel('Intensity (a.u.)');
title(sprintf('RAW mean spectra - shakeflask%d', wantedFlask));

%% PREPROCESS ALL INDIVIDUAL SPECTRA FOR ALL PIPELINES

Results = struct([]);

for p = 1:numel(pipelines)
    pipe = pipelines{p};

    Results(p).name = pipe.name;
    Results(p).time = struct([]);

    allMetrics = [];
    allProc = [];
    allTimeId = [];

    for t = 1:nTime
        Xraw = D(t).Xraw;

        [Xproc, Xbasecorr, Xbaseline] = preprocessMatrix(Xraw, x, pipe);

        M = computeMetrics(Xraw, Xproc, Xbasecorr, x, noiseWin, analWins);

        Results(p).time(t).timeLabel    = D(t).timeLabel;
        Results(p).time(t).nSpectra     = size(Xraw,1);
        Results(p).time(t).Xraw         = Xraw;
        Results(p).time(t).Xproc        = Xproc;
        Results(p).time(t).Xbasecorr    = Xbasecorr;
        Results(p).time(t).Xbaseline    = Xbaseline;
        Results(p).time(t).metrics      = M;
        Results(p).time(t).meanRaw      = mean(Xraw, 1);
        Results(p).time(t).meanProc     = mean(Xproc, 1);
        Results(p).time(t).stdProc      = std(Xproc, 0, 1);

        allMetrics = [allMetrics; M]; %#ok<AGROW>
        allProc    = [allProc; Xproc]; %#ok<AGROW>
        allTimeId  = [allTimeId; t*ones(size(Xproc,1),1)]; %#ok<AGROW>
    end

    Results(p).allMetrics = allMetrics;
    Results(p).allProc    = allProc;
    Results(p).allTimeId  = allTimeId;
    Results(p).summary    = summarizeMetrics(allMetrics);
end

%% PRINT OVERALL SUMMARY TABLE

fprintf('\n============================================================\n');
fprintf('OVERALL PIPELINE SUMMARY - shakeflask%d\n', wantedFlask);
fprintf('============================================================\n');

summaryCell = cell(numel(pipelines), 8);

for p = 1:numel(pipelines)
    S = Results(p).summary;
    summaryCell{p,1} = Results(p).name;
    summaryCell{p,2} = S.meanNoiseStd;
    summaryCell{p,3} = S.meanPeakToNoise;
    summaryCell{p,4} = S.meanBaselineMean;
    summaryCell{p,5} = S.meanFracNeg;
    summaryCell{p,6} = S.meanRoughness;
    summaryCell{p,7} = S.meanProcDynamicRange;
    summaryCell{p,8} = S.meanReplicateSpread;
end

summaryTable = cell2table(summaryCell, 'VariableNames', ...
    {'Pipeline','NoiseStd','PeakToNoise','BaselineMean','FracNegative','Roughness','DynamicRange','ReplicateSpread'});

disp(summaryTable);

%% PLOT PREPROCESSED MEAN SPECTRA FOR EACH PIPELINE

figure('Name', sprintf('Mean spectra after preprocessing - shakeflask%d', wantedFlask), ...
    'Position', [100 100 1400 900]);

tl = tiledlayout(ceil(numel(pipelines)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for p = 1:numel(pipelines)
    nexttile;
    hold on;

    for t = 1:nTime
        mu = Results(p).time(t).meanProc;
        plot(x, mu, 'Color', C(t,:), 'LineWidth', 1.0);

        [~, ix] = min(abs(x - labelXpos));
        text(x(ix)+5, mu(ix), D(t).timeLabel, ...
            'Color', C(t,:), 'FontSize', 8, 'Interpreter', 'none');
    end

    grid on;
    xlabel('Raman shift (cm^{-1})');
    ylabel('Processed intensity');
    title(Results(p).name, 'Interpreter', 'none');
end

title(tl, sprintf('Preprocessed mean spectra - shakeflask%d', wantedFlask));

%% BOXPLOTS OF METRICS ACROSS PIPELINES

metricNames = {'noiseStd', 'peakToNoise', 'baselineMean', 'fracNegative', 'roughness'};
prettyNames = {'Noise std', 'Peak-to-noise', 'Baseline mean', 'Fraction negative', 'Roughness'};

figure('Name', sprintf('Pipeline metric comparison - shakeflask%d', wantedFlask), ...
    'Position', [100 100 1400 700]);

tl2 = tiledlayout(2,3, 'TileSpacing', 'compact', 'Padding', 'compact');

for m = 1:numel(metricNames)
    nexttile;
    vals = [];
    groups = [];

    for p = 1:numel(pipelines)
        v = Results(p).allMetrics.(metricNames{m});
        vals = [vals; v]; %#ok<AGROW>
        groups = [groups; p*ones(numel(v),1)]; %#ok<AGROW>
    end

    boxplot(vals, groups, 'Labels', strrep({Results.name}, '_', '\_'));
    xtickangle(25);
    ylabel(prettyNames{m});
    title(prettyNames{m});
    grid on;
end

title(tl2, sprintf('Metrics from individual spectra - shakeflask%d', wantedFlask));

%% PER-TIMEPOINT SUMMARY TABLE FOR EACH PIPELINE

fprintf('\n============================================================\n');
fprintf('PER-TIMEPOINT SUMMARIES\n');
fprintf('============================================================\n');

for p = 1:numel(pipelines)
    fprintf('\n--- %s ---\n', Results(p).name);

    timeCell = cell(nTime, 6);
    for t = 1:nTime
        M = Results(p).time(t).metrics;
        timeCell{t,1} = Results(p).time(t).timeLabel;
        timeCell{t,2} = mean(M.noiseStd);
        timeCell{t,3} = mean(M.peakToNoise);
        timeCell{t,4} = mean(M.baselineMean);
        timeCell{t,5} = mean(M.fracNegative);
        timeCell{t,6} = mean(M.roughness);
    end

    Ttime = cell2table(timeCell, 'VariableNames', ...
        {'Timepoint','NoiseStd','PeakToNoise','BaselineMean','FracNegative','Roughness'});
    disp(Ttime);
end

%% BASELINE EXAMPLE PLOT FOR ONE REPRESENTATIVE SINGLE SPECTRUM

% Choose one timepoint roughly in the middle
tMid = ceil(nTime/2);

% Choose one representative spectrum = median total intensity
XrawMid = D(tMid).Xraw;
areaVals = trapz(x, XrawMid, 2);
[~, iRep] = min(abs(areaVals - median(areaVals)));

yRaw = XrawMid(iRep,:);

pipe = pipelines{examplePipelineIdx};
[yProc, yBaseCorr, yBaseline] = preprocessSingle(yRaw, x, pipe);

figure('Name', sprintf('Baseline example - %s - %s', pipe.name, D(tMid).timeLabel), ...
    'Position', [100 100 1200 700]);

subplot(2,1,1);
plot(x, yRaw, 'k', 'LineWidth', 1.0); hold on;
plot(x, yBaseline, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('Raman shift (cm^{-1})');
ylabel('Intensity (a.u.)');
title(sprintf('Raw spectrum + estimated baseline | %s | %s', ...
    D(tMid).timeLabel, pipe.name), 'Interpreter', 'none');
legend({'Raw','Estimated baseline'}, 'Location', 'best');

subplot(2,1,2);
plot(x, yBaseCorr, 'b', 'LineWidth', 1.0); hold on;
plot(x, yProc, 'm', 'LineWidth', 1.0);
grid on;
xlabel('Raman shift (cm^{-1})');
ylabel('Intensity');
title('Baseline-corrected and final preprocessed spectrum');
legend({'Baseline-corrected','Final processed'}, 'Location', 'best');

%% OPTIONAL: PLOT RAW VS PROCESSED FOR ALL TIMEPOINTS USING BEST PIPELINE
% Change "bestPipelineIdx" manually after inspecting summaries
bestPipelineIdx = examplePipelineIdx;

figure('Name', sprintf('Raw vs processed means - %s', Results(bestPipelineIdx).name), ...
    'Position', [100 100 1000 200+120*nTime]);

% tiled layout: nTime rows, 2 columns (left = raw, right = processed)
tl3 = tiledlayout(nTime, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% compute consistent y-limits for each column
rawAll = vertcat(Results(bestPipelineIdx).time.meanRaw);
procAll = vertcat(Results(bestPipelineIdx).time.meanProc);
ylimRaw  = [min(rawAll, [], 'all'),  max(rawAll, [], 'all')];
ylimProc = [min(procAll, [], 'all'),  max(procAll, [], 'all')];

for t = 1:nTime
    % left: raw mean
    nexttile((t-1)*2 + 1);
    plot(x, Results(bestPipelineIdx).time(t).meanRaw, 'k', 'LineWidth', 1.0);
    grid on;
    if t == 1
        title('Raw mean');
    end
    if t == nTime
        xlabel('Raman shift (cm^{-1})');
    end
    ylabel(D(t).timeLabel, 'Interpreter', 'none');
    ylim(ylimRaw);

    % right: processed mean
    nexttile((t-1)*2 + 2);
    plot(x, Results(bestPipelineIdx).time(t).meanProc, 'r', 'LineWidth', 1.0);
    grid on;
    if t == 1
        title('Processed mean');
    end
    if t == nTime
        xlabel('Raman shift (cm^{-1})');
    end
    ylim(ylimProc);
end

% link x-axes across all tiles for synchronized zoom/pan
ax = findall(gcf,'Type','axes');
linkaxes(ax, 'x');
sgtitle(sprintf('Raw vs processed mean spectra - %s', Results(bestPipelineIdx).name), 'Interpreter', 'none');

%% ========================= LOCAL FUNCTIONS =========================

function D = collectShakeflaskSpectra(baseDir, wantedFlask, power_mw, t_ms, shiftWin)

    dateDirs = dir(baseDir);
    dateDirs = dateDirs([dateDirs.isdir]);
    dateDirs = dateDirs(~ismember({dateDirs.name}, {'.','..'}));

    D = struct([]);
    count = 0;

    for i = 1:numel(dateDirs)
        dateName = dateDirs(i).name;
        sfDir = fullfile(baseDir, dateName, sprintf('shakeflask%d', wantedFlask));

        if ~isfolder(sfDir)
            continue;
        end

        files = dir(fullfile(sfDir, '*.spc'));
        if isempty(files)
            continue;
        end

        keep = false(numel(files),1);
        for k = 1:numel(files)
            nm = files(k).name;
            hasPower = contains(lower(nm), sprintf('%dmw', power_mw));
            hasTime  = contains(lower(nm), sprintf('%dms', t_ms));
            keep(k) = hasPower && hasTime;
        end
        files = files(keep);

        if isempty(files)
            continue;
        end

        Xraw = [];
        xRef = [];

        for k = 1:numel(files)
            f = fullfile(files(k).folder, files(k).name);
            s = tgspcread(f);

            x = s.X(:);
            y = s.Y(:);

            mask = x >= shiftWin(1) & x <= shiftWin(2);
            x = x(mask);
            y = y(mask);

            % mild spike suppression before comparison
            y = hampel(y, 11, 3);

            if isempty(xRef)
                xRef = x;
                Xraw = zeros(numel(files), numel(xRef));
            else
                if numel(x) ~= numel(xRef) || any(abs(x - xRef) > 1e-6)
                    error('Spectra do not share identical Raman shift axis in %s', sfDir);
                end
            end

            Xraw(k,:) = y(:).';
        end

        count = count + 1;
        D(count).timeLabel   = dateName;
        D(count).ramanShift  = xRef;
        D(count).Xraw        = Xraw;
        D(count).nFiles      = size(Xraw,1);
    end

    % sort by date label if parseable as yy.mm_dd or similar text order
    if ~isempty(D)
        [~, ord] = sort({D.timeLabel});
        D = D(ord);
    end
end

function [Xproc, Xbasecorr, Xbaseline] = preprocessMatrix(Xraw, x, pipe)

    n = size(Xraw,1);
    Xbasecorr = zeros(size(Xraw));
    Xproc     = zeros(size(Xraw));
    Xbaseline = zeros(size(Xraw));

    for i = 1:n
        [Xproc(i,:), Xbasecorr(i,:), Xbaseline(i,:)] = preprocessSingle(Xraw(i,:), x, pipe);
    end
end

function [yProc, yBaseCorr, yBaseline] = preprocessSingle(yRaw, x, pipe)

    y = yRaw(:);

    switch lower(pipe.baseline)
        case 'msbackadj'
            yAdj = msbackadj(x(:), y, ...
                'WindowSize', 160, ...
                'StepSize', 70, ...
                'EstimationMethod', 'quantile', ...
                'QuantileValue', 0.10);
            yBaseline = y - yAdj;
            yBaseCorr = yAdj;

        case 'asls'
            yBaseline = aslsBaseline(y, pipe.aslsLambda, pipe.aslsP, pipe.aslsNIter);
            yBaseCorr = y - yBaseline;

        otherwise
            error('Unknown baseline method: %s', pipe.baseline);
    end

    yWork = yBaseCorr;

    % --- Savitzky-Golay smoothing (use sgolayfilt derivative signature if supported) ---
    
% --- SG smoothing / derivative using sgolay coefficients (robust) ---
frame = round(pipe.sgolayFrame);
npts  = numel(yWork);

% clamp frame to signal length and ensure odd and > order
frame = min(frame, npts);
if mod(frame,2) == 0, frame = frame - 1; end
if frame <= pipe.sgolayOrder
    frame = pipe.sgolayOrder + 1;
    if mod(frame,2) == 0, frame = frame + 1; end
end

if frame > npts
    % too short: skip
else
    % compute SG matrix/coefficients
    k = pipe.sgolayOrder;
    % sgolay returns matrix g of size frame x (k+1)
    g = sgolay(k, frame);

    if pipe.sgolayDeriv == 0
        % smoothing only: use central convolution vector (column 1)
        coeff = g(:,1) / sum(g(:,1)); % normalize (sgolay returns polynomial basis; normalization keeps scale)
        yWork = conv(yWork, coeff(end:-1:1), 'same'); % flip for conv alignment
    else
        der = max(0, min(round(pipe.sgolayDeriv), frame-1));
        dx = mean(diff(x));
        if dx <= 0 || ~isfinite(dx)
            error('Invalid dx computed from x: %g', dx);
        end
        % SG derivative coefficients: factorial scaling / dx^der
        % g(:, der+1) gives polynomial coefficients for derivative order 'der'
        coeff = g(:, der+1) * factorial(der) / (dx^der);
        % convolution (flip for alignment)
        yWork = conv(yWork, coeff(end:-1:1), 'same');
    end
end

    % Normalization
    if pipe.useSNV
        mu = mean(yWork);
        sd = std(yWork);
        yWork = (yWork - mu) ./ max(sd, eps);
    elseif pipe.useVectorNorm
        nv = norm(yWork);
        yWork = yWork ./ max(nv, eps);
    end

    yProc = yWork(:).';
    yBaseCorr = yBaseCorr(:).';
    yBaseline = yBaseline(:).';
end

function z = aslsBaseline(y, lambda, p, nIter)
    % Asymmetric least squares baseline
    y = y(:);
    n = numel(y);

    D = diff(speye(n), 2);
    w = ones(n,1);

    for it = 1:nIter
        W = spdiags(w, 0, n, n);
        z = (W + lambda*(D'*D)) \ (w .* y);
        w = p * (y > z) + (1-p) * (y < z);
    end

    z = full(z);
end

function M = computeMetrics(Xraw, Xproc, Xbasecorr, x, noiseWin, analWins)

    n = size(Xraw,1);

    maskN = x >= noiseWin(1) & x <= noiseWin(2);
    if ~any(maskN)
        error('noiseWin is outside the Raman shift range.');
    end

    maskA = false(size(x));
    for j = 1:size(analWins,1)
        maskA = maskA | (x >= analWins(j,1) & x <= analWins(j,2));
    end

    M.noiseStd         = zeros(n,1);
    M.peakToNoise      = zeros(n,1);
    M.baselineMean     = zeros(n,1);
    M.fracNegative     = zeros(n,1);
    M.roughness        = zeros(n,1);
    M.procDynamicRange = zeros(n,1);
    M.replicateSpread  = zeros(n,1);

    % replicate spread = median std across wavelengths within same timepoint
    replicateSpreadTime = median(std(Xproc, 0, 1));
    M.replicateSpread(:) = replicateSpreadTime;

    for i = 1:n
        yProc = Xproc(i,:);
        yBase = Xbasecorr(i,:);

        noiseStd = std(yBase(maskN));
        peakSig  = max(abs(yProc(maskA)));

        M.noiseStd(i)         = noiseStd;
        M.peakToNoise(i)      = peakSig / max(noiseStd, eps);
        M.baselineMean(i)     = mean(yBase(maskN));
        M.fracNegative(i)     = mean(yBase(maskN) < 0);
        M.roughness(i)        = std(diff(yProc));
        M.procDynamicRange(i) = max(yProc) - min(yProc);
    end
end


function S = summarizeMetrics(M)
    % Helper: collect numeric values from struct-array field and ignore NaNs.
    getFieldVec = @(s, fname) reshape([s.(fname)], 1, []); % row vector

    safeMean = @(v) ...
        ( isempty(v(~isnan(v))) * NaN + ~isempty(v(~isnan(v))) * mean(v(~isnan(v))) );

    v = getFieldVec(M, 'noiseStd');
    S.meanNoiseStd = safeMean(v);

    v = getFieldVec(M, 'peakToNoise');
    S.meanPeakToNoise = safeMean(v);

    v = getFieldVec(M, 'baselineMean');
    S.meanBaselineMean = safeMean(v);

    v = getFieldVec(M, 'fracNegative');
    S.meanFracNeg = safeMean(v);

    v = getFieldVec(M, 'roughness');
    S.meanRoughness = safeMean(v);

    v = getFieldVec(M, 'procDynamicRange');
    S.meanProcDynamicRange = safeMean(v);

    v = getFieldVec(M, 'replicateSpread');
    S.meanReplicateSpread = safeMean(v);
end