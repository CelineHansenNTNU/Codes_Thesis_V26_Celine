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
% Celine Hansen + ChatGPT, 2026

clear; clc; close all;

%% USER INPUT

baseDir      = fullfile(pwd, 'Shakeflask/');
wantedFlask  = 1;              % <<< choose shakeflask here
%wantedGroup = [2 3 4 5 6 7 8 9 10 11 12]; % use this for intermediate
power_mw     = 450;
t_ms         = 1000;
shiftWin     = [0 3400];

% IMPORTANT:
% Old noise window [2400 2600] does NOT fit when shiftWin = [400 1800].
% Choose a quieter region inside the cropped range.
noiseWin     = [680 710];

% Analysis windows where real Raman peaks are expected
analWins     = [ ...
    1050 1400;   % sugar / fingerprint region
    900 1050;  % strong peak region
    1550 1700  % broad high-wavenumber feature in your data
    ];

% Which pipeline to use for detailed baseline example plot
exampleMethodIdx = 1;

% Plot labels for mean spectra
labelXpos = 1180;

% Color gradient for spectra over time
% Two-color gradient: first timepoint -> last timepoint.
% Examples:
%   blue to red:     [0.15 0.45 0.90] -> [0.90 0.15 0.10]
%   light to dark:   [0.75 0.86 1.00] -> [0.00 0.20 0.70]
%   purple to green: [0.45 0.20 0.75] -> [0.10 0.65 0.35]
gradientStartColor = [0.45 0.20 0.75];  % first timepoint, RGB in [0,1]
gradientEndColor   = [0.10 0.65 0.35];  % last timepoint, RGB in [0,1]
showAllColorbarTicks = true;           % true = tick for every timepoint; false = only first/last

%% Choose preprocessing method to plot
% raw = no preprocessing
% asls = use PreprocessingSpectras_AsLS.m
% msbackadj = use PreprocessingSpectras2.m
preprocMethods = ["asls", "msbackadj"]; % this plots all
useSNV = false ; % or true


%% LOAD INDIVIDUAL SPECTRA

D = collectShakeflaskSpectra(baseDir, wantedFlask, power_mw, t_ms,shiftWin); % for SF data
%D = collectIntermediateSpectra(baseDir, wantedFlask, power_mw, t_ms, shiftWin); % for intermediate data

if isempty(D)
    error('No spectra found for shakeflask%d.', wantedFlask);
end

fprintf('Found %d timepoints for shakeflask%d:\n', numel(D), wantedFlask);
disp(string({D.timeLabel})');

x = D(1).ramanShift(:);
nTime = numel(D);

%% PLOT RAW MEAN SPECTRA

figure('Name', sprintf('RAW mean spectra - Exp 1 shakeflask%d', wantedFlask));
ax = axes; hold(ax,'on');
% Make one-color gradient from first to last timepoint
C = timeGradientColors_local(nTime, gradientStartColor, gradientEndColor);

hRaw = gobjects(nTime,1);

for t = 1:nTime
    mu = mean(D(t).Xraw, 1);
    hRaw(t) = plot(x, mu, 'Color', C(t,:), 'LineWidth', 1.2);
end

applyTimeColorbar_local(ax, C, string({D.timeLabel}), showAllColorbarTicks);

grid on;
xlabel('Raman shift (cm^{-1})');
ylabel('Intensity (a.u.)');
title(sprintf('RAW mean spectra - %d', wantedFlask));

%% Preprocess all individual spectra for selected methods
Results = struct([]);
for p = 1:numel(preprocMethods)
    method = lower(string(preprocMethods(p)));
    Results(p).name = char(method);
    Results(p).time = struct([]);
    allMetrics = [];
    allProc = [];
    allTimeId = [];
    for t = 1:nTime
        Xraw = D(t).Xraw;
        [Xproc, Xbasecorr, Xbaseline] = preprocessMatrix_usingOwnScripts(Xraw, x, method, useSNV);
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
fprintf('OVERALL PIPELINE SUMMARY - %d\n', wantedFlask);
fprintf('============================================================\n');

summaryCell = cell(numel(Results), 8);

for p = 1:numel(Results)
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

figure('Name', sprintf('Mean spectra after preprocessing - Exp 1 shakeflask%d', wantedFlask), ...
    'Position', [100 100 1400 900]);

tl = tiledlayout(ceil(numel(Results)/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for p = 1:numel(Results)
    axp = nexttile;
    hold(axp, 'on');

    hProc = gobjects(nTime,1);

    for t = 1:nTime
        mu = Results(p).time(t).meanProc;
        hProc(t) = plot(x, mu, 'Color', C(t,:), 'LineWidth', 1.0);
    end

    grid on;
    xlabel('Raman shift (cm^{-1})');
    ylabel('Processed intensity');
    title(Results(p).name, 'Interpreter', 'none');

    applyTimeColorbar_local(axp, C, string({D.timeLabel}), showAllColorbarTicks);
end

%title(tl, sprintf('Preprocessed mean spectra - Exp 1 shakeflask%d', wantedFlask));

%% PER-TIMEPOINT SUMMARY TABLE FOR EACH PIPELINE

fprintf('\n============================================================\n');
fprintf('PER-TIMEPOINT SUMMARIES\n');
fprintf('============================================================\n');

for p = 1:numel(Results)
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


%% ========================= LOCAL FUNCTIONS =========================


function C = timeGradientColors_local(nTime, startColor, endColor)
    % Create a two-color gradient from first to last timepoint.
    startColor = reshape(startColor, 1, 3);
    endColor   = reshape(endColor,   1, 3);

    if nTime == 1
        C = startColor;
        return;
    end

    a = linspace(0, 1, nTime).';
    C = (1-a).*startColor + a.*endColor;
    C = max(0, min(1, C));
end

function applyTimeColorbar_local(ax, C, timeLabels, showAllTicks)
    % Attach a colorbar that maps shade to timepoint order.
    % The plotted lines use explicit RGB values, so the colorbar is set manually
    % with the same colormap and CLim values.
    colormap(ax, C);
    clim(ax, [1 size(C,1)]);
    cb = colorbar(ax);
    cb.Label.String = 'Timepoint';
    cb.Label.Interpreter = 'none';

    nTime = numel(timeLabels);
    if showAllTicks || nTime <= 8
        cb.Ticks = 1:nTime;
        cb.TickLabels = cellstr(timeLabels);
    else
        cb.Ticks = [1 nTime];
        cb.TickLabels = {char(timeLabels(1)), char(timeLabels(end))};
    end
    cb.TickLabelInterpreter = 'none';
end

function D = collectShakeflaskSpectra(baseDir, wantedGroup, power_mw, t_ms, shiftWin)

    dateDirs = dir(baseDir);
    dateDirs = dateDirs([dateDirs.isdir]);
    dateDirs = dateDirs(~ismember({dateDirs.name}, {'.','..'}));

    D = struct([]);
    count = 0;

    for i = 1:numel(dateDirs)
        dateName = dateDirs(i).name;
        sfDir = fullfile(baseDir, dateName, sprintf('shakeflask%d', wantedGroup));

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

    % sort timepoints numerically
    if ~isempty(D)
        timeNums = nan(numel(D),1);

        for ii = 1:numel(D)
            lbl = string(D(ii).timeLabel);
    
            % Finds number after T, e.g. T24 -> 24
            tok = regexp(lbl, 'T(\d+)', 'tokens', 'once');
    
            if ~isempty(tok)
                timeNums(ii) = str2double(tok{1});
            else
                % fallback if folder name is only a number or contains another number
                tok2 = regexp(lbl, '\d+', 'match', 'once');
                if ~isempty(tok2)
                    timeNums(ii) = str2double(tok2);
                end
            end
        end
        [~, ord] = sort(timeNums, 'ascend');
        D = D(ord);
    end
end

function [Xproc, Xbasecorr, Xbaseline] = preprocessMatrix_usingOwnScripts(Xraw, x, method, useSNV)

    method = lower(string(method));

    switch method

        case "raw"
            Xproc = Xraw;
            Xbasecorr = Xraw;
            Xbaseline = zeros(size(Xraw));

        case "asls"
            Xproc = PreprocessingSpectras_AsLS(Xraw, x, useSNV);
            Xbasecorr = Xproc;
            Xbaseline = Xraw - Xproc;

        case "msbackadj"
            Xproc = PreprocessingSpectras2(Xraw, x, useSNV);
            Xbasecorr = Xproc;
            Xbaseline = Xraw - Xproc;

        otherwise
            error('Unknown preprocessing method: %s', method);
    end

    if isvector(Xproc)
        Xproc = Xproc(:).';
    end

    if isvector(Xbasecorr)
        Xbasecorr = Xbasecorr(:).';
    end

    if isvector(Xbaseline)
        Xbaseline = Xbaseline(:).';
    end
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

% function D = collectIntermediateSpectra(baseDir, wantedGroup, power_mw, t_ms, shiftWin)
% 
%     [files, G] = collect_spc_files_by_group_local(baseDir, power_mw, t_ms);
% 
%     keep = G == wantedGroup;
%     files = files(keep);
% 
%     if isempty(files)
%         D = struct([]);
%         return;
%     end
% 
%     Xraw = [];
%     xRef = [];
% 
%     for k = 1:numel(files)
%         f = files{k};
%         s = tgspcread(f);
% 
%         x = s.X(:);
%         y = s.Y(:);
% 
%         mask = x >= shiftWin(1) & x <= shiftWin(2);
%         x = x(mask);
%         y = y(mask);
% 
%         if isempty(xRef)
%             xRef = x;
%             Xraw = zeros(numel(files), numel(xRef));
%         else
%             if numel(x) ~= numel(xRef) || any(abs(x - xRef) > 1e-6)
%                 error('Spectra do not share identical Raman shift axis in group %d.', wantedGroup);
%             end
%         end
% 
%         Xraw(k,:) = y(:).';
%     end
% 
%     D = struct([]);
%     D(1).timeLabel  = sprintf('Group %d', wantedGroup);
%     D(1).ramanShift = xRef;
%     D(1).Xraw       = Xraw;
%     D(1).nFiles     = size(Xraw,1);
% end

function [files, G] = collect_spc_files_by_group_local(baseDir, power_mw, t_ms)

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