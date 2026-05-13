% Script for checking differences between original Raman spectra and
% preprocessed spectra for ONE chosen shakeflask across all timepoints
% Adapted for new folder structure:
%   Shakeflask/dateFolder/shakeflaskX/*.spc
%
% Needs access to 
%   - collectMeanSpectra2.m
%   - PreprocessingSpectras2.m
%   - PreprocessingSpectras_AsLS.m
%
% Celine Hansen + Copilot 

clear; clc; close all;

%% USER INPUT

baseDirs = {fullfile(pwd, 'Raman_mixed_new_samples/88/')};
sugars   = {'88'};   % kept only for compatibility with collectMeanSpectra

power_mw = 450;              % laser power in mW
t_ms     = 1000;             % integration time in ms
shiftWin = [400 1800];         % Raman shift window

wantedFlask = 9;             % CHOOSE SHAKEFLASK HERE

% Evaluation windows
noiseWin = [680 710];      % region assumed to contain no Raman peaks
analWin  = [1000 1300];        % region where peaks of interest are present

%% COLLECT ALL MEAN SPECTRA FROM ALL DATES / SHAKEFLASKS

S = collectMeanSpectra2(baseDirs, sugars, power_mw, t_ms, shiftWin);

if isempty(S)
    error('No spectra were collected. Check folder structure or filename patterns.');
end

%% EXTRACT DATE + SHAKEFLASK NUMBER FROM TAG

flaskNum  = nan(numel(S),1);
timeNum   = nan(numel(S),1);
timeLabel = strings(numel(S),1);

for k = 1:numel(S)
    tok = regexp(S(k).tag, '(T\d+)\s*\|\s*shakeflask(\d+)', 'tokens', 'once', 'ignorecase');

    if ~isempty(tok)
        timeLabel(k) = string(tok{1});                 % e.g. T0
        flaskNum(k)  = str2double(tok{2});            % e.g. 1

        % convert 'T48' -> 48
        tclean = regexprep(tok{1}, '[Tt]', '');
        timeNum(k) = str2double(tclean);
    else
        warning('Could not parse tag: %s', S(k).tag);
    end
end

%% KEEP ONLY THE CHOSEN SHAKEFLASK

idx = find(flaskNum == wantedFlask);

if isempty(idx)
    error('No entries found for shakeflask%d.', wantedFlask);
end

[~, ord] = sort(timeNum(idx));
idx = idx(ord);

Ssel       = S(idx);
timeSel    = timeLabel(idx);
timeNumSel = timeNum(idx);

fprintf('Found %d timepoints for shakeflask%d:\n', numel(Ssel), wantedFlask);
disp(timeSel);

%% PLOT ORIGINAL MEAN SPECTRA IN ONE FIGURE

figure('Name', sprintf('Original mean spectra - shakeflask%d', wantedFlask));
ax1 = axes;
hold(ax1, 'on');

colors = lines(numel(Ssel));

for k = 1:numel(Ssel)
    x  = Ssel(k).ramanShift(:);
    mu = Ssel(k).meanSpectrum(:);

    plot(ax1, x, mu, ...
        'Color', colors(k,:), ...
        'LineWidth', 1.0);

    % Put date label near the right end of the curve
    xTarget = 800;
    [~, idxText] = min(abs(x - xTarget));
    xText = x(idxText);
    yText = mu(idxText);

    text(ax1, xText + 8, yText, char(timeSel(k)), ...
        'Color', colors(k,:), ...
        'FontSize', 9, ...
        'Interpreter', 'none', ...
        'VerticalAlignment', 'middle');
end

grid(ax1, 'on');
xlabel(ax1, 'Raman shift (cm^{-1})');
ylabel(ax1, 'Intensity (a.u.)');
title(ax1, sprintf('Original mean spectra - shakeflask%d', wantedFlask));

% Expand x-limits a bit so text fits on the right
xlim(ax1, [min(Ssel(1).ramanShift) max(Ssel(1).ramanShift)+80]);

hold(ax1, 'off');

%% PLOT PREPROCESSED MEAN SPECTRA IN ONE FIGURE

figure('Name', sprintf('Preprocessed mean spectra - shakeflask%d', wantedFlask));
ax2 = axes;
hold(ax2, 'on');

for k = 1:numel(Ssel)
    x  = Ssel(k).ramanShift(:);
    mu = Ssel(k).meanSpectrum(:).';

    % Applying full preprocessing
    muProc = PreprocessingSpectras_AsLS(mu, x, false);

    plot(ax2, x, muProc(:), ...
        'Color', colors(k,:), ...
        'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s (n=%d)', timeSel(k), Ssel(k).nFiles));
end

grid(ax2, 'on');
xlabel(ax2, 'Raman shift (cm^{-1})');
ylabel(ax2, 'Intensity (a.u.)');
title(ax2, sprintf('Preprocessed mean spectra - shakeflask%d', wantedFlask));
legend(ax2, 'show', 'Location', 'eastoutside', 'Interpreter', 'none');
hold(ax2, 'off');

%% QUALITY METRICS FOR EACH TIMEPOINT

fprintf('\n============================================================\n');
fprintf('Quality metrics for shakeflask%d\n', wantedFlask);
fprintf('============================================================\n');

for k = 1:numel(Ssel)

    x  = Ssel(k).ramanShift(:);
    mu = Ssel(k).meanSpectrum(:).';

    % Full preprocessing with SNV
    muProc = PreprocessingSpectras_AsLS(mu, x, false);

    % Baseline-only version (no SNV) for evaluating baseline/noise
    muBase = PreprocessingSpectras_AsLS(mu, x, false);

    % Masks for windows
    maskN = x >= noiseWin(1) & x <= noiseWin(2);
    maskA = x >= analWin(1)  & x <= analWin(2);

    if ~any(maskN)
        warning('Noise window is outside Raman shift range for %s', timeSel(k));
        continue;
    end

    if ~any(maskA)
        warning('Analysis window is outside Raman shift range for %s', timeSel(k));
        continue;
    end

    % Baseline statistics in noise window
    baselineMean = mean(muBase(maskN));
    fracNeg      = mean(muBase(maskN) < 0);

    % Noise reduction factor
    NRF = std(mu(maskN)) / max(std(muBase(maskN)), eps);

    % Peak prominence / peak-to-noise
    rawDet  = mu     - movmedian(mu,     301, 'omitnan');
    baseDet = muBase - movmedian(muBase, 301, 'omitnan');

    [~,~,~,promRaw]  = findpeaks(rawDet(maskA),  x(maskA));
    [~,~,~,promBase] = findpeaks(baseDet(maskA), x(maskA));

    if isempty(promRaw)
        PNR_raw = NaN;
    else
        PNR_raw = max(promRaw) / max(std(mu(maskN)), eps);
    end

    if isempty(promBase)
        PNR_pre = NaN;
    else
        PNR_pre = max(promBase) / max(std(muBase(maskN)), eps);
    end

    gainPNR = PNR_pre ./ max(PNR_raw, eps);

    % Band power and noise power
    xA = x(maskA);
    sig_raw  = mu(maskA);
    sig_base = muBase(maskA);

    noise_raw  = mu(maskN);
    noise_base = muBase(maskN);

    bandPower_raw_time  = trapz(xA, abs(sig_raw));
    bandPower_base_time = trapz(xA, abs(sig_base));

    noisePow_raw  = var(noise_raw);
    noisePow_base = var(noise_base);

    SNR_lin_raw  = bandPower_raw_time  / (noisePow_raw  + eps);
    SNR_lin_base = bandPower_base_time / (noisePow_base + eps);

    SNR_dB_raw  = 10 * log10(SNR_lin_raw);
    SNR_dB_base = 10 * log10(SNR_lin_base);

    % Peak-based SNR from baseline-corrected spectrum
    sigPeak     = max(muBase(maskA));
    noiseStd    = std(muBase(maskN));
    SNR_peak    = sigPeak / (noiseStd + eps);
    SNR_peak_dB = 20 * log10(SNR_peak);

    % Print results
    fprintf('\n[%s]\n', timeSel(k));
    fprintf('  NRF(noise)       = %.3f\n', NRF);
    fprintf('  PNR raw          = %.3f\n', PNR_raw);
    fprintf('  PNR pre          = %.3f\n', PNR_pre);
    fprintf('  PNR gain         = %.3f\n', gainPNR);
    fprintf('  Baseline mean    = %.3f\n', baselineMean);
    fprintf('  Fraction negative= %.3f\n', fracNeg);
    fprintf('  BandPower raw    = %.3g\n', bandPower_raw_time);
    fprintf('  BandPower base   = %.3g\n', bandPower_base_time);
    fprintf('  NoisePow raw     = %.3g\n', noisePow_raw);
    fprintf('  NoisePow base    = %.3g\n', noisePow_base);
    fprintf('  SNR_lin raw      = %.3g\n', SNR_lin_raw);
    fprintf('  SNR_lin base     = %.3g\n', SNR_lin_base);
    fprintf('  SNR_dB raw       = %.3f dB\n', SNR_dB_raw);
    fprintf('  SNR_dB base      = %.3f dB\n', SNR_dB_base);
    fprintf('  SNR_peak         = %.3f\n', SNR_peak);
    fprintf('  SNR_peak_dB      = %.3f dB\n', SNR_peak_dB);
end

%% OPTIONAL: TILE PLOT WITH ONE PANEL PER TIMEPOINT (ORIGINAL)

figure('Name', sprintf('Original mean spectra by timepoint - shakeflask%d', wantedFlask));
tiledlayout('flow');

for k = 1:numel(Ssel)
    x  = Ssel(k).ramanShift(:);
    mu = Ssel(k).meanSpectrum(:).';

    nexttile;
    plot(x, mu, 'LineWidth', 1.0);
    grid on;
    title(sprintf('%s (n=%d)', timeSel(k), Ssel(k).nFiles), 'Interpreter', 'none');
    xlabel('Raman shift (cm^{-1})');
    ylabel('Intensity (a.u.)');
end

%% OPTIONAL: TILE PLOT WITH ONE PANEL PER TIMEPOINT (PREPROCESSED)

figure('Name', sprintf('Preprocessed mean spectra by timepoint - shakeflask%d', wantedFlask));
tiledlayout('flow');

for k = 1:numel(Ssel)
    x  = Ssel(k).ramanShift(:);
    mu = Ssel(k).meanSpectrum(:).';

    muProc = PreprocessingSpectras_AsLS(mu, x, false);

    nexttile;
    plot(x, muProc(:), 'LineWidth', 1.0);
    grid on;
    title(sprintf('%s (n=%d)', timeSel(k), Ssel(k).nFiles), 'Interpreter', 'none');
    xlabel('Raman shift (cm^{-1})');
    ylabel('Intensity (a.u.)');
end