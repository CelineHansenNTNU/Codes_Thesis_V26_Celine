function [X_processed, ramanShift, X_raw] = PreprocessingSpectras2(inputData, ramanShiftOrRange, useSNV)
    % Flexible preprocessing of Raman spectra
    %
    % Option 1:
    %   [X_processed, ramanShift, X_raw] = PreprocessingSpectras(spectrumFiles, shiftRange, useSNV)
    %   -> reads .spc files, native grid, Hampel filter, preprocessing
    %
    % Option 2:
    %   [X_processed, ramanShift, X_raw] = PreprocessingSpectras(Xin, ramanShift, useSNV)
    %   -> preprocesses already imported spectra matrix
    %
    % Celine Hansen + Copilot, modified 09.03.26

    if nargin < 3 || isempty(useSNV)
        useSNV = true;
    end

    % Case 1: input is file list
    if iscell(inputData)
        spectrumFiles = inputData;
        shiftRange = ramanShiftOrRange;

        S = cellfun(@tgspcread, spectrumFiles, 'uni', 0);
        Xs = cellfun(@(s) s.X(:), S, 'uni', 0);
        Ys = cellfun(@(s) s.Y(:), S, 'uni', 0);

        nSpec = numel(spectrumFiles);

        xmin = max(cellfun(@(x) min(x), Xs));
        xmax = min(cellfun(@(x) max(x), Xs));
        xmin = max(xmin, shiftRange(1));
        xmax = min(xmax, shiftRange(2));

        tol = 1e-6;
        x0 = Xs{1};
        sameGrid = true;
        for i = 2:nSpec
            if numel(Xs{i}) ~= numel(x0) || any(abs(Xs{i} - x0) > tol)
                sameGrid = false;
                break;
            end
        end

        if ~sameGrid
            error(['Native preprocessing requires identical X grids across files. ' ...
                   'Your spectra do not all share the same Raman shift axis.']);
        end

        mask = (x0 >= xmin) & (x0 <= xmax);
        ramanShift = x0(mask);

        X_raw = zeros(nSpec, numel(ramanShift));
        for i = 1:nSpec
            y = hampel(Ys{i}, 15, 3);
            X_raw(i, :) = y(mask).';
        end

        X = X_raw;

    % Case 2: input is already a matrix
    else
        X = inputData;
        ramanShift = ramanShiftOrRange(:);
        X_raw = X;
    end

    x = ramanShift(:);

    % Baseline correction
    for i = 1:size(X, 1)
        yi = X(i, :).';
        yi = msbackadj(x, yi, 'WindowSize', 160, 'StepSize', 75, 'EstimationMethod', 'quantile', 'QuantileValue', 0.10);
        X(i, :) = yi.';
    end

    % SNV
    if useSNV
        mu = mean(X, 2);
        sd = std(X, 0, 2);
        X = (X - mu) ./ max(sd, eps);
    end

    % Savitzky-Golay smoothing
    X = sgolayfilt(X, 2, 9, [], 2);

    X_processed = X;
end