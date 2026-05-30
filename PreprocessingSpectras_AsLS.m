
function [X_processed, ramanShift, X_raw] = PreprocessingSpectras_AsLS(inputData, ramanShiftOrRange, useSNV)
% Flexible preprocessing of Raman spectra: AsLS baseline -> SNV (opt) or mean-centering -> SG
% Usage:
%  [X_processed, ramanShift, X_raw] = PreprocessingSpectras2_AsLS(spectrumFiles, shiftRange, useSNV)
%  [X_processed, ramanShift, X_raw] = PreprocessingSpectras2_AsLS(Xin, ramanShift, useSNV)
%
% Default: useSNV = true

    if nargin < 3 || isempty(useSNV)
        useSNV = true;
    end

    % Input handling 
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

    else
        X = inputData;
        ramanShift = ramanShiftOrRange(:);
        X_raw = X;
    end

    x = ramanShift(:);
    nSpec = size(X,1);

    % AsLS baseline correction parameters (change as needed)
    lambda = 1e4;    % smoothness (increase -> smoother baseline)
    p      = 1e-4;  % asymmetry (small -> baseline stays below peaks)
    nIter  = 10;     % IRLS iterations

    % Apply AsLS sample-by-sample
    for i = 1:nSpec
        yi = X(i, :).';
        % if there are NaNs, interpolate linearly before baseline
        if any(isnan(yi))
            xi = (1:numel(yi)).';
            yi = fillmissing(yi, 'linear');
        end
        z = asls_baseline(yi, lambda, p, nIter);
        X(i, :) = (yi - z).';
    end

    % Savitzky-Golay smoothing
    % SG parameters (change as needed)
    sgOrder = 2;
    sgFrame = 9; % must be odd and > sgOrder
    % apply filter along dimension 2 (rows are spectra)
    X = sgolayfilt(X, sgOrder, sgFrame, [], 2);

    % --- SNV (optional) ---
    if useSNV
        mu = mean(X, 2);
        sd = std(X, 0, 2);
        X = (X - mu) ./ max(sd, eps);

    else
        X = X - mean(X, 2);
    end

    X_processed = X;
end

% -------------------------------------------------------------------------
function z = asls_baseline(y, lambda, p, nIter)
% ASLS_BASELINE   Asymmetric least squares baseline estimation
%  z = asls_baseline(y, lambda, p, nIter)
%  y: column vector, lambda>0, p in (0,1), nIter positive integer

    y = y(:);
    n = numel(y);
    if n < 3
        z = y;
        return;
    end
    if nargin < 4, nIter = 10; end

    % second-difference operator D (sparse)
    e = ones(n,1);
    D = spdiags([e -2*e e], 0:2, n-2, n);
    B = D' * D; % (n x n) sparse

    w = ones(n,1);
    for it = 1:nIter
        W = spdiags(w, 0, n, n);
        A = W + lambda * B;
        % solve A * z = W * y  (sparse backslash)
        z = A \ (w .* y);
        % update weights: smaller weight for positive residuals (peaks)
        r = y - z;
        w = p * (r > 0) + (1 - p) * (r <= 0);
        % ensure positive weights
        w = max(w, eps);
    end
end
