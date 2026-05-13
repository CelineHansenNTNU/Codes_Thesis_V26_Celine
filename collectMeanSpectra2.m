function S = collectMeanSpectra2(baseDirs, sugars, power_mw, t_ms, shiftRange)
% Function for collecting mean Raman spectrum of replicates
% Adapted for folder structure:
%   Shakeflask/dateFolder/shakeflaskX/*.spc
%
% Returns one mean spectrum per shakeflaskX per date folder.
%
% Inputs kept compatible with old script:
%   baseDirs   - cell array, should now typically contain {'.../Shakeflask'}
%   sugars     - not really used anymore, kept for compatibility
%   power_mw   - laser power used in filenames
%   t_ms       - acquisition time used in filenames
%   shiftRange - [xmin xmax]
%
% Output:
%   S(k) with fields:
%       tag
%       baseDir
%       sub
%       ramanShift
%       meanSpectrum
%       nFiles

    S = struct('tag',{},'baseDir',{},'sub',{},'ramanShift',{},'meanSpectrum',{},'nFiles',{});

    % Robust filename patterns
    patList = { ...
        sprintf('*_%dmw_%dms_*.spc',  power_mw, t_ms), ...
        sprintf('*_%dmW_%dms_*.spc',  power_mw, t_ms), ...
        sprintf('*_%dmw_%ds_*.spc',   power_mw, t_ms), ...
        sprintf('*_%dmW_%ds_*.spc',   power_mw, t_ms), ...
        sprintf('*_%dmw_%05dms_*.spc',power_mw, t_ms), ...
        sprintf('*_%dmW_%05dms_*.spc',power_mw, t_ms), ...
        sprintf('*_%dmw_%05ds_*.spc', power_mw, t_ms), ...
        sprintf('*_%dmW_%05ds_*.spc', power_mw, t_ms) ...
    };

    kOut = 0;

    for b = 1:numel(baseDirs)
        rootDir = baseDirs{b};

        % Find first-level folders, expected to be date folders
        D1 = dir(rootDir);
        D1 = D1([D1.isdir]);
        D1 = D1(~ismember({D1.name},{'.','..'}));

        for iDate = 1:numel(D1)
            dateFolder = D1(iDate).name;
            datePath   = fullfile(rootDir, dateFolder);

            % Look for shakeflask folders inside each date folder
            D2 = dir(datePath);
            D2 = D2([D2.isdir]);
            D2 = D2(~ismember({D2.name},{'.','..'}));

            subNames = {D2.name};
            use = startsWith(lower(subNames), 'shakeflask');
            flaskSubs = subNames(use);

            if isempty(flaskSubs)
                continue;
            end

            for iSub = 1:numel(flaskSubs)
                sub = flaskSubs{iSub};
                subPath = fullfile(datePath, sub);

                files = {};

                for pi = 1:numel(patList)
                    tmp = dir(fullfile(subPath, patList{pi}));
                    if ~isempty(tmp)
                        files = [files, fullfile({tmp.folder}, {tmp.name})]; %#ok<AGROW>
                    end
                end

                files = unique(files);

                if isempty(files)
                    warning('Found no files for %s / %s (power=%dmW, t=%dms)', ...
                        dateFolder, sub, power_mw, t_ms);
                    continue;
                end

                % Existing routine for reading, deglitching and interpolating
                % Srep = plotReplicateStatsNew(files, shiftRange, false);
                Srep = plotReplicateStatsNew(files, shiftRange, false, "native");
                x   = Srep.ramanShift(:);
                YI  = Srep.YI;
                mu  = mean(YI, 2);

                kOut = kOut + 1;
                S(kOut).tag          = sprintf('%s | %s', dateFolder, sub);
                S(kOut).baseDir      = datePath;
                S(kOut).sub          = sprintf('%s_%s', dateFolder, sub);
                S(kOut).ramanShift   = x;
                S(kOut).meanSpectrum = mu(:).';
                S(kOut).nFiles       = size(YI, 2);
            end
        end
    end

    if isempty(S)
        warning('No spectra were collected. Check folder structure and filename patterns.');
    end
end