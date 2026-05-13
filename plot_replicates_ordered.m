function out = plot_replicates_ordered(folderPath, power_mw, t_ms, shiftRange, preprocMethod, useSNV)
% Plot all replicate spectra in acquisition order for one folder
% Marks first and last measurements and shows mean +- SD

    if nargin < 5 || isempty(preprocMethod), preprocMethod = 'msbackadj'; end
    if nargin < 6 || isempty(useSNV), useSNV = false; end

    % Find matching files 
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

    files = {};
    for pi = 1:numel(patList)
        tmp = dir(fullfile(folderPath, patList{pi}));
        if ~isempty(tmp)
            files = [files, fullfile({tmp.folder}, {tmp.name})];
        end
    end
    files = unique(files);

    if isempty(files)
        error('No matching .spc files found in %s', folderPath);
    end

    % Sort by acquisition number in filename
    orderNum = nan(numel(files),1);
    for i = 1:numel(files)
        [~, name, ~] = fileparts(files{i});
        % expecting ending like ..._00091
        tok = regexp(name, '_(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            orderNum(i) = str2double(tok{1});
        else
            tok2 = regexp(name, '(\d+)', 'tokens');
            if ~isempty(tok2)
                orderNum(i) = str2double(tok2{end}{1});
            end
        end
    end

    [orderNum, ord] = sort(orderNum);
    files = files(ord);

    % Read and preprocess all spectra
    switch lower(string(preprocMethod))
        case "msbackadj"
            [Xproc, ramanShift, Xraw] = PreprocessingSpectras2(files, shiftRange, useSNV);
        case "asls"
            [Xproc, ramanShift, Xraw] = PreprocessingSpectras_AsLS(files, shiftRange, useSNV);
        case "peaknorm"
            [Xproc, ramanShift, Xraw] = PreprocessingSpectrasNew(files, shiftRange, useSNV);
        otherwise
            error('Unknown preprocMethod: %s', preprocMethod);
    end

    % Xproc and Xraw are nRep x nWavenumbers
    mu_raw = mean(Xraw, 1, 'omitnan');
    sd_raw = std(Xraw, 0, 1, 'omitnan');
    mu_proc = mean(Xproc, 1, 'omitnan');
    sd_proc = std(Xproc, 0, 1, 'omitnan');

    first_raw = Xraw(1,:);
    last_raw = Xraw(end,:);
    first_proc = Xproc(1,:);
    last_proc = Xproc(end,:);

    % distance from each replicate to mean spectrum
    d_raw = sqrt(sum((Xraw - mu_raw).^2, 2));
    d_proc = sqrt(sum((Xproc - mu_proc).^2, 2));

    % print worst replicate outliers
    nShow = min(15, numel(files));
    [~, idxProc] = sort(d_proc, 'descend');
    [~, idxRaw] = sort(d_raw, 'descend');

    fprintf('\nWorst %d replicates by preprocessed distance-to-mean:\n', nShow);
    for k = 1:nShow
        ii = idxProc(k);
        [~, name, ext] = fileparts(files{ii});
        fprintf('%2d) idx = %3d | order = %5d | distProc = %8.2f | %s%s\n', k, ii, orderNum(ii), d_proc(ii), name, ext);
    end

    fprintf('\nWorst %d replicates by raw distance-to-mean:\n', nShow);
    for k = 1:nShow
        ii = idxRaw(k);
        [~, name, ext] = fileparts(files{ii});
        fprintf('%2d) idx = %3d | order = %5d | distProc = %8.2f | %s%s\n', k, ii, orderNum(ii), d_raw(ii), name, ext);
    end

    % Plot raw Raman spectra
    figure('Name', 'Replicates in acquisition order: raw', 'Units', 'normalized', 'Position', [0.08 0.10 0.80 0.70]);
    tiledlayout(2,1,"TileSpacing","compact","Padding","compact");

    nexttile;
    plot(ramanShift, Xraw.', 'Color', [0.75 0.75 0.75]); hold on;
    plot(ramanShift, first_raw, 'b-', 'LineWidth', 1.8, 'DisplayName', 'First');
    plot(ramanShift, last_raw, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Last');
    plot(ramanShift, mu_raw, 'k-', 'LineWidth', 2.0, 'DisplayName', 'Mean');
    ylim(prctile(Xraw(:), [1 99]));
    grid on;
    xlabel('Raman shift (cm^{-1})');
    ylabel('Intensity');
    title(sprintf('Raw replicates | n = %d | %s', size(Xraw,1), folderPath), 'Interpreter','none');
    legend('Location','best');

    nexttile;
    fill([ramanShift; flipud(ramanShift)], [mu_raw(:)-sd_raw(:); flipud(mu_raw(:)+sd_raw(:))], [0.85 0.85 0.85], 'EdgeColor','none'); hold on;
    plot(ramanShift, mu_raw, 'k-', 'LineWidth',2);
    plot(ramanShift, first_raw, 'b-', 'LineWidth', 1.2);
    plot(ramanShift, last_raw, 'r-', 'LineWidth', 1.2);
    grid on;
    xlabel('Raman shift (cm^{-1})');
    ylabel('Intensity');
    title('Raw mean \pm SD, with first and last');
    legend('Mean \pm SD', 'Mean', 'First', 'Last', 'Location','best');

    % Plot preprocessed
    figure('Name', 'Replicates in acquisition order: preprocessed', 'Units','normalized', 'Position',[0.10 0.12 0.80 0.70]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(ramanShift, Xproc.', 'Color', [0.75 0.75 0.75]); hold on;
    plot(ramanShift, first_proc, 'b-', 'LineWidth',1.8, 'DisplayName','First');
    plot(ramanShift, last_proc, 'r-', 'LineWidth', 1.8, 'DisplayName','Last');
    plot(ramanShift, mu_proc, 'k-', 'LineWidth',2.0, 'DisplayName','Mean');
    grid on;
    xlabel('Raman shift (cm^{-1})');
    ylabel('Intensity');
    title(sprintf('Preprocessed replicates | %s | SNV = %d', preprocMethod, useSNV));
    legend('Location','best');

    nexttile;
    fill([ramanShift; flipud(ramanShift)], [mu_proc(:)-sd_proc(:); flipud(mu_proc(:)+sd_proc(:))], [0.85 0.85 0.85], 'EdgeColor','none'); hold on;
    plot(ramanShift, mu_proc, 'k-', 'LineWidth',2);
    plot(ramanShift, first_proc, 'b-', 'LineWidth',1.2);
    plot(ramanShift, last_proc, 'r-', 'LineWidth',1.2);
    grid on;
    xlabel('Raman shift (cm^{-1})');
    ylabel('Intensity');
    title('Preprocessed mean \pm SD, with first and last');
    legend('Mean \pm SD', 'Mean', 'First', 'Last', 'Location','best');

    % Order diagnostics
    figure('Name', 'Replicate order diagnostics', 'Units','normalized', 'Position', [0.18 0.18 0.65 0.45]);
    tiledlayout(1,2, 'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(orderNum, d_raw, '-o'); hold on;
    worstIdx_raw = idxRaw(1:min(5,numel(idxRaw)));
    plot(orderNum(worstIdx_raw), d_raw(worstIdx_raw), 'ro', 'MarkerSize',8, 'LineWidth',1.5); grid on;
    xlabel('Measurement order number');
    ylabel('Distance to mean spectrum');
    title('Raw: distance to mean');

    nexttile;
    plot(orderNum, d_proc, '-o'); hold on;
    worstIdx = idxProc(1:min(5,numel(idxProc)));
    plot(orderNum(worstIdx), d_proc(worstIdx), 'ro', 'MarkerSize',8, 'LineWidth',1.5); grid on;
    xlabel('Measurement order number');
    ylabel('Distance to mean spectrum');
    title('Preprocessed: distance to mean');

    T = table((1:numel(files))', orderNum(:), string(files(:)), d_raw(:), d_proc(:), 'VariableNames', {'IndexInSeries','MeasurementOrder','File','DistRaw','DistProc'});

    T_sorted_proc = sortrows(T, 'DistProc', 'descend');
    T_sorted_raw  = sortrows(T, 'DistRaw',  'descend');

    fprintf('\nTop 10 worst replicates by preprocessed distance:\n');
    disp(T_sorted_proc(1:min(10,height(T_sorted_proc)), :));

    fprintf('\nTop 10 worst replicates by raw distance:\n');
    disp(T_sorted_raw(1:min(10,height(T_sorted_raw)), :));

    % output
    out.files = files;
    out.orderNum = orderNum;
    out.ramanShift = ramanShift;
    out.Xraw = Xraw;
    out.Xproc = Xproc;
    out.meanRaw = mu_raw;
    out.sdRaw = sd_raw;
    out.meanProc = mu_proc;
    out.sdProc = sd_proc;
    out.distRaw = d_raw;
    out.distProc = d_proc;
    out.replicateTable = T;
    out.repTableSortedProc = T_sorted_proc;
    out.repTableSortedRaw = T_sorted_raw;
end
