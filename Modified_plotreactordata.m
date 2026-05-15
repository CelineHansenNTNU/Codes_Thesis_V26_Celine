superfile="/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/bioreactor_package/ramanspectrarun1_perfectrunnoproblems";
%reactors=struct2cell(dir(superfile));
%reactors=reactors(:,3:end);

reactorStruct = dir(superfile);
reactorStruct = reactorStruct([reactorStruct.isdir]);
reactorStruct = reactorStruct(~ismember({reactorStruct.name}, {'.','..'}));

% Keep only reactor folders, e.g. reactor1/reactor2 or reactorA/reactorB
reactorNames = string({reactorStruct.name});
keepReactors = startsWith(lower(reactorNames), "reactor");

reactorStruct = reactorStruct(keepReactors);

% Optional: display which folders are used
disp("Reactor folders used:")
disp(string({reactorStruct.name})')

reactors = struct2cell(reactorStruct);

%read data
data={size(reactors,2),1};
%data = cell(size(reactors,2),1);
for i=1:size(reactors,2)
    %samples=struct2cell(dir(fullfile(reactors{2,i},reactors{1,i})));
    %samples=samples(:,3:end);

    sampleStruct = dir(fullfile(reactors{2,i}, reactors{1,i}));
    sampleStruct = sampleStruct([sampleStruct.isdir]);
    sampleStruct = sampleStruct(~ismember({sampleStruct.name}, {'.','..'}));
    samples = struct2cell(sampleStruct);
    data{i,1}=cell(size(samples,2),3);
    for j=1:size(samples,2)
        %temp=struct2cell(dir(fullfile(samples{2,j},samples{1,j})));
        %temp=temp(:,3:end);

        % ensuring the datetime is a text converted to a string
        %traw = temp{3,1}; 
        %tstr = string(traw);

        tempStruct = dir(fullfile(samples{2,j}, samples{1,j}));
        tempStruct = tempStruct(~ismember({tempStruct.name}, {'.','..'}));
        tempStruct = tempStruct(~[tempStruct.isdir]);
        tempStruct = tempStruct(endsWith({tempStruct.name}, '.spc', 'IgnoreCase', true));
        if isempty(tempStruct)
            warning('No .spc files found in: %s', fullfile(samples{2,j}, samples{1,j}));
            continue;
        end
        temp = struct2cell(tempStruct);
        traw = temp{3,1}; 
        tstr = string(traw);

        % removing trailing dots
        tstr = regexprep(tstr, '([A-Za-z]+)\.', '$1');
        dt = datetime(tstr, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss', 'Locale','nb_NO');
        data{i,1}{j,2}=posixtime(dt);

        fileList = fullfile(temp(2,:), temp(1,:));
        data{i,1}{j,1}=snr_spectrum_from_files_4(fileList,[-400,18000]);
        data{i,1}{j,3}=snr_spectrum_from_files_4(fileList,[400,1800]);
    end
    
    mask=cellfun(@(x) x.mu(35)>2000,data{i,1}(:,1));
    data{i,1}=data{i,1}(mask,:);
    %fprintf('Reactor %d: %d rows before mask\n', i, size(data{i,1},1));
    %disp(cellfun(@isempty, data{i,1}(:,1))')
end


superfile2="/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/bioreactor_package/ramanspectrarun1_perfectrunnoproblemssecondfeeding";
%reactors2=struct2cell(dir(superfile2));
%reactors2=reactors2(:,3:end);

reactorStruct2 = dir(superfile2);
reactorStruct2 = reactorStruct2([reactorStruct2.isdir]);
reactorStruct2 = reactorStruct2(~ismember({reactorStruct2.name}, {'.','..'}));

reactorNames2 = string({reactorStruct2.name});
keepReactors2 = startsWith(lower(reactorNames2), "reactor");

reactorStruct2 = reactorStruct2(keepReactors2);

disp("Reactor folders used second feeding:")
disp(string({reactorStruct2.name})')

reactors2 = struct2cell(reactorStruct2);

%read data
data2={size(reactors2,2),1};
%data2 = cell(size(reactors2,2),1);
for i=1:size(reactors2,2)
    %samples=struct2cell(dir(fullfile(reactors2{2,i},reactors2{1,i})));
    %samples=samples(:,3:end);

    sampleStruct = dir(fullfile(reactors2{2,i}, reactors2{1,i}));
    sampleStruct = sampleStruct([sampleStruct.isdir]);
    sampleStruct = sampleStruct(~ismember({sampleStruct.name}, {'.','..'}));
    samples = struct2cell(sampleStruct);
    data2{i,1}=cell(size(samples,2),3);
    for j=1:size(samples,2)

        tempStruct = dir(fullfile(samples{2,j}, samples{1,j}));
        tempStruct = tempStruct(~ismember({tempStruct.name}, {'.','..'}));
        tempStruct = tempStruct(~[tempStruct.isdir]);
        tempStruct = tempStruct(endsWith({tempStruct.name}, '.spc', 'IgnoreCase', true));
        if isempty(tempStruct)
            warning('No .spc files found in: %s', fullfile(samples{2,j}, samples{1,j}));
            continue;
        end
        temp = struct2cell(tempStruct);
        traw = temp{3,1}; 
        tstr = string(traw);

        % removing trailing dots
        tstr = regexprep(tstr, '([A-Za-z]+)\.', '$1');
        dt = datetime(tstr, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss', 'Locale','nb_NO');
        data2{i,1}{j,2}=posixtime(dt);

        %temp=struct2cell(dir(fullfile(samples{2,j},samples{1,j})));
        %temp=temp(:,3:end);
        %data2{i,1}{j,2}=posixtime(datetime(temp{3,1}));
        data2{i,1}{j,1}=snr_spectrum_from_files_4(fullfile(temp(2,:),temp(1,:)),[-400,18000]);
        data2{i,1}{j,3}=snr_spectrum_from_files_4(fullfile(temp(2,:),temp(1,:)),[400,1800]);
    end
    
    mask=cellfun(@(x) x.mu(35)>2000,data2{i,1}(:,1));
    data2{i,1}=data2{i,1}(mask,:);

    %fprintf('Reactor %d: %d rows before mask\n', i, size(data2{i,1},1));
    %disp(cellfun(@isempty, data2{i,1}(:,1))')
end

data{1,1}=[data{1,1};data2{1,1}];
data{2,1}=[data{2,1};data2{2,1}];
%for i=1:2
%    [~,b]=sort(cell2mat(data{i,1}(:,2)));
%    data{i,1}=data{i,1}(b,:);
%end 
for i = 1:2
    validRows = ~cellfun(@isempty, data{i,1}(:,2));
    data{i,1} = data{i,1}(validRows,:);
    if isempty(data{i,1})
        warning('data{%d,1} is empty after filtering.',i);
        continue;
    end
    [~, b] = sort(cell2mat(data{i,1}(:,2)));
    data{i,1} = data{i,1}(b,:);
end
sugarorder=["glucose", "xylose", "mannitol"];
timeofinnoculation=posixtime(datetime(2026,4,22,12,28,0));

%%
MQdir = "/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/bioreactor_package/MQ";
MQfiles=cell(1,100);
for i=1:100
    MQfiles{i} = fullfile(MQdir, sprintf('MQ_450mw_1000ms_%05d.spc', i-1));
end

readmq=snr_spectrum_from_files_4(MQfiles,[-100,4000]);



%%
model=load("Multi_VB_res_msbackadj_SNVon.mat");
%model=load("Multi_VB_res_msbackadj_noSNV.mat");
%model=load("Multi_VB_res_asls_SNVon.mat");
%model = load("Multi_VB_res_asls_noSNV.mat");
%model = load("SF_VB_70_30_MultiPLSR_SF9_MS_withSNV_maxLV20.mat");
%model = load("SF_VB_70_30_MultiPLSR_SF9_asls_withSNV_maxLV20.mat");
%predictionvector=model.out.BETA_locked;
predictionvector=model.results.BETA;
reactorA=load("/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/bioreactor_package/20260422_10_49_39_W17_RAMAN_Glucose_A2.mat");
reactorB=load("/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/bioreactor_package/20260422_11_20_20_W17_RAMAN_Glucose_B2.mat");
reactorAbiomassntime=[reactorA.VarMem.measurements(:,2),reactorA.VarMem.ProcTime'+posixtime(datetime(reactorA.clk_0))];
reactorBbiomassntime=[reactorB.VarMem.measurements(:,2),reactorB.VarMem.ProcTime'+posixtime(datetime(reactorB.clk_0))];
reactorANIRntime=[reactorA.VarMem.Instrument.NIR',reactorA.VarMem.ProcTime'+posixtime(datetime(reactorA.clk_0))];
reactorBNIRntime=[reactorB.VarMem.Instrument.NIR',reactorB.VarMem.ProcTime'+posixtime(datetime(reactorB.clk_0))];
%%
hplc=struct();
hplc.results=readcell("/Users/celinehansen/Documents/NMPC/Matlab/Raman_HPLC_master/bioreactor_package/RAMAN week 17.xlsx",sheet="Results");
hplc.A=readcell("RAMAN week 17.xlsx",sheet="2026 Raman A2");
hplc.B=readcell("RAMAN week 17.xlsx",sheet="2026 Raman B2");
hplc_sugarsncode=[hplc.results(5:end,7:9),hplc.results(5:end,2)];
hplc_timencode_A=[hplc.A(28:end,7),hplc.A(28:end,14)];
hplc_timencode_B=[hplc.B(35:end,7),hplc.B(35:end,14)];
[~,idxA]=ismember(hplc_timencode_A(:,2),hplc_sugarsncode(:,4));
[~,idxB]=ismember(hplc_timencode_B(:,2),hplc_sugarsncode(:,4));
hplc_timencode_A=hplc_timencode_A(idxA~=0,:); hplc_timencode_B=hplc_timencode_B(idxB~=0,:);
idxA=idxA(idxA~=0);idxB=idxB(idxB~=0);
hplc_sugarntime_A=[hplc_sugarsncode(idxA,1:3), num2cell(posixtime(datetime(hplc_timencode_A(:,1))))];
hplc_sugarntime_B=[hplc_sugarsncode(idxB,1:3), num2cell(posixtime(datetime(hplc_timencode_B(:,1))))];
hplc_sugarntime_A=cellfun(@(x) num2cell(x.*isnumeric(x)),hplc_sugarntime_A ,'UniformOutput',false);
hplc_sugarntime_A=cellfun(@(x) x(1,1),hplc_sugarntime_A );
hplc_sugarntime_B=cellfun(@(x) num2cell(x.*isnumeric(x)),hplc_sugarntime_B ,'UniformOutput',false);
hplc_sugarntime_B=cellfun(@(x) x(1,1),hplc_sugarntime_B );




%%


colouring1=data{1,1}(:,2);
colouring2=data{2,1}(:,2);
colouring1=cell2mat(colouring1);
colouring2=cell2mat(colouring2);
colouring1map=interp1([colouring1(1),colouring1(end)],[[0 0 1]; [0 1 0]],colouring1);
colouring2map=interp1([colouring2(1),colouring2(end)],[[0 0 1]; [0 1 0]],colouring2);
timebar1=(interp1([1,5],[colouring1(1),colouring1(end)],1:5)-timeofinnoculation)/3600;
timebar2=(interp1([1,5],[colouring2(1),colouring2(end)],1:5)-timeofinnoculation)/3600;

figure();
for i=1:size(data{1,1},1)
    plot(data{1,1}{i,1}.ramanShift,data{1,1}{i,1}.mu,Color=colouring1map(i,:))
    hold on
    grid on
    colormap(colouring1map)
    cbar=colorbar;
    cbar.Ticks=linspace(cbar.Limits(1), cbar.Limits(2), 5);
    cbar.TickLabels={timebar1};
    cbar.Label.String="Hours since Innoculation";
    ylabel("Intensity (a.u.)")
    xlabel("Stokes Shift (cm-1)")
end
title("Raw  spectra  reactor A")

figure();
for i=1:size(data{1,1},1)
    plot(data{1,1}{i,3}.ramanShift,data{1,1}{i,3}.mu3,Color=colouring1map(i,:))
    hold on
    grid on
    colormap(colouring1map)
    cbar=colorbar;
    cbar.Ticks=linspace(cbar.Limits(1), cbar.Limits(2), 5);
    cbar.TickLabels={timebar1};
    cbar.Label.String="Hours since Innoculation";
    ylabel("Intensity (a.u.)")
    xlabel("Stokes Shift (cm-1)")
end
title("Processed spectra reactor A")

figure();
for i=1:size(data{2,1},1)
    plot(data{2,1}{i,1}.ramanShift,data{2,1}{i,1}.mu,Color=colouring2map(i,:))
    hold on
    grid on
    colormap(colouring2map)
    cbar=colorbar;
    cbar.Ticks=linspace(cbar.Limits(1), cbar.Limits(2), 5);
    cbar.TickLabels={timebar2};
    cbar.Label.String="Hours since Innoculation";
    ylabel("Intensity (a.u.)")
    xlabel("Stokes Shift (cm-1)")
end
title("Raw  spectra  reactor B")

figure();
for i=1:size(data{2,1},1)
    plot(data{2,1}{i,3}.ramanShift,data{2,1}{i,3}.mu3,Color=colouring2map(i,:))
    hold on
    grid on
    colormap(colouring2map)
    cbar=colorbar;
    cbar.Ticks=linspace(cbar.Limits(1), cbar.Limits(2), 5);
    cbar.TickLabels={timebar2};
    cbar.Label.String="Hours since Innoculation";
    ylabel("Intensity (a.u.)")
    xlabel("Stokes Shift (cm-1)")
end
title("Processed spectra reactor B")

%%

figure()
tiledlayout(3,2)

nexttile;
plot((reactorAbiomassntime(:,2)-timeofinnoculation)/3600,reactorAbiomassntime(:,1))
xlim([1,34])%([(data{1,1}{1,2}-timeofinnoculation)/3600,34])
grid on
ylabel("Estimated Biomass (g/L)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlabel("Time since innoculation (h)")
title("Estimated Biomass via NIR OD probe, reactor A")


nexttile;
plot((reactorBbiomassntime(:,2)-timeofinnoculation)/3600,reactorBbiomassntime(:,1))
xlim([1,34])%([(data{2,1}{1,2}-timeofinnoculation)/3600,34])
grid on
ylabel("Estimated Biomass (g/L)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlabel("Time since innoculation (h)")
title("Estimated Biomass via NIR OD probe, reactor B")

nexttile;
for i=1:3
    plot((cellfun(@(x) (x-timeofinnoculation)/3600,hplc_sugarntime_A(:,4))),cell2mat(hplc_sugarntime_A(:,i)).*(4/5),'o')
    grid on
    hold on
end
ylabel("Estimeated Sugar (mM)")
xlabel("Time since innoculation (h)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlim([1,34])%([(data{1,1}{1,2}-timeofinnoculation)/3600,34])
legend(sugarorder)
title("HPLC measured sugars over time, reactor A")

nexttile;
for i=1:3
    plot((cellfun(@(x) (x-timeofinnoculation)/3600,hplc_sugarntime_B(:,4))),cell2mat(hplc_sugarntime_B(:,i)).*(4/5),'o')
    grid on
    hold on
end
ylabel("Estimeated Sugar (mM)")
xlabel("Time since innoculation (h)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlim([1,34])%([(data{2,1}{1,2}-timeofinnoculation)/3600,34])
legend(sugarorder)
title("HPLC measured sugars over time, reactor B")

nexttile;
sugarpreds=cellfun(@(x) [1;x.mu3]'*predictionvector ,data{1,1}(:,3) ,'UniformOutput' ,false);
sugarpreds=cell2mat(sugarpreds);
for i=1:3
    plot((cellfun(@(x) (x-timeofinnoculation)/3600,data{1,1}(:,2))),sugarpreds(:,i)-sugarpreds(1,i))
    grid on
    hold on
end
ylabel("Estimeated Sugar (mM)")
xlabel("Time since innoculation (h)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlim([1,34])%([(data{1,1}{1,2}-timeofinnoculation)/3600,34])
legend(sugarorder)
title("Change in raman Predicted sugars over time, reactor A")

nexttile;
sugarpreds2=cellfun(@(x) [1;x.mu3]'*predictionvector ,data{2,1}(:,3) ,'UniformOutput' ,false);
sugarpreds2=cell2mat(sugarpreds2);
for i=1:3
    plot((cellfun(@(x) (x-timeofinnoculation)/3600,data{2,1}(:,2))),sugarpreds2(:,i)-sugarpreds2(1,i))
    grid on
    hold on
end
ylabel("Estimeated Sugar (mM)")
xlabel("Time since innoculation (h)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlim([1,34])%([(data{2,1}{1,2}-timeofinnoculation)/3600,34])
legend(sugarorder)
title("Change in Raman predicted sugars over time, reactor B")

%% Raman predictions on absolute concentrations
tRamanA = cellfun(@(x) (x-timeofinnoculation)/3600, data{1,1}(:,2));
tRamanB = cellfun(@(x) (x-timeofinnoculation)/3600, data{2,1}(:,2));
tHPLCA = cellfun(@(x) (x-timeofinnoculation)/3600, hplc_sugarntime_A(:,4));
tHPLCB = cellfun(@(x) (x-timeofinnoculation)/3600, hplc_sugarntime_B(:,4));

sugarpreds = cellfun(@(x) [1; x.mu3]' * predictionvector, ...
    data{1,1}(:,3), 'UniformOutput', false);
sugarpreds = cell2mat(sugarpreds);

sugarpreds2 = cellfun(@(x) [1; x.mu3]' * predictionvector, ...
    data{2,1}(:,3), 'UniformOutput', false);
sugarpreds2 = cell2mat(sugarpreds2);

% Optional: remove negative predicted concentrations
%sugarpreds  = max(sugarpreds, 0);
%sugarpreds2 = max(sugarpreds2, 0);

%% Offset-correct Raman predictions to first HPLC value

hplcA = cell2mat(hplc_sugarntime_A(:,1:3)).*(4/5);
hplcB = cell2mat(hplc_sugarntime_B(:,1:3)).*(4/5);

% First nonzero / valid HPLC value for each sugar
hplcA0 = zeros(1,3);
hplcB0 = zeros(1,3);

for i = 1:3
    idxA0 = find(~isnan(hplcA(:,i)) & hplcA(:,i) > 0, 1, 'first');
    idxB0 = find(~isnan(hplcB(:,i)) & hplcB(:,i) > 0, 1, 'first');

    hplcA0(i) = hplcA(idxA0,i);
    hplcB0(i) = hplcB(idxB0,i);
end

sugarpreds_corr  = sugarpreds  - sugarpreds(1,:)  + hplcA0;
sugarpreds2_corr = sugarpreds2 - sugarpreds2(1,:) + hplcB0;

% Optional: only for visualisation
sugarpreds_corr  = max(sugarpreds_corr, 0);
sugarpreds2_corr = max(sugarpreds2_corr, 0);

%% Separate comparison figure: Raman vs HPLC for each sugar

figure()
tiledlayout(3,2, 'TileSpacing','compact', 'Padding','compact')

for i = 1:3

    % Reactor A
    nexttile;
    plot(tRamanA, sugarpreds_corr(:,i), '-', 'LineWidth', 1.2, ...
        'DisplayName','Raman PLSR prediction')
    hold on
    plot(tHPLCA, cell2mat(hplc_sugarntime_A(:,i)).*(4/5), 'ko', ...
        'MarkerSize', 5, 'DisplayName','HPLC measured')
    grid on
    xline(11.5,'--b', 'feeding 1')
    xline(22,'--b', 'feeding 2')
    xlim([1,34])
    ylim([0 inf])
    ylabel("Concentration (mM)")
    xlabel("Time since inoculation (h)")
    title("Reactor A: " + sugarorder(i))
    legend('Location','best')

    % Reactor B
    nexttile;
    plot(tRamanB, sugarpreds2_corr(:,i), '-', 'LineWidth', 1.2, ...
        'DisplayName','Raman PLSR prediction')
    hold on
    plot(tHPLCB, cell2mat(hplc_sugarntime_B(:,i)).*(4/5), 'ko', ...
        'MarkerSize', 5, 'DisplayName','HPLC measured')
    grid on
    xline(11.5,'--b', 'feeding 1')
    xline(22,'--b', 'feeding 2')
    xlim([1,34])
    ylim([0 inf])
    ylabel("Concentration (mM)")
    xlabel("Time since inoculation (h)")
    title("Reactor B: " + sugarorder(i))
    legend('Location','best')

end

sgtitle("Raman PLSR predicted and HPLC measured sugar concentrations in bioreactors")

%%

areas1=cellfun(@(x) sum(x.mu(1772:end)) ,data{1,1}(:,1) ,'UniformOutput' ,false);
areas1=cell2mat(areas1);
areas2=cellfun(@(x) sum(x.mu(1772:end)) ,data{2,1}(:,1) ,'UniformOutput' ,false);
areas2=cell2mat(areas2);



figure()
tiledlayout(2,2)

nexttile;
plot((reactorANIRntime(:,2)-timeofinnoculation)/3600,reactorANIRntime(:,1))
xlim([1,34])%([(data{1,1}{1,2}-timeofinnoculation)/3600,34])
grid on
ylabel("TCD (t.b.d)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlabel("Time since innoculation (h)")
title("NIR OD probe")


nexttile;
plot((reactorBNIRntime(:,2)-timeofinnoculation)/3600,reactorBNIRntime(:,1))
xlim([1,34])%([(data{2,1}{1,2}-timeofinnoculation)/3600,34])
grid on
ylabel("TCD (t.b.d)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlabel("Time since innoculation (h)")
title("NIR OD probe")



nexttile
plot((cellfun(@(x) (x-timeofinnoculation)/3600,data{1,1}(:,2))),sum(readmq.mu(1772:end))./areas1-1)
grid on
ylabel("Change in Attenuation (Water attenuation)")
xlabel("time since innoculation (h)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlim([1,34])
title("change in OD over time, reactor A")



nexttile
plot((cellfun(@(x) (x-timeofinnoculation)/3600,data{2,1}(:,2))),sum(readmq.mu(1772:end))./areas2-1)
grid on
ylabel("Change in Attenuation (Water attenuation)")
xlabel("time since innoculation (h)")
xline(11.5,'--b', 'feeding 1')
xline(22,'--b', 'feeding 2')
xlim([1,34])
title("change in OD over time, reactor B")

%%

figure();
for i=1:size(data{2,1},1)
    plot(data{2,1}{i,1}.ramanShift,data{1,1}{i,1}.sd3,Color=colouring2map(i,:))
    hold on
    grid on
    colormap(colouring2map)
    cbar=colorbar;
    cbar.Ticks=linspace(cbar.Limits(1), cbar.Limits(2), 5);
    cbar.TickLabels={timebar2};
    cbar.Label.String="Hours since Innoculation";
    ylabel("Intensity (a.u.^2)")
    xlabel("Stokes Shift (cm-1)")
    ylim([0,0.1])
end
title("noise processed reactor A")

figure();
for i=1:size(data{2,1},1)
    plot(data{2,1}{i,1}.ramanShift,data{2,1}{i,1}.sd3,Color=colouring2map(i,:))
    hold on
    grid on
    colormap(colouring2map)
    cbar=colorbar;
    cbar.Ticks=linspace(cbar.Limits(1), cbar.Limits(2), 5);
    cbar.TickLabels={timebar2};
    cbar.Label.String="Hours since Innoculation";
    ylabel("Intensity (a.u.^2)")
    xlabel("Stokes Shift (cm-1)")
    ylim([0,0.1])
end
title("noise processed reactor B")