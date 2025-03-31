clear;clc
%% input
TracerID = 'glc_12';
MappingID = 'mapping1';
DataSet = {'data1_1','data1_2','data1_3'};
Conditions = {'CM','LA+','LA'};
%% read data
f1 = figure;
f1.Position = [1971,9,1134,956];
set(groot,'DefaultAxesFontSize',15);
SSR_MID = zeros(numel(Conditions),1);
for i = 1:numel(Conditions)
    % obtain the VariablesRefernce, BoundaryVariablesRefernce, STable, ToSimulateMS and ToSimulateFlux
    load(strcat('Tracer_',TracerID,'/FluxBalance&EMUBalance/IntermediateVariable_',MappingID,'_',DataSet{i},'.mat'));
    % VariablesRefernce
    VariablesRefernce = AllVariables.VariablesRefernce;
    NumVariables = height(VariablesRefernce);
    SizeVariables = VariablesRefernce.emu_size;
    XVar = VariablesRefernce.emu_var;
    dXdVVar = regexprep(VariablesRefernce.emu_var,'X','dXdV');
    dAdV = AllVariables.dAdV;
    % BoundaryVariablesRefernce
    BoundaryVariablesRefernce = AllVariables.BoundaryVariablesRefernce;
    NumBoundaryVariables  = height(BoundaryVariablesRefernce);
    SizeBoundaryVariables = BoundaryVariablesRefernce.emu_size;
    YVar = BoundaryVariablesRefernce.emu_var;
    dYdVVar = S2_NonlinearGradient(NumBoundaryVariables,BoundaryVariablesRefernce);
    dBdV = AllVariables.dBdV;
    % STable
    STable = AllVariables.STable;
    SMatrix = STable.Variables;
    NumFluxes = width(SMatrix);
    % ToSimulateMS
    ToSimulateMS = AllVariables.ToSimulateMS;
    ToSimulateMS_MID = AllVariables.ToSimulateMS_MID;
    % ToSimulateFlux
    ToSimulateFlux = AllVariables.ToSimulateFlux;
    % Load the raw flux values
    FluxValues_raw = readtable(strcat('Tracer_',TracerID,'/Data_raw/FluxValues_estimation_raw.xlsx'),"VariableNamingRule","preserve");
    MappingAtom = readtable(strcat('Tracer_',TracerID,'/MappingAtom/',regexprep(MappingID,'mapping','mapping_atom_'),'.xlsx'),"VariableNamingRule","preserve");
    Vhat = zeros(NumFluxes,1);
    for j = 1:height(FluxValues_raw)
        if FluxValues_raw.Reversibility(j) == 0
            index = ismember(MappingAtom.("SubstrateIDs(atoms)"),FluxValues_raw.("SubstrateIDs(atoms)"){j}) & strcmp(MappingAtom.("ProductIDs(atoms)"),FluxValues_raw.("ProductIDs(atoms)"){j});
            Vhat(index) = FluxValues_raw{j,5+2*i-1};
        else
            index_f = ismember(MappingAtom.("SubstrateIDs(atoms)"),FluxValues_raw.("SubstrateIDs(atoms)"){j}) & strcmp(MappingAtom.("ProductIDs(atoms)"),FluxValues_raw.("ProductIDs(atoms)"){j});
            index_b = ismember(MappingAtom.("SubstrateIDs(atoms)"),FluxValues_raw.("ProductIDs(atoms)"){j}) & strcmp(MappingAtom.("ProductIDs(atoms)"),FluxValues_raw.("SubstrateIDs(atoms)"){j});
            Vhat(index_f) = FluxValues_raw{j,5+2*i-1};
            Vhat(index_b) = FluxValues_raw{j,5+2*i};
        end
    end
    V = Vhat;
    % MID
    [X, dXdV] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,V,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
    NumToSimulateMS_MID = height(ToSimulateMS_MID);
    emuSimulatedMID = zeros(NumToSimulateMS_MID,1);
    for j = 1:NumToSimulateMS_MID
        condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{j});
        emuSimulatedMID(j) = X(condition);
    end
    SSR_MID(i) = (str2double(ToSimulateMS_MID.Measured)-emuSimulatedMID)'*(AllVariables.SigmaInv(1:NumToSimulateMS_MID,1:NumToSimulateMS_MID)).^2*(str2double(ToSimulateMS_MID.Measured)-emuSimulatedMID);
    %% Plot
    for ii = 1:height(ToSimulateMS)
        subplot(numel(Conditions),height(ToSimulateMS),numel(Conditions)*(i-1)+ii);
        Size = ToSimulateMS.Size(ii);
        emu = ToSimulateMS.EMUIDs{ii};
        mid = strcat('M',cellfun(@num2str, num2cell(0:Size),'UniformOutput',false));
        indexEmu = startsWith(ToSimulateMS_MID.EMUIDs_mid,strcat(emu,'_'));
        data = [cell2mat(cellfun(@str2num,ToSimulateMS_MID.Measured(indexEmu),'UniformOutput',false)),emuSimulatedMID(indexEmu)];
        databar=bar(data,'grouped','EdgeColor','none');
        databar(1).FaceColor =  [0.00, 0.68, 0.92];
        databar(2).FaceColor =  [0.9961    0.7608    0.5373];
        xticks(1:1:Size+1);
        xticklabels(mid);
        ylim([0,1])
        set(gca,'XGrid','off','YGrid','on')
        ylabel("Relative fraction");
        title(regexprep(emu,'_',' '),'Units', 'Normalized','Position',[0.5,1],'FontSize',15,'FontWeight','bold')
    end
end
%%
hl = legend('Experimental MID', 'Simulated MID','Orientation', 'horizontal');
set(hl,'Units', 'normalized', 'Position', [0.23, 0.96, 0.05, 0.03],'Box','off','FontWeight','bold');
annotation('textbox', [0.78, 0.85, 0.1, 0.1], 'String', ['SSR = ',num2str(SSR_MID(1))],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'none','FontSize',13,'FontWeight','bold');
annotation('textbox', [0.78, 0.55, 0.1, 0.1], 'String', ['SSR = ',num2str(SSR_MID(2))],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'none','FontSize',13,'FontWeight','bold');
annotation('textbox', [0.78, 0.25, 0.1, 0.1], 'String', ['SSR = ',num2str(SSR_MID(3))],'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FitBoxToText', 'on', 'BackgroundColor', 'white', 'EdgeColor', 'none','FontSize',13,'FontWeight','bold');
saveas(f1,strcat('Tracer_',TracerID,'/VerifyMyModel/TestRawFluxes/TestRawFluxes_',MappingID),'png')