%==========================================================================
%============================ Goodness of fit =============================
%==========================================================================
clear;clc
%% input
TracerID = 'glc_12';
MappingID = 'mapping9';
DataSet = {'data1_1','data1_2','data1_3'};
Conditions = {'CM','LA+','LA'};
%% read data
f1 = figure;
f1.Position = [1921          49        1920         956];
set(groot,'DefaultAxesFontSize',15);
SSR_raw = zeros(numel(Conditions),1);
SSR_mine = zeros(numel(Conditions),1);
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
    % read Vhat
    Vhat_table = readtable(strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/Vhat_',MappingID,'.xlsx'),"VariableNamingRule","preserve");
    Vhat_total = [Vhat_table.CM,Vhat_table.("LA+"),Vhat_table.LA];
    %% Solve the model
    % MID
    V = Vhat_total(:,i);
    [X, dXdV] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,V,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
    NumToSimulateMS_MID = height(ToSimulateMS_MID);
    emuSimulatedMID = zeros(NumToSimulateMS_MID,1);
    for j = 1:NumToSimulateMS_MID
        condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{j});
        emuSimulatedMID(j) = X(condition);
    end
    % Flux
    SimulatedFlux = zeros(height(ToSimulateFlux),1);
    for k = 1:height(ToSimulateFlux)
        switch ToSimulateFlux.Reversibility(k)
            case 0
                SimulatedFlux(k) = V(ToSimulateFlux.FluxIndex_f(k));
            case 1
                SimulatedFlux(k) = V(ToSimulateFlux.FluxIndex_f(k))-V(ToSimulateFlux.FluxIndex_b(k));
        end
    end
    % Compute raw SSR 
    ToSimulateData = [str2double(ToSimulateMS_MID.Measured);ToSimulateFlux.Measured];
    SimulatedData_raw = [ToSimulateMS_MID.Simulated;ToSimulateFlux.Simulated]; 
    SSR_raw(i) = (SimulatedData_raw-ToSimulateData)'*(AllVariables.SigmaInv.^2)*(SimulatedData_raw-ToSimulateData);
    % Compute mine SSR
    SimulatedData_mine = [emuSimulatedMID;SimulatedFlux];
    SSR_mine(i) = (SimulatedData_mine-ToSimulateData)'*(AllVariables.SigmaInv.^2)*(SimulatedData_mine-ToSimulateData);
    % Plot
    for j = 1:height(ToSimulateMS)+4
        if j<=height(ToSimulateMS)
            subplot(numel(Conditions),height(ToSimulateMS)+4,7*(i-1)+j);
            Size = ToSimulateMS.Size(j);
            emu = ToSimulateMS.EMUIDs{j};
            mid = strcat('M',cellfun(@num2str, num2cell(0:Size),'UniformOutput',false));
            indexEmu = startsWith(ToSimulateMS_MID.EMUIDs_mid,strcat(emu,'_'));
            data = [cell2mat(cellfun(@str2num,ToSimulateMS_MID.Measured(indexEmu),'UniformOutput',false)),ToSimulateMS_MID.Simulated(indexEmu),emuSimulatedMID(indexEmu)];
            databar=bar(data,'grouped','EdgeColor','none');
            databar(1).FaceColor =  [0.00, 0.68, 0.92];
            databar(2).FaceColor =  [1, 0.49, 0.43];
            databar(3).FaceColor =  [0.9961    0.7608    0.5373];
            xticks(1:1:Size+1);
            xticklabels(mid);
            ylim([0,1])
            set(gca, 'XGrid', 'off', 'YGrid', 'on');
            ylabel("Relative fraction");
            title(strcat(Conditions{i},': ',regexprep(emu,'_',' ')),'Units', 'Normalized','Position',[0.5,1],'FontSize',15,'FontWeight','bold')
        else
            subplot(numel(Conditions),height(ToSimulateMS)+4,[7*(i-1)+3+1,7*(i-1)+7]);
            factor = 1*ones(height(ToSimulateFlux),1);
            for k = 1:height(ToSimulateFlux)
                if strcmp(ToSimulateFlux.ReactionID{k},'Biomass') 
                    factor(k) = 1;
                elseif strcmp(ToSimulateFlux.ReactionID{k},'mAb')
                    factor(k) = 1;
                end
            end
            data = [ToSimulateFlux.Measured,ToSimulateFlux.Simulated,SimulatedFlux].*factor;
            databar=bar(data,'grouped','EdgeColor','none');
            databar(1).FaceColor =  [0.00, 0.68, 0.92];
            databar(2).FaceColor =  [1, 0.49, 0.43];
            databar(3).FaceColor =  [0.9961    0.7608    0.5373];
            xticks(1:1:height(ToSimulateFlux));
            xticklabels(ToSimulateFlux.ReactionID);
            ylabel("Net Rates");
            set(gca, 'XGrid', 'off', 'YGrid', 'on');
            title(strcat(Conditions{i},': ','Flux'),'Units', 'Normalized','Position',[0.5,1],'FontSize',15,'FontWeight','bold')
        end
    end
end
%%
hl = legend('Experimental MID', 'Simulated MID (Raw)','Simulated MID (New)','Orientation', 'horizontal');
set(hl,'Units', 'normalized', 'Position', [0.28, 0.96, 0.05, 0.03],'Box','off','FontWeight','bold');
annotation('textbox', [0.78, 0.83, 0.1, 0.1], 'String', {['Raw SSR = ',num2str(SSR_raw(1))],['New SSR = ',num2str(SSR_mine(1))]},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FitBoxToText', 'on', 'BackgroundColor', 'none', 'EdgeColor', 'none','FontSize',13,'FontWeight','bold');
annotation('textbox', [0.78, 0.53, 0.1, 0.1], 'String', {['Raw SSR = ',num2str(SSR_raw(2))],['New SSR = ',num2str(SSR_mine(2))]},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FitBoxToText', 'on', 'BackgroundColor', 'none', 'EdgeColor', 'none','FontSize',13,'FontWeight','bold');
annotation('textbox', [0.78, 0.23, 0.1, 0.1], 'String', {['Raw SSR = ',num2str(SSR_raw(3))],['New SSR = ',num2str(SSR_mine(3))]},'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FitBoxToText', 'on', 'BackgroundColor', 'none', 'EdgeColor', 'none','FontSize',13,'FontWeight','bold');
saveas(f1,strcat('Tracer_',TracerID,'/VerifyMyModel/GoodnessOfFit/GoodnessOfFit_',MappingID),'png')