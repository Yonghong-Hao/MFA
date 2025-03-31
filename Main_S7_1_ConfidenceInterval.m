%==========================================================================
%===================== Confidence interval estimate =======================
%==========================================================================
clear;clc
tic
%% input
TracerID = 'glc_12';
MappingID = 'mapping1';
DataSet = {'data1_1','data1_2','data1_3'};
Conditions = {'CM','LA+','LA'};
%% read data
for i = 1%:numel(Conditions)
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
    NumMetabolites = height(SMatrix);
    % ToSimulateMS
    ToSimulateMS = AllVariables.ToSimulateMS;
    ToSimulateMS_MID = AllVariables.ToSimulateMS_MID;
    NumToSimulateMS_MID = height(ToSimulateMS_MID);
    % ToSimulateFlux
    ToSimulateFlux = AllVariables.ToSimulateFlux;
    % read Vhat
    Vhat_table = readtable(strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/Vhat_',MappingID,'.xlsx'),"VariableNamingRule","preserve");
    Vhat_total = [Vhat_table.CM,Vhat_table.("LA+"),Vhat_table.LA];
    Vhat = Vhat_total(:,i);
    %% Grid search
    % optimal solution
    [X, ~] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,Vhat,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
    X_Simulated = zeros(NumToSimulateMS_MID,1);
    for ii = 1:NumToSimulateMS_MID
        condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{ii});
        X_Simulated(ii) = X(condition);
    end
    V_Simulated = zeros(height(ToSimulateFlux),1);
    for ii = 1:height(ToSimulateFlux)
        switch ToSimulateFlux.Reversibility(ii)
            case 0
                V_Simulated(ii) = Vhat(ToSimulateFlux.FluxIndex_f(ii));
            case 1
                V_Simulated(ii) = Vhat(ToSimulateFlux.FluxIndex_f(ii))-Vhat(ToSimulateFlux.FluxIndex_b(ii));
        end
    end
    SimulatedData = [X_Simulated;V_Simulated];
    ToSimulateData = [str2double(ToSimulateMS_MID.Measured);ToSimulateFlux.Measured];
    Phi_Vhat = (SimulatedData-ToSimulateData)'*(AllVariables.SigmaInv.^2)*(SimulatedData-ToSimulateData);
    %% Backward search
    SearchedData = struct('BackwardSearch',{},'OptimalSolution',{},'ForwardSearch',{});
    parfor j = 1:NumFluxes
        BackwardData = zeros(1000,3);
        if log10(Vhat(j))>1
            StepSize = Vhat(j)/10000
        else
            StepSize = Vhat(j)/1000;
        end
        OldV = Vhat;
        kb = 1;
        while true
            [NewVBackward, exitflag,fval] = S6_BackwardMinNewVOldV(SMatrix,OldV,NumFluxes,NumMetabolites,j,StepSize);
            [X, ~] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,NewVBackward,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
            X_Simulated = zeros(NumToSimulateMS_MID,1);
            for ii = 1:NumToSimulateMS_MID
                condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{ii});
                X_Simulated(ii) = X(condition);
            end
            V_Simulated = zeros(height(ToSimulateFlux),1);
            for ii = 1:height(ToSimulateFlux)
                switch ToSimulateFlux.Reversibility(ii)
                    case 0
                        V_Simulated(ii) = NewVBackward(ToSimulateFlux.FluxIndex_f(ii));
                    case 1
                        V_Simulated(ii) = NewVBackward(ToSimulateFlux.FluxIndex_f(ii))-NewVBackward(ToSimulateFlux.FluxIndex_b(ii));
                end
            end
            SimulatedData = [X_Simulated;V_Simulated];
            ToSimulateData = [str2double(ToSimulateMS_MID.Measured);ToSimulateFlux.Measured];
            Phi_V = (SimulatedData-ToSimulateData)'*(AllVariables.SigmaInv.^2)*(SimulatedData-ToSimulateData);
            disp([kb,NewVBackward(j),Phi_V-Phi_Vhat])
            BackwardData(kb,:) = [kb,NewVBackward(j),Phi_V];
            if Phi_V-Phi_Vhat>chi2inv(1-0.05,1) || NewVBackward(j)<1e-20 || kb == 1000
                break
            end
            OldV = NewVBackward;
            kb = kb+1;
        end
        SearchedData(j).OptimalSolution = array2table([0,Vhat(j),Phi_Vhat],"VariableNames",{'iteration','Vhat','Phi_Vhat'});
        SearchedData(j).BackwardSearch = array2table(BackwardData(1:kb,:),"VariableNames",{'iteration','NewV','Phi_NewV'});
    end
    %% Forward search
    parfor j = 1:NumFluxes
        ForwardData = zeros(1000,3);
        if log10(Vhat(j))>1
            StepSize = Vhat(j)/10000
        else
            StepSize = Vhat(j)/1000;
        end
        OldV = Vhat;
        kf = 1;
        while true
            [NewVForward, exitflag,fval] = S6_ForwardMinNewVOldV(SMatrix,OldV,NumFluxes,NumMetabolites,j,StepSize);
            [X, ~] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,NewVForward,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
            X_Simulated = zeros(NumToSimulateMS_MID,1);
            for ii = 1:NumToSimulateMS_MID
                condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{ii});
                X_Simulated(ii) = X(condition);
            end
            V_Simulated = zeros(height(ToSimulateFlux),1);
            for ii = 1:height(ToSimulateFlux)
                switch ToSimulateFlux.Reversibility(ii)
                    case 0
                        V_Simulated(ii) = NewVForward(ToSimulateFlux.FluxIndex_f(ii));
                    case 1
                        V_Simulated(ii) = NewVForward(ToSimulateFlux.FluxIndex_f(ii))-NewVForward(ToSimulateFlux.FluxIndex_b(ii));
                end
            end
            SimulatedData = [X_Simulated;V_Simulated];
            ToSimulateData = [str2double(ToSimulateMS_MID.Measured);ToSimulateFlux.Measured];
            Phi_V = (SimulatedData-ToSimulateData)'*(AllVariables.SigmaInv.^2)*(SimulatedData-ToSimulateData);
            disp([kf,NewVForward(j),Phi_V-Phi_Vhat])
            ForwardData(kf,:) = [kf,NewVForward(j),Phi_V];
            if Phi_V-Phi_Vhat>chi2inv(1-0.05,1) || NewVForward(j)>100 || kf == 1000
                break
            end
            OldV = NewVForward;
            kf = kf+1;
        end
        SearchedData(j).ForwardSearch = array2table(ForwardData(1:kf,:),"VariableNames",{'iteration','NewV','Phi_NewV'});
    end
    save(strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/ConfidenceInterval_',MappingID,'_',DataSet{i},'.mat'),'SearchedData')
end
toc