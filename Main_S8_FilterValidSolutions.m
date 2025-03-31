clear;clc
%% input
TracerID = 'Tracer_glc_12';
MappingID = 'mapping1';
DataID_1 = 'data1_1';
DataID_2 = 'data1_2';
DataID_3 = 'data1_3';
threshold_FB = 1e-9;
threshold_EB = 0.01;

%% 合并单向流为可逆流
% Mapping without reversibility
Mapping_Irreversibility = readtable(strcat(TracerID,"/MappingAtom/",regexprep(MappingID,'mapping','mapping_atom_'),".xlsx"),"VariableNamingRule","preserve");
SubstrateTotal = Mapping_Irreversibility.("SubstrateIDs(atoms)");
ProductTotal = Mapping_Irreversibility.("ProductIDs(atoms)");

% Mapping with reversibility
Mapping_Reversibility = readtable(strcat(TracerID,"/Data_raw/",regexprep(MappingID,'mapping','MappingAtom_raw_'),".xlsx"),"VariableNamingRule","preserve");
NumReactions = height(Mapping_Reversibility);
ReactionIDs = Mapping_Reversibility.ReactionIDs;
Reaction_f_index = zeros(numel(ReactionIDs),1);
Reaction_b_index = Reaction_f_index;
Reversibility = Mapping_Reversibility.Reversibility;
for i = 1:NumReactions
    Reversibility_i = Reversibility(i);
    Substrate_i = Mapping_Reversibility.("SubstrateIDs(atoms)"){i};
    Product_i = Mapping_Reversibility.("ProductIDs(atoms)"){i};
    switch Reversibility_i
        case 0
            Reaction_f_index(i) = find((ismember(SubstrateTotal,Substrate_i)+ismember(ProductTotal,Product_i))==2);
            Reaction_b_index(i) = inf;
        case 1
            Reaction_f_index(i) = find((ismember(SubstrateTotal,Substrate_i)+ismember(ProductTotal,Product_i))==2);
            Reaction_b_index(i) = find((ismember(SubstrateTotal,Product_i)+ismember(ProductTotal,Substrate_i))==2);
    end
end
ReactionsTotal = table(ReactionIDs,Reaction_f_index,Reaction_b_index,Reversibility);

%% 创建结构体，将数据读取存储在该结构体 (only flux values)
Condition = {'Parental','FH diminished'};
NumCondition = numel(Condition);
%Result_Parental = struct('FileName',{},'FluxValue',{},'LossFunctionValue',{});
%Result_FH_diminished = struct('FileName',{},'FluxValue',{},'LossFunctionValue',{});
Result_Parental = struct('FileName',{},'FluxValue',{});
Result_FH_diminished = struct('FileName',{},'FluxValue',{});
for i = 1:NumCondition
    Condition_i = Condition{i};
    FolderPath = strcat(TracerID,"/FluxValues_LossFunctionValues/",MappingID,'/',Condition_i);
    FilePattern_V_i = fullfile(FolderPath, '*FluxValues*.xlsx');
    Files_V_i = dir(FilePattern_V_i);

    switch Condition_i
        case Condition{1}
            for i_V = 1:numel(Files_V_i)
                % read the VStore and fvalStore table -- Monoalyer
                Result_Parental(i_V).FileName = Files_V_i(i_V).name;
                Result_Parental(i_V).FluxValue = readtable(fullfile(FolderPath, Files_V_i(i_V).name),"VariableNamingRule","preserve","ReadRowNames",true);
            end
        case Condition{2}
            for i_V = 1:numel(Files_V_i)
                % read the VStore and fvalStore table -- Spheroid
                Result_FH_diminished(i_V).FileName = Files_V_i(i_V).name;
                Result_FH_diminished(i_V).FluxValue = readtable(fullfile(FolderPath, Files_V_i(i_V).name),"VariableNamingRule","preserve","ReadRowNames",true);
            end
    end
end

%%%%%%%%%%%%%%% Condition 1 %%%%%%%%%%%%%%%
%% combine all flux values
NumIteration = numel(Result_Parental);
VTotal_condition1 = [];
for i = 1:NumIteration
    VTotal_i = Result_Parental(i).FluxValue.Variables;
    VTotal_condition1 = [VTotal_condition1,VTotal_i];
end
NumSamples_condition_1 = width(VTotal_condition1);

%% obtain the VariablesRefernce, BoundaryVariablesRefernce, STable and ToSimulateMS
load(strcat(TracerID,'/FluxBalance&EMUBalance/IntermediateVariable_',MappingID,'_',DataID_1,'.mat'));
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
ToSimulateMS_MID = table('Size',[0 5],'VariableTypes',{'string','string','double','double','double'},'VariableNames',{'EMUIDs_mid','Measured','Simulated','Weight','Size'});
for i = 1:height(ToSimulateMS)
    ToSimulate = strcat(ToSimulateMS.EMUIDs{i},'_m',cellfun(@num2str, num2cell(0:ToSimulateMS.Size(i)),'UniformOutput',false));
    MIDDataMean = eval(ToSimulateMS.Measured{i});
    emuSimulatedMID = eval(ToSimulateMS.Simulated{i});
    weightMID = eval(ToSimulateMS.Std{i});
    Size = repmat(ToSimulateMS.Size(i),ToSimulateMS.Size(i)+1,1);
    ToSimulateSet_new = table(ToSimulate',MIDDataMean',emuSimulatedMID',weightMID',Size,'VariableNames',{'EMUIDs_mid','Measured','Simulated','Weight','Size'});
    ToSimulateMS_MID = vertcat(ToSimulateMS_MID,ToSimulateSet_new);
end
%% Traverse each set of fluxes
% test S balance
SBalance = zeros(1,NumSamples_condition_1);
for s = 1:NumSamples_condition_1
    V_s = VTotal_condition1(:,s);
    balance_right = SMatrix*V_s;
    if all(abs(balance_right)<threshold_FB)
        SBalance(s) = 1;
    end
end
VTotal_condition1_new = VTotal_condition1(:,logical(SBalance));
NumSamples_condition_1_new = width(VTotal_condition1_new);

% test EMU balance
EMUBalance = zeros(1,NumSamples_condition_1_new);
RMSEStore_condition1 = zeros(1,NumSamples_condition_1_new);
for s = 1:NumSamples_condition_1_new
    V_s = VTotal_condition1_new(:,s);
    % 计算Vhat下的emu mid
    [X, ~] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,V_s,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
    T = S4_ConstraintForX0(VariablesRefernce);
    constriantX = T*X;
    if all((abs(constriantX-1)<1e-5 | abs(constriantX)<1e-5)==true)
        EMUBalance(s) = 1;
    end
    % 选取MID模拟值
    NumToSimulateMS_MID = height(ToSimulateMS_MID);
    emuSimulatedMID = zeros(NumToSimulateMS_MID,1);
    for j = 1:NumToSimulateMS_MID
        condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{j});
        emuSimulatedMID(j) = X(condition);
    end
    % RMSE of MID of EMUs
    ToSimulateData = str2double(ToSimulateMS_MID.Measured);
    RMSE = sqrt(1/NumToSimulateMS_MID*((emuSimulatedMID-ToSimulateData)'*(emuSimulatedMID-ToSimulateData)));
    if RMSE<threshold_EB && EMUBalance(s)==1
        RMSEStore_condition1(s) = RMSE;
    end
end

% combine the VTotal_condition1_new and RMSEStore
data_combine = unique([RMSEStore_condition1;VTotal_condition1_new]','rows','stable');
RMSEStore_condition1_unique = data_combine(:,1);
VTotal_condition1_new_unique = data_combine(:,2:end);

% choose the best solution
[RMSEStore_condition1_unique_hat, min_index] = min(RMSEStore_condition1_unique);
VTotal_condition1_new_unique_hat = VTotal_condition1_new_unique(min_index,:);

%%%%%%%%%%%%%%% Condition 2 %%%%%%%%%%%%%%%
%% combine all flux values
NumIteration = numel(Result_FH_diminished);
VTotal_condition2 = [];
for i = 1:NumIteration
    VTotal_i = Result_FH_diminished(i).FluxValue.Variables;
    VTotal_condition2 = [VTotal_condition2,VTotal_i];
end
NumSamples_condition_2 = width(VTotal_condition2);

%% obtain the VariablesRefernce, BoundaryVariablesRefernce, STable and ToSimulateMS
load(strcat(TracerID,'/FluxBalance&EMUBalance/IntermediateVariable_',MappingID,'_',DataID_2,'.mat'));
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
ToSimulateMS_MID = table('Size',[0 5],'VariableTypes',{'string','string','double','double','double'},'VariableNames',{'EMUIDs_mid','Measured','Simulated','Weight','Size'});
for i = 1:height(ToSimulateMS)
    ToSimulate = strcat(ToSimulateMS.EMUIDs{i},'_m',cellfun(@num2str, num2cell(0:ToSimulateMS.Size(i)),'UniformOutput',false));
    MIDDataMean = eval(ToSimulateMS.Measured{i});
    emuSimulatedMID = eval(ToSimulateMS.Simulated{i});
    weightMID = eval(ToSimulateMS.Weight{i});
    Size = repmat(ToSimulateMS.Size(i),ToSimulateMS.Size(i)+1,1);
    ToSimulateSet_new = table(ToSimulate',MIDDataMean',emuSimulatedMID',weightMID',Size,'VariableNames',{'EMUIDs_mid','Measured','Simulated','Weight','Size'});
    ToSimulateMS_MID = vertcat(ToSimulateMS_MID,ToSimulateSet_new);
end
%% Traverse each set of fluxes
% test S balance
SBalance = zeros(1,NumSamples_condition_2);
for s = 1:NumSamples_condition_2
    V_s = VTotal_condition2(:,s);
    balance_right = SMatrix*V_s;
    if all(abs(balance_right)<threshold_FB)
        SBalance(s) = 1;
    end
end
VTotal_condition2_new = VTotal_condition2(:,logical(SBalance));
NumSamples_condition_2_new = width(VTotal_condition2_new);

% test EMU balance
EMUBalance = zeros(1,NumSamples_condition_2_new);
RMSEStore_condition2 = zeros(1,NumSamples_condition_2_new);
for s = 1:NumSamples_condition_2_new
    V_s = VTotal_condition2_new(:,s);
    % 计算Vhat下的emu mid
    [X, ~] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,V_s,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);
    T = S4_ConstraintForX0(VariablesRefernce);
    constriantX = T*X;
    if all((abs(constriantX-1)<1e-5 | abs(constriantX)<1e-5)==true)
        EMUBalance(s) = 1;
    end
    % 选取MID模拟值
    NumToSimulateMS_MID = height(ToSimulateMS_MID);
    emuSimulatedMID = zeros(NumToSimulateMS_MID,1);
    for j = 1:NumToSimulateMS_MID
        condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{j});
        emuSimulatedMID(j) = X(condition);
    end
    % RMSE of MID of EMUs
    ToSimulateData = str2double(ToSimulateMS_MID.Measured);
    RMSE = sqrt(1/NumToSimulateMS_MID*((emuSimulatedMID-ToSimulateData)'*(emuSimulatedMID-ToSimulateData)));
    if RMSE<threshold_EB && EMUBalance(s)==1
        RMSEStore_condition2(s) = RMSE;
    end
end

% combine the VTotal_condition1_new and RMSEStore
data_combine = unique([RMSEStore_condition2;VTotal_condition2_new]','rows','stable');
RMSEStore_condition2_unique = data_combine(:,1);
VTotal_condition2_new_unique = data_combine(:,2:end);

% choose the best solution
[RMSEStore_condition2_unique_hat, min_index] = min(RMSEStore_condition2_unique);
VTotal_condition2_new_unique_hat = VTotal_condition2_new_unique(min_index,:);

%%%%%%%%%%%%%%% Variations/(Controls+Variation) %%%%%%%%%%%%%%%
%% Variations/Controls
NumSolution_1 = height(VTotal_condition1_new_unique);
NumSolution_2 = height(VTotal_condition2_new_unique);
VTotal_condition1_new_unique_CombineReverse = zeros(NumSolution_1,NumReactions);
VTotal_condition2_new_unique_CombineReverse = zeros(NumSolution_2,NumReactions);
for i = 1:NumReactions
    ReactionID_i = ReactionsTotal.ReactionIDs{i};
    Reversibility_i = ReactionsTotal.Reversibility(i);
    Reaction_f_index_i = ReactionsTotal.Reaction_f_index(i);
    Reaction_b_index_i = ReactionsTotal.Reaction_b_index(i);
    switch Reversibility_i
        case 0
            VTotal_condition1_new_unique_CombineReverse(:,i) = VTotal_condition1_new_unique(:,Reaction_f_index_i);
            VTotal_condition2_new_unique_CombineReverse(:,i) = VTotal_condition2_new_unique(:,Reaction_f_index_i);
        case 1
            VTotal_condition1_new_unique_CombineReverse(:,i) = VTotal_condition1_new_unique(:,Reaction_f_index_i)-VTotal_condition1_new_unique(:,Reaction_b_index_i);
            VTotal_condition2_new_unique_CombineReverse(:,i) = VTotal_condition2_new_unique(:,Reaction_f_index_i)-VTotal_condition2_new_unique(:,Reaction_b_index_i);
    end

end

% compute RatioMatrix and regulate
RatioMatrix = zeros(NumSolution_1*NumSolution_2,NumReactions);
nonzeros_index = zeros(NumSolution_1*NumSolution_2,1);
for i = 1:NumSolution_1
    for j = 1:NumSolution_2
        V_condition1_i = VTotal_condition1_new_unique_CombineReverse(i,:);
        V_condition2_j = VTotal_condition2_new_unique_CombineReverse(j,:);
        for k = 1:NumReactions
            V_condition1_i_k = V_condition1_i(k);
            V_condition2_j_k = V_condition2_j(k);
            if V_condition1_i_k*V_condition2_j_k>=0
            RatioMatrix(i*j,k) = V_condition2_j_k./(V_condition1_i_k+V_condition2_j_k);
            else
                RatioMatrix(i*j,k) = -V_condition2_j_k./(V_condition2_j_k-V_condition1_i_k);
            end
        end
        nonzeros_index(i*j) = 1;
    end
end
RatioMatrix = RatioMatrix(logical(nonzeros_index),:);
RatioMatrixTable = array2table(RatioMatrix);
RatioMatrixTable.Properties.VariableNames = ReactionIDs;
writetable(RatioMatrixTable,strcat(TracerID,'/FluxAnalysis/Filtering of estimates/RatioMatrix_',MappingID,'.xlsx'))

