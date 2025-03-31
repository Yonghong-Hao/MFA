clear;clc
%% input
TracerID = 'Tracer_glc_12';
MappingID = 'mapping9';
%DataID = 'data1_1'; % CM
%DataID = 'data1_2'; % LA+
DataID = 'data1_3'; % LA
Conditions = {'CM','LA+','LA'};
% modify 默认空集
AddingBoundaryNodes = {'cd'};

%% Decomposition of the network of EMU reactions
[SimplifiedEMUReactions, UnsimplifiedEMUReactions, STable, BoundaryProperty, MappingAtomTable, ToSimulateMS, ToSimulateMS_raw, UnvisitedEMUList, EMUsTotal] = S2_EMUDecomposition(TracerID,MappingID,DataID,AddingBoundaryNodes);

%% Stoichiometric matrix and boundary variables
ModelID = strcat(MappingID,'_',DataID);
% Stoichiometric matrix
MetabolitesinS = STable.Row;
SMatrix = STable.Variables;
Fluxes = STable.Properties.VariableNames;
NumFluxes = numel(Fluxes);

% the emu set to simulate
NumToSimulate = height(ToSimulateMS);
ToSimulateSplit = regexp(ToSimulateMS.EMUIDs,'_','split');
ToSimulateMets = cell(NumToSimulate,1);
for i = 1:NumToSimulate
    ToSimulateSplit_i = ToSimulateSplit{i};
    ToSimulateSplit_i(end)=[];
    ToSimulateMets{i} = strjoin(ToSimulateSplit_i,'_');
end
ToSimulateMets = unique(ToSimulateMets);

% Tracer
TracerFile = strcat(TracerID,'/Data_raw/LabelledPattern.xlsx');
TracerTable = readtable(TracerFile,"VariableNamingRule","preserve");

%% Get all EMU reactions
EMUReactionTable = UnsimplifiedEMUReactions;

% store equations (products) and variables (reactants) to the refernce
[VariablesRefernce, BoundaryReference, emuReactant, emuProduct] = S2_VariableReference(EMUReactionTable,TracerTable,UnvisitedEMUList); %table 这里还可以继续添加常变量的输入，例如被标记代谢物底物emu
%writetable(VariablesRefernce,strcat(TracerID,'/FluxBalance&EMUBalance/VariablesRefernce_',MappingID,'_',DataID,'.xlsx'));
% 将每个反应中emu扩展为mid of emu
EMUMIDTable = table('Size',[0 4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'ReactantMid','ProductMid','Flux','Size'});
for r = 1:height(EMUReactionTable)
    reactant_r = emuReactant{r};
    product_r = emuProduct{r};
    flux_r = EMUReactionTable.Flux(r);
    Size_r = EMUReactionTable.Size(r);
    % 判断是否属于condensation reaction
    NumReactant = numel(emuReactant{r});
    if NumReactant>1 % this is a condensation reaction
        midTerms = S2_CondensationReaction(reactant_r,product_r,flux_r,Size_r);
    elseif NumReactant==1 % this is a cleavage reaction or unimolecular reaction
        midTerms = S2_NonCondensationReaction(reactant_r,product_r,flux_r,Size_r);
    else
        disp('others');
    end
    EMUMIDTable = vertcat(EMUMIDTable,midTerms);
end

%% test
if ~isempty(setdiff(EMUMIDTable.ProductMid,VariablesRefernce.emu_mid))
    disp('WARNING: The following EMUs are in EMUMIDTable.ProductMid, VariablesRefernce.emu_mid')
    diff1 = setdiff(EMUMIDTable.ProductMid,VariablesRefernce.emu_mid);
    disp(diff1)
elseif ~isempty(setdiff(VariablesRefernce.emu_mid,EMUMIDTable.ProductMid))
    disp('WARNING: The following EMUs are in VariablesRefernce.emu_mid, not in EMUMIDTable.ProductMid')
    diff2 = setdiff(VariablesRefernce.emu_mid,EMUMIDTable.ProductMid);
    disp(diff2)
end

%% 按emu size分类写线性方程组
% 变量替换：v->V()
EMUMIDTable.Flux = cellfun(@(x) regexprep(x, 'v(\d+)', 'V($1)'),EMUMIDTable.Flux,'UniformOutput',false);
Fluxes = cellfun(@(x) regexprep(x, 'v(\d+)', 'V($1)'),Fluxes,'UniformOutput',false);

% 边界变量（包括卷积类型变量）
BoundaryVariablesRefernce = S2_BoundaryVariablesRefernce(EMUMIDTable,VariablesRefernce,BoundaryReference);
%writetable(BoundaryVariablesRefernce,strcat(TracerID,'/FluxBalance&EMUBalance/BoundaryVariablesRefernce_',MappingID,'_',DataID,'.xlsx'));
% 构建一个大的A, B矩阵
NumVariables = height(VariablesRefernce);
NumBoundaryVariables  = height(BoundaryVariablesRefernce);

% 将AMatrix, BMatrix系数存储为数值矩阵，实则dAdV,dBdV
dAdV = cell(NumVariables,NumVariables);
dBdV = cell(NumVariables,NumBoundaryVariables);

% 遍历EMUMIDTable中每个反应
% 生成项
for r = 1:height(EMUMIDTable)
    ReactantMid_r = EMUMIDTable.ReactantMid{r};
    ProductMid_r = EMUMIDTable.ProductMid{r};
    EqIndex = find(strcmp(VariablesRefernce.emu_mid,ProductMid_r));
    ProductEMU_r = VariablesRefernce.emu_emu{EqIndex};
    ProductMid_r_met = EMUsTotal.Mets{strcmp(EMUsTotal.EMUs,ProductEMU_r)};
    SimulatedVarIndex = find(strcmp(VariablesRefernce.emu_mid,ReactantMid_r));
    BoundaryVarIndex = find(strcmp(BoundaryVariablesRefernce.emu_mid,string(expand(str2sym(ReactantMid_r)))));
    if ~isempty(SimulatedVarIndex)
        FluxTerms = EMUMIDTable.Flux{r};
        SplitFluxTerms = regexp(FluxTerms,' \+ ','split');
        v = zeros(1,NumFluxes);
        for i = 1:numel(SplitFluxTerms)
            SplitCoefficient = regexp(SplitFluxTerms{i},'\d*\.*\d*(?=\*)','match');
            SplitFlux = regexp(SplitFluxTerms{i},'(?<=\*)V\(\d+\)','match');
            indexFlux = ismember(Fluxes,SplitFlux);
            StoichCoefficient = SMatrix(ismember(STable.Row,ProductMid_r_met),indexFlux);
            v(indexFlux) = v(indexFlux)+StoichCoefficient*str2double(SplitCoefficient);
        end
        dAdV{EqIndex,SimulatedVarIndex} = v;
    end
    if ~isempty(BoundaryVarIndex)
        FluxTerms = EMUMIDTable.Flux{r};
        SplitFluxTerms = regexp(FluxTerms,' \+ ','split');
        v = zeros(1,NumFluxes);
        for i = 1:numel(SplitFluxTerms)
            SplitCoefficient = regexp(SplitFluxTerms{i},'\d*\.*\d*(?=\*)','match');
            SplitFlux = regexp(SplitFluxTerms{i},'(?<=\*)V\(\d+\)','match');
            indexFlux = ismember(Fluxes,SplitFlux);
            StoichCoefficient = SMatrix(ismember(STable.Row,ProductMid_r_met),indexFlux);
            v(indexFlux) = v(indexFlux)+StoichCoefficient*str2double(SplitCoefficient);
        end
        dBdV{EqIndex,BoundaryVarIndex} = -v;
    end
end

% 消耗项（直接对相应行生成项取负）
for ii = 1:height(VariablesRefernce)
    fromAMatrix = dAdV(ii,~cellfun('isempty',dAdV(ii,:)));
    fromBMatrix = dBdV(ii,~cellfun('isempty',dBdV(ii,:)));
    if ~isempty(fromAMatrix)
        v = zeros(1,NumFluxes);
        for jj = 1:numel(fromAMatrix)
            v = v-fromAMatrix{jj};
        end
        fromAMatrix = v;
    else
        fromAMatrix = zeros(1,NumFluxes);
    end
    if ~isempty(fromBMatrix)
        v = zeros(1,NumFluxes);
        for jj = 1:numel(fromBMatrix)
            v = v+fromBMatrix{jj};
        end
        fromBMatrix = v;
    else
        fromBMatrix = zeros(1,NumFluxes);
    end
    dAdV{ii,ii} = fromAMatrix+fromBMatrix;
end

%% 变量名称
% Size
SizeVariables = VariablesRefernce.emu_size;
SizeBoundaryVariables = BoundaryVariablesRefernce.emu_size;
%X, Y变量名称
XVar = VariablesRefernce.emu_var;
YVar = BoundaryVariablesRefernce.emu_var;
%dXdV, dYdV的变量名称
dXdVVar = regexprep(VariablesRefernce.emu_var,'X','dXdV');
dYdVVar = S2_NonlinearGradient(NumBoundaryVariables,BoundaryVariablesRefernce);

%% save the vatiables
AllVariables.STable = STable;
AllVariables.ToSimulateMS = ToSimulateMS;
AllVariables.VariablesRefernce = VariablesRefernce;
AllVariables.BoundaryVariablesRefernce = BoundaryVariablesRefernce;
AllVariables.dAdV = dAdV;
AllVariables.dBdV = dBdV;
AllVariables.dXdVVar = dXdVVar;
AllVariables.dYdVVar = dYdVVar;