function [SimplifiedEMUReactions, UnsimplifiedEMUReactions, STable, BoundaryProperty, MappingAtomTable, ToSimulateMS, ToSimulateMS_raw, UnvisitedEMUList, EMUsTotal] = S2_EMUDecomposition(TracerID,MappingID,DataID,AddingBoundaryNodes)
%% Read and organize all raw information of models
% tracer (metabolite) -- Tracer
Tracer = strcat(regexp(TracerID,'(?<=_)\w+(?=_)','match'),'_tracer');

% property of boundary nodes -- BoundaryProperty
BoundaryProperty = table('Size',[0 2],'VariableTypes',{'string','string'},'VariableNames',{'Boundary nodes','Property'});
%BoundaryProperty = [BoundaryProperty;{Tracer,'tracer'}];
BoundaryProperty = [BoundaryProperty;Tracer,cellstr(repmat('tracer',numel(Tracer),1))];

% metabolic network & atom transformation -- MappingAtomTable
MappingAtomFile = strcat(TracerID,'/MappingAtom/',replace(MappingID,'mapping','mapping_atom_'),'.xlsx');
MappingAtomTable = readtable(MappingAtomFile,"VariableNamingRule","preserve");
MappingAtomTable = removevars(MappingAtomTable,"ReactionIDs");

% metabolic network (excluding atom transform) -- MappingTable
MappingAtomCell = [MappingAtomTable.FluxIDs,MappingAtomTable.("SubstrateIDs(atoms)"),MappingAtomTable.("ProductIDs(atoms)")];
MappingCell = cellfun(@(x) regexprep(x,'\([a-z]+\)|\([a-z]+\,[a-z]+\)',''),MappingAtomCell,'UniformOutput',false);
MappingTable = cell2table(MappingCell,"VariableNames",{'flux_ID','substrate_IDs','product_IDs'});

% total metabolites in the metabolic network -- MetabolitesTotal
Metabolites = regexp({MappingTable.substrate_IDs{:};MappingTable.product_IDs{:}},'+','split');
MetabolitesTotal = unique(regexprep([Metabolites{:}]','\d*\.*\d*\*',''));
NumMetabolites = numel(MetabolitesTotal);
NumFluxes = height(MappingAtomTable);

%% Expand metabolic network to EMU reaction network -- EMUsTotal
% extract all nodes in atom transformation -- MetaboliteAtomsTotal
MetaboliteAtoms = regexp({MappingAtomTable.("SubstrateIDs(atoms)"){:};MappingAtomTable.("ProductIDs(atoms)"){:}},'+','split');
MetaboliteAtomsTotal = unique(regexprep([MetaboliteAtoms{:}]','\d*\.*\d*\*',''));
indexDel = zeros(numel(MetaboliteAtomsTotal),1);
for i = 1:numel(MetaboliteAtomsTotal)
    if ~ismember('(',MetaboliteAtomsTotal{i})
        indexDel(i) = 1;
    end
end
MetaboliteAtomsTotal(logical(indexDel)) = [];

% organize the serial number of every atom in each EMU -- EMUsTotal (unassembled)
Mets = cellfun(@(x) regexp(x,'\w+(?=\()','match'),MetaboliteAtomsTotal,'UniformOutput',false);
Atoms = cellfun(@(y) regexprep(y,'\)',''),cellfun(@(x) regexp(x,'(?<=\().+','match'),MetaboliteAtomsTotal,'UniformOutput',false),'UniformOutput',false);
AtomsIndex = cell(height(Atoms),1); % the index of every atom in a EMU
AtomsNumber = AtomsIndex; % the serial numbers of every atom in a EMU
for n = 1:height(Atoms)
    size_n = length(char(regexprep(Atoms{n},',[a-z]+','')));
    AtomIndex_n = cell(size_n,1);
    AtomNumber_n = cell(size_n,1);
    for k = 1:size_n
        AtomIndex_nk = nchoosek(1:size_n,k);
        AtomIndex_n{k} = AtomIndex_nk;
        AtomNumber_nk = char([]);
        for i = 1:k
            AtomNumber_nk(:,i) = mydec2hex(AtomIndex_nk(:,i));
        end
        AtomNumber_n{k} = AtomNumber_nk;
    end
    AtomsIndex{n} = AtomIndex_n;
    AtomsNumber{n} = AtomNumber_n;
end
EMUsTotal = table(MetaboliteAtomsTotal,Mets,Atoms,AtomsIndex,AtomsNumber);

% assembel EMUsTotal -- EMUsTotal
% EMUsTotal.EMUs -- EMUs
EMUs = {};
for n = 1:height(EMUsTotal)
    AtomNumber_0 = {};
    for k = 1:height(EMUsTotal.AtomsNumber{n})
        AtomNumber_0 = [AtomNumber_0;EMUsTotal.AtomsNumber{n}(k)];
        EMUs_new = strcat(EMUsTotal.Mets{n},'_',AtomNumber_0{k});
        EMUs = unique([EMUs;EMUs_new]);
    end
end
% EMUsTotal.Mets -- Mets
Mets = cellfun(@(x) regexp(x,'\w+(?=\_[0-9])|\w+(?=\_[A-G])','match'),EMUs,'UniformOutput',false);
Mets1 = cellfun(@(x) strcat('(?<=',x,'_',').+'),Mets,'UniformOutput',false);
% EMUsTotal.Atoms
Atoms = cellfun(@(emu,meta) regexp(emu,meta,'match'),EMUs,Mets1,'UniformOutput',false);
Atoms = [Atoms{:}]';
% EMUsTotal.Size
Size = [];
for i = 1:height(EMUs)
    Size(i) = length(char(Atoms{i}));
end
Size = Size';
% EMUsTotal (assembled)
EMUsTotal = table(Mets,EMUs,Atoms,Size);

%% Obtain S matrix
% assembel the raw S matrix -- STable_raw
STable_raw = table('Size',[NumMetabolites NumFluxes],'VariableTypes',string(repmat('double',NumFluxes,1)),'VariableNames',MappingAtomTable.FluxIDs,'RowNames',MetabolitesTotal);
STable_raw_variables = STable_raw.Variables;
for i = 1:NumFluxes
    % reactant
    ReactantIterms = Metabolites{1,i};
    for j = 1:numel(ReactantIterms)
        ReactantTerms_j = regexp(ReactantIterms{j},'(\d*\.*\d*)\*(\w+)','tokens');
        if isempty(ReactantTerms_j)
            RIndex = find(strcmp(MetabolitesTotal,ReactantIterms{j}));
            STable_raw_variables(RIndex,i) = STable_raw_variables(RIndex,i)-1;
        else
            RCoefficient = ReactantTerms_j{1}{1};
            RIndex = find(strcmp(MetabolitesTotal,ReactantTerms_j{1}{2}));
            STable_raw_variables(RIndex,i) = STable_raw_variables(RIndex,i)-str2double(RCoefficient);
        end
    end
    % product
    ProductIterms = Metabolites{2,i};
    for j = 1:numel(ProductIterms)
        ProductTerms_j = regexp(ProductIterms{j},'(\d*\.*\d*)\*(\w+)','tokens');
        if isempty(ProductTerms_j)
            PIndex = find(strcmp(MetabolitesTotal,ProductIterms{j}));
            STable_raw_variables(PIndex,i) = STable_raw_variables(PIndex,i)+1;
        else
            PCoefficient = ProductTerms_j{1}{1};
            PIndex = find(strcmp(MetabolitesTotal,ProductTerms_j{1}{2}));
            STable_raw_variables(PIndex,i) = STable_raw_variables(PIndex,i)+str2double(PCoefficient);
        end
    end
end
STable_raw.Variables = STable_raw_variables;

%% Analyze the properties of all nodes in the metabolic network according the raw S matrix
% only output (no limit to the number of arrows output, only the type) -- Source
Source = {};
% only input (no limit to the number of arrows input, only the type) -- Sink
Sink = {};
% single-in & single-out -- Source_Sink
Source_Sink = {};
% no fluxes
NoFluxes = {};
for i = 1:NumMetabolites
    Row_i = STable_raw{i,:};
    count1 = sum(Row_i>0);
    count2 = sum(Row_i<0);
    if count1==0 && count2>=1
        Source = [Source;STable_raw.Row{i}];
    elseif count1>=1 && count2==0
        Sink = [Sink;STable_raw.Row{i}];
    elseif count1==1 && count2==1
        Source_Sink = [Source_Sink;STable_raw.Row{i}];
    elseif count1==0 && count2==0
        NoFluxes = [NoFluxes;STable_raw.Row{i}];
    end
end
% add source-nodes into BoundaryProperty
BoundaryProperty = [BoundaryProperty;Source,cellstr(repmat('source',numel(Source),1))];
% add sink-nodes into BoundaryProperty
BoundaryProperty = [BoundaryProperty;Sink,cellstr(repmat('sink',numel(Sink),1))];
% add nodes that have the same one-way input and output -- Source_Sink_pairs
Source_Sink_pair = {};
DeleteNodes = table('Size',[0 3],'VariableTypes',{'string','string','string'},'VariableNames',{'InputFlux','DelNode','OutputFlux'});
NodesMets = STable_raw.Row;
Fluxes = STable_raw.Properties.VariableNames;
Variables = STable_raw.Variables;
for i = 1:numel(Source_Sink)
    DelNode = Source_Sink{i};
    indexNodes = logical(strcmp(NodesMets,DelNode));
    indexOutputFlux = Variables(indexNodes,:)<0;
    indexToNode = Variables(:,indexOutputFlux)>0;
    ToNode = NodesMets(indexToNode);
    indexInputFlux = Variables(indexNodes,:)>0;
    indexFromNode = Variables(:,indexInputFlux)<0;
    FromNode = NodesMets(indexFromNode);
    if ~isempty(setdiff(ToNode,FromNode)) && ~isempty(setdiff(FromNode,ToNode))
        DeleteNodes = [DeleteNodes;{Fluxes{indexInputFlux},DelNode,Fluxes{indexOutputFlux}}];
    else
        Source_Sink_pair = [Source_Sink_pair;{DelNode,ToNode}];
    end
end
if ~isempty(Source_Sink_pair)
    BoundaryProperty = [BoundaryProperty;Source_Sink_pair(:,1),cellstr(repmat('source-sink pair',height(Source_Sink_pair),1))];
end
% add AddingBoundaryNodes -- Source_Sink_pairs
if ~isempty(AddingBoundaryNodes)
BoundaryProperty = [BoundaryProperty;AddingBoundaryNodes,cellstr(repmat('addition',numel(AddingBoundaryNodes),1))];
end
% add NoFluxes nodes into BoundaryProperty
if ~isempty(NoFluxes)
    BoundaryProperty = [BoundaryProperty;NoFluxes,cellstr(repmat('no fluxes',height(NoFluxes),1))];
end

%% Remove the metabolites that are not included in the equilibrium constraint -- STable & SMatrix
BoundaryNodes = unique(BoundaryProperty.("Boundary nodes"));
[BalancedNodes,BalanceIndex] = setdiff(MetabolitesTotal,BoundaryNodes);
STable = STable_raw(BalanceIndex,:);
SMatrix = table2array(STable);

%% Set metabolite symmetry
% atomic equivalence relation -- EquivalentTable
SplitMetAtom = cellfun(@(x) regexprep(x,'\)',''),regexp(MetaboliteAtomsTotal,'\(|,','split'),'UniformOutput',false);
EquivalentTable = {};
for i = 1:numel(MetaboliteAtomsTotal)
    MetAtom = SplitMetAtom{i};
    Mets = MetAtom(1);
    Atoms = SplitMetAtom{i}{2};
    Mets = repmat(Mets,length(char(Atoms)),1);
    AtomsID_1 = cellstr(mydec2hex(1:length(char(Atoms))));
    if numel(MetAtom)==2
        AtomsID_2 = AtomsID_1;
    elseif numel(MetAtom)==3
        Atoms = SplitMetAtom{i}{3};
        AtomsID_2 = cellstr(mydec2hex(length(char(Atoms)):-1:1));
    end
    EquivalentAtoms = table(Mets,AtomsID_1,AtomsID_2,'VariableNames',{'Mets','AtomsID_1','AtomsID_2'});
    EquivalentTable = [EquivalentTable;EquivalentAtoms];
end
EquivalentTable = unique(EquivalentTable);

% equivalent EMUs -- EMUsTotal
EquivalentAtoms = cell(height(EMUsTotal),1);
for i = 1:height(MetabolitesTotal)
    if ismember(MetabolitesTotal(i),[EMUsTotal.Mets{:}])
        IndexInEMUsTotal = ismember([EMUsTotal.Mets{:}]',MetabolitesTotal(i));
        IndexInEquivalentTable = ismember(EquivalentTable.Mets,MetabolitesTotal(i));
        EquivalentAtoms(IndexInEMUsTotal) = replace([EMUsTotal.Atoms{IndexInEMUsTotal}]',EquivalentTable.AtomsID_1(IndexInEquivalentTable),EquivalentTable.AtomsID_2(IndexInEquivalentTable));
    end
end
EquivalentAtoms = cellfun(@(x) sort(x),EquivalentAtoms,'UniformOutput',false);
EquivalentEMUs = cellfun(@(x,y) strcat(x,'_',y),EMUsTotal.Mets,EquivalentAtoms);
EMUsTotal = addvars(EMUsTotal,EquivalentEMUs,EquivalentAtoms,'Before','Size');
EquivalentTerms = sort(string([EMUsTotal.EMUs,EMUsTotal.EquivalentEMUs]),2);
[~,iEquivalentTerms,~] = unique(EquivalentTerms,'rows');
EMUsTotal = EMUsTotal(iEquivalentTerms,:);

%% Set data set（net flux measurements, MS measurements (to simulate set and MID data), tracer labeled patterns）
% MS measurements
ID = regexprep(DataID,'data','');
ToSimulateFile = strcat(TracerID,'/Data_raw/ToSimulateMS_',ID,'_',MappingID,'.xlsx');
ToSimulateMS_raw = readtable(ToSimulateFile, "VariableNamingRule","preserve");
ToSimulateMet = regexp(ToSimulateMS_raw.EMUIDs,'\w+(?=_\d+|_[A-Z]+)','match');
NumToSimulate_raw = height(ToSimulateMS_raw);
ToSimulateAtom = cell(NumToSimulate_raw,1);
for i = 1:NumToSimulate_raw
ToSimulateAtom{i} = regexprep(regexprep(ToSimulateMS_raw.EMUIDs{i},ToSimulateMet{i},''),'_','');
end
ToSimulateMS_raw = addvars(ToSimulateMS_raw,ToSimulateMet,ToSimulateAtom,'Before','Size');

% take the intersection of ToSimulate set and emu set in the model
ia = ismember(cellfun(@(x) char(x),ToSimulateMS_raw.ToSimulateMet,'UniformOutput',false),STable_raw.Row);
ToSimulateMS = ToSimulateMS_raw(ia,:);

% delete the emus which can not be simulated（Source, AddedBoundary, Tracer, Source_Sink_pair）
if ~isempty(Source_Sink_pair)
    ib = ~ismember(cellfun(@(x) char(x),ToSimulateMS.ToSimulateMet,'UniformOutput',false),[Source;AddingBoundaryNodes;Tracer;Source_Sink_pair(:,1)]);
else
    ib = ~ismember(cellfun(@(x) char(x),ToSimulateMS.ToSimulateMet,'UniformOutput',false),[Source;AddingBoundaryNodes;Tracer]);
end
ToSimulateMS = ToSimulateMS(ib,:);
%ia = ismember(cellfun(@(x) char(x),ToSimulateTable_raw.ToSimulateMet,'UniformOutput',false),BalancedNodes);

NumToSimulate = height(ToSimulateMS);

%% Set the start set of EMU decomposition（initial set is ToSimulateTable）and the end set（initial set is BoundaryNodes）
% the start set
UnvisitedEMUList = table('Size',[0 4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'EMUs','Mets','Atoms','Size'});
ToSimulateTable_copy = removevars(ToSimulateMS,{'Measured','Simulated','Std'});
%ToSimulateTable_copy = removevars(ToSimulateMS,{'Measured','Simulated'});
ToSimulateTable_copy.Properties.VariableNames{'EMUIDs'} = 'EMUs';
ToSimulateTable_copy.Properties.VariableNames{'ToSimulateMet'} = 'Mets';
ToSimulateTable_copy.Properties.VariableNames{'ToSimulateAtom'} = 'Atoms';
UnvisitedEMUList = vertcat(UnvisitedEMUList,ToSimulateTable_copy);

% the end set
VisitedEMUList = {};
StoppedNodes = setdiff(BoundaryNodes,Sink);
for i = 1:height(StoppedNodes)
    index = strcmp([EMUsTotal.Mets{:}]',StoppedNodes{i});
    VisitedEMUList = [VisitedEMUList;table(EMUsTotal.EMUs(index,:),EMUsTotal.Mets(index,:),EMUsTotal.Atoms(index,:),EMUsTotal.Size(index,:),'VariableNames',{'EMUs','Mets','Atoms','Size'})];
end


%% EMU Decomposition
EMUReactionTable = table('Size',[0 3],'VariableTypes',{'string','string','double'},'VariableNames',{'Reaction','Flux','Size'});
while true
    % 从UnvisitedEMUList中选出size最大的EMU -- EMU
    [MaximumSize,IndexMax] = max(UnvisitedEMUList.Size);
    %%%%%%%%%%%%%%%%%%%%%%%%%% Break condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if MaximumSize<0
        break
    end
    EMU = UnvisitedEMUList.EMUs{IndexMax};

    % 找出EMU的equivalent EMU in EMUsTotal -- EquivalentEMU
    IndexEquivalent = find(strcmp(EMUsTotal.EMUs,EMU));
    if isempty(IndexEquivalent)
        IndexEquivalent = find(strcmp(EMUsTotal.EquivalentEMUs,EMU));
        EMU = EMUsTotal.EMUs{IndexEquivalent};
        EquivalentEMU = EMUsTotal.EquivalentEMUs{IndexEquivalent};
    else
        EquivalentEMU = EMUsTotal.EquivalentEMUs{IndexEquivalent};
    end

    % 若EMU不在VisitedEMUList中，则开始模拟该EMU
    if ~ismember(EMU,VisitedEMUList.EMUs)
        VisitedEMUList = unique([VisitedEMUList;UnvisitedEMUList(IndexMax,:)]);
        if strcmp(EMU,EquivalentEMU)
            %%%%%%%%%%%%%%%%若EMU的等价EMU是自己，则开始不含等价EMU的模拟流程%%%%%%%%%%%%%%%%
            % 得到ToSimulate的信息（Met, Atom） -- ToSimulateTerms
            ToSimulateTerms = UnvisitedEMUList(IndexMax,:);
            ToSimulateEMU = ToSimulateTerms.EMUs{:};
            ToSimulateMet = ToSimulateTerms.Mets{:};
            ToSimulateAtomID = ToSimulateTerms.Atoms{:};
            ToSimulateSize = ToSimulateTerms.Size;

            % 根据生成物，在MappingCell中，溯源其所在反应 -- IndexProducts (vector)
            SplitMappingCell = cellfun(@(x) regexp(x,'+','split'),MappingCell,'UniformOutput',false); % 含反应系数
            IndexProducts = logical(cell2mat(cellfun(@(x) sum(x),cellfun(@(x) strcmp(cellfun(@(y) regexprep(y,'\d*\.*\d*\*',''),x,'UniformOutput',false),ToSimulateMet),SplitMappingCell(:,3),'UniformOutput',false),'UniformOutput',false)));

            % 根据索引，在MappingAtomCell中，找到对应项 -- ReactionTerms, ReactantTerms, ProductTerms
            ReactionIterms = MappingAtomCell(IndexProducts,1);
            ReactantIterms = MappingAtomCell(IndexProducts,2);
            ProductIterms = MappingAtomCell(IndexProducts,3);

            % 遍历上述每个反应
            for i = 1:sum(IndexProducts)
                % 若EMU含等价EMU，删除等价原子字符 -- ReactantEMU_i
                ReactantEMU_i = ReactantIterms{i};
                if ismember(',',ReactantEMU_i)
                    ReactantEMU_i = regexprep(ReactantEMU_i,',[a-z]+','');
                end

                % 生成物项，及其反应系数-- ProductEMU_i, StoichiometricNumber_P
                ProductEMU_i = ProductIterms{i};
                patternStoichiometricNumber = strcat('\d*\.*\d*(?=*',ToSimulateMet,'\()');
                StoichiometricNumber_P = str2double(regexp(ProductEMU_i,patternStoichiometricNumber,'match'));
                if isempty(StoichiometricNumber_P)
                    StoichiometricNumber_P = 1;
                end

                % 从ProductEMU_i中获取要模拟的EMU对应代谢物的所有原子 -- ProductAtomAll
                % 从ToSimulateAtomID中获取要模拟的EMU的原子位置索引 -- IndexProductAtom
                patternAtom = strcat('(?<=',ToSimulateMet,'\()[a-z]+\)|(?<=',ToSimulateMet,'\()[a-z]+,[a-z]+\)');
                ProductAtomAll = char(regexprep(regexp(ProductEMU_i,patternAtom,'match'),'\)',''));
                NumProductAtomAll = width(ProductAtomAll);
                NumSameProduct = height(ProductAtomAll);
                IndexProductAtom = zeros(NumProductAtomAll,1);
                for j = 1:NumProductAtomAll
                    IndexProductAtom(j) = contains(ToSimulateAtomID,mydec2hex(j));
                end

                % NumSameProduct考虑到了像'glc(abcdef)->pyr(abc)+pyr(def)'的情况
                for jj = 1:NumSameProduct
                    % 根据要模拟的EMU的原子索引IndexProductAtom，从ProductAtomAll中提取相应原子 -- ProductAtoms
                    ProductAtomAll_jj = ProductAtomAll(jj,:);
                    ProductAtoms = ProductAtomAll_jj(logical(IndexProductAtom));

                    % 对于要模拟的EMU中的每个原子，在相应反应物EMU中找到对应原子 -- ReactantEMUTable
                    ReactantEMUTable = table('Size',[0 3],'VariableTypes',{'double','string','string'},'VariableNames',{'StoichiometricNumber','Met','AtomID'});
                    for k = 1:ToSimulateSize
                        % 遍历ProductAtoms中每个原子，根据该原子找到对应反应物 -- Met
                        patternMet = strcat('\w+(?=\([a-z]*',ProductAtoms(k),'[a-z]*\))');
                        Met = regexp(ReactantEMU_i,patternMet,'match');

                        % 根据Met，找到反应物对应的化学计量数 -- StoichiometricNumber_R
                        patternStoichiometricNumber = strcat('\d*\.*\d*(?=*',Met{:},')');
                        StoichiometricNumber_R = str2double(regexp(ReactantEMU_i,patternStoichiometricNumber,'match'));
                        if isempty(StoichiometricNumber_R) || StoichiometricNumber_R==0
                            StoichiometricNumber_R = 1;
                        end

                        % 根据Met，找到对应反应物的原子 -- Atoms
                        patternAtoms = strcat('(?<=',Met,'\()[a-z]+');
                        Atoms = regexp(ReactantEMU_i,patternAtoms,'match');

                        % 据'numel(Atoms{1})>1'判断反应物是否存在像'ru5p(abcde)+ru5p(fghij)->'的情况
                        % 使ProductAtoms(k)对应到正确的项，并提取对应原子次序编号 -- AtomID
                        if numel(Atoms{1})>1
                            for kk = 1:numel(Atoms{1})
                                atoms_kk = strfind(char(Atoms{1}{kk}),ProductAtoms(k));
                                if ~isempty(atoms_kk)
                                    AtomID = mydec2hex(atoms_kk);
                                    Met = strcat('(',num2str(kk),')',Met);
                                end
                            end
                        else
                            AtomID = mydec2hex(strfind(char(Atoms{:}),ProductAtoms(k)));
                        end
                        % 遍历ProductAtoms中每个原子，溯源到Reactant EMU, 并将StoichiometricNumber_R, Met, AtomID存在ReactantEMUTable
                        ReactantEMUTable = [ReactantEMUTable;{StoichiometricNumber_R,Met,AtomID}];
                    end

                    % 合并ReactantEMUTable中同类项 -- ReactantEMUTableGroup
                    ReactantEMUTableGroup = groupsummary(ReactantEMUTable,["Met", "StoichiometricNumber"],@(x) strjoin(x,''),"AtomID");
                    ReactantEMUTableGroup.fun1_AtomID = cellfun(@(x) sort(x),ReactantEMUTableGroup.fun1_AtomID,'UniformOutput',false);
                    ReactantEMUTableGroup.Met = regexprep(ReactantEMUTableGroup.Met,'\(.+\)','');
                    % copy of Met -- Met_copy
                    Met_copy = regexprep(ReactantEMUTableGroup.Met,'\(.+\)','');

                    % % 遍历反应物中每一项（当项数大于1时为聚合反应）
                    % for m = 1:height(ReactantEMUTableGroup)
                    %     % 若反应物化学计量数不为1，则用*连接 (e.g. 2*glu)
                    %     if ReactantEMUTableGroup.StoichiometricNumber(m)~=1
                    %         ReactantEMUTableGroup.Met(m) = strcat(num2str(ReactantEMUTableGroup.StoichiometricNumber(m)),'*',ReactantEMUTableGroup.Met);
                    %     end
                    % end

                    % 根据ReactantEMUTableGroup中每一项，合并成相应EMU -- ReactantEMU
                    ReactantEMU = strcat(ReactantEMUTableGroup.Met,'_',ReactantEMUTableGroup.fun1_AtomID);
                    ReactantAtomID = ReactantEMUTableGroup.fun1_AtomID;
                    % copy of ReactantEMU -- ReactantEMU_copy
                    ReactantEMU_copy = strcat(Met_copy,'_',ReactantEMUTableGroup.fun1_AtomID);

                    % 遍历ReactantEMU中每个EMU
                    for ii = 1:numel(ReactantEMU)
                        % 若该EMU有等价EMU，只保留EMUsTotal.EMU项
                        if ~ismember(ReactantEMU{ii},EMUsTotal.EMUs)
                            iReactantEMU = strcmp(EMUsTotal.EquivalentEMUs,ReactantEMU{ii});
                            ReactantEMU{ii} = EMUsTotal.EMUs{iReactantEMU};
                            ReactantAtomID{ii} = EMUsTotal.Atoms{iReactantEMU};
                            ReactantEMU_copy(ii) = strcat(Met_copy,'_',EMUsTotal.Atoms{iReactantEMU}); %不含化学计量数
                        end
                        % 若反应物化学计量数不为1，则用*连接 (e.g. 2*glu)
                        if ReactantEMUTableGroup.StoichiometricNumber(ii)~=1
                            ReactantEMU(ii) = strcat(num2str(ReactantEMUTableGroup.StoichiometricNumber(ii)),'*',ReactantEMU);
                        end
                    end

                    % 若该反应为聚合反应，将多个反应物串联成一项，用逗号连接 -- ReactantEMU
                    if numel(ReactantEMU)>1
                        ReactantEMU = strcat('(',strjoin(ReactantEMU,', '),')'); %含化学计量数
                    end

                    % 将反应物与生成物串联，用@连接 -- ReactionEMU
                    ReactionEMU = strcat(ReactantEMU,'@',ToSimulateEMU);
                    Flux = strcat(num2str(1/StoichiometricNumber_P/NumSameProduct),'*',ReactionIterms{i});
                    Size = ToSimulateSize;

                    % 将新的EMU反应加入EMUReactionTable
                    EMUReactionTable = [EMUReactionTable;{ReactionEMU,Flux,Size}];
                    % 若UnvisitedEMUList中没有新溯源的EMU，则将新的EMU反应物加入UnvisitedEMUList
                    UnvisitedEMU_new = table(ReactantEMU_copy,Met_copy,ReactantAtomID,ReactantEMUTableGroup.GroupCount,'VariableNames',{'EMUs','Mets','Atoms','Size'});
                    for ii = 1:height(UnvisitedEMU_new)
                        UnvisitedEMU_new_ii = UnvisitedEMU_new(ii,:);
                        if ~ismember(UnvisitedEMU_new_ii.EMUs,UnvisitedEMUList.EMUs)
                            UnvisitedEMUList = vertcat(UnvisitedEMUList,UnvisitedEMU_new_ii);
                        end
                    end
                end
            end
        else
            %%%%%%%%%%%%%%%%若EMU的等价EMU不是自己，则开始含等价EMU的模拟流程%%%%%%%%%%%%%%%%
            % % 从EMUsTotal提取等价EMU项 -- EquivalentTerm
            % EquivalentTerm = EMUsTotal(IndexEquivalent,:);
            % EquivalentTerm = table(EquivalentTerm.EquivalentEMUs,EquivalentTerm.Mets,EquivalentTerm.EquivalentAtoms,EquivalentTerm.Size,'VariableNames',{'EMUs','Mets','Atoms','Size'});

            % 得到ToSimulate的信息（Met, Atom） -- ToSimulateTerms
            ToSimulateTerms = UnvisitedEMUList(IndexMax,:);
            ToSimulateEMU = ToSimulateTerms.EMUs{:};
            ToSimulateMet = ToSimulateTerms.Mets{:};
            ToSimulateAtomID = ToSimulateTerms.Atoms{:};
            ToSimulateSize = ToSimulateTerms.Size;

            % 根据生成物，在MappingCell中，溯源其所在反应 -- IndexProducts (vector)
            SplitMappingCell = cellfun(@(x) regexp(x,'+','split'),MappingCell,'UniformOutput',false);
            IndexProducts = logical(cell2mat(cellfun(@(x) sum(x),cellfun(@(x) strcmp(cellfun(@(y) regexprep(y,'\d*\.*\d*\*',''),x,'UniformOutput',false),ToSimulateMet),SplitMappingCell(:,3),'UniformOutput',false),'UniformOutput',false)));

            % 根据索引，在MappingAtomCell中，找到对应项 -- ReactionTerms, ReactantTerms, ProductTerms
            ReactionIterms = MappingAtomCell(IndexProducts,1);
            ReactantIterms = MappingAtomCell(IndexProducts,2);
            ProductIterms = MappingAtomCell(IndexProducts,3);

            % 遍历上述每个反应
            for i = 1:sum(IndexProducts)
                % 若反应物EMU含等价EMU，删除等价原子字符 -- ReactantEMU_i
                ReactantEMU_i = ReactantIterms{i};
                if ismember(',',ReactantEMU_i)
                    ReactantEMU_i = regexprep(ReactantEMU_i,',[a-z]+','');
                end

                % 生成物项，及其反应系数-- ProductEMU_i, StoichiometricNumber_P
                ProductEMU_i = ProductIterms{i};
                patternStoichiometricNumber = strcat('\d*\.*\d*(?=*',ToSimulateMet,'\()');
                StoichiometricNumber_P = str2double(regexp(ProductEMU_i,patternStoichiometricNumber,'match'));
                if isempty(StoichiometricNumber_P)
                    StoichiometricNumber_P = 1;
                end

                % 从ProductEMU_i中获取要模拟的EMU对应代谢物的所有原子，并用','将对称的原子分隔开 -- ProductAtomAll
                % 从ToSimulateAtomID中获取要模拟的EMU的原子位置索引 -- IndexProductAtom
                patternAtom = strcat('(?<=',ToSimulateMet,'\()[a-z]+\)|(?<=',ToSimulateMet,'\()[a-z]+,[a-z]+\)');
                ProductAtomAll = regexp(char(regexprep(regexp(ProductEMU_i,patternAtom,'match'),'\)','')),',','split');
                NumProductAtomAll = numel(ProductAtomAll{1});
                IndexProductAtom = zeros(NumProductAtomAll,1);
                % *******这里没有考虑'glc(abcdef)->pyr(abc)+pyr(def)'的情况*******
                for j = 1:NumProductAtomAll
                    IndexProductAtom(j) = contains(ToSimulateAtomID,mydec2hex(j));
                end

                % 根据要模拟的EMU的原子索引IndexProductAtom，从ProductAtomAll中提取相应原子 -- ProductAtoms
                % 含等价原子，且这里只考虑旋转180度的对称（jj=1:2）
                ProductAtoms = cell(1,2);
                for jj = 1:2
                    IndexProductAtom_jj = ProductAtomAll{jj};
                    ProductAtoms{jj} = IndexProductAtom_jj(logical(IndexProductAtom));
                end

                % 依次考虑原原子序列和对称原子序列
                for jj = 1:2
                    ProductAtoms_jj = ProductAtoms{jj};
                    % 对于要模拟的EMU中的每个原子，在相应反应物EMU中找到对应原子 -- ReactantEMUTable
                    ReactantEMUTable = table('Size',[0 3],'VariableTypes',{'double','string','string'},'VariableNames',{'StoichiometricNumber','Met','AtomID'});
                    for k = 1:ToSimulateSize
                        % 遍历ProductAtoms中每个原子，根据该原子找到对应反应物 -- Met
                        patternMet = strcat('\w+(?=\([a-z]*',ProductAtoms_jj(k),'[a-z]*\))');
                        Met = regexp(ReactantEMU_i,patternMet,'match');

                        % 根据Met，找到反应物对应的化学计量数 -- StoichiometricNumber_R
                        patternStoichiometricNumber = strcat('\d*\.*\d*(?=*',Met{:},')');
                        StoichiometricNumber_R = str2double(regexp(ReactantEMU_i,patternStoichiometricNumber,'match'));
                        if isempty(StoichiometricNumber_R) || StoichiometricNumber_R==0
                            StoichiometricNumber_R = 1;
                        end

                        % 根据Met，找到对应反应物的原子 -- Atoms
                        patternAtoms = strcat('(?<=',Met,'\()[a-z]+');
                        Atoms = regexp(ReactantEMU_i,patternAtoms,'match');

                        % *******这里没有考虑否存在像'ru5p(abcde)+ru5p(fghij)->'的情况*******
                        % 使ProductAtoms(k)对应到正确的项，并提取对应原子次序编号 -- AtomID
                        AtomID = mydec2hex(strfind(char(Atoms{:}),ProductAtoms_jj(k)));

                        % 遍历ProductAtoms中每个原子，溯源到Reactant EMU, 并将StoichiometricNumber_R, Met, AtomID存在ReactantEMUTable
                        ReactantEMUTable = [ReactantEMUTable;{StoichiometricNumber_R,Met,AtomID}];
                    end

                    % 合并ReactantEMUTable中同类项 -- ReactantEMUTableGroup
                    ReactantEMUTableGroup = groupsummary(ReactantEMUTable,["Met", "StoichiometricNumber"],@(x) strjoin(x,''),"AtomID");
                    ReactantEMUTableGroup.fun1_AtomID = cellfun(@(x) sort(x),ReactantEMUTableGroup.fun1_AtomID,'UniformOutput',false);
                    % copy of Met -- Met_copy
                    Met_copy = ReactantEMUTableGroup.Met;

                    % % 遍历反应物中每一项（当项数大于1时为聚合反应）
                    % for m = 1:height(ReactantEMUTableGroup)
                    %     %若反应物化学计量数不为1，则用*连接 (e.g. 2*glu)
                    %     if ReactantEMUTableGroup.StoichiometricNumber(m)~=1
                    %         ReactantEMUTableGroup.Met(m) = strcat(num2str(ReactantEMUTableGroup.StoichiometricNumber(m)),'*',ReactantEMUTableGroup.Met);
                    %     end
                    % end

                    % 根据ReactantEMUTableGroup中每一项，合并成相应EMU -- ReactantEMU
                    ReactantEMU = strcat(ReactantEMUTableGroup.Met,'_',ReactantEMUTableGroup.fun1_AtomID);
                    ReactantAtomID = ReactantEMUTableGroup.fun1_AtomID;
                    % copy of ReactantEMU -- ReactantEMU_copy
                    ReactantEMU_copy = strcat(Met_copy,'_',ReactantEMUTableGroup.fun1_AtomID);

                    % 遍历ReactantEMU中每个EMU
                    for ii = 1:numel(ReactantEMU)
                        % 若该EMU有等价EMU，只保留EMUsTotal.EMU项
                        if ~ismember(ReactantEMU{ii},EMUsTotal.EMUs)
                            iReactantEMU = strcmp(EMUsTotal.EquivalentEMUs,ReactantEMU{ii});
                            ReactantEMU{ii} = EMUsTotal.EMUs{iReactantEMU};
                            ReactantAtomID{ii} = EMUsTotal.Atoms{iReactantEMU};
                            ReactantEMU_copy(ii) = strcat(Met_copy,'_',EMUsTotal.Atoms{iReactantEMU}); %不含化学计量数
                        end
                        % 若反应物化学计量数不为1，则用*连接 (e.g. 2*glu)
                        if ReactantEMUTableGroup.StoichiometricNumber(ii)~=1
                            ReactantEMU(ii) = strcat(num2str(ReactantEMUTableGroup.StoichiometricNumber(ii)),'*',ReactantEMU);
                        end
                    end

                    % 若该反应为聚合反应，将多个反应物串联成一项，用逗号连接 -- ReactantEMU
                    if numel(ReactantEMU)>1
                        ReactantEMU = strcat('(',strjoin(ReactantEMU,', '),')'); %含化学计量数
                    end

                    % 将反应物与生成物串联，用@连接 -- ReactionEMU
                    ReactionEMU = strcat(ReactantEMU,'@',ToSimulateEMU);
                    Flux = strcat(num2str(1/2/StoichiometricNumber_P),'*',ReactionIterms{i});
                    Size = ToSimulateSize;

                    % 将新的EMU反应加入EMUReactionTable
                    EMUReactionTable = [EMUReactionTable;{ReactionEMU,Flux,Size}];
                    % 若UnvisitedEMUList中没有新溯源的EMU，则将新的EMU反应物加入UnvisitedEMUList
                    UnvisitedEMU_new = table(ReactantEMU_copy,Met_copy,ReactantAtomID,ReactantEMUTableGroup.GroupCount,'VariableNames',{'EMUs','Mets','Atoms','Size'});
                    for ii = 1:height(UnvisitedEMU_new)
                        UnvisitedEMU_new_ii = UnvisitedEMU_new(ii,:);
                        if ~ismember(UnvisitedEMU_new_ii.EMUs,UnvisitedEMUList.EMUs)
                            UnvisitedEMUList = vertcat(UnvisitedEMUList,UnvisitedEMU_new_ii);
                        end
                    end
                end
            end
            UnvisitedEMUList.Size(ismember(UnvisitedEMUList.EMUs,EquivalentEMU))=-inf; %以防万一
        end
        UnvisitedEMUList.Size(IndexMax)=-1; %finish the simulation
    else
        UnvisitedEMUList.Size(IndexMax)=-2; %have been simulated
    end
end
% 就此得到的所有模拟EMU反应，作为未简化的EMU Reactions输出 -- UnsimplifiedEMUReactions
UnsimplifiedEMUReactions = EMUReactionTable;

%% 进一步对EMUReactionTable简化，得到SimplifiedEMUReactions
% 根据DeleteNodes，从EMUReactionTable中所有相关流 -- FluxesAll
FluxesAll = regexp(EMUReactionTable.Flux,'(?<=\*).+','match');

% 根据FluxesAll，从EMUReactionTable中提取相关行 -- DeleteEMUReactions
% 记录所属行在EMUReactionTable中的索引 -- indexDel1
DeleteEMUReactions = {};
indexDel1 = zeros(height(EMUReactionTable),1);
for i = 1:height(DeleteNodes)
    indexInputFluxAll = cell2mat(cellfun(@(x) strcmp(x,DeleteNodes.InputFlux{i}),FluxesAll,'UniformOutput',false));
    indexOutputFluxAll = cell2mat(cellfun(@(x) strcmp(x,DeleteNodes.OutputFlux{i}),FluxesAll,'UniformOutput',false));
    if sum(indexInputFluxAll)~=0 && sum(indexOutputFluxAll)~=0
        DeleteEMUReactions = [DeleteEMUReactions;EMUReactionTable(logical(indexInputFluxAll+indexOutputFluxAll),:)];
        indexDel1(logical(indexInputFluxAll),:) = 1;
        indexDel1(logical(indexOutputFluxAll),:) = 1;
    end
end

% 若可删除节点非空（sum(indexDel1)~=0）
if sum(indexDel1)~=0
    % copy of DeleteEMUReactions -- AddEMUReactions (在循环中更新)
    AddEMUReactions = unique(DeleteEMUReactions);
    while true
        % 合并相同节点，首尾串联反应，存储新的反应 -- NewEMUReactions
        NewEMUReactions = table('Size',[0 3],'VariableTypes',{'string','string','double'},'VariableNames',{'Reaction','Flux','Size'});
        % 标记AddEMUReactions中可串联的节点 -- LabelNodes
        LabelDelNodes = zeros(height(AddEMUReactions),1);

        % 从AddEMUReactions中提取反应项 -- DeleteReactants, DeleteReactantsSplit
        DeleteReactants = regexp(AddEMUReactions.Reaction,'.+(?=@)','match');
        DeleteReactants = [DeleteReactants{:}]';
        DeleteReactantsSplit = regexp(regexprep(DeleteReactants,'\(|\)',''),', ','split');

        % 从AddEMUReactions中提取生成项
        DeleteProducts = regexp(AddEMUReactions.Reaction,'(?<=@).+','match');
        DeleteProducts = [DeleteProducts{:}]';

        % 记录可能是串联节点的位置索引 -- Locb
        NumDel = numel(DeleteProducts);
        Locb = zeros(NumDel,1);
        % 遍历DeleteProducts中所有生成物
        for iD = 1:NumDel
            % 判断第iD个生成物是否为其他反应中的反应物
            conditon = ismember(DeleteProducts{iD},[DeleteReactantsSplit{:}]);
            if conditon
                Locb(iD) = iD;
            end
        end

        % 从Locb中选取最小(~=0)的开始遍历
        index1 = min(Locb(Locb>0));

        %%%%%%%%%%%%%%%%%%%%%%%% Break condition %%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(index1)
            break
        end

        % 对应DeleteProducts{index1}的代谢物 -- metaDel
        metaDel = EMUsTotal.Mets{strcmp(EMUsTotal.EMUs,DeleteProducts{index1})};
        % 若该代谢物不属于DeleteNodes，则进入以下循环，直到取到DeleteNodes
        while ~ismember(metaDel,DeleteNodes.DelNode)
            Locb(index1) = 0;
            index1 = min(Locb(Locb>0));
            if ~isempty(index1)
                metaDel = EMUsTotal.Mets{strcmp(EMUsTotal.EMUs,DeleteProducts{index1})};
            else
                break
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%% Break condition %%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(index1)
            break
        end

        % 标记可串联节点（生成物） -- LabelDelNodes
        LabelDelNodes(index1) = 1;

        % 根据DeleteProducts(index1)在DeleteReactantsSplit找到相同节点的位置索引 -- index2
        index2 = find(cell2mat(cellfun(@(x) ismember(DeleteProducts(index1),x),DeleteReactantsSplit,'UniformOutput',false)));
        % 遍历index2中的所有位置索引
        for i2 = 1:numel(index2)
            % 标记可串联节点（反应物） -- LabelDelNodes
            LabelDelNodes(index2(i2)) = 1;
            if numel(DeleteReactantsSplit{index2(i2)})==1
                % 若索引到的反应是单分子反应，直接连接DeleteReactants{index1}和DeleteProducts{index2(i2)}
                Reaction = strcat(DeleteReactants{index1},'@',DeleteProducts{index2(i2)});
            else
                % 若索引到的反应是聚合反应，则替换DeleteReactants{index2(i2)}中的对应EMU，再进行连接
                DeleteReactants{index2(i2)} = replace(DeleteReactants{index2(i2)},DeleteProducts{index1},DeleteReactants{index1});
                Reaction = strcat(DeleteReactants{index2(i2)},'@',DeleteProducts{index2(i2)});
            end

            % 从SMatrix中提取出metaDel所在行 -- metarow
            metarow = SMatrix(strcmp(BalancedNodes,metaDel),:);
            % 从metarow中提取出首反应流及流系数 -- Flux1, CoeFlux1
            Flux1 = STable.Properties.VariableNames{metarow>0};
            CoeFlux1 = abs(metarow(metarow>0));
            % 从metarow中提取出尾反应流及流系数 -- Flux2, CoeFlux2
            Flux2 = STable.Properties.VariableNames{metarow<0};
            CoeFlux2 = abs(metarow(metarow<0));

            Flux1 = strcat(num2str(CoeFlux1/CoeFlux2),'*',Flux1);
            Flux = replace(AddEMUReactions.Flux{index2(i2)},Flux2,Flux1);
            Size = AddEMUReactions.Size(index2(i2));
            NewEMUReactions = [NewEMUReactions;{Reaction,Flux,Size}];
        end
        AddEMUReactions(logical(LabelDelNodes),:) = [];
        AddEMUReactions = unique([AddEMUReactions;NewEMUReactions]);
    end
    AddEMUReactions.Flux = cellfun(@(x) strcat(num2str(eval(string(regexp(x,'.+(?=\*v)','match')))),'*v',regexp(x,'(?<=\*v).+','match')),AddEMUReactions.Flux,'UniformOutput',false);
    EMUReactionTable(logical(indexDel1),:) = [];
    EMUReactionTable = vertcat(EMUReactionTable,AddEMUReactions);
end

%% 合并EMUReactionTable中同类项
CombinedEMUReactionTable1 = groupsummary(UnsimplifiedEMUReactions,["Reaction", "Size"],@(x) strjoin(x,' + '),"Flux");
UnsimplifiedEMUReactions = table(CombinedEMUReactionTable1.Reaction,CombinedEMUReactionTable1.fun1_Flux,CombinedEMUReactionTable1.Size,'VariableNames',{'Reaction','Flux','Size'});
CombinedEMUReactionTable2 = groupsummary(EMUReactionTable,["Reaction", "Size"],@(x) strjoin(x,' + '),"Flux");
SimplifiedEMUReactions = table(CombinedEMUReactionTable2.Reaction,CombinedEMUReactionTable2.fun1_Flux,CombinedEMUReactionTable2.Size,'VariableNames',{'Reaction','Flux','Size'});

end