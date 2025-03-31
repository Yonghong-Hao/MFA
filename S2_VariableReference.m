function [VariablesRefernce, BoundaryReference, emuReactant, emuProduct] = S2_VariableReference(EMUReactionTable,TracerTable,UnvisitedEMUList)
%% emu reactant
emuReactant = regexp(regexprep(cellstr(regexp(EMUReactionTable.Reaction,'.+(?=@)','match')),'\(|\)|',''),', ','split');
% select the emus in condensation reactions
Condensation = cell2mat(cellfun(@numel,emuReactant,'UniformOutput',false))>1;
emuCondensation = emuReactant(Condensation);
emuCondensationTotal = unique([emuCondensation{:}]);
emuCondensationSize = zeros(numel(emuCondensationTotal),1);
for i = 1:numel(emuCondensationTotal)
    emu = emuCondensationTotal{i};
    emuSplit = regexp(emu,'_','split');
    atom = emuSplit(end);
    emuCondensationSize(i) = length(atom{:});
end
emuCondensation_emu = cell(numel(emuCondensationTotal),1);
emuCondensation_mid = cell(numel(emuCondensationTotal),1);
emuCondensationSize_mid = zeros(0,1);
for i = 1:numel(emuCondensationTotal)
    emu = emuCondensationTotal(i);
    emuCondensation_emu(i) = {repmat(emu,1,emuCondensationSize(i)+1)};
    emuCondensation_mid(i) = {strcat(emu, '_m', cellfun(@num2str, num2cell(0:emuCondensationSize(i)),'UniformOutput',false))};
    emuCondensationSize_mid = [emuCondensationSize_mid;repmat(emuCondensationSize(i),emuCondensationSize(i)+1,1)];
end
emuCondensation_mid = table([emuCondensation_emu{:}]',[emuCondensation_mid{:}]',emuCondensationSize_mid,'VariableNames',{'emu_emu','emu_mid','emu_size'});

% treat the emus not in condensation reactions
NonCondensation = ~Condensation;
emuNonCondensation = emuReactant(NonCondensation);
emuNonCondensation_Size = table([emuNonCondensation{:}]',EMUReactionTable.Size(NonCondensation),'VariableNames',{'emuReactant','Size'});
emuNonCondensationTotal = unique(emuNonCondensation_Size);

emuNonCondensation_emu = cell(numel(emuNonCondensationTotal.emuReactant),1);
emuNonCondensation_mid = cell(numel(emuNonCondensationTotal.emuReactant),1);
emuNonCondensationSize_mid = zeros(0,1);
for i = 1:numel(emuNonCondensationTotal.emuReactant)
    emu = emuNonCondensationTotal.emuReactant(i);
    emuNonCondensation_emu(i) = {repmat(emu,1,emuNonCondensationTotal.Size(i)+1)};
    emuNonCondensation_mid(i) = {strcat(emu, '_m', cellfun(@num2str, num2cell(0:emuNonCondensationTotal.Size(i)),'UniformOutput',false))};
    emuNonCondensationSize_mid = [emuNonCondensationSize_mid;repmat(emuNonCondensationTotal.Size(i),emuNonCondensationTotal.Size(i)+1,1)];
end
emuNonCondensation_mid = table([emuNonCondensation_emu{:}]',[emuNonCondensation_mid{:}]',emuNonCondensationSize_mid,'VariableNames',{'emu_emu','emu_mid','emu_size'});
emuReactant_mid = union(emuCondensation_mid,emuNonCondensation_mid,'stable');

%% emu product
emuProduct = regexp(EMUReactionTable.Reaction,'(?<=@).+','match');
emuProduct_Size = table([emuProduct{:}]',EMUReactionTable.Size,'VariableNames',{'emuProduct','Size'});
ProductTotal = unique(emuProduct_Size);

emuProduct_emu = cell(numel(ProductTotal.emuProduct),1);
emuProduct_mid = cell(numel(ProductTotal.emuProduct),1);
emuProductSize_mid = zeros(0,1);
for i = 1:numel(ProductTotal.emuProduct)
    emu = ProductTotal.emuProduct(i);
    emuProduct_emu(i) = {repmat(emu,1,ProductTotal.Size(i)+1)};
    emuProduct_mid(i) = {strcat(emu, '_m', cellfun(@num2str, num2cell(0:ProductTotal.Size(i)),'UniformOutput',false))};
    emuProductSize_mid = [emuProductSize_mid;repmat(ProductTotal.Size(i),ProductTotal.Size(i)+1,1)];
end
emuProduct_mid = table([emuProduct_emu{:}]',[emuProduct_mid{:}]',emuProductSize_mid,'VariableNames',{'emu_emu','emu_mid','emu_size'});
emu_midTotal = union(emuReactant_mid,emuProduct_mid,'stable');
% add emu_atom to emu_midTotal
emuTotal = regexprep(emu_midTotal.emu_mid,'_m\d+','');
meta_atom = regexp(emuTotal,'_','split');
emu_atom = cellfun(@(x) x(end),meta_atom,'UniformOutput',false);
% add emu_meta to emu_midTotal
emu_meta = cellfun(@(x,y) regexp(x,strcat('\w+(?=_',y,'_)'),"match"),emu_midTotal.emu_mid,emu_atom,'UniformOutput',false);
emu_midTotal = addvars(emu_midTotal,emu_meta,emu_atom,Before='emu_size');

%% Total variables (from EMUReactionTable)
VariablesTotal = UnvisitedEMUList.EMUs;
%Variables = VariablesTotal(UnvisitedEMUList.Size==-1);
BoundaryVariables = VariablesTotal(UnvisitedEMUList.Size==-2);

%% boundary variables
emuBoundary_emu = {};
emuBoundary_mid = {};
emuBoundary_atom = {};
emuBoundary_size = zeros(0,1);
for bv = 1:numel(BoundaryVariables)
    indexBoundary = ismember(emu_midTotal.emu_emu,BoundaryVariables(bv));
    if sum(indexBoundary)~=0
        emuBoundary_emu = [emuBoundary_emu;emu_midTotal.emu_emu(indexBoundary)];
        emuBoundary_mid = [emuBoundary_mid;emu_midTotal.emu_mid(indexBoundary)];
        emuBoundary_atom = [emuBoundary_atom;emu_midTotal.emu_atom(indexBoundary)];
        emuBoundary_size = [emuBoundary_size;emu_midTotal.emu_size(indexBoundary)];
    end
end

emuBoundary_var = cell(numel(emuBoundary_mid),1);
for bv = 1:numel(emuBoundary_mid)
    emu_bv = emuBoundary_mid{bv};
    emu_bv_split = regexp(emu_bv,'\_','split');
    if strcmp(emu_bv_split{end},'m0')
        emuBoundary_var{bv} = '1';
    else
        emuBoundary_var{bv} = '0';
    end
end

for t = 1:height(TracerTable)
    emuBoundary_var(strcmp(emuBoundary_mid,TracerTable.emu_tracer_mid{t})) = {num2str(TracerTable.emu_tracer_var(t))};
end

BoundaryReference = table(emuBoundary_emu,emuBoundary_atom,emuBoundary_mid,emuBoundary_var,emuBoundary_size,'VariableNames',{'emuBoundary_emu','emuBoundary_atom','emuBoundary_mid','emuBoundary_var','emuBoundary_size'});

%% combine
[emu_mid, ia] = setdiff(emu_midTotal.emu_mid,emuBoundary_mid,'stable');
emu_emu = emu_midTotal.emu_emu(ia);
emu_atom = emu_midTotal.emu_atom(ia);
emu_var = strcat('X_', cellfun(@num2str, num2cell(1:numel(emu_mid)),'UniformOutput',false))';
emu_size = emu_midTotal.emu_size(ia);
emu_val = cell(numel(emu_mid),1);
VariablesRefernce = table(emu_emu,emu_atom,emu_mid,emu_var,emu_val,emu_size);
end