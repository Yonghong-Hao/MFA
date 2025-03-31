%==========================================================================
%===================== Confidence interval estimate =======================
%==========================================================================
clear;clc
%% input
TracerID = 'glc_12';
MappingID = 'mapping1';
DataSet = {'data1_1','data1_2','data1_3'};
Conditions = {'CM','LA+','LA'};
%% Merge reversible reactions
% Mapping without reversibility
Mapping_Irreversible = readtable(strcat('Tracer_',TracerID,"/MappingAtom/",regexprep(MappingID,'mapping','mapping_atom_'),".xlsx"),"VariableNamingRule","preserve");
SubstrateTotal = Mapping_Irreversible.("SubstrateIDs(atoms)");
ProductTotal = Mapping_Irreversible.("ProductIDs(atoms)");

% Mapping with reversibility
Mapping_Reversible = readtable(strcat('Tracer_',TracerID,"/Data_raw/",regexprep(MappingID,'mapping','MappingAtom_raw_'),".xlsx"),"VariableNamingRule","preserve");
NumReactions = height(Mapping_Reversible);
ReactionIDs = Mapping_Reversible.ReactionIDs;
Reaction_f_index = zeros(numel(ReactionIDs),1);
Reaction_b_index = Reaction_f_index;
Reversibility = Mapping_Reversible.Reversibility;
for i = 1:NumReactions
    Reversibility_i = Reversibility(i);
    Substrate_i = Mapping_Reversible.("SubstrateIDs(atoms)"){i};
    Product_i = Mapping_Reversible.("ProductIDs(atoms)"){i};
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
%% read data
SearchedData_CM = load(strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/ConfidenceInterval_',MappingID,'_',DataSet{1},'.mat'));
SearchedData_LA_high = load(strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/ConfidenceInterval_',MappingID,'_',DataSet{2},'.mat'));
SearchedData_LA = load(strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/ConfidenceInterval_',MappingID,'_',DataSet{3},'.mat'));
ConfidenceInterval_CM = table('Size',[height(Mapping_Irreversible),3],'VariableTypes',{'double','double','double'},'VariableNames',{'Vhat_CM','LB_CM','UB_CM'});
for i = 1:height(Mapping_Irreversible)
    ConfidenceInterval_CM.LB_CM(i) = SearchedData_CM.SearchedData(i).BackwardSearch.NewV(end);
    ConfidenceInterval_CM.UB_CM(i) = SearchedData_CM.SearchedData(i).ForwardSearch.NewV(end);
    ConfidenceInterval_CM.Vhat_CM(i) = SearchedData_CM.SearchedData(i).OptimalSolution.Vhat;
end
ConfidenceInterval_LA_high = table('Size',[height(Mapping_Irreversible),3],'VariableTypes',{'double','double','double'},'VariableNames',{'Vhat_LA+','LB_LA+','UB_LA+'});
for i = 1:height(Mapping_Irreversible)
    ConfidenceInterval_LA_high.("LB_LA+")(i) = SearchedData_CM.SearchedData(i).BackwardSearch.NewV(end);
    ConfidenceInterval_LA_high.("UB_LA+")(i) = SearchedData_CM.SearchedData(i).ForwardSearch.NewV(end);
    ConfidenceInterval_LA_high.("Vhat_LA+")(i) = SearchedData_CM.SearchedData(i).OptimalSolution.Vhat;
end
ConfidenceInterval_LA = table('Size',[height(Mapping_Irreversible),3],'VariableTypes',{'double','double','double'},'VariableNames',{'Vhat_LA','LB_LA','UB_LA'});
for i = 1:height(Mapping_Irreversible)
    ConfidenceInterval_LA.LB_LA(i) = SearchedData_CM.SearchedData(i).BackwardSearch.NewV(end);
    ConfidenceInterval_LA.UB_LA(i) = SearchedData_CM.SearchedData(i).ForwardSearch.NewV(end);
    ConfidenceInterval_LA.Vhat_LA(i) = SearchedData_CM.SearchedData(i).OptimalSolution.Vhat;
end
ConfidenceInterval_Irreversible = [ConfidenceInterval_CM,ConfidenceInterval_LA_high,ConfidenceInterval_LA];
ConfidenceInterval_Reversible = zeros(height(ReactionsTotal),3*numel(Conditions));
for i = 1:height(Mapping_Reversible)
    if ReactionsTotal.Reversibility(i)==0
        ConfidenceInterval_Reversible(i,:) = ConfidenceInterval_Irreversible{ReactionsTotal.Reaction_f_index(i),:};
    else
        %ConfidenceInterval_Reversible(i,:) = ConfidenceInterval_Irreversible{ReactionsTotal.Reaction_f_index(i),:}-ConfidenceInterval_Irreversible{ReactionsTotal.Reaction_b_index(i),:};
        ConfidenceInterval_Reversible(i,1:3) = [ConfidenceInterval_CM{ReactionsTotal.Reaction_f_index(i),1}-ConfidenceInterval_CM{ReactionsTotal.Reaction_b_index(i),1},ConfidenceInterval_CM{ReactionsTotal.Reaction_f_index(i),2:3}-flip(ConfidenceInterval_CM{ReactionsTotal.Reaction_b_index(i),2:3})];
        ConfidenceInterval_Reversible(i,4:6) = [ConfidenceInterval_CM{ReactionsTotal.Reaction_f_index(i),1}-ConfidenceInterval_LA_high{ReactionsTotal.Reaction_b_index(i),1},ConfidenceInterval_LA_high{ReactionsTotal.Reaction_f_index(i),2:3}-flip(ConfidenceInterval_LA_high{ReactionsTotal.Reaction_b_index(i),2:3})];
        ConfidenceInterval_Reversible(i,7:9) = [ConfidenceInterval_CM{ReactionsTotal.Reaction_f_index(i),1}-ConfidenceInterval_LA{ReactionsTotal.Reaction_b_index(i),1},ConfidenceInterval_LA{ReactionsTotal.Reaction_f_index(i),2:3}-flip(ConfidenceInterval_LA{ReactionsTotal.Reaction_b_index(i),2:3})];
    end
end
ConfidenceInterval = array2table(ConfidenceInterval_Reversible,"VariableNames",{'Vhat_CM','LB_CM','UB_CM','Vhat_LA+','LB_LA+','UB_LA+','Vhat_LA','LB_LA','UB_LA'});
Mapping_Confidence = [Mapping_Reversible,ConfidenceInterval];
% delete some fluxes that have no practical significance
switch MappingID
    case 'mapping1'
        DeleteFluxSet = [];
    case 'mapping8'
        DeleteFluxSet = {'v11','v12','v13','v60','v61','v62','v63','v64','v65','v66'};
end
if ~isempty(DeleteFluxSet)
    ConfidenceInterval_new = ConfidenceInterval(~ismember(Mapping_Confidence.FluxIDs,DeleteFluxSet),:);
    Mapping_Confidence_new = removevars(Mapping_Confidence(~ismember(Mapping_Confidence.FluxIDs,DeleteFluxSet),2:end),'Note');
else
    ConfidenceInterval_new = ConfidenceInterval;
    Mapping_Confidence_new = Mapping_Confidence;
end
writetable(Mapping_Confidence_new,strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/ConfidenceInterval_total_',MappingID,'.xlsx'))
%% plot
f = figure;
set(groot,'DefaultAxesFontSize',10);
f.Position = [1921,49,1920,954];
for i = 1:numel(Conditions)
    subplot(numel(Conditions),1,i)
    factor = ones(height(Mapping_Confidence_new),1);
    for j = 1:height(Mapping_Confidence_new)
        if strcmp(Mapping_Confidence_new.ReactionIDs{j},'Biomass')
            factor(j) = 1;
        elseif strcmp(Mapping_Confidence_new.ReactionIDs{j},'mAb')
            factor(j) = 1;
        end
    end
    b = bar(ConfidenceInterval_new{:,3*(i-1)+1}.*factor,'EdgeColor','none','FaceColor',[0.6157    0.7804    0.8941]);
    hold on
    CI = [ConfidenceInterval_new{:,3*(i-1)+2}.*factor,ConfidenceInterval_new{:,3*i}.*factor];
    x = 1:1:height(CI);
    hold on
    scatter(x,CI(:,1),35,[0    0.6000    0.8196],"_")
    hold on
    scatter(x,CI(:,2),35,[0    0.6000    0.8196],"_")
    hold on
    plot([x;x],CI','Color',[0    0.6000    0.8196])
    xticks(1:1:height(Mapping_Confidence_new));
    xticklabels(Mapping_Confidence_new.ReactionIDs);
    ylabel("Net Flux");
    %ylim([-1,3])
    title(strcat(Conditions{i}),'Units', 'Normalized','Position',[0.5,1],'FontSize',15,'FontWeight','bold')
end
saveas(f,strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/ConfidenceInterval_bar_',MappingID),'png')