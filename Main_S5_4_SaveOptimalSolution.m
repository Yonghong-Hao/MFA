clear;clc
%% input
TracerID = 'glc_12';
MappingID = 'mapping27';
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
Vhat_Irreversible = readtable(strcat('Tracer_',TracerID,"/VerifyMyModel/ConfidenceInterval/Vhat_",MappingID,".xlsx"),"VariableNamingRule","preserve");
Vhat_Reversible = zeros(height(ReactionsTotal),numel(Conditions));
for i = 1:height(Mapping_Reversible)
    if ReactionsTotal.Reversibility(i)==0
        Vhat_Reversible(i,:) = Vhat_Irreversible{ReactionsTotal.Reaction_f_index(i),6:8};
    else
        Vhat_Reversible(i,:) = Vhat_Irreversible{ReactionsTotal.Reaction_f_index(i),6:8}-Vhat_Irreversible{ReactionsTotal.Reaction_b_index(i),6:8};
    end
end
VhatTable = array2table(Vhat_Reversible,"VariableNames",{'Vhat_CM','Vhat_LA+','Vhat_LA'});
Mapping_Vhat = [Mapping_Reversible,VhatTable];
% delete some fluxes that have no practical significance
switch MappingID
    case 'mapping1'
        DeleteFluxSet = [];
    case 'mapping8'
        DeleteFluxSet = {'Cit.c','Cit.m','Cit.mix','Pyr.c','Pyr.m','Pyr.m2','Pyr.mix','Akg.c','Akg.m','Akg.mix'};
    case 'mapping9'
        DeleteFluxSet = {'Cit.c','Cit.m','Cit.mix','Mal.c','Mal.m','Mal.mix','Fum.c','Fum.m','Fum.mix','Suc Dilution.m','Suc Dilution.u','Suc Sink','Fum Dilution.e','Fum Dilution.u','Fum Sink.d','Pyr.c','Pyr.m','Pyr.mix'};
    case 'mapping20'
        DeleteFluxSet = {'Cit.c','Cit.m','Cit.mix','Mal.c','Mal.m','Mal.mix','Fum.c','Fum.m','Fum.mix'};
    case 'mapping25'
        DeleteFluxSet = {'Cit.c','Cit.m','Cit.mix','Pyr.c','Pyr.m','Pyr.mix','Akg.c','Akg.m','Akg.mix','Glu.c','Glu.m','Glu.mix','Mal.c','Mal.m','Mal.mix','Fum.c','Fum.m','Fum.mix'};
    case 'mapping26'
        DeleteFluxSet = {'Cit.c','Cit.m','Cit.mix','Mal.c','Mal.m','Mal.mix','Pyr.c','Pyr.m','Pyr.mix'};
    case 'mapping27'
        DeleteFluxSet = [];
end
if ~isempty(DeleteFluxSet)
    Mapping_Vhat_new = Mapping_Vhat(~ismember(Mapping_Vhat.ReactionIDs,DeleteFluxSet),:);
else
    Mapping_Vhat_new = Mapping_Vhat;
end
writetable(removevars(Mapping_Vhat_new,"Note"),strcat('Tracer_',TracerID,'/VhatTotal/Vhat_total_',MappingID,'.xlsx'))
