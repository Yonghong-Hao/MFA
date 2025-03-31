clear;clc
load('Data Source/biot201700518-sup-0001-suppagm_cm-s1');
%% Mapping and Flux estimation
Reactions = m.rates;
ReactionIDs = cell(width(Reactions),1);
SubstrateIDs = cell(width(Reactions),1);
ProductIDs = SubstrateIDs;
Reversibility = SubstrateIDs;
FluxValue_f = SubstrateIDs;
FluxValue_b = SubstrateIDs;
for i = 1:width(Reactions)
    reaction_i = Reactions(i).flx;
    if numel(reaction_i) == 1
        ReactionIDs{i} = Reactions(i).id;
        Reversibility{i} = 0;
        SubstrateIDs{i} = reaction_i.sub.id;
        ProductIDs{i} = reaction_i.prod.id;
        FluxValue_f{i} = reaction_i.val;
        FluxValue_b{i} = 0;
    else
        Reversibility{i} = 1;
        ReactionIDs{i} = Reactions(i).id;
        SubstrateIDs{i} = reaction_i(1).sub.id;
        ProductIDs{i} = reaction_i(1).prod.id;
        FluxValue_f{i} = reaction_i(1).val;
        FluxValue_b{i} = reaction_i(2).val;
    end
end
ReactionTable = table(ReactionIDs,SubstrateIDs,ProductIDs,Reversibility,FluxValue_f,FluxValue_b);

%% MID tracer glc labelled pyr_123,lac_123,lac_23
clear;clc
load('Data Source/biot201700518-sup-0001-suppagm_cm-s1');
MeasuredData = f.mnt;
ToSimulateIDs = cell(3,1);
MeasuredVal = ToSimulateIDs;
SimulatedVal_Raw = ToSimulateIDs;
for i = 1:width(MeasuredData)
    if strcmp(MeasuredData(i).id,'Pyr 174') || strcmp(MeasuredData(i).id,'Lac 233') || strcmp(MeasuredData(i).id,'Lac 261')
        ToSimulateIDs{i} = MeasuredData(i).id;
        MeasuredVal{i} = mat2str([MeasuredData(i).res.data]);
        SimulatedVal_Raw{i} = mat2str([MeasuredData(i).res.fit]);
    elseif strcmp(MeasuredData(i).type,'Flux')
        ToSimulateIDs{i} = MeasuredData(i).id;
        MeasuredVal{i} = MeasuredData(i).res.data;
        SimulatedVal_Raw{i} = MeasuredData(i).res.fit;
    end
end
MeasuredTable = table(ToSimulateIDs,MeasuredVal,SimulatedVal_Raw);
writetable(MeasuredTable,"Tracer_glc_12/Data_raw/MeasredData_1.xlsx");












