%% Input
clear;clc
MappingID = 'mapping27';
TracerID = 'glc_12';
%% Align all reaction IDs
RawMapping = readtable("Tracer_glc_12/Data_raw/MappingAtom_raw_1.xlsx","VariableNamingRule","preserve");
AlignedMapping = readtable(strcat('Tracer_',TracerID,'/Data_raw/',regexprep(MappingID,'mapping','MappingAtom_raw_'),'.xlsx'),"VariableNamingRule","preserve");
for i = 1:height(AlignedMapping)
    f_index = ismember(RawMapping.("SubstrateIDs(atoms)"),AlignedMapping.("SubstrateIDs(atoms)"){i}) & ismember(RawMapping.("ProductIDs(atoms)"),AlignedMapping.("ProductIDs(atoms)"){i});
    b_index = ismember(RawMapping.("SubstrateIDs(atoms)"),AlignedMapping.("ProductIDs(atoms)"){i}) & ismember(RawMapping.("ProductIDs(atoms)"),AlignedMapping.("SubstrateIDs(atoms)"){i});
    index = f_index | b_index;
    if sum(index) == 1
        AlignedMapping.ReactionIDs{i} = RawMapping.ReactionIDs{index};
    end
end
writetable(AlignedMapping,strcat('Tracer_',TracerID,'/Data_raw/',regexprep(MappingID,'mapping','MappingAtom_raw_'),'.xlsx'));

%% Align reactions to simulate
ToSimulateFlux1 = readtable("Tracer_glc_12/Data_raw/ToSimulateFlux_1_1_mapping1.xlsx","VariableNamingRule","preserve");
for i = 1:height(ToSimulateFlux1)
    reaction_index = ismember(AlignedMapping.ReactionIDs,ToSimulateFlux1.ReactionID{i});
    if sum(reaction_index) == 1
        ToSimulateFlux1.("SubstrateIDs(atom)"){i} = AlignedMapping.("SubstrateIDs(atoms)"){reaction_index};
        ToSimulateFlux1.("ProductIDs(atom)"){i} = AlignedMapping.("ProductIDs(atoms)"){reaction_index};
        ToSimulateFlux1.Reversibility(i) = AlignedMapping.Reversibility(reaction_index);
    end
end
ToSimulateFlux2 = readtable("Tracer_glc_12/Data_raw/ToSimulateFlux_1_2_mapping1.xlsx","VariableNamingRule","preserve");
for i = 1:height(ToSimulateFlux2)
    reaction_index = ismember(AlignedMapping.ReactionIDs,ToSimulateFlux2.ReactionID{i});
    if sum(reaction_index) == 1
        ToSimulateFlux2.("SubstrateIDs(atom)"){i} = AlignedMapping.("SubstrateIDs(atoms)"){reaction_index};
        ToSimulateFlux2.("ProductIDs(atom)"){i} = AlignedMapping.("ProductIDs(atoms)"){reaction_index};
        ToSimulateFlux2.Reversibility(i) = AlignedMapping.Reversibility(reaction_index);
    end
end
ToSimulateFlux3 = readtable("Tracer_glc_12/Data_raw/ToSimulateFlux_1_3_mapping1.xlsx","VariableNamingRule","preserve");
for i = 1:height(ToSimulateFlux3)
    reaction_index = ismember(AlignedMapping.ReactionIDs,ToSimulateFlux3.ReactionID{i});
    if sum(reaction_index) == 1
        ToSimulateFlux3.("SubstrateIDs(atom)"){i} = AlignedMapping.("SubstrateIDs(atoms)"){reaction_index};
        ToSimulateFlux3.("ProductIDs(atom)"){i} = AlignedMapping.("ProductIDs(atoms)"){reaction_index};
        ToSimulateFlux3.Reversibility(i) = AlignedMapping.Reversibility(reaction_index);
    end
end
writetable(ToSimulateFlux1,strcat('Tracer_',TracerID,'/Data_raw/ToSimulateFlux_1_1_',MappingID,'.xlsx'))
writetable(ToSimulateFlux2,strcat('Tracer_',TracerID,'/Data_raw/ToSimulateFlux_1_2_',MappingID,'.xlsx'))
writetable(ToSimulateFlux3,strcat('Tracer_',TracerID,'/Data_raw/ToSimulateFlux_1_3_',MappingID,'.xlsx'))

