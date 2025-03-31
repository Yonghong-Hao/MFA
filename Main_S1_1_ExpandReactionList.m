clear; clc;
%%将带有可逆性的反应列表转换成不带可逆性的反应列表
%% input
TracerID = 'Tracer_glc_12';
MappingID = 'mapping_27';
%%
MappingAtomFile_raw = strcat(TracerID,'/Data_raw/',regexprep(MappingID,'mapping','MappingAtom_raw'),'.xlsx');
MappingAtomTable_raw = readtable(MappingAtomFile_raw,"VariableNamingRule","preserve");
NumReactions_raw = height(MappingAtomTable_raw);
MappingAtomTable_reversibility = MappingAtomTable_raw;
MappingAtomTable_reversibility.("SubstrateIDs(atoms)") = MappingAtomTable_raw.("ProductIDs(atoms)");
MappingAtomTable_reversibility.("ProductIDs(atoms)") = MappingAtomTable_raw.("SubstrateIDs(atoms)");
MappingAtomTable_reversibility(MappingAtomTable_reversibility.Reversibility==0,:)=[];
NumReactions = height(MappingAtomTable_reversibility)+NumReactions_raw;
MappingAtomTable_reversibility.FluxIDs = strcat('v',cellfun(@num2str, num2cell(1+NumReactions_raw:NumReactions),'UniformOutput',false))';
MappingAtomTable = vertcat(MappingAtomTable_raw,MappingAtomTable_reversibility);
MappingAtomTable = removevars(MappingAtomTable,"Reversibility");

% 保存到'FreFluxinput'
writetable(MappingAtomTable,strcat(TracerID,'/MappingAtom/',regexprep(MappingID,'mapping','mapping_atom'),'.xlsx'))

