clear;clc
%% input
TracerID = 'glc_12';
MappingID = 'mapping9';
DataSet = {'data1_1','data1_2','data1_3'};
Conditions = {'CM','LA+','LA'};
%% Choose the optimal solution
% read data
Result_CM = struct('FileName',{},'FluxValue',{},'LossFunctionValue',{});
Result_LA_high = struct('FileName',{},'FluxValue',{},'LossFunctionValue',{});
Result_LA = struct('FileName',{},'FluxValue',{},'LossFunctionValue',{});
for i = 1:numel(Conditions)
    FolderPath = strcat('Tracer_',TracerID,"/FluxValues_LossFunctionValues/",MappingID,'/',Conditions{i});
    FilePattern_V_i = fullfile(FolderPath, '*FluxValues*.xlsx');
    Files_V_i = dir(FilePattern_V_i);
    FilePattern_LF = fullfile(FolderPath, '*LossFunctionValues*.xlsx');
    Files_LF_i = dir(FilePattern_LF);
    switch Conditions{i}
        case Conditions{1}
            for i_V = 1:numel(Files_V_i)
                Result_CM(i_V).FileName = Files_V_i(i_V).name;
                Result_CM(i_V).FluxValue = readtable(fullfile(FolderPath, Files_V_i(i_V).name),"VariableNamingRule","preserve","ReadRowNames",true);
                Result_CM(i_V).LossFunctionValue = readtable(fullfile(FolderPath, Files_LF_i(i_V).name),"VariableNamingRule","preserve");
            end
        case Conditions{2}
            for i_V = 1:numel(Files_V_i)
                Result_LA_high(i_V).FileName = Files_V_i(i_V).name;
                Result_LA_high(i_V).FluxValue = readtable(fullfile(FolderPath, Files_V_i(i_V).name),"VariableNamingRule","preserve","ReadRowNames",true);
                Result_LA_high(i_V).LossFunctionValue = readtable(fullfile(FolderPath, Files_LF_i(i_V).name),"VariableNamingRule","preserve");
            end
        case Conditions{3}
            for i_V = 1:numel(Files_V_i)
                Result_LA(i_V).FileName = Files_V_i(i_V).name;
                Result_LA(i_V).FluxValue = readtable(fullfile(FolderPath, Files_V_i(i_V).name),"VariableNamingRule","preserve","ReadRowNames",true);
                Result_LA(i_V).LossFunctionValue = readtable(fullfile(FolderPath, Files_LF_i(i_V).name),"VariableNamingRule","preserve");
            end
    end
end
% CM
MinLFSet = zeros(width(Result_CM),2);
for i = 1:width(Result_CM)
    [MinLFSet(i,2),MinLFSet(i,1)] = min(Result_CM(i).LossFunctionValue.Variables);
end
[~,MinLF_iteration] = min(MinLFSet(:,2));
Vhat_1 = Result_CM(MinLF_iteration).FluxValue{:,MinLFSet(MinLF_iteration)};
% LA_high
MinLFSet = zeros(width(Result_LA_high),2);
for i = 1:width(Result_LA_high)
    [MinLFSet(i,2),MinLFSet(i,1)] = min(Result_LA_high(i).LossFunctionValue.Variables);
end
[~,MinLF_iteration] = min(MinLFSet(:,2));
Vhat_2 = Result_LA_high(MinLF_iteration).FluxValue{:,MinLFSet(MinLF_iteration)};
% LA
MinLFSet = zeros(width(Result_LA),2);
for i = 1:width(Result_LA)
    [MinLFSet(i,2),MinLFSet(i,1)] = min(Result_LA(i).LossFunctionValue.Variables);
end
[~,MinLF_iteration] = min(MinLFSet(:,2));
Vhat_3 = Result_LA(MinLF_iteration).FluxValue{:,MinLFSet(MinLF_iteration)};
% save
MappingAtomTable = readtable(strcat('Tracer_',TracerID,'/MappingAtom/',regexprep(MappingID,'mapping','mapping_atom_'),'.xlsx'),"VariableNamingRule","preserve");
MappingAtomTable_Vhat = addvars(MappingAtomTable,Vhat_1,Vhat_2,Vhat_3,'NewVariableNames',Conditions);
writetable(MappingAtomTable_Vhat,strcat('Tracer_',TracerID,'/VerifyMyModel/ConfidenceInterval/Vhat_',MappingID,'.xlsx'));
