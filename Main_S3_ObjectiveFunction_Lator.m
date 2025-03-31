%% Processing of data simulating extracellular fluxes
% the raw extralcellular flux data
ToSimulateFlux_raw = readtable(strcat(TracerID,'/Data_raw/ToSimulateFlux_',regexprep(DataID,'data',''),'_',MappingID,".xlsx"),"VariableNamingRule","preserve");
% the raw MappingAtom data
fluxTotal = MappingAtomTable.FluxIDs;
reactantTotal = MappingAtomTable.("SubstrateIDs(atoms)");
productTotal = MappingAtomTable.("ProductIDs(atoms)");
% 取交集
NumToSimulateFlux_raw = height(ToSimulateFlux_raw);
ToSimulateFluxIndex = zeros(NumToSimulateFlux_raw,1);
ToSimulateFlux_raw = addvars(ToSimulateFlux_raw,zeros(NumToSimulateFlux_raw,1),zeros(NumToSimulateFlux_raw,1),'After','Simulated','NewVariableNames',{'FluxIndex_f','FluxIndex_b'});
for i = 1:NumToSimulateFlux_raw
    reactant_i = ToSimulateFlux_raw.("SubstrateIDs(atom)"){i};
    product_i = ToSimulateFlux_raw.("ProductIDs(atom)"){i};
    reversibility_i = ToSimulateFlux_raw.Reversibility(i);
    fluxVal_i = ToSimulateFlux_raw.Measured(i);
    switch reversibility_i
        case 0
            fluxIndex = ismember(reactantTotal,reactant_i) & ismember(productTotal,product_i);
            if sum(fluxIndex) == 1
                ToSimulateFluxIndex(i) = 1;
                ToSimulateFlux_raw.FluxIndex_f(i) = str2double(regexprep(fluxTotal{fluxIndex},'v',''));
            end
        case 1
            fluxIndex_f = ismember(reactantTotal,reactant_i) & ismember(productTotal,product_i);
            fluxIndex_b = ismember(reactantTotal,product_i) & ismember(productTotal,reactant_i);
            if sum(fluxIndex_f) == 1 && sum(fluxIndex_b) == 1
                ToSimulateFluxIndex(i) = 1;
                ToSimulateFlux_raw.FluxIndex_f(i) = str2double(regexprep(fluxTotal{fluxIndex_f},'v',''));
                ToSimulateFlux_raw.FluxIndex_b(i) = str2double(regexprep(fluxTotal{fluxIndex_b},'v',''));
            end
    end
end
ToSimulateFlux_raw_overlap = ToSimulateFlux_raw(logical(ToSimulateFluxIndex),:);
GlcIN_index = ismember(ToSimulateFlux_raw_overlap.("SubstrateIDs(atom)"),'glc_tracer(abcdef)');
GlcIN = ToSimulateFlux_raw_overlap.Measured(GlcIN_index);
% normalization to GlcIN
ToSimulateFlux_raw_overlap.Measured = ToSimulateFlux_raw_overlap.Measured/GlcIN;
if sum(~isnan(ToSimulateFlux_raw_overlap.Simulated))==height(ToSimulateFlux_raw_overlap)
    ToSimulateFlux_raw_overlap.Simulated = ToSimulateFlux_raw_overlap.Simulated/GlcIN;
end
if sum(~isnan(ToSimulateFlux_raw_overlap.Std))==height(ToSimulateFlux_raw_overlap)
    ToSimulateFlux_raw_overlap.Std = ToSimulateFlux_raw_overlap.Std/GlcIN;
end
ToSimulateFlux = ToSimulateFlux_raw_overlap;

%% Processing of data simulating MID
% 扩展成mid of emus--include weight
ToSimulateMS_MID = table('Size',[0 5],'VariableTypes',{'string','string','double','double','double'},'VariableNames',{'EMUIDs_mid','Measured','Simulated','Std','Size'});
%NumDataSamples = width(ToSimulateMS)-2;
for i = 1:height(ToSimulateMS)
    ToSimulate = strcat(ToSimulateMS.EMUIDs{i},'_m',cellfun(@num2str, num2cell(0:ToSimulateMS.Size(i)),'UniformOutput',false));
    MIDDataMean = eval(ToSimulateMS.Measured{i});
    % test MID
    disp(sum(MIDDataMean));
    emuSimulatedMID = eval(ToSimulateMS.Simulated{i});
    Size = repmat(ToSimulateMS.Size(i),ToSimulateMS.Size(i)+1,1);
    if isnan(ToSimulateMS.Std(i))
        weightMID = repmat(ToSimulateMS.Std(i),ToSimulateMS.Size(i)+1,1);
    else
        weightMID = eval(ToSimulateMS.Std{i});
    end
    ToSimulateSet_new = table(ToSimulate',MIDDataMean',emuSimulatedMID',weightMID,Size,'VariableNames',{'EMUIDs_mid','Measured','Simulated','Std','Size'});
    ToSimulateMS_MID = vertcat(ToSimulateMS_MID,ToSimulateSet_new);
end

% 标准误差逆--SigmaInv
NumToSimulateMS_MID = height(ToSimulateMS_MID);
NumToSimulateFlux = height(ToSimulateFlux);
if sum(isnan(ToSimulateMS_MID.Std))==NumToSimulateMS_MID && sum(isnan(ToSimulateFlux.Std))==NumToSimulateFlux
    Sigma_MS = ones(NumToSimulateMS_MID,1)*0.02;
    Sigma_Flux = (ToSimulateFlux.Measured*0.050844+0.0027995);
    Sigma = [Sigma_MS;Sigma_Flux];
    Sigma(Sigma==0) = 1e-8;
else
    Sigma_MS = ToSimulateMS_MID.Std;
    Sigma_Flux = ToSimulateFlux.Std;
    Sigma = [Sigma_MS;Sigma_Flux];
end
SigmaInv = diag(1./Sigma);
AllVariables.ToSimulateMS_MID = ToSimulateMS_MID;
AllVariables.ToSimulateFlux = ToSimulateFlux;
AllVariables.SigmaInv = SigmaInv;
save(strcat(TracerID,'/FluxBalance&EMUBalance/IntermediateVariable_',MappingID,'_',DataID,'.mat'),'AllVariables')
%% write the the objective function file
FileContent = [
    "function [f, g] = ObjectiveFunction(U,ToSimulateMS_MID,ToSimulateFlux,VariablesRefernce,SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar,SigmaInv)"
    "% 得到V"
    "V = exp(U');"
    "%tic"
    "% 计算代谢物MID"
    "[X, dXdV] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,V,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar);"
    "%time = toc"
    ""
    "NumToSimulateMS_MID = height(ToSimulateMS_MID);"
    "% 选取优化变量的模拟值和优化变量的梯度值"
    "emuSimulatedMID = zeros(NumToSimulateMS_MID,1);"
    "emuSimulatedGradient = zeros(NumToSimulateMS_MID,NumFluxes);"
    "for j = 1:NumToSimulateMS_MID"
    "    condition = strcmp(VariablesRefernce.emu_mid,ToSimulateMS_MID.EMUIDs_mid{j});"
    "    emuSimulatedMID(j) = X(condition);"
    "    emuSimulatedGradient(j,:) = dXdV(condition,:).*exp(U);"
    "end"
    ""
    "% 计算外部流"
    "NumToSimulateFlux = height(ToSimulateFlux);"
    "SimulatedFlux = zeros(NumToSimulateFlux,1);"
    "FluxSimulatedGradient = zeros(NumToSimulateFlux,NumFluxes);"
    "for i = 1:NumToSimulateFlux"
    "    reversibility_i = ToSimulateFlux.Reversibility(i);"
    "    FluxIndex_f_i = ToSimulateFlux.FluxIndex_f(i);"
    "    FluxIndex_b_i = ToSimulateFlux.FluxIndex_b(i);"
    "    switch reversibility_i"
    "        case 0"
    "            SimulatedFlux(i) = V(FluxIndex_f_i);"
    "            FluxSimulatedGradient(i,FluxIndex_f_i) = 1;"
    "        case 1"
    "            SimulatedFlux(i) = V(FluxIndex_f_i)-V(FluxIndex_b_i);"
    "            FluxSimulatedGradient(i,FluxIndex_f_i) = 1;"
    "            FluxSimulatedGradient(i,FluxIndex_b_i) = -1;"
    "    end"
    "    FluxSimulatedGradient(i,:) = FluxSimulatedGradient(i,:).*exp(U);"
    "end"
    ""
    "% objective function"
    "SimulatedData = [emuSimulatedMID;SimulatedFlux];"
    "ToSimulateData = [str2double(ToSimulateMS_MID.Measured);ToSimulateFlux.Measured];"
    "f = (SimulatedData-ToSimulateData)'*(SigmaInv.^2)*(SimulatedData-ToSimulateData);"
    "%if f<0.1"
    "%disp(f)"
    "%end"
    "% gradient of the objective function"
    "SimulatedGradient = [emuSimulatedGradient;FluxSimulatedGradient];"
    "g = 2*SimulatedGradient'*(SigmaInv.^2)*(SimulatedData-ToSimulateData);"
    "end"
    ];

% 指定要生成的.m文件的文件名
FileName = 'ObjectiveFunction.m';

% 创建.m文件并写入内容
fid = fopen(FileName, 'w'); % open a file
fprintf(fid, '%s\n', FileContent); % write
fclose(fid); % close the file

% 显示生成的.m文件的路径
disp(['生成的.m文件保存在：' pwd '/' FileName]);

%% write the the variance-weighted sum of squared residuals function file
FileContent = [
    "function [phi, XSimulated,dXdVSimulated] = Phi(X,dXdV,ToSimulateSet,NumToSimulate,NumFluxes,VariablesRefernce,SigmaInv)"
    "% 选取优化变量的模拟值和优化变量的梯度值"
    "XSimulated = zeros(NumToSimulate,1);"
    "dXdVSimulated = zeros(NumToSimulate,NumFluxes);"
    "for j = 1:NumToSimulate"
    "    condition = strcmp(VariablesRefernce.emu_mid,ToSimulateSet.emuToSimulate{j});"
    "    XSimulated(j) = X(condition);"
    "    dXdVSimulated(j,:) = dXdV(condition,:);"
    "end"
    ""
    "% phi function"
    "phi = (XSimulated-str2double(ToSimulateSet.emuMID))'*SigmaInv*(XSimulated-str2double(ToSimulateSet.emuMID));"
    "end"
    ];

% 指定要生成的.m文件的文件名
FileName = 'Phi.m';

% 创建.m文件并写入内容
fid = fopen(FileName, 'w'); % open a file
fprintf(fid, '%s\n', FileContent); % write
fclose(fid); % close the file

% 显示生成的.m文件的路径
disp(['生成的.m文件保存在：' pwd '/' FileName]);