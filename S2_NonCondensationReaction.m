function midTerms = S2_NonCondensationReaction(reactant_r,product_r,flux_r,Size_r)
%% emu reactant
emu_midR = strcat(reactant_r,'_m',cellfun(@num2str, num2cell(0:Size_r),'UniformOutput',false));

%% emu product
emu_midP = strcat(product_r,'_m',cellfun(@num2str, num2cell(0:Size_r),'UniformOutput',false));

%% emu flux
emu_flux = repmat(flux_r,Size_r+1,1);
% 编辑flux格式
for k = 1:numel(emu_flux)
    str = emu_flux{k};
    if ismember('+',str)
        str = strcat('(',str,')');
    end
    emu_flux{k} = str;
end

%% emu size
emu_size = repmat(Size_r,Size_r+1,1);

%% combine emu reactant and product
midTerms = table(emu_midR',emu_midP',emu_flux,emu_size,'VariableNames',{'ReactantMid','ProductMid','Flux','Size'});
end