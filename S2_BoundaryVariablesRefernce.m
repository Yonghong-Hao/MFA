function BoundaryVariablesRefernce = S2_BoundaryVariablesRefernce(EMUMIDTable,VariablesRefernce,BoundaryReference)
[emu_mid, ia] = setdiff(EMUMIDTable.ReactantMid,EMUMIDTable.ProductMid,'stable');
emu_mid = string(expand(str2sym(emu_mid)));
emu_size = EMUMIDTable.Size(ia);
emu_var = cell(numel(emu_mid),1);
ReferenceTotal = [VariablesRefernce.emu_mid,VariablesRefernce.emu_var;BoundaryReference.emuBoundary_mid,BoundaryReference.emuBoundary_var];
for i = 1:numel(emu_mid)
    Formula_i = emu_mid{i};
    sumIterms = regexp(Formula_i,' \+ ','split');
    sumItem_num = '0';
    for j = 1:numel(sumIterms)
        multiIterms = regexp(sumIterms{j},'\*','split');
        multiIterms_num = '1';
        for k = 1:numel(multiIterms)
            multiIterms_k = multiIterms{k};
            if isempty(regexprep(multiIterms_k,'\d+',''))
                multiper_new = multiIterms_k;
            elseif ismember('^',multiIterms_k)
                multiper_split = regexp(multiIterms_k,'\^','split');
                index = ismember(ReferenceTotal(:,1),multiper_split(1));
                multiper_new = strcat(ReferenceTotal(index,2),'^',multiper_split(2));
            else
                index = ismember(ReferenceTotal(:,1),multiIterms_k);
                multiper_new = ReferenceTotal(index,2);
            end
            multiIterms_num = strcat(multiIterms_num,'*',multiper_new);
        end
        sumItem_num = strcat(sumItem_num,' + ',multiIterms_num);
    end
    emu_var{i} = string(expand(str2sym(sumItem_num)));
end
BoundaryVariablesRefernce = table(emu_mid,emu_var,emu_size);
end
