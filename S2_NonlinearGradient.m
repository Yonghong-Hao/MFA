function dYdVVar = S2_NonlinearGradient(NumBoundaryVariables,BoundaryVariablesRefernce)
emu_var = BoundaryVariablesRefernce.emu_var;
dYdVVar = cell(NumBoundaryVariables,1);
for i = 1:NumBoundaryVariables
    termsTotal_i = regexp(emu_var{i},'\s\+\s','split');
    termsSum_i = cell(1,numel(termsTotal_i));
    for j = 1:numel(termsTotal_i)
        termsTotal_j = regexp(termsTotal_i{j},'\*','split');
        termsSum_j = cell(1,numel(termsTotal_j));
        for k = 1:numel(termsTotal_j)
            termsTotal_j_k = termsTotal_j{k};
            pattern = regexp(termsTotal_j_k,'X_\d+','match');
            if ~isempty(pattern)
                termsTotal_j_copy = termsTotal_j;
                if ismember('^',termsTotal_j_k)
                    termsTotal_j_k_split = regexp(termsTotal_j_k,'\^','split');
                    termsTotal_j_copy{k} = strcat(termsTotal_j_k_split{2},'*',termsTotal_j_k_split{1},'^(',termsTotal_j_k_split{2},'-1',')','*',regexprep(termsTotal_j_k_split{1},'X_(\d+)', 'dXdV_$1'));
                else
                    termsTotal_j_copy{k} = regexprep(termsTotal_j_k,'X_(\d+)', 'dXdV_$1');
                end
                termsSum_j{k} = strjoin(termsTotal_j_copy,'*');
            else
                termsTotal_j_copy = termsTotal_j;
                termsTotal_j_copy{k} = '0';
                termsSum_j{k} = strjoin(termsTotal_j_copy,'*');
            end

        end
        termsSum_i{j} = strjoin(termsSum_j,' + ');
    end
    dYdVVar{i} = string(expand(str2sym(strjoin(termsSum_i,' + '))));
end
end
