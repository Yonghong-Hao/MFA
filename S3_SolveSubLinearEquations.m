function VariablesTotal = S3_SolveSubLinearEquations(Result_inSize,BoundaryReference,V)
%% 按emu size分类求解线性方程组
VariablesTotal = table('Size',[0 3],'VariableTypes',{'string','double','string'},'VariableNames',{'emu_mid','emu_var','emu_size'});
Result_inSizeNew = Result_inSize;
for s = 1:height(Result_inSize)
    % 计算右端项
    BoundaryReference_s = Result_inSizeNew(s).BoundaryReference_sub;
    b_s = Result_inSizeNew(s).b_sub;
    for r = 1:height(b_s)
        b_sr = b_s{r};
        if isempty(b_sr)
            b_sr = 0;
        else
            b_sr = eval(replace(b_sr,BoundaryReference_s.emuBoundary_mid,BoundaryReference_s.emuBoundary_var));
        end
        Result_inSizeNew(s).b_sub{r} = b_sr;
    end
    % 计算系数矩阵
    A_s = Result_inSizeNew(s).A_sub;
    for i = 1:height(A_s)
        for j = 1:width(A_s)
            A_sij = A_s{i,j};
            if isempty(A_sij)
                A_sij = 0;
            else
                A_sij = eval(char(A_sij));
            end
            Result_inSizeNew(s).A_sub{i,j} = A_sij;
        end
    end
    % 求解x_s
    x_s = cell2mat(Result_inSizeNew(s).A_sub)\cell2mat(Result_inSizeNew(s).b_sub);
    Result_inSizeNew(s).VariablesRefernce_sub.emu_var = x_s;
    VariablesRefernce_s = Result_inSizeNew(s).VariablesRefernce_sub;
    VariablesTotal = vertcat(VariablesTotal,VariablesRefernce_s);
    % 用x_s中的部分替换下一个Size的Boundary变量
    if s+1 <=height(Result_inSize)
        BoundaryReference_sNew = Result_inSizeNew(s+1).BoundaryReference_sub;
        for t = 1:height(BoundaryReference_sNew)
            BoundaryReference_stNew = ReplaceVariableToValue(BoundaryReference_sNew.emuBoundary_mid{t},VariablesTotal.emu_mid,VariablesTotal.emu_var);
            Result_inSizeNew(s+1).BoundaryReference_sub.emuBoundary_var{t} = num2str(eval(replace(BoundaryReference_stNew,BoundaryReference.emuBoundary_mid,BoundaryReference.emuBoundary_var)));
        end
    end
end
end