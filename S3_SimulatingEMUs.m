function [X, dXdV] = S3_SimulatingEMUs(SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,V,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar)
SizeSet = unique(SizeVariables);
X = zeros(NumVariables,1);
dXdV = zeros(NumVariables,NumFluxes);
%% 计算Size1中相关变量
s = 1;
indexAVar_1 = SizeVariables==SizeSet(s);
dAdV_1 = dAdV(indexAVar_1,indexAVar_1);
indexBVar_1 = SizeBoundaryVariables==SizeSet(s);
dBdV_1 = dBdV(indexAVar_1,indexBVar_1);
YVar_1 = YVar(indexBVar_1);
[NumAVariables_1, NumBVariables_1]  = size(dBdV_1);

% 计算矩阵AVar_1-->Aval_1（代入V）
AVal_1 = zeros(NumAVariables_1,NumAVariables_1);
for i = 1:NumAVariables_1
    for j = 1:NumAVariables_1
        dAdV_1ij = dAdV_1{i,j};
        if ~isempty(dAdV_1ij)
            AVal_1(i,j) = dAdV_1ij*V;
        end
    end
end

% 计算矩阵BVar_1-->BVal_1（代入V）
BVal_1 = zeros(NumAVariables_1,NumBVariables_1);
for i = 1:NumAVariables_1
    for j = 1:NumBVariables_1
        dBdV_1ij = dBdV_1{i,j};
        if ~isempty(dBdV_1ij)
            BVal_1(i,j) = dBdV_1ij*V;
        end
    end
end

% 计算向量YVar_1-->YVal_1（代入X）
Formula = YVar_1;
YVal_1 = zeros(NumBVariables_1,1);
for i = 1:NumBVariables_1
    Formula_i = Formula{i};
    if isempty(regexp(Formula_i,'X','match'))
        % if strcmp(Formula_i,'0')
        %     YVal_s(i) = 0;
        % elseif strcmp(Formula_i,'1')
        %     YVal_s(i) = 1;
        YVal_1(i) = str2double(Formula_i);
    else
        sumIterms = regexp(Formula_i,' \+ ','split');
        sumItem_num = 0;
        for j = 1:numel(sumIterms)
            multiIterms = regexp(sumIterms{j},'\*','split');
            multiIterms_num = 1;
            for k = 1:numel(multiIterms)
                index = ismember(XVar,multiIterms(k));
                multiper_new = X(index);
                multiIterms_num = multiIterms_num*multiper_new;
            end
            sumItem_num = sumItem_num+multiIterms_num;
        end
        YVal_1(i) = sumItem_num;
    end
end

% 计算向量X，并更新X
XVal_1 = S3_SolveLinearEquations(AVal_1,BVal_1*YVal_1);
X(indexAVar_1) = XVal_1;

% 计算向量dXdV（代入AVal, dAdV, dBdV, YVal）
iterm_1 = zeros(NumAVariables_1,NumFluxes);
for i = 1:NumAVariables_1
    iterm_1i = zeros(1,NumFluxes);
    for j = 1:NumAVariables_1
        dAdV_1ij = dAdV_1{i,j};
        if ~isempty(dAdV_1ij)
            iterm_1i = iterm_1i+dAdV_1ij*XVal_1(j);
        end
    end
    iterm_1(i,:) = iterm_1i;
end
iterm_2 = zeros(NumAVariables_1,NumFluxes);
for i = 1:NumAVariables_1
    iterm_2i = zeros(1,NumFluxes);
    for j = 1:NumBVariables_1
        dBdV_1ij = dBdV_1{i,j};
        if ~isempty(dBdV_1ij)
            iterm_2i = iterm_2i+dBdV_1ij*YVal_1(j);
        end
    end
    iterm_2(i,:) = iterm_2i;
end
dXdV_1 = pinv(AVal_1)*(iterm_2-iterm_1);
dXdV(indexAVar_1,:) = dXdV_1;

%% 计算Size 1+中相关变量
for s = 2:numel(SizeSet)
    indexAVar_s = SizeVariables==SizeSet(s);
    dAdV_s = dAdV(indexAVar_s,indexAVar_s);
    indexBVar_s = SizeBoundaryVariables==SizeSet(s);
    dBdV_s = dBdV(indexAVar_s,indexBVar_s);
    YVar_s = YVar(indexBVar_s);
    dYdVVar_s = dYdVVar(indexBVar_s);
    [NumAVariables_s, NumBVariables_s]  = size(dBdV_s);
    % 计算矩阵AVar_s-->Aval_s（代入V）
    AVal_s = zeros(NumAVariables_s,NumAVariables_s);
    for i = 1:NumAVariables_s
        for j = 1:NumAVariables_s
            dAdV_sij = dAdV_s{i,j};
            if ~isempty(dAdV_sij)
                AVal_s(i,j) = dAdV_sij*V;
            end
        end
    end

    % 计算矩阵BVar_s--BVal_s（代入V）
    BVal_s = zeros(NumAVariables_s,NumBVariables_s);
    for i = 1:NumAVariables_s
        for j = 1:NumBVariables_s
            dBdV_sij = dBdV_s{i,j};
            if ~isempty(dBdV_sij)
                BVal_s(i,j) = dBdV_sij*V;
            end
        end
    end

    % 计算向量YVar_s-->YVal_s（代入X）
    Formula = YVar_s;
    YVal_s = zeros(NumBVariables_s,1);
    for i = 1:NumBVariables_s
        Formula_i = Formula{i};
        if isempty(regexp(Formula_i,'X','match'))
            YVal_s(i) = str2double(Formula_i);
        else
            sumIterms = regexp(Formula_i,' \+ ','split');
            sumItem_num = 0;
            for j = 1:numel(sumIterms)
                multiIterms = regexp(sumIterms{j},'\*','split');
                multiIterms_num = 1;
                for k = 1:numel(multiIterms)
                    multiIterms_k = multiIterms{k};
                    pattern = regexp(multiIterms_k,'X_\d+','match');
                    if isempty(pattern)
                        multiper_new = str2double(multiIterms_k);
                    elseif ismember('^',multiIterms_k)
                        multiIterms_k_split = regexp(multiIterms_k,'\^','split');
                        index = ismember(XVar,multiIterms_k_split(1));
                        multiper_new = X(index)^str2double(multiIterms_k_split(2));
                    else
                        index = ismember(XVar,multiIterms_k);
                        multiper_new = X(index);
                    end
                    multiIterms_num = multiIterms_num*multiper_new;
                end
                sumItem_num = sumItem_num+multiIterms_num;
            end
            YVal_s(i) = sumItem_num;
        end
    end

    % 计算向量X，并更新X
    XVal_s = S3_SolveLinearEquations(AVal_s,BVal_s*YVal_s);
    X(indexAVar_s) = XVal_s;

    % 计算矩阵dYdV_s（代入dXdV, XVal）
    Formula = dYdVVar_s;
    dYdV_s = zeros(NumBVariables_s,NumFluxes);
    for i = 1:NumBVariables_s
        Formula_i = Formula{i};
        if strcmp(Formula_i,'0')
            dYdV_s(i) = 0;
        else
            sumIterms = regexp(Formula_i,' \+ ','split');
            sumItem_num = zeros(1,NumFluxes);
            for j = 1:numel(sumIterms)
                multiIterms = regexp(sumIterms{j},'\*','split');
                multiIterms_num = 1;
                for k = 1:numel(multiIterms)
                    multiIterms_k = multiIterms{k};
                    if ismember(multiIterms_k,dXdVVar) % dXdV iterm -- vector
                        index = ismember(dXdVVar,multiIterms_k);
                        multiper_new = dXdV(index,:);                        
                    elseif ismember('^',multiIterms_k) % X^m iterm
                        multiIterms_k_split = regexp(multiIterms_k,'\^','split');
                        index = ismember(XVar,multiIterms_k_split{1});
                        multiper_new = X(index)^str2double(multiIterms_k_split{2});
                    elseif ismember(multiIterms_k,XVar) % X iterm
                        index = ismember(XVar,multiIterms_k);
                        multiper_new = X(index);
                    else % number iterm
                        multiper_new = str2double(multiIterms_k);
                    end
                    multiIterms_num = multiIterms_num*multiper_new;
                end
                sumItem_num = sumItem_num+multiIterms_num;
            end
            dYdV_s(i,:) = sumItem_num;
        end
    end

    % 计算矩阵dXdV（代入AVal, dAdV, dBdV, YVal）
    iterm_1 = zeros(NumAVariables_s,NumFluxes);
    for i = 1:NumAVariables_s
        iterm_1i = zeros(1,NumFluxes);
        for j = 1:NumAVariables_s
            dAdV_sij = dAdV_s{i,j};
            if ~isempty(dAdV_sij)
                iterm_1i = iterm_1i+dAdV_sij*XVal_s(j);
            end
        end
        iterm_1(i,:) = iterm_1i;
    end
    iterm_2 = zeros(NumAVariables_s,NumFluxes);
    for i = 1:NumAVariables_s
        iterm_2i = zeros(1,NumFluxes);
        for j = 1:NumBVariables_s
            dBdV_sij = dBdV_s{i,j};
            if ~isempty(dBdV_sij)
                iterm_2i = iterm_2i+dBdV_sij*YVal_s(j);
            end
        end
        iterm_2(i,:) = iterm_2i;
    end
    iterm_3 = BVal_s*dYdV_s;
    dXdV_s = pinv(AVal_s)*(iterm_2+iterm_3-iterm_1);
    dXdV(indexAVar_s,:) = dXdV_s;
end
end