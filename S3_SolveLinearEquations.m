%% 线性方程的最小范数最小二乘解
function X = S3_SolveLinearEquations(A,b)
X = lsqminnorm(A,b,1e-8);

% %% 带约束最小二乘问题
%function X = S3_SolveLinearEquations(C,d,s,NumVariables_s)
% %约束
% A = [];
% b = [];
% NumsEMU = NumVariables_s/(s+1);
% Aeq = zeros(NumsEMU,NumVariables_s);
% for i = 1:NumsEMU
% Aeq(i,(i-1)*(s+1)+1:i*(s+1)) = ones(1,s+1);
% end
% beq = ones(NumsEMU,1);
% lb = zeros(NumVariables_s,1);
% ub = ones(NumVariables_s,1);
% %求解
% X = lsqlin(C,d,A,b,Aeq,beq,lb,ub);
end