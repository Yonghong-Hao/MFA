% === 最小化OldV和NewV
function [NewVForward, exitflag, fval] = S6_ForwardMinNewVOldV(SMatrix,OldV,NumFluxes,NumMetabolites,i,StepSize)
V0 = OldV;
A = [];
b = [];
Aeq = [SMatrix;zeros(1,NumFluxes)];
Aeq(NumMetabolites+1,i) = 1;
beq = [zeros(NumMetabolites,1)];
beq(NumMetabolites+1) = OldV(i)+StepSize;
lb = 1e-10*ones(NumFluxes,1);
ub = 100*ones(NumFluxes,1);
myfun = @(NewVForward) (NewVForward-OldV)'*(NewVForward-OldV);
nonlcon = [];
options = optimoptions("fmincon","Algorithm","interior-point","MaxFunctionEvaluations",1000000,"MaxIterations",10000,"EnableFeasibilityMode",true,"SubproblemAlgorithm","cg");
%options = optimoptions("fmincon","MaxFunctionEvaluations",1000000,"MaxIterations",10000,"ConstraintTolerance",1e-20);
[NewVForward, fval, exitflag] = fmincon(myfun,V0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end