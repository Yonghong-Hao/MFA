% % === 最小化OldV和NewV
% function [NewVBackward, exitflag,fval] = S6_BackwardMinNewVOldV(SMatrix,OldV,NumFluxes,NumMetabolites,i,StepSize)
% V0 = OldV*(1-0.001);
% A = [];
% b = [];
% Aeq = [SMatrix;zeros(1,NumFluxes)];
% Aeq(NumMetabolites+1,i) = 1;
% beq = [zeros(NumMetabolites,1)];
% beq(NumMetabolites+1) = OldV(i)-StepSize;
% lb = 1e-10*ones(NumFluxes,1);
% ub = 100*ones(NumFluxes,1);
% myfun = @(NewVBackward) (NewVBackward-OldV)'*(NewVBackward-OldV);
% nonlcon = [];
% options = optimoptions("fmincon","Algorithm","interior-point","MaxFunctionEvaluations",1000000,"MaxIterations",10000,"EnableFeasibilityMode",true,"SubproblemAlgorithm","cg");
% %options = [];
% %options = optimoptions("fmincon","MaxFunctionEvaluations",1000000,"MaxIterations",10000,"ConstraintTolerance",1e-20);
% [NewVBackward, fval, exitflag] = fmincon(myfun,V0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% end
% === 最小化OldV和NewV
function [NewVBackward, exitflag,fval] = S6_BackwardMinNewVOldV(MappingID,SMatrix,OldV,j,StepSize,ToSimulateMS_MID,ToSimulateFlux,VariablesRefernce,SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar,SigmaInv)
V0 = OldV'*(1-0.0001);
lb = log(ones(NumFluxes,1)*1e-10); %这里有些模型是提供流上下界的
ub = log(ones(NumFluxes,1)*100);
switch MappingID
    case 'mapping1'
        nonlcon = @(U) S4_mycon_mapping1_CI(U,SMatrix,j,StepSize,OldV);
    case 'mapping8'
        nonlcon = @(U) S4_mycon_mapping8(U,S);
    case 'mapping9'
        nonlcon = @(U) S4_mycon_mapping9(U,S);
    case 'mapping20'
        nonlcon = @(U) S4_mycon_mapping20(U,S);
    case 'mapping25'
        nonlcon = @(U) S4_mycon_mapping25(U,S);
    case 'mapping26'
        nonlcon = @(U) S4_mycon_mapping26(U,S);
    case 'mapping27'
        nonlcon = @(U) S4_mycon_mapping27(U,S);
end
options = optimoptions('fmincon','MaxIterations',50,'MaxFunctionEvaluations',10000,'SpecifyObjectiveGradient',true,'Algorithm','interior-point');
[U, fval, exitflag] = fmincon(@(U) ObjectiveFunction(U,ToSimulateMS_MID,ToSimulateFlux,VariablesRefernce,SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar,SigmaInv),log(V0),[],[],[],[],lb,ub,nonlcon,options);
NewVBackward = exp(U)';
end