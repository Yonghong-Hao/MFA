%==========================================================================
%============================ Flux Estimation =============================
%==========================================================================
%clear;clc
profile on
% input
CalculationTimes = 100;
Flux_ub = 10000;
MaxIterations = 30;
S = SMatrix;
T = S4_ConstraintForX0(VariablesRefernce);
[Tsize1,Tsize2] = size(T);

%% initial parameters
%============================= LogTransform ===============================
% some numbers
NumVariables = height(VariablesRefernce);
[NumMetabolites,NumFluxes] = size(S);

% some initial values
rng default %保证可复现
V0Matrix = lhsdesign(CalculationTimes,NumFluxes)*Flux_ub;

%  some storage matrix
VStore = zeros(NumFluxes,CalculationTimes);
fvalStore = zeros(1,CalculationTimes);

% optimal parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb = log(ones(NumFluxes,1)*1e-10); %这里有些模型是提供流上下界的
ub = log(ones(NumFluxes,1)*Flux_ub);

%mAb_f = ToSimulateFlux.FluxData(strcmp(ToSimulateFlux.ReactionID,'mAb'));
switch MappingID
    case 'mapping1'
        nonlcon = @(U) S4_mycon_mapping1(U,S);
    case 'mapping9'
        nonlcon = @(U) S4_mycon_mapping9(U,S);
    case 'mapping20'
        nonlcon = @(U) S4_mycon_mapping20(U,S);
end
Iterations = 50;
options = optimoptions('fmincon','MaxIterations',Iterations,'MaxFunctionEvaluations',10000,'SpecifyObjectiveGradient',true);

%% solve the optimal problem
ConvergenceRate = [];
k = 0;
while true
        %nonlcon = @(U) S4_mycon(U,S,addingConS,addingConb);
        Delresults = zeros(CalculationTimes,1);
        exitflagStore = Delresults;
        parfor i=1:CalculationTimes
            V0 = V0Matrix(i,:);
            U0 = log(V0);
            % solve
            [U, fval, exitflag, output,~,grad,~] = fmincon(@(U) ObjectiveFunction(U,ToSimulateMS_MID,ToSimulateFlux,VariablesRefernce,SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar,SigmaInv),U0,A,b,Aeq,beq,lb,ub,nonlcon,options);
            %disp(grad)
            fvalStore(i) = fval;
            V = exp(U);
            VStore(:,i) = V;
            exitflagStore(i) = exitflag;
            if exitflag<0
                Delresults(i) = 1;
            end
        end
   
    %% choose the best solution from NumSamples solutions in each tissue
    % delete bad solutions
    fvalStore(:,logical(Delresults)) = [];
    VStore(:,logical(Delresults)) = [];
    % choose the best solution
    MinRes = min(fvalStore);
    ConvergenceRate = [ConvergenceRate;[(k+1)*Iterations,MinRes]];
    disp(MinRes)
    CalculationTimes = width(VStore);
    V0Matrix = VStore';
    k = k+1;
    if k >MaxIterations
        break
    end
end
index_Phihat = find(fvalStore==MinRes);
Vhat = VStore(:,index_Phihat);

% congvergence rate
scatter(ConvergenceRate(:,1),log10(ConvergenceRate(:,2)),'filled');
%writetable(addvars(MappingAtomTable,Vhat,'After','product_IDs(atom)'),strcat('Tracer_gln_12345/Results/Vhat_',ModelID,'.xlsx'));
%writematrix(VStore,strcat('Tracer_gln_12345/Results/VStore_',ModelID,'.xlsx'));

profile off
profile viewer