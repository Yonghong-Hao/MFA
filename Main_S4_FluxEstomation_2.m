%==========================================================================
%============================ Flux Estimation =============================
%==========================================================================
%clear;clc
profile on
% input
CalculationTimes = 500;
Iterations = 50;
IntervalCount = 60;
Flux_ub = 100;
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
%V0Matrix = VStore';

%  some storage matrix
VStore = zeros(NumFluxes,CalculationTimes);
fvalStore = zeros(1,CalculationTimes);
iterationStore = zeros(1,CalculationTimes);
sampleTotal = regexprep(strcat('sample',cellstr(num2str((1:CalculationTimes)'))),' ','');

% optimal parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb = log(ones(NumFluxes,1)*1e-10); %这里有些模型是提供流上下界的
ub = log(ones(NumFluxes,1)*Flux_ub);

switch MappingID
    case 'mapping1'
        nonlcon = @(U) S4_mycon_mapping1(U,S);
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
options = optimoptions('fmincon','MaxIterations',Iterations,'MaxFunctionEvaluations',10000,'SpecifyObjectiveGradient',true);

%% solve the optimal problem
ConvergenceRate = zeros(IntervalCount,2);
IterationCount0 = 0;
for j = 1:IntervalCount
    Delresults = zeros(CalculationTimes,1);
    exitflagStore = Delresults;
    %sampleTotal = sampleTotal0;
    parfor i=1:CalculationTimes
        V0 = V0Matrix(i,:);
        U0 = log(V0);
        % solve
        [U, fval, exitflag, output,~,grad,~] = fmincon(@(U) ObjectiveFunction(U,ToSimulateMS_MID,ToSimulateFlux,VariablesRefernce,SizeVariables,SizeBoundaryVariables,NumVariables,NumFluxes,dAdV,dBdV,XVar,YVar,dXdVVar,dYdVVar,SigmaInv),U0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        %disp(output)
        iterationStore(i) = output.iterations;
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
    sampleTotal(logical(Delresults)) = [];
    % choose the best solution
    [MinRes,index_Phihat_j] = min(fvalStore);
    IterationCount = iterationStore(index_Phihat_j);
    ConvergenceRate(j,:) = [IterationCount0+IterationCount,MinRes];
    IterationCount0 = IterationCount0+IterationCount;
    disp(MinRes)
    % save the loss function values and the flux values
    VStoreTable = array2table(VStore, 'RowNames', fluxTotal','VariableNames',sampleTotal);
    fvalStoreTable = array2table(fvalStore,'VariableNames',sampleTotal);
    % 保存每次流估计和损失函数值的结果
    writetable(VStoreTable,strcat(TracerID,'/FluxValues_LossFunctionValues/',MappingID,'/',Conditions{str2double(regexprep(DataID,'data1_',''))},'/FluxValues_',ModelID,'_iteration',num2str(IntervalCount),'_',num2str(j),'.xlsx'),'WriteRowNames',true);
    writetable(fvalStoreTable,strcat(TracerID,'/FluxValues_LossFunctionValues/',MappingID,'/',Conditions{str2double(regexprep(DataID,'data1_',''))},'/LossFunctionValues_',ModelID,'_iteration',num2str(IntervalCount),'_',num2str(j),'.xlsx'));

    CalculationTimes = width(VStore);
    V0Matrix = VStore';
end
% congvergence rate
scatter(ConvergenceRate(:,1),log10(ConvergenceRate(:,2)),"filled");
writetable(table(ConvergenceRate(:,1),ConvergenceRate(:,2),'VariableNames',{'Iterations','LossfunctionValue'}),strcat(TracerID,'/VerifyMyModel/ConvergenceRate/ConvergenceRate_',ModelID,'.xlsx'));

profile off
profile viewer