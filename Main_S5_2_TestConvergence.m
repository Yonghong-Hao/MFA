clear;clc
%% input
TracerID = 'glc_12';
MappingID = 'mapping9';
DataSet = {'data1_1','data1_2','data1_3'};
Conditions = {'CM','LA+','LA'};
%% Test Convergence
% read convergence data
Convergence_1 = readtable(strcat('Tracer_',TracerID,'/VerifyMyModel/ConvergenceRate/ConvergenceRate_',MappingID,'_',DataSet{1},'.xlsx'),"VariableNamingRule","preserve");
Convergence_2 = readtable(strcat('Tracer_',TracerID,'/VerifyMyModel/ConvergenceRate/ConvergenceRate_',MappingID,'_',DataSet{2},'.xlsx'),"VariableNamingRule","preserve");
Convergence_3 = readtable(strcat('Tracer_',TracerID,'/VerifyMyModel/ConvergenceRate/ConvergenceRate_',MappingID,'_',DataSet{3},'.xlsx'),"VariableNamingRule","preserve");
Convergence_total = [(50:50:50*60)',Convergence_1.LossfunctionValue,Convergence_2.LossfunctionValue,Convergence_3.LossfunctionValue];

% plot
f2 = figure;
set(groot,'DefaultAxesFontSize',15);
f2.Position = [1945         376        1849         489];
for i = 1:3
    subplot(1,3,i)
    scatter(Convergence_total(:,1),log10(Convergence_total(:,1+i)),25,[0.2588    0.5725    0.7765],"filled")
    %scatter(Convergence_total(:,1),Convergence_total(:,1+i),25,[0.2588    0.5725    0.7765],"filled")
    xlabel('Iterations');
    ylabel('log10(loss function values)')
    %ylim([0,4])
    title(Conditions{i},"FontSize",15,"FontWeight","bold");
end
saveas(f2,strcat('Tracer_',TracerID,'/VerifyMyModel/ConvergenceRate/ConvergenceRate_',MappingID),'png')