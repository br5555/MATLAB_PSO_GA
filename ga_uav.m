clc
clear all
num_vars = 1;
Fitnessfunc = @uav_error_function_yaw;%@uav_error_function;%@uav_error_function_lin;
upper = 1e11*ones(1,num_vars);
lower = 0*ones(1,num_vars);
options = gaoptimset;
%options = gaoptimset(options,'MutationFcn', @mutationgaussian);
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'TolFun',1e-15);
options = gaoptimset(options,'PopInitRange',[lower ;upper]);
options = gaoptimset(options,'Generations',100);
options = gaoptimset(options,'PopulationSize',1000);
options = gaoptimset(options,'UseParallel', 'always');

%options = gaoptimset(options,'SelectionFcn',{@selectionroulette});                     
%options = gaoptimset(options,'MutationFcn',{@mutationuniform, 0.1});                      
%options = gaoptimset(options,'CrossoverFcn', {@crossoversinglepoint});                      
options = gaoptimset(options,'EliteCount',5);                      
options = gaoptimset(options,'StallGenLimit',100);                      
options = gaoptimset(options,'PlotFcns',{@gaplotbestf, @gaplotbestindiv, @gaplotexpectation, @gaplotgenealogy});                       
%options = gaoptimset(options,'UseParallel', true);                      
options = gaoptimset(options,'Vectorized', 'off');                      
options = gaoptimset(options,'Display', 'iter'); 
%options = gaoptimset(options,'OutputFcn',@gaoutfun);
IntCon = 1;
%%

[x,fval,exitflag,output,population,score] = ...
ga(Fitnessfunc,num_vars,[],[],[],[],lower, upper,[],IntCon,options);