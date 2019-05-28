clc
clear all
close all

num_vars = 9;
Fitnessfunc = @uav_error_function;%@uav_error_function_lin;
upper = 1e11*ones(1,num_vars);
lower = 0*ones(1,num_vars);

options = optimoptions('particleswarm', 'SwarmSize', 100,'Display','iter','PlotFcn','pswplotbestf', 'HybridFcn', {@fmincon},'MaxIterations',200,  'UseParallel' ,true);

tic
[kOptimized,fval,exitflag,output] = particleswarm(Fitnessfunc, num_vars, lower, upper, options);
toc