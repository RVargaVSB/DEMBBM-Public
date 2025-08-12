%restoredefaultpath;  X
%addpath(pwd); 

clc
clear
close all

addpath(genpath('ProgramFiles'))
addpath(genpath('Problems'))

profile clear
profile on

% Set problem

%%% Tasks SimplySupportedBeam([x,y,z],Element Radius, "middle/top",
%%% "time/steps", [Force dir start end], [unload time (-1 = no unload) unload speed],stop on crack)

% problem = SimplySupportedBeam([5,1,1],0.1,"middle","time",[-4e5 3 0 0.01],[-1,0],0);
% problem = Konzolka();
problem = TensionProblem([2,1,1],0.1,0.8,1500000);

problem.material = {Material()};
problem.lambda = [0.705399 0.705399];

problem.elements = SetElements(problem);

problem.solver = Solver();
problem.drawer = Drawer();
problem.currentState = CurrentState(problem);

problem = problem.solver.Solve(problem);

profile off
p = profile('info');