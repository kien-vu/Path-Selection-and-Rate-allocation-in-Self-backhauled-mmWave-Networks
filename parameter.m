
% This function is used to load the optimization solver, and declare the
% simulation paramters.
loadsolvers();
global N_Flow % number of flows 
N_Flow = 2;
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
N_SubF = 4;
global N_BSs % the number of SC BSs
N_BSs = 6;
global N_Actions % The number of actions : Each user is assumed to choose two paths among 4 paths, then 
% the number of actions is the combination of choosing 2 from 4.
N_Actions = 6; % number of actions can be large, which depends on your choise
global Iterations % the number of iterations before eating
Iterations = 10000;
global alpha1 % amount of data, user 1 divides for subflow1 and 2
global alpha2 % amount of data, user 2 divides for subflow3 and 4
alpha1 = log(1 + 4)*8 ; % We assume that each flow is divided equally log(1 + 4)
alpha2 = alpha1; % 
global ops
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0,'debug',1,'beeponproblem',1,'saveduals',1); % set the interal solver to be mosek

