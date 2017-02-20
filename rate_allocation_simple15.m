function[served_rate, transmit_power] = rate_allocation_simple15(vq, selectedaction, delay_r) 
global N_Flow % number of flows 
N_Flow = 2;
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
N_SubF = 4;
global N_BSs % the number of SC BSs
N_BSs = 6;
global N_Actions % The number of actions : Each user is assumed to choose two paths among 4 paths, then 
% the number of actions is the combination of choosing 2 from 4.
N_Actions = 6;
global ops
% selectedaction = [1 1]';
% vq = 10*rand( N_SubF, 1 + N_BSs);
[nSF, nBS] = size(vq); % number of subflows and BSs
Delay_const = zeros(nBS, nSF);
indicator_BSs = routingtable(selectedaction);
arrival  =  zeros(nSF, 1);
for nbs = 1:nBS
    % Transmit power constraints
    for nsf = 1:nSF
        Delay_const(nbs, nsf) = min(delay_r(nsf, nbs), log(1+5) );
    end
end
% Routing table
RT.action(1).route1 = [1 2 3];
RT.action(1).route2 = [1 4 5];
RT.action(2).route1 = [1 2 3];
RT.action(2).route2 = [1 6];
RT.action(3).route1 = [1 2 3];
RT.action(3).route2 = [1 7];
RT.action(4).route1 = [1 4 5];
RT.action(4).route2 = [1 7];
RT.action(5).route1 = [1 4 5];
RT.action(5).route2 = [1 6];
RT.action(6).route1 = [1 7];
RT.action(6).route2 = [1 6];
%% 2 Solving Rate Allocation Problem, which uses SOCP method, Yalmip toolbox
% ops = sdpsettings('solver','mosek'); % fmincon set the interal solver to be mosek
  ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0,'debug',1,'beeponproblem',1,'saveduals',1); % set the interal solver to be mosek
%   lowerbound = 0.5 * rand(nBS, nSF);
%   lowerbound(1, :) = 0.2 *rand(1, nSF);
% Assign the lower bound
  lowerbound = Delay_const;

% Declare the slack variables
  u_next = ones(nBS, nSF);
  tx_next = zeros(nBS, nSF);
  n = 15; % ensuring that error accurancy is less than 10^-5 
  ObjFun = 1e-8;
  lowerbound = lowerbound.* indicator_BSs;
while(1)% Use to find the adaptive control parameter by using convex concave procedure
% Define variable
% rate        = sdpvar(nBS, nSF);
    tx          = sdpvar(nBS, nSF);
    x_max       = sdpvar(nSF,1);
    % Nonconvex variables
    u         = sdpvar(nBS, nSF);
    kappa     = sdpvar(n+4,nSF,'full'); % for convex approximation of exponential function
% Define constraints and Objectives
    constraints = []; % contain all the constraints
    constraints = [constraints, x_max >= 0];
    constraints = [constraints, tx >= 0];
    constraints = [constraints, u >= 1];
    for nbs = 1:nBS
        % Transmit power constraints
        for nsf = 1:nSF
            if indicator_BSs(nbs,nsf) == 0
                constraints = [constraints, tx(nbs, nsf) == 0];
            end
            constraints = [constraints, tx(nbs, nsf) <= 5];
        end
        if nbs >= 2
            constraints = [constraints, sum(tx(nbs, :)) <= 10]; % maximum power constraint
        end
    end
        constraints = [constraints, sum(tx(1, :)) <= 20]; % maximum power constraint
        
% Here is to deal with the approximated convex constraints
for player = 1:N_Flow % rate_allocated = zeros(N_BSs+1, N_SubF);
        if player == 1
            UE = 8;
            [~, sizea] = size(RT.action(selectedaction(player)).route1);
            bs_list1 = RT.action(selectedaction(player)).route1;
            for n1 = 1:sizea
                if n1 == 1 % if MBS
%                     constraints = [constraints, (1+ tx(bs_list1(n1), 1))
%                     >= exp(5+x_max(bs_list2(n1), 1))  ]; which is
%                     replaced by the below constraints, 20170213
                    k = 1;
                    constraints = [constraints, (1+ tx(bs_list1(n1), k)) >= exp(lowerbound(bs_list1(n1), k))   ];
                    constraints = [constraints,  kappa(:,k) >= 0];
                    constraints = [constraints, (1+ tx(bs_list1(n1), k)) >= kappa(1,k)  ];
                    constraints = [constraints, cone([2 + (  x_max(k, bs_list1(n1)) )/(2^(n-1)); 1-kappa(2,k)], kappa(2,k) +1)];
                    constraints = [constraints, cone([5/3 + (  x_max(k, bs_list1(n1)))/(2^(n)); 1-kappa(3,k)], kappa(3,k) +1)];
                    constraints = [constraints, cone([2*kappa(2,k); 1-kappa(4,k)],kappa(4,k) +1)];
                    constraints = [constraints, 19/72 + kappa(3,k) + 1/24*kappa(4,k) <= kappa(5,k)];
                    for  mVar = 6:n+4
                        constraints = [constraints, cone([2*kappa(mVar-1,k);1-kappa(mVar,k)], kappa(mVar,k) +1)];
                    end
                    constraints = [constraints, cone([2*kappa(n+4,k); 1-kappa(1,k)], 1+kappa(1,k))];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list1(n1), 1); tx(bs_list1(n1), 1)./2], (2+ tx(bs_list1(n1), 1))./2  ) ];
                    constraints = [constraints, 2 * u(bs_list1(n1), 1)*u_next(bs_list1(n1), 1)./(1 + tx_next(bs_list1(n1-1), 1))  - (1+ tx(bs_list1(n1-1), 1)) * (u_next(bs_list1(n1), 1)./(1 + tx_next(bs_list1(n1-1), 1))).^2 >= exp(lowerbound(bs_list1(n1), 1)) ];
                end
            end
            % 20170210 stop here                  
            [~, sizeb] = size(RT.action(selectedaction(player)).route2);
            bs_list2 = RT.action(selectedaction(player)).route2;
            for n1 = 1:sizeb
                if n1 == 1 % if MBS
                    k = 2;
                    constraints = [constraints, (1+ tx(bs_list2(n1), k)) >= exp(lowerbound(bs_list2(n1), k))   ];
                    constraints = [constraints,  kappa(:,k) >= 0];
                    constraints = [constraints, (1+ tx(bs_list2(n1), k)) >= kappa(1,k)  ];
                    constraints = [constraints, cone([2 + (  x_max(k, bs_list2(n1)) )/(2^(n-1)); 1-kappa(2,k)], kappa(2,k) +1)];
                    constraints = [constraints, cone([5/3 + (  x_max(k, bs_list2(n1)))/(2^(n)); 1-kappa(3,k)], kappa(3,k) +1)];
                    constraints = [constraints, cone([2*kappa(2,k); 1-kappa(4,k)],kappa(4,k) +1)];
                    constraints = [constraints, 19/72 + kappa(3,k) + 1/24*kappa(4,k) <= kappa(5,k)];
                    for  mVar = 6:n+4
                        constraints = [constraints, cone([2*kappa(mVar-1,k);1-kappa(mVar,k)], kappa(mVar,k) +1)];
                    end
                    constraints = [constraints, cone([2*kappa(n+4,k); 1-kappa(1,k)], 1+kappa(1,k))];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list2(n1), 2); tx(bs_list2(n1), 2)./2], (2+ tx(bs_list2(n1), 2))./2) ];
                    constraints = [constraints, 2* u(bs_list2(n1), 2) * u_next(bs_list2(n1), 2)./(1 + tx_next(bs_list2(n1-1), 2))  - (1+ tx(bs_list2(n1-1), 2)) * (u_next(bs_list2(n1), 2)./(1+tx_next(bs_list2(n1-1), 2))).^2 >= exp(lowerbound(bs_list2(n1), 2)) ];
                end
            end        
            else % if player 2, then update subflow 3 & 4
            UE = 9;
            [~, sizea] = size(RT.action(selectedaction(player)).route1);
            bs_list1 = RT.action(selectedaction(player)).route1;
            for n1 = 1:sizea
                if n1 == 1 % if MBS
                    k = 3;
                    constraints = [constraints, x_max(3, bs_list1(n1))>= (lowerbound(bs_list1(n1), k))   ];
                    constraints = [constraints,  kappa(:,k) >= 0];
                    constraints = [constraints, (1+ tx(bs_list1(n1), 3)) >= kappa(1,k)  ];
                    constraints = [constraints, cone([2 + ( x_max(3, bs_list1(n1)) )/(2^(n-1)); 1-kappa(2,k)], kappa(2,k) +1)];
                    constraints = [constraints, cone([5/3 + ( x_max(3, bs_list1(n1)))/(2^(n)); 1-kappa(3,k)], kappa(3,k) +1)];
                    constraints = [constraints, cone([2*kappa(2,k); 1-kappa(4,k)],kappa(4,k) +1)];
                    constraints = [constraints, 19/72 + kappa(3,k) + 1/24*kappa(4,k) <= kappa(5,k)];
                    for  mVar = 6:n+4
                        constraints = [constraints, cone([2*kappa(mVar-1,k);1-kappa(mVar,k)], kappa(mVar,k) +1)];
                    end
                    constraints = [constraints, cone([2*kappa(n+4,k); 1-kappa(1,k)], 1+kappa(1,k))];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list1(n1), 3); tx(bs_list1(n1), 3)./2], (2+ tx(bs_list1(n1), 3))./2) ];
                    constraints = [constraints, 2* u(bs_list1(n1), 3) * u_next(bs_list1(n1), 3)./(1 + tx_next(bs_list1(n1-1), 3))  - (1+ tx(bs_list1(n1-1), 3)) * (u_next(bs_list1(n1), 3)./(1+tx_next(bs_list1(n1-1), 3))).^2 >= exp(lowerbound(bs_list1(n1), 3)) ];
                end
            end
            % 20170210 stop here                   
            [~, sizeb] = size(RT.action(selectedaction(player)).route2);
            bs_list2 = RT.action(selectedaction(player)).route2;
            for n1 = 1:sizeb
                if n1 == 1 % if MBS
                    k = 4;
                    constraints = [constraints, x_max(4, bs_list2(n1)) >= (lowerbound(bs_list2(n1), k))   ];
                    constraints = [constraints,  kappa(:,k) >= 0];
                    constraints = [constraints, (1+ tx(bs_list2(n1), 4)) >= kappa(1,k)  ];
                    constraints = [constraints, cone([2 + (  x_max(4, bs_list2(n1)) )/(2^(n-1)); 1-kappa(2,k)], kappa(2,k) +1)];
                    constraints = [constraints, cone([5/3 + ( x_max(4, bs_list2(n1)))/(2^(n)); 1-kappa(3,k)], kappa(3,k) +1)];
                    constraints = [constraints, cone([2*kappa(2,k); 1-kappa(4,k)],kappa(4,k) +1)];
                    constraints = [constraints, 19/72 + kappa(3,k) + 1/24*kappa(4,k) <= kappa(5,k)];
                    for  mVar = 6:n+4
                        constraints = [constraints, cone([2*kappa(mVar-1,k);1-kappa(mVar,k)], kappa(mVar,k) +1)];
                    end
                    constraints = [constraints, cone([2*kappa(n+4,k); 1-kappa(1,k)], 1+kappa(1,k))];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list2(n1), 4); tx(bs_list2(n1), 4)./2], (2+ tx(bs_list2(n1), 4))./2 ) ];
                    constraints = [constraints, 2* u(bs_list2(n1), 4) * u_next(bs_list2(n1), 4)./(1 + tx_next(bs_list2(n1-1), 4))  - (1+ tx(bs_list2(n1-1), 4)) * (u_next(bs_list2(n1), 4)./(1+tx_next(bs_list2(n1-1), 4))).^2 >= exp(lowerbound(bs_list2(n1), 4)) ];
                end
            end
        end   
end

% Objective function
     obj =  -sum(vq(:,1)' * x_max);
% Solve the problem
     sol = optimize(constraints, obj, ops); % solve the problem optimize replaced sdpsolve
     x_max = value(x_max);
     tx = value(tx);
     u  = value(u);
     obj = value(obj);
% Check the results
     sol.info;
     if sol.problem == 0
         tx_next  = tx;
         u_next = u;
         temp = - obj;
         if temp <= ObjFun + 1e-2   % abs(temp - ObjFun) <= 1e-2 
            transmit_power = tx;
            % calculate the serving rate
            for nbs = 1:nBS
                for nsf = 1:nSF
                    served_rate(nsf,nbs) = log(1 + tx(nbs,nsf));
                end
            end          
            break;
         end
         ObjFun = temp;
     else
         disp('wrong!');
         sol.info;
         yalmiperror(sol.problem)
     end
     
end % end while
    
end % end function
