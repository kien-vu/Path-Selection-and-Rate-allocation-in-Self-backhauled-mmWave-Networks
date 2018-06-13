% function [ served_rate ] = rate_allocation( ) %vq, selectedaction, out_data, t
clc
global alpha1 % fraction of data, user 1 divides for subflow1
global alpha2 % fraction of data, user 2 divides for subflow3
global N_Flow % number of flows 
N_Flow = 2;
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
N_SubF = 4;
global N_BSs % the number of SC BSs
N_BSs = 6;
global N_Actions % The number of actions : Each user is assumed to choose two paths among 4 paths, then 
% the number of actions is the combination of choosing 2 from 4.
N_Actions = 6;

selectedaction = [1 2]';
vq = 10*rand( N_SubF, 1 + N_BSs);
[nSF, nBS] = size(vq); % number of subflows and BSs
served_rate = zeros(nSF, nBS);
Delay_const = zeros(nBS, nSF);
indicator_BSs = routingtable(selectedaction);
%     for nbs = 1:nBS
%         for nsf = 1:nSF
%             Delay_const(nbs, nsf) = 0.5 * (t-2*0.05) - out_data(nsf, nbs);
%         end
%     end
%     for bs = 1:nBS
%         served_rate(1, bs) = indicator_BSs(bs,1) *  max( vq(1, bs) +  (alpha1   ), out_data(1, bs));
%         served_rate(2, bs) = indicator_BSs(bs,2) *  max( vq(2, bs) + ((1-alpha1)), out_data(2, bs));
%         served_rate(3, bs) = indicator_BSs(bs,3) *  max( vq(3, bs) +  (alpha2   ), out_data(3, bs));
%         served_rate(4, bs) = indicator_BSs(bs,4) *  max( vq(4, bs) +  ((1-alpha2)), out_data(4, bs));
%     end

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
%     
% LoadingSolvers();
  ops = sdpsettings('solver','mosek'); % set the interal solver to be mosek
%   ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0,'debug',1,'beeponproblem',1,'saveduals',1); % set the interal solver to be mosek
  lowerbound = rand(nBS, nSF);
  b = 5*rand(nBS, nSF) +1;
  
  u_next = ones(nBS, nSF);
  tx_next = ones(nBS, nSF);
    
while(1)
% Define variable
%     rate        = sdpvar(nBS, nSF);
    tx          = sdpvar(nBS, nSF);
    x_max       = sdpvar(nBS, nSF);
    % Nonconvex variables
    u = sdpvar(nBS, nSF);

% Define constraints and Objectives
    constraints = []; % contain all the constraints
    constraints = [constraints, x_max >= 0];
    constraints = [constraints, tx >= 0];
    constraints = [constraints, u >= 0];

    for nbs = 1:nBS
        for nsf = 1:nSF
            if indicator_BSs(nbs,nsf) == 0
                constraints = [constraints, tx(nbs, nsf) == 0];
                constraints = [constraints, x_max(nbs, nsf) == 0];
            end
            % Delay constraint 
            constraints = [constraints, x_max(nbs, nsf) >= lowerbound(nbs, nsf)];
            constraints = [constraints, x_max(nbs, nsf) <= b(nbs, nsf)];
        end
        % Transmit power constraints
        constraints = [constraints, sum(tx(nbs, :)) <= 10]; % maximum power constraint
        
    end

    % Here is f***ing difficult - 
for player = 1:N_Flow % rate_allocated = zeros(N_BSs+1, N_SubF);
        if player == 1
            UE = 8;
            [~, sizea] = size(RT.action(selectedaction(player)).route1);
            bs_list1 = RT.action(selectedaction(player)).route1;
            for n1 = 1:sizea
                if n1 == 1 % if MBS
                    constraints = [constraints, (1+ tx(bs_list1(n1), 1)) <= x_max(bs_list1(n1), 1) ];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list1(n1-1), 1); tx(bs_list1(n1-1), 1)./2], (2+ tx(bs_list1(n1-1), 1))./2  ) ];
                    constraints = [constraints, 2 * u(bs_list1(n1-1))*u_next(bs_list1(n1-1), 1)./(1 + tx_next(bs_list1(n1), 1))  - (1+ tx(bs_list1(n1), 1)) * (u_next(bs_list1(n1-1), 1)./(1+tx_next(bs_list1(n1), 1))).^2 >= exp(lowerbound(bs_list1(n1), 1)) ];
                end
            end
            % 20170210 stop here
                      
            [~, sizeb] = size(RT.action(selectedaction(player)).route2);
            bs_list2 = RT.action(selectedaction(player)).route2;
            for n1 = 1:sizeb
                if n1 == 1 % if MBS
                    constraints = [constraints, (1+ tx(bs_list2(n1), 2)) <= exp(5+x_max(bs_list2(n1), 2))  ];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list2(n1-1), 2); tx(bs_list2(n1-1), 2)./2], (2+ tx(bs_list2(n1-1), 2))./2) ];
                    constraints = [constraints, 2* u(bs_list2(n1-1)) * u_next(bs_list2(n1-1), 2)./(1 + tx_next(bs_list2(n1), 2))  - (1+ tx(bs_list2(n1), 2)) * (u_next(bs_list2(n1-1), 2)./(1+tx_next(bs_list2(n1), 2))).^2 >= exp(x_max(bs_list2(n1), 2)) ];
                end
            end
            
            else % if player 2, then update subflow 3 & 4
            UE = 9;
            [~, sizea] = size(RT.action(selectedaction(player)).route1);
            bs_list1 = RT.action(selectedaction(player)).route1;
            for n1 = 1:sizea
                if n1 == 1 % if MBS
                    constraints = [constraints, (1+ tx(bs_list1(n1), 3)) <= x_max(bs_list1(n1), 3) ];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list1(n1-1), 3); tx(bs_list1(n1-1), 3)./2], (2+ tx(bs_list1(n1-1), 3))./2) ];
                    constraints = [constraints, 2* u(bs_list1(n1-1)) * u_next(bs_list1(n1-1), 3)./(1 + tx_next(bs_list1(n1), 3))  - (1+ tx(bs_list1(n1), 3)) * (u_next(bs_list1(n1-1), 3)./(1+tx_next(bs_list1(n1), 3))).^2 >= exp(x_max(bs_list1(n1), 3)) ];
                end
            end
            % 20170210 stop here
                      
            [~, sizeb] = size(RT.action(selectedaction(player)).route2);
            bs_list2 = RT.action(selectedaction(player)).route2;
            for n1 = 1:sizeb
                if n1 == 1 % if MBS
                    constraints = [constraints, (1+ tx(bs_list2(n1), 4)) <= x_max(bs_list2(n1), 4) ];
                else
                    % Cone constraint for SC BSs, and convex approximation
                    constraints = [constraints, cone( [ u(bs_list2(n1-1), 4); tx(bs_list2(n1-1), 4)./2], (2+ tx(bs_list2(n1-1), 4))./2 ) ];
                    constraints = [constraints, 2* u(bs_list2(n1-1)) * u_next(bs_list2(n1-1), 4)./(1 + tx_next(bs_list2(n1), 4))  - (1+ tx(bs_list2(n1), 4)) * (u_next(bs_list2(n1-1), 4)./(1+tx_next(bs_list2(n1), 4))).^2 >= exp(x_max(bs_list2(n1), 4)) ];
                end
            end
        end   
end


    % Objective function
     obj =  - vq(1,:) * x_max(:,1) - vq(2,:) * x_max(:,2) ...
            - vq(3,:) * x_max(:,3) - vq(4,:) * x_max(:,4);
% Solve the problem
     sol = optimize(constraints, obj, ops); % solve the problem optimize replaced sdpsolve
     x_max = value(x_max);
     tx = value(tx);
     u  = value(u);
% Check the results
     sol.info;
     if sol.problem == 0
         tx_next  = tx;
         u_next = u;
         temp = - value(obj);
         if (temp - ObjFun) <= 1e-2 
            served_rate = x_max;
            break;
         end
     else
         disp('PA, something went wrong!');
         sol.info;
         yalmiperror(sol.problem)
     end
     
end % end while
    

    
% end

