function [indicator_bs] = routingtable(selectedaction)
% This function is going to show the routing table for each route: Based on
% selected action and the topology (given)
% That helps to flood the traffic
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
global N_BSs % the number of SC BSs
global N_Flow % number of flows or players


% Access to the data
indicator_bs = zeros(N_BSs+1,N_SubF);
% selectedaction = [1 2]';
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


for player = 1:N_Flow
    [~, sizea] = size(RT.action(selectedaction(player)).route1); % 
    bs_list1 = RT.action(selectedaction(player)).route1;
    [~, sizeb] = size(RT.action(selectedaction(player)).route2);
    bs_list2 = RT.action(selectedaction(player)).route2;
    if player == 1
        for n1 = 1:sizea
               indicator_bs(bs_list1(n1),1) = 1; % subflow 1
        end

        for n2 = 1:sizeb
               indicator_bs(bs_list2(n2),2) = 1; % subflow 2
        end
    else
        for n1 = 1:sizea
               indicator_bs(bs_list1(n1),3) = 1; % subflow 1
        end
        
        for n2 = 1:sizeb
               indicator_bs(bs_list2(n2),4) = 1; % subflow 2
        end
    end
    % Last case for MBS
    indicator_bs(1,:) = [1 1 1 1]; % traffic always goes through the MBS

end


end

