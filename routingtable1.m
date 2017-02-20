function [ rate_allocated ] = routingtable1( selectedaction )
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
global N_BSs % the number of SC BSs
global N_Flow % number of flows 

global alpha1 % fraction of data, user 1 divides for subflow1
global alpha2 % fraction of data, user 2 divides for subflow3

max_rate = 30;

% Access to the data
% indicator_bs = zeros(N_Flow,N_BSs+1, N_SubF);
rate_allocated = zeros(N_BSs+1, N_SubF);
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


for player = 1:N_Flow % rate_allocated = zeros(N_BSs+1, N_SubF);
        if player == 1
        UE = 8;
        [~, sizea] = size(RT.action(selectedaction(player)).route1);
        bs_list1 = RT.action(selectedaction(player)).route1;
        for n1 = 1:(sizea-1)
            if n1 == 1 
               rate_allocated(bs_list1(n1),1) = max_rate*alpha1 - rate(bs_list1(n1),bs_list1(n1+1));
            else
               rate_allocated(bs_list1(n1),1) = rate(bs_list1(n1-1),bs_list1(n1)) - rate(bs_list1(n1),bs_list1(n1+1));
            end
        end
        % 
        rate_allocated(bs_list1(sizea),1) = rate(bs_list1(sizea-1),bs_list1(sizea)) - rate(bs_list1(sizea), UE);
        
        [~, sizeb] = size(RT.action(selectedaction(player)).route2);
        bs_list2 = RT.action(selectedaction(player)).route2;
        for n2 = 1:(sizeb-1)
            if n2 == 1 
               rate_allocated(bs_list2(n2),2) = max_rate*(1-alpha1) - rate(bs_list2(n2),bs_list2(n2+1));
            else
               rate_allocated(bs_list2(n2),2) = rate(bs_list2(n2-1),bs_list2(n2)) - rate(bs_list2(n1),bs_list2(n2+1));
            end
        end
        % 
        rate_allocated(bs_list2(sizeb),2) = rate(bs_list2(sizeb-1),bs_list2(sizeb)) - rate(bs_list2(sizeb), UE);
        
        else % if player 2, then update subflow 3 & 4
        UE = 9;
        [~, sizea] = size(RT.action(selectedaction(player)).route1);
        bs_list1 = RT.action(selectedaction(player)).route1;
        for n1 = 1:(sizea-1)
            if n1 == 1 
               rate_allocated(bs_list1(n1),3) = max_rate*alpha2 - rate(bs_list1(n1),bs_list1(n1+1));
            else
               rate_allocated(bs_list1(n1),3) = rate(bs_list1(n1-1),bs_list1(n1)) - rate(bs_list1(n1),bs_list1(n1+1));
            end
        end
        % 
        rate_allocated(bs_list1(sizea),3) = rate(bs_list1(sizea-1),bs_list1(sizea)) - rate(bs_list1(sizea), UE);
        
        [~, sizeb] = size(RT.action(selectedaction(player)).route2);
        bs_list2 = RT.action(selectedaction(player)).route2;
        for n2 = 1:(sizeb-1)
            if n2 == 1 
               rate_allocated(bs_list2(n2),4) = max_rate*(1-alpha2) - rate(bs_list2(n2),bs_list2(n2+1));
            else
               rate_allocated(bs_list2(n2),4) = rate(bs_list2(n2-1),bs_list2(n2)) - rate(bs_list2(n1),bs_list2(n2+1));
            end
        end
        % 
        rate_allocated(bs_list2(sizeb),4) = rate(bs_list2(sizeb-1),bs_list2(sizeb)) - rate(bs_list2(sizeb), UE);
        end   
end






end

