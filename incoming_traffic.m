function [ incoming_rate ] = incoming_traffic( selectedaction, tx )
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
global N_BSs % the number of SC BSs
global N_Flow % number of flows 
global alpha1

% Access to the data
% indicator_bs = zeros(N_Flow,N_BSs+1, N_SubF);
incoming_rate = zeros(N_BSs+1, N_SubF);
% selectedaction = [1 2]';
% Routing table
% [distance_cal, h]  = topology();
% channel_gain( distance_cal(a,b), h(a,b) )
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
            k = 1;
            if n1 == 1 
               incoming_rate(bs_list1(n1),1) = alpha1;
            else
               incoming_rate(bs_list1(n1),1) =  log( 1 + tx(bs_list1(n1-1) , k));
            end
        end
        % 
        incoming_rate(bs_list1(sizea),1) =  log( 1 + tx(bs_list1(sizea-1) , k));
        
        [~, sizeb] = size(RT.action(selectedaction(player)).route2);
        bs_list2 = RT.action(selectedaction(player)).route2;
        k = 2;
        for n2 = 1:(sizeb-1)
            if n2 == 1 
               incoming_rate(bs_list2(n2),2) = alpha1;
            else
               incoming_rate(bs_list2(n2),2) = log( 1 + tx(bs_list2(n2-1) , k));
            end
        end
        % 
        incoming_rate(bs_list2(sizeb),2) = log( 1 + tx(bs_list2(sizeb-1) , k)); %rate(bs_list2(sizeb-1),bs_list2(sizeb)) - rate(bs_list2(sizeb), UE);
        
        else % if player 2, then update subflow 3 & 4
        UE = 9;
        [~, sizea] = size(RT.action(selectedaction(player)).route1);
        bs_list1 = RT.action(selectedaction(player)).route1;
        k = 3;
        for n1 = 1:(sizea-1)
            if n1 == 1 
               incoming_rate(bs_list1(n1),3) = alpha1;
            else
               incoming_rate(bs_list1(n1),3) = log( 1 + tx(bs_list1(n1-1) , k));
            end
        end
        % 
        incoming_rate(bs_list1(sizea),3) = log( 1 + tx(bs_list1(sizea-1) , k));
        
        [~, sizeb] = size(RT.action(selectedaction(player)).route2);
        bs_list2 = RT.action(selectedaction(player)).route2;
        k = 4;
        for n2 = 1:(sizeb-1)
            if n2 == 1 
               incoming_rate(bs_list2(n2),4) = alpha1;
            else
               incoming_rate(bs_list2(n2),4) = log( 1 + tx(bs_list2(n2-1) , k));
            end
        end
        % 
        incoming_rate(bs_list2(sizeb),4) = log( 1 + tx(bs_list2(sizeb-1) , k));
        end   
end
incoming_rate = incoming_rate';
end

