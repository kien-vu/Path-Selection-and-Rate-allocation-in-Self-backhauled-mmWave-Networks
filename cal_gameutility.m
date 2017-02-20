function [ gameutility_observe ] = cal_gameutility( networkqueue )
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
global N_BSs % the number of SC BSs
global N_Flow % number of flows 
global N_Actions % The number of actions : Each user is assumed to choose two paths among 4 paths, then 

% Access to the data
% indicator_bs = zeros(N_Flow,N_BSs+1, N_SubF);
% 
pre_queueBS = ones(N_BSs+1,1);
pre_servedrate = 1000*ones(N_Flow,N_SubF);
pre_queueRoute = ones(N_SubF,1);

% 
gameutility_observe = ones(N_Flow,N_Actions);

gameutility_by_route = ones(N_Flow,N_SubF); % the gain for each route

% 4 main routes
Route(1).R = [1 2 3];
Route(2).R = [1 4 5];
Route(3).R = [1 7];
Route(4).R = [1 6];



%%
for n = 1:N_BSs+1
    pre_queueBS(n,1) = sum(networkqueue(:,n));
end
for r = 1:N_SubF
    [~, sizea] = size(Route(r).R );
    bs_list1 = Route(r).R;
    for n1 = 1:sizea
    	pre_queueRoute(r,1) = max(pre_queueRoute(r,1), pre_queueBS(bs_list1(n1),1));
    end
end

for player = 1:N_Flow
    for r = 1:N_SubF
        if player == 1
            UE = 8;
            [~, sizea] = size(Route(r).R );
            bs_list1 = Route(r).R;
            for n1 = 1:(sizea-1)
                if n1 == 1 
                   pre_servedrate(player,r) = min(pre_servedrate(player,r), rate(bs_list1(n1),bs_list1(n1+1)));
                else
                   pre_servedrate(player,r) = min(pre_servedrate(player,r), rate(bs_list1(n1),bs_list1(n1+1)));
                end
            end
            %
            pre_servedrate(player,r) = min(pre_servedrate(player,r), rate(bs_list1(sizea), UE));
        
        else % if player 2, 
            UE = 9;         
            [~, sizea] = size(Route(r).R );
            bs_list1 = Route(r).R;
            for n1 = 1:(sizea-1)
                if n1 == 1 
                   pre_servedrate(player,r) = min(pre_servedrate(player,r), rate(bs_list1(n1),bs_list1(n1+1)));
                else
                   pre_servedrate(player,r) = min(pre_servedrate(player,r), rate(bs_list1(n1),bs_list1(n1+1)));
                end
            end
            %
            pre_servedrate(player,r) = min(pre_servedrate(player,r), rate(bs_list1(sizea), UE));
        end
    end
end


%%

for player = 1:N_Flow
        gameutility_by_route(player,:) = ...
            [ (pre_servedrate(player,1)/pre_queueRoute(1))
              (pre_servedrate(player,2)/pre_queueRoute(2))
              (pre_servedrate(player,3)/pre_queueRoute(3))
              (pre_servedrate(player,4)/pre_queueRoute(4)) ]';
end


for player = 1:N_Flow
    if player == 1
            gameutility_observe(player,:) = ...
            [   gameutility_by_route(player,1) + gameutility_by_route(player,2) 
                gameutility_by_route(player,1) + gameutility_by_route(player,3) 
                gameutility_by_route(player,1) + gameutility_by_route(player,4) 
                gameutility_by_route(player,2) + gameutility_by_route(player,3) 
                gameutility_by_route(player,2) + gameutility_by_route(player,4) 
                gameutility_by_route(player,3) + gameutility_by_route(player,4) ]';
    else
            gameutility_observe(player,:) = ...
            [   gameutility_by_route(player,1) + gameutility_by_route(player,2) 
                gameutility_by_route(player,1) + gameutility_by_route(player,3) 
                gameutility_by_route(player,1) + gameutility_by_route(player,4) 
                gameutility_by_route(player,2) + gameutility_by_route(player,3) 
                gameutility_by_route(player,2) + gameutility_by_route(player,4) 
                gameutility_by_route(player,3) + gameutility_by_route(player,4) ]';
    end
end





end