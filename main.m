
%%
% Authors: Trung Kien Vu, 
% Date Update on 2015 May 06
% Date Update on 2015 Set 11
% MATLAB version: 8.4 (R2014b) License host:
%   Username: vkien
%   Host ID: 842b2bc09a54 (eth3)
% CVX Version 2.1, Build 1103 (9714d49)


clc;
clear all;
close all;
rng('shuffle'); %  the random number generator based on the current time so that rand, randi, and randn produce a different sequence of numbers after each time you call rng.
random_seed = rng;
rng(random_seed);
% load the parameter
parameter;

global Iterations
global N_Flow % number of flows 
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
global N_BSs % the number of SC BSs
global N_Actions % The number of actions : Each user is assumed to choose two paths among 4 paths, then 
global Congested_BSs
Iterations = Iterations + 1; 
Congested_Route = zeros(N_Flow,N_SubF,N_Actions,Iterations); % The congested level of each route
Congested_BSs   = rand(N_BSs+N_Flow+1,Iterations); % The congested level of each BS
% 1./Congested_BSs
% Learning parameters
Utility_Estimate   = zeros(N_Flow,N_Actions,Iterations); % The utility function
Utility  = zeros(N_Flow,Iterations); % The utility function

Regret          = zeros(N_Flow,N_Actions,Iterations); % The regret
Probability     = 1/N_Actions*ones(N_Flow,N_Actions,Iterations); % The regret
Utility_Observe   = zeros(N_Flow,N_Actions,Iterations); % The utility function
BG_probability   = 1/N_Actions*ones(N_Flow,N_Actions,Iterations); % The utility function
selected_action = zeros(N_Flow,Iterations);


% learning rate
kappa = [20 5 100 200]; % Boltzmann temperature
X = [1 2 3 4 5 6]; % Action index
Congested_BSs1   = rand(N_BSs+N_Flow+1); % The congested level of each BS
Utility_Estimate(:,:,1) = calculate_utility( Congested_BSs1);
calculate_utility( Congested_BSs1)
% Probability(1,:,1) = [0.5 0.1 0.1 0.1 0.1 0.1];
% Probability(1,:,1) = [0.1 0.1 0.1 0.1 0.5 0.1];
same_action = 0;
for iter = 2:Iterations
     Congested_BSs1 = rand(N_BSs+N_Flow+1); % The congested level of each BS
    % The player selects the action based on the probability distribution
    % from previous state t-1
     selected_action(1,iter) = get_x_from_pmf(X,Probability(1,:,iter-1));
     selected_action(2,iter) = get_x_from_pmf(X,Probability(2,:,iter-1));
    % Receive the feeback if two players play the same action or not
    if selected_action(1,iter)  ~= selected_action(2,iter)         
         Utility_Observe(:,:,iter) = calculate_utility( Congested_BSs1); % Congested_BSs(:,iter)
    else % should punnish when playing same action
        same_action = same_action + 1;
%         fprintf('\n Action is played by both players is %d - %d\n', selected_action(1,iter) , selected_action(2,iter) )
        Utility_Observe(:,:,iter) = (1/10)*calculate_utility( Congested_BSs1);
    end
    for player = 1:N_Flow
        for action = 1:N_Actions    
            % Check actions and Update rules for learning procedure
            if  player == 1 
                if action == selected_action(1,iter)
                    Utility_Estimate(player,action,iter) = Utility_Estimate(player,action,iter-1) + 1/(1+iter).^(0.5) ...
                        *(Utility_Observe(player,action,iter) - Utility_Estimate(player,action,iter-1));
                    Utility(player,iter) = Utility_Estimate(player,action,iter);
                else
                    Utility_Estimate(player,action,iter) = Utility_Estimate(player,action,iter-1);
                end
            else
                if action == selected_action(2,iter)
                    Utility_Estimate(player,action,iter) = Utility_Estimate(player,action,iter-1) + 1/(1+iter).^(0.5) ...
                        *(Utility_Observe(player,action,iter) - Utility_Estimate(player,action,iter-1));
                    Utility(player,iter) = Utility_Estimate(player,action,iter);
                else
                    Utility_Estimate(player,action,iter) = Utility_Estimate(player,action,iter-1);
                end
            end
            
            % Regret
            Regret(player,action,iter) = max(0,Regret(player,action,iter-1) + 1/(1+iter).^(0.55) ...
                *(Utility_Observe(player,action,iter) - Utility_Estimate(player,action,iter) - Regret(player,action,iter-1)));
            
        end
        for action = 1:N_Actions
            % Calculate the BG probability distributions
            BG_probability(player,action,iter) = exp(1/kappa(1)*max(Regret(player,action,iter),0))/ ...
                                                    sum(exp(1/kappa(1)*max(Regret(player,:,iter),0))); % do I need to extract the current action
            % Probability 
            Probability(player,action,iter) = Probability(player,action,iter-1) + 1/(1+iter).^(0.6) ...
                *(BG_probability(player,action,iter) - Probability(player,action,iter-1));
        end
        %% end update rules check if the probability distribution for each action does not change, then convegence.
        
    end
end

fprintf('Number of times both players playing the same action is %d ~ %f%%\n', same_action, 100*same_action/Iterations)

disp('Finnished');
    for player = 1:N_Flow
        figure
        for action = 1:N_Actions  
            cdfplot(BG_probability(player,action,:)); hold on;
            %ecdf(BG_probability(player,action,:)); hold on;
        end
        ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
        xlabel('Probability that each player plays for each action','FontSize',12.6,'FontName','Times New Roman')
        
%         figure
%         for action = 1:N_Actions  
%             cdfplot(Regret(player,action,:)); hold on;
%             %ecdf(BG_probability(player,action,:)); hold on;
%         end
%         ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
%         xlabel('Probability of Regret value','FontSize',12.6,'FontName','Times New Roman')
        
    end

        [f,xi] = ksdensity(selected_action(1,:));
        figure
        plot(xi,f);
        ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
        xlabel('Probability of action 1 value','FontSize',12.6,'FontName','Times New Roman')       
        [f,xi] = ksdensity(selected_action(1,:));
        figure
        plot(xi,f);
        ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
        xlabel('Probability of action 2 value','FontSize',12.6,'FontName','Times New Roman')
%     figure
%     hist(X,sort(action1)); % qqplot normplot ksdensity
%     figure
%     hist(X,sort(action2));

% serie = rand;
% file = strcat('2017Feb_','rand_',num2str(serie),'_','On_',num2str(date),'.mat');
% save(file);
% 
% % figure
% % x = normrnd(5,1,100,1);
% % y = wblrnd(2,0.5,100,1);
% % qqplot(x,y);



    
    