
%%
% Authors: Trung Kien Vu
%    Date Update on 2017 Feb 14
%    MATLAB version: 9.1 (R2016b)
%    OS: Windows 7 amd64 version 6.1
%    Java version: 1.7.0_60
%    Yalmip lastest version
%    Name     Status             Version    Location
% ------------------------------------------------------------
%    Mosek    selected,default   7.1.0.12   {cvx}\mosek\w64
%    SDPT3                       4.0        {cvx}\sdpt3
%    SeDuMi                      1.34       {cvx}\sedumi
%%

clc;
%  clear all;
%  close all;
% rng('shuffle'); %  the random number generator based on the current time so that rand, randi, and randn produce a different sequence of numbers after each time you call rng.
% random_seed = rng;
rng(random_seed);
% load the parameter
parameter;

global Iterations
global N_Flow % number of flows 
global N_SubF % number of each subflows for each flow, assume that each flow is divided equally
global N_BSs % the number of SC BSs
global N_Actions % The number of actions : Each user is assumed to choose two paths among 4 paths, then 
global alpha1 % fraction of data, user 1 divides for subflow1
global alpha2 % fraction of data, user 1 divides for subflow1
Iterations = Iterations + 1; 

%% Learning parameters
Learning_Utility_Estimate   = zeros(N_Flow,N_Actions,Iterations); % The utility function
Utility                     = zeros(N_Flow,Iterations); % The utility function
Regret                      = zeros(N_Flow,N_Actions,Iterations); % The regret
Probability                 = 1/N_Actions*ones(N_Flow,N_Actions,Iterations); % The regret
Learning_Utility_Observe    = zeros(N_Flow,N_Actions,Iterations); % The utility function
BG_probability              = 1/N_Actions*ones(N_Flow,N_Actions,Iterations); % The utility function
selected_action             = zeros(N_Flow,Iterations);
% classical learning
class_Probability                 = 1/N_Actions*ones(N_Flow,N_Actions,Iterations); % The classical regret learning
class_selected_action             = zeros(N_Flow,Iterations);

% learning rate and initial setting
kappa = [2 5 100 200]; % Boltzmann temperature
X = [1 2 3 4 5 6]; % Action index
%% Lyapunov Parameters
network_queue   = zeros( N_SubF, 1 + N_BSs, Iterations); % Number of base stations, with number of subflows (4)
virtual_queue   = 1/Iterations*ones( N_SubF, 1 + N_BSs, Iterations); % Number of base stations, with number of subflows (4)
serving_rate    = zeros( N_SubF, 1 + N_BSs, Iterations); % serving rate
remaining_data  = zeros( N_SubF, 1 + N_BSs, Iterations); % reamining data
incoming_data   = zeros( N_SubF, 1 + N_BSs, Iterations); % incomming data
auxiliary_var   = zeros( N_SubF, 1 + N_BSs, Iterations); % auxililary variables
indicator_bs    = zeros(N_BSs+1,N_SubF,Iterations);      % which BS turns on
delay_require   = zeros( N_SubF, 1 + N_BSs, Iterations); % Data required for delay
transmit_power  = zeros( 1 + N_BSs, N_SubF, Iterations); % transmit power
% Optimization setup
beta1 = 3; % lets say we can guarantee up to 3
epsilon1 = 0.05; % Probabilistic delay constraint
%% Simulation starts
% Learning part runs for a long term period, let says after 10 slots, while
% Rate allocation runs for short term period, at each time slot.
same_action = 0;
for iter = 2:Iterations
    % Note to perform learning part in long term period, a frame, let say 10 slots,
    % while the rate allocation is done at each time slot basic, subframe.
    if iter == 2 || mod(iter-2,10) == 0
     if iter == 2 % At the beginning assume that the network is not congested at all
        network_queue(:,:,iter-1) = ones(N_SubF, 1 + N_BSs); % The congested level of each BS
     end
    % The player selects the action based on the probability distribution
    % from previous state t-1
     selected_action(1,iter) = get_x_from_pmf(X,Probability(1,:,iter-1));
     selected_action(2,iter) = get_x_from_pmf(X,Probability(2,:,iter-1));
    % Receive the feeback if two players play the same action or not, now
    % just check number of playing same actions
    if check_action(selected_action(1,iter), selected_action(2,iter)) == 1       
         Learning_Utility_Observe(:,:,iter) = cal_gameutility( network_queue(:,:,iter-1)); % Congested_BSs(:,iter)
    else % should punnish when playing same action
        same_action = same_action + 1;
        Learning_Utility_Observe(:,:,iter) = 1/2 * cal_gameutility( network_queue(:,:,iter-1) );
    end
    % Select for the baseline - classical learning
     class_selected_action(1,iter) = get_x_from_pmf(X,class_Probability(1,:,iter-1));
     class_selected_action(2,iter) = get_x_from_pmf(X,class_Probability(2,:,iter-1));
    for player = 1:N_Flow
        for action = 1:N_Actions    
            % Check actions and Update rules for learning procedure
            if  player == 1 
                if action == selected_action(1,iter)
                    Learning_Utility_Estimate(player,action,iter) = Learning_Utility_Estimate(player,action,iter-1) + 1/(1+iter).^(0.5) ...
                        *(Learning_Utility_Observe(player,action,iter) - Learning_Utility_Estimate(player,action,iter-1));
                    Utility(player,iter) = Learning_Utility_Estimate(player,action,iter);
                else
                    Learning_Utility_Estimate(player,action,iter) = Learning_Utility_Estimate(player,action,iter-1);
                end
            else
                if action == selected_action(2,iter)
                    Learning_Utility_Estimate(player,action,iter) = Learning_Utility_Estimate(player,action,iter-1) + 1/(1+iter).^(0.5) ...
                        *(Learning_Utility_Observe(player,action,iter) - Learning_Utility_Estimate(player,action,iter-1));
                    Utility(player,iter) = Learning_Utility_Estimate(player,action,iter);
                else
                    Learning_Utility_Estimate(player,action,iter) = Learning_Utility_Estimate(player,action,iter-1);
                end
            end  
            % Regret update
            Regret(player,action,iter) = max(0,Regret(player,action,iter-1) + 1/(1+iter).^(0.55) ...
                *(Learning_Utility_Observe(player,action,iter) - Learning_Utility_Estimate(player,action,iter) - Regret(player,action,iter-1)));
            
        end
        for action = 1:N_Actions
            % Calculate the BG probability distributions
            BG_probability(player,action,iter) = exp(1/kappa(1)*max(Regret(player,action,iter), 0))/ ...
                                                    sum( exp(1/kappa(1)*max(Regret(player,:,iter), 0)) ); % do I need to extract the current action
            % Probability 
            Probability(player,action,iter) = Probability(player,action,iter-1) + 1/(1+iter).^(0.6) ...
                *(BG_probability(player,action,iter) - Probability(player,action,iter-1));
        end
        % end update rules check if the probability distribution for each action does not change, then convegence.
    end
    
% Lyapunov starts from here at begining of a frame
        % Auxiliary Variable Selection
        auxiliary_var(:,:,iter) = auxiliary_var_selection(virtual_queue(:,:,iter));
        indicator_BSs = routingtable (selected_action(:,iter));
        auxiliary_var(:,:,iter) = auxiliary_var(:,:,iter).*indicator_BSs';
        % Calculate the remain traffic by adding the incoming and subtract the outcoming
        % Calculate the delay required
        for nbs = 1:N_BSs+1
            for nsf = 1:N_SubF
                remaining_data(nsf,nbs,iter)  = max(remaining_data(nsf,nbs,iter) - serving_rate(nsf,nbs,iter-1), 0) + incoming_data(nsf,nbs,iter-1); % outcommingdata X
                delay_require(nsf,nbs,iter)   = indicator_BSs(nbs, nsf) * max(0, alpha1*(0 - beta1*epsilon1) + remaining_data(nsf,nbs,iter) ); % this is a simplify constraint
            end
        end
        % Rate Allocation
        [serving_rate(:,:,iter), transmit_power(:,:,iter) ] = rate_allocation_simple(virtual_queue(:,:,iter), ...
            selected_action(:,iter), delay_require(:,:,iter) );
        incoming_data(:,:,iter) =   incoming_traffic(selected_action(:,iter), transmit_power(:,:,iter));
        % Network queue update 
        network_queue(:,:,iter + 1) = network_queue_update(selected_action(:,iter),...
            network_queue(:,:,iter), serving_rate(:,:,iter), incoming_data(:,:,iter)' );
        % Virtual Queue Update 
        virtual_queue(:,:,iter + 1) = virtual_queue_update(virtual_queue(:,:,iter),...
            auxiliary_var(:,:,iter), serving_rate(:,:,iter));
    
    %% Update other mod(iter+2,10) not equal 0
    else % the action from 2 to 11 will be the same, playing the same action as before
    for player = 1:N_Flow
        for action = 1:N_Actions 
            selected_action(1,iter) = selected_action(1,iter-1);
            selected_action(2,iter) = selected_action(2,iter-1);  
            Probability(player,action,iter) = Probability(player,action,iter-1);
            BG_probability(player,action,iter) = BG_probability(player,action,iter-1);
            Regret(player,action,iter) = Regret(player,action,iter-1);
            Learning_Utility_Estimate(player,action,iter) = Learning_Utility_Estimate(player,action,iter-1);
            Utility(player,iter) = Utility(player, iter-1);
            Learning_Utility_Observe(player,action,iter) = Learning_Utility_Observe(player,action,iter-1);
        end
    end
% During the short TTI, 1 ms, there is no incoming data for MBS, while the
% SCs try to push all flow data
        auxiliary_var(:,:,iter) = auxiliary_var_selection(virtual_queue(:,:,iter));
        indicator_BSs = routingtable (selected_action(:,iter));
        auxiliary_var(:,:,iter) = auxiliary_var(:,:,iter).*indicator_BSs';
        % Calculate the remain traffic by adding the incoming and subtract the outcoming
        % Calculate the delay required
        for nbs = 1:N_BSs+1
            for nsf = 1:N_SubF
                remaining_data(nsf,nbs,iter)  = max(remaining_data(nsf,nbs,iter) - serving_rate(nsf,nbs,iter-1), 0) + incoming_data(nsf,nbs,iter-1); % outcommingdata X
                delay_require(nsf,nbs,iter)   = indicator_BSs(nbs, nsf) * max(0, -alpha1*beta1*epsilon1 + remaining_data(nsf,nbs,iter) ); % this is a simplify constraint
            end
        end
        % Rate Allocation
        serving_rate(:,:,iter)  =   serving_rate(:,:,iter-1);
        transmit_power(:,:,iter)=   transmit_power(:,:,iter-1) ; % possible transmit power
        % need to determind any remaining traffic in the queue at previous time slot
        incoming_data(:,:,iter) =   incoming_traffic_0(selected_action(:,iter), transmit_power(:,:,iter), network_queue(:,:,iter-1)); % network_queue(:,:,iter-1)
        % Network queue update, 
        alpha1 = 0; % no data arrive this point
        network_queue(:,:,iter + 1) = network_queue_update(selected_action(:,iter),...
            network_queue(:,:,iter), serving_rate(:,:,iter), incoming_data(:,:,iter)' );
        % Virtual Queue Update 
        virtual_queue(:,:,iter + 1) = virtual_queue_update(virtual_queue(:,:,iter),...
            auxiliary_var(:,:,iter), serving_rate(:,:,iter));
       alpha1 = alpha2 ; % restore
    end
%% baselines
        
end

fprintf('\nNumber of times both players choosing the path with common BS is %d ~ %f%%\n', same_action*9, 100*same_action/Iterations*9)

fprintf('\nFinished\n');

serie = rand;
file = strcat('multihop','No_Iterations_',num2str(Iterations),'_Arrival_',num2str(alpha1),'_','rand',num2str(serie),'_','On',num2str(date),'.mat');
cd Results
save(file);
cd ..

plotdata;





