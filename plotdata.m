
% plot the probability that each action is taken
     for player = 1:N_Flow
        figure
        h = histogram(sort(selected_action(player,:)),'Normalization','probability'); % qqplot normplot ksdensity
        ylabel('Probability','FontSize',12.6,'FontName','Times New Roman')
        xlabel('Actions','FontSize',12.6,'FontName','Times New Roman')
     end
     
     
     
     

     
     


%    for player = 1:N_Flow
%         figure
%         for action = 1:N_Actions  
%             cdfplot(BG_probability(player,action,:)); hold on;
%             %ecdf(BG_probability(player,action,:)); hold on;
%         end
%         ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
%         xlabel('Probability that each player plays for each action','FontSize',12.6,'FontName','Times New Roman')
        
%         figure
%         for action = 1:N_Actions  
%             cdfplot(Regret(player,action,:)); hold on;
%             %ecdf(BG_probability(player,action,:)); hold on;
%         end
%         ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
%         xlabel('Probability of Regret value','FontSize',12.6,'FontName','Times New Roman')
        
%     end

%         [f,xi] = ksdensity(selected_action(1,:));
%         figure
%         plot(xi,f);
%         ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
%         xlabel('Probability of action 1 value','FontSize',12.6,'FontName','Times New Roman')       
%         [f,xi] = ksdensity(selected_action(1,:));
%         figure
%         plot(xi,f);
%         ylabel('CDF','FontSize',12.6,'FontName','Times New Roman')
%         xlabel('Probability of action 2 value','FontSize',12.6,'FontName','Times New Roman')
        
  

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



    
    



