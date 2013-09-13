% Shortcut to RBoltzmann testing & graphing
% Apply after training weights

   % After weights are trained, try this. 
   
   % Initialize with one step of activity
   prevStep = placeCells(2,:);
   
   % Cycle placeCell timesteps
   for j = 1:Npats
    % prevStep = placeCells(j,:);
        
    %% Apply RBoltzmannMachine
        
    
    %% Apply temporal RBoltzmann connectivity (hidden to hidden)    
        [tweights2, ~, hidden_real] = RBoltzmann_real(tweights, prevStep, Nhidden, Nhidden);
        % discard tweights2 (no learning!)
        prevStep = hidden_real;
        
        pathRecon(j,1:Nhidden) = hidden_real;
        pathRecon(j,Nhidden+1) = patterns(j,Nvisible+1); % coordinates
        pathRecon(j,Nhidden+2) = patterns(j,Nvisible+2);
        
        disp(j)
        
   end % End cycling patterns
    
   
%    for t = 0:9 % Cycle each place cell (test)
%      figure
%      % t = 0
%       for q = 1:10
%         %% Graphing %%
%         spikeskip = 1; % No skipping!
%         subplot(3,4,q) % Plot each cell individually
%     
%         % plot(placeCells(1:length(placeCells), Nhidden+1),placeCells(1:length(placeCells), Nhidden+2),'.','MarkerSize',8);
%         hold on; plot(placeCells(spikeskip*find(placeCells(1:spikeskip:end, (q+(10*t)))>0.05), Nhidden+1),placeCells(spikeskip*find(placeCells(1:spikeskip:end, (q+(10*t)))>0.05), Nhidden+2),'R.','MarkerSize',8);
%         title('Abstract model','fontsize',12)
%         xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12) 
%         set(gca,'xlim',[min(placeCells(:,Nhidden+1)) max(placeCells(:,Nhidden+1))],'ylim',[min(placeCells(:,Nhidden+2)) max(placeCells(:,Nhidden+2))])
%       end
%    end
   
figure 
for t = 1:99 % Cycle each place cell (test)
     % t = 0
        %% Graphing %%
        spikeskip = 1; % No skipping!
    
        plot(placeCells(1:end, Nhidden+1),placeCells(1:end, Nhidden+2),'B.','MarkerSize',4);
        hold on; plot(placeCells(spikeskip*find(placeCells(1:spikeskip:end, (t))>0.05), Nhidden+1),placeCells(spikeskip*find(placeCells(1:spikeskip:end, (t))>0.05), Nhidden+2),'R.','MarkerSize',8);
        title('Desired Path','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12) 
        set(gca,'xlim',[min(placeCells(:,Nhidden+1)) max(placeCells(:,Nhidden+1))],'ylim',[min(placeCells(:,Nhidden+2)) max(placeCells(:,Nhidden+2))])
end

figure
for t = 1:99 % Cycle pathReconstruction
     % t = 0
        %% Graphing %%
        spikeskip = 1; % No skipping!
    
        % plot(placeCells(1:length(placeCells), Nhidden+1),placeCells(1:length(placeCells), Nhidden+2),'.','MarkerSize',8);
        hold on; plot(pathRecon(spikeskip*find(pathRecon(1:spikeskip:end, (t))>0.05), Nhidden+1),pathRecon(spikeskip*find(pathRecon(1:spikeskip:end, (t))>0.05), Nhidden+2),'R.','MarkerSize',8);
        title('Recreated Path','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12) 
        set(gca,'xlim',[min(placeCells(:,Nhidden+1)) max(placeCells(:,Nhidden+1))],'ylim',[min(placeCells(:,Nhidden+2)) max(placeCells(:,Nhidden+2))])
      
end
