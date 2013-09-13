% Shortcut to RBoltzmann testing & graphing
% Apply after training weights

   % After weights are trained, try this. 
    
   % Cycle patterns
   for j = 1:Npats
        input = patterns(j,:); % Select one of the patterns
    
        
    %% Apply RBoltzmannMachine
       [weights2, tweights2, partialcost1, hidden_real] = RBoltzmann_temporal_real( weights1, tweights, input, prevStep, Nvisible, Nhidden);

      % [weights2, partialcost1, hidden_real] = RBoltzmann_real( weights1, input, Nvisible, Nhidden);
        % discard weights2 (no learning!)
        
        placeCells(j,1:Nhidden) = hidden_real;
        placeCells(j,Nhidden+1) = patterns(j,Nvisible+1); % coordinates
        placeCells(j,Nhidden+2) = patterns(j,Nvisible+2);
        
        disp(j)
        
   end % End cycling patterns
    
   for t = 0:9 % Cycle each place cell
      figure
     % t = 0
      for q = 1:10
        %% Graphing %%
        spikeskip = 1;
        subplot(3,4,q) % Plot each cell individually
    
        plot(placeCells(1:length(placeCells), Nhidden+1),placeCells(1:length(placeCells), Nhidden+2),'.','MarkerSize',8);
        hold on; plot(placeCells(spikeskip*find(placeCells(1:spikeskip:end, (q+(10*t)))>0.05), Nhidden+1),placeCells(spikeskip*find(placeCells(1:spikeskip:end, (q+(10*t)))>0.05), Nhidden+2),'R.','MarkerSize',8);
        title('Abstract model','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12) 
        set(gca,'xlim',[min(placeCells(:,Nhidden+1)) max(placeCells(:,Nhidden+1))],'ylim',[min(placeCells(:,Nhidden+2)) max(placeCells(:,Nhidden+2))])
      end
   end
