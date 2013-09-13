% Shortcut to RBoltzmann testing & graphing
% Apply after training weights

   % After weights are trained, try this. 
    
   % Cycle patterns
   for j = 1:Npats
        input = patterns(j,:); % Select one of the patterns
    
        
    %% Apply RBoltzmannMachine
        [weights2, partialcost1, hidden] = RBoltzmann( weights1, input, Nvisible, Nhidden);
        % discard weights2 (no learning!)
        
        placeCells(j,1:Nhidden) = hidden;
        placeCells(j,Nhidden+1) = patterns(j,Nvisible+1); % coordinates
        placeCells(j,Nhidden+2) = patterns(j,Nvisible+2);
        
        disp(j)
        
   end % End cycling patterns
    
   for t = 1:10 % Cycle each place cell
    %% Graphing %%
    spikeskip = 1;
    subplot(3,4,t) % Plot each cell individually
    
    plot(placeCells(1:length(placeCells), Nvisible+1),placeCells(1:length(placeCells), Nvisible+2),'.','MarkerSize',8);
    hold on; plot(placeCells(spikeskip*find(placeCells(1:spikeskip:end, t)>0), Nvisible+1),placeCells(spikeskip*find(placeCells(1:spikeskip:end, t)>0), Nvisible+2),'R.','MarkerSize',8);
    title('Abstract model','fontsize',12)
    xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12) 
    set(gca,'xlim',[min(placeCells(:,Nvisible+1)) max(placeCells(:,Nvisible+1))],'ylim',[min(placeCells(:,Nvisible+2)) max(placeCells(:,Nvisible+2))])
   end
