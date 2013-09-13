%% This script plots some of the place cell reconstructions. It's currently set to plot 10 of them on a 3 x 4 grid. You have to modify numbers if you want to plot different place cells outputs.    
    for t = 0:9 % Cycle each pathRecon (test)
     figure
     % t = 0
      for q = 1:10
        %% Graphing %%
        spikeskip = 1; % No skipping!
        subplot(3,4,q) % Plot each cell individually
    
        % plot(placeCells(1:length(placeCells), Nhidden+1),placeCells(1:length(placeCells), Nhidden+2),'.','MarkerSize',8);
        hold on; plot(pathRecon(spikeskip*find(pathRecon(1:spikeskip:end, (q+(10*t)))>0.40), Nhidden+1),pathRecon(spikeskip*find(pathRecon(1:spikeskip:end, (q+(10*t)))>0.40), Nhidden+2),'R.','MarkerSize',8);
        title('Abstract model','fontsize',12)
        xlabel('Position (m)','fontsize',12), ylabel('Position (m)','fontsize',12) 
        set(gca,'xlim',[min(placeCells(:,Nhidden+1)) max(placeCells(:,Nhidden+1))],'ylim',[min(placeCells(:,Nhidden+2)) max(placeCells(:,Nhidden+2))])
      end
   end