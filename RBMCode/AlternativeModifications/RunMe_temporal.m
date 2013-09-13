% RBoltzmann script which calls RBoltmann.m and sets the parameters of the
% number of layers and how they interact


% Originally written by Mark Shaw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nvisible is the dimensionality of the visible layer
% Nhidden is the dimensionality of the hidden layer

Nvisible = 100; % 100 grid cell visible units
Nhidden = 99; % 100 place cell hidden units
Ntimestep = 2000;
Npats = 31; % number of points in rat's route

% for skip = 100:100:360000
%     % This is the pattern matrix
%         patterns(skip/100,:) = savedAbstract(skip,:);
%     % This is the input data from the grid cell firing. 
% end

%%% Present each pattern to the network, and calculate the activity
%%% using repeated firing on each other (but not on the self)

% % Setup weights matrix (random between -0.5 and 0.5)
% for i = 1:Nvisible,
%     for j = 1:Nhidden,
%         weights1(i,j) = rand - 0.5;
%     end
% end
% 
%  % Add bias weights
%  for i = 1:Nvisible
%     weights1(i,Nhidden+1) = 0;
%  end
 
 % Setup weights for temporal connectivity (random between -0.5 and 0.5)
for y = 1:Nhidden
    for j = 1:Nhidden
        tweights(y,j) = rand - 0.5;
    end
end

 % Add bias weights
 for y = 1:Nhidden
    tweights(y,Nhidden+1) = 1;
 end

%  % Initialize prevStep
%  for j = 1:Nhidden
%      prevStep(j) = 0;
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Firing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:Ntimestep % Repeat Ntimestep times
   oldcost = 0; % Reset summed cost
   cost = 0; % Reset cost
   disp(t)
   % Cycle patterns
   for j = 1:Npats
        input = placeCells(j,:); % Select a timestep of placeCell activity
        
    
        
    %% Apply RBoltzmannMachine
        % Call temporal learning & RBM
        % Hidden-to-hidden weights
        [tweights, partialcost1, hidden] = RBoltzmann_real(tweights, input, Nhidden, Nhidden);
        % disp(weights1) % Check and balance - what are the weights now
        prevStep = hidden;
        
    %% COST %%
        % Cost is the sum for all patterns of the difference in Energy in
        % the (+) phase and (-) phase.
        % Cost only calculated for first RBM (visible & hidden layer)
        cost = cost + partialcost1;
        
   end % End cycling patterns
    
    %% Graphing %%
        timestepmatrix(t) = t; % Setup the x axis to be t, as in current timestep
        costmatrix(t) = cost; % Setup the y axis to be cost
        
        plot(timestepmatrix,costmatrix) % Plot timestep vs weighted summed cost
        ylabel('Cost')
        xlabel('Timestep');
% End cycling time steps    
end 