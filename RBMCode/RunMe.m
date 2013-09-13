% RBoltzmann script which calls RBoltmann.m and sets the parameters of the
% number of layers and how they interact

% This version includes temporal blurring of alpha = 0.1 whereby (0.1 * the
% previous grid cell layer) is added to the current grid cell layer at every
% time point. 


% Originally written by Mark Shaw

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('100GridCellsFiringInputData.mat')

% Nvisible is the dimensionality of the visible layer
% Nhidden is the dimensionality of the hidden layer

Nvisible = 100; % 100 grid cell visible units
Nhidden = 100; % 100 place cell hidden units
Ntimestep = 4000;
Npats = 3600; % number of points in rat's route

for skip = 100:100:360000
    % This is the pattern matrix
        patterns(skip/100,:) = savedAbstract(skip,:);
    % This is the input data from the grid cell firing. 
end

%%% Present each pattern to the network, and calculate the activity
%%% using repeated firing on each other (but not on the self)

% Setup weights matrix (random between -0.5 and 0.5)
for i = 1:Nvisible,
    for j = 1:Nhidden,
        weights1(i,j) = rand - 0.5;
    end
end

 % Add bias weights
 for i = 1:Nvisible
    weights1(i,Nhidden+1) = 0;
 end
 

 % Initialize prevStep
 for j = 1:(Nhidden+2)
     prevGridCellLayer(j) = 0;
 end
 
 alpha = 0.1; % Amount of temporal blurring

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Firing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:Ntimestep % Repeat Ntimestep times
   oldcost = 0; % Reset summed cost
   cost = 0; % Reset cost
   disp(t)
   % Cycle patterns
   for j = 1:Npats
        gridCellLayer = patterns(j,:); % Select one of the patterns
        input = gridCellLayer * alpha + prevGridCellLayer * (1-alpha);
    
        
    %% Apply RBoltzmannMachine
        [weights1, partialcost1, hidden] = RBoltzmann_real( weights1, input, Nvisible, Nhidden);
        % disp(weights1) % Check and balance - what are the weights now
        
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