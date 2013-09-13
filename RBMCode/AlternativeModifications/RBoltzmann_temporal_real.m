function [ weights, temporalWeights, partialcost, hidden_real ] = RBoltzmann_temporal_real( weights, temporalWeights, visible, prevStep_real, Nvisible, Nhidden)
%%% RBoltzmann.m is a Restricted Boltzmann machine demonstration
%%% RBoltzmann demonstrates a probabilistic learning algorithm. The system
%%% should learn to replicate the patterns they are presented with.

%% Requires RunMe.m to run

           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        energyMinus = 0; % Reset the energy of the minus phase
        energyPlus = 0; % Reset the energy of the plus phase


       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% (+) Phase
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % In the (+) phase, clamp the visible units to the pattern and update the hidden
       % units once.
       
       %Select positions from hidden layer one at a time to update
            for j = 1:(Nhidden)  % j cycles hidden units
                oldsum = 0;
                for i = 1:Nvisible % i cycles visible units
                       %% sum = Energy (weighted summed inputs)
                       sum = oldsum + visible(i) * weights(i,j);
                       oldsum = sum;
                       temp = visible(i)*weights(i,j);
                       % sum is the weighted summed input
                       % This is the same as the energy function
                       energyPlus = energyPlus + temp;
                    % 
                end
                
              %  sum = oldsum + temporalWeights(j) * prevStep_real(j);
                
                %%  ADDON TEMPORAL INPUT from previous step (t-1)
                 for y = 1:Nhidden % y cycles prev step hidden units
                     % Add prevStep hidden layer (t-1) firing onto hidden layer
                     sum = oldsum + prevStep_real(y) * temporalWeights(y,j);
                     oldsum = sum;
                 end
                
            %% ACTIVATION FN
                probability = (1 / (1 + exp(-sum))); % Boltzmann activation function
                % Roll a random number, if it's less than the probability
                % fire, if greater than don't fire
                if rand >= probability
                   hidden(j) = 0; 
                else
                   hidden(j) = 1;
                end
                
                hidden_real(j) = probability;
                 
            end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% (-) Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Generate the reconstruction of the visible layer
            
        %Select positions from hidden layer one at a time and generate
        %reconstruction
            for i = 1:Nvisible  % i cycles reconstruction units
                oldsum = 0;
                for j = 1:Nhidden % j cycles hidden units
                       sum = oldsum + hidden(j) * weights(i,j);
                       oldsum = sum;
                       % sum is the weighted summed input
                       % This is the same as the energy function
                end
                
            %% ACTIVATION FN (2nd)
                probability = (1 / (1 + exp(-sum))); % Boltzmann activation function
                % Roll a random number, if it's less than the probability
                % fire, if greater than don't fire
                if rand >= probability
                   reconstruction(i) = 0;
                else
                   reconstruction(i) = 1;
                end
                
                reconstruction_real(i) = probability;
                 
            end
            
        %% Update hidden layer
        
        %Select positions from visible layer (reconstruction) one at a time
        %and generate hidden layer
            for j = 1:Nhidden % j cycles hidden units
                oldsum = 0;
                for i = 1:Nvisible % i cycles visible(reconstruction) units
                       sum = oldsum + reconstruction(i) * weights(i,j);
                       oldsum = sum;
                       temp = reconstruction(i)*weights(i,j);
                       % sum is the weighted summed input
                       % This is the same as the energy function
                    % 
                    energyMinus = energyMinus + temp;
                end
                
               %%  ADDON TEMPORAL INPUT from previous step (t-1)
                 for y = 1:Nhidden % y cycles prev step hidden units
                     % Add prevStep hidden layer (t-1) firing onto hidden layer
                     sum = oldsum + prevStep_real(y) * temporalWeights(y,j);
                     oldsum = sum;
                 end
                
            %% ACTIVATION FN (3rd)
                probability = (1 / (1 + exp(-sum))); % Boltzmann activation function
                % Roll a random number, if it's less than the probability
                % fire, if greater than don't fire
                if rand >= probability
                   hidden2(j) = 0;
                else
                   hidden2(j) = 1;
                end
                
                hidden2_real(j) = probability;
            end
                    
           
        % Alter the weights based on similarity of reconstruction and
        % visible layer
        
        g = 0.03; % Set the learning rate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% LEARNING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            %% Update weights
            for i = 1:Nvisible
                for j = 1:(Nhidden+1)
                    if j == (Nhidden+1)
                        % Leave bias weights alone
                        deltaweights(i,j) = 0;                      
                    else
                        % Change weights by difference score
                        deltaweights(i,j) = g * (visible(i)*hidden_real(j) - hidden2_real(j)*reconstruction_real(i));
                    end
                end
            end
            
             
         %  Update weights using deltaweights
            oldweights = weights;
            weights = oldweights + deltaweights;
            
            
            %% Update temporal weights
            for y = 1:Nhidden
                for j = 1:(Nhidden+1)
                    if j == (Nhidden+1)
                        % Leave bias weights alone
                        deltaweights2(y,j) = 0;                      
                    else
                        % Change weights by difference score
                        deltaweights2(y,j) = g * (prevStep_real(y)*hidden_real(j) - hidden2_real(j)*prevStep_real(y));
                    end
                end
            end
            
        
            
         %  Update weights using deltaweights
            oldTemporalweights = temporalWeights;
            temporalWeights = oldTemporalweights + deltaweights2;
            
%% Calculate cost for this pattern
% Cost needs to be summed over all patterns in the script
partialcost = energyPlus - energyMinus;
        
        