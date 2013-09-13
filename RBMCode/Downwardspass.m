function [ reconstruction ] = Downwardspass( weights, hidden, Nvisible, Nhidden)
%%% RBoltzmann.m is a Restricted Boltzmann machine demonstration
%%% RBoltzmann demonstrates a probabilistic learning algorithm. The system
%%% should learn to replicate the patterns they are presented with.

        
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
                
            %% ACTIVATION FN
                probability = (1 / (1 + exp(-sum))); % Boltzmann activation function
                % Roll a random number, if it's less than the probability
                % fire, if greater than don't fire
                if rand >= probability
                   reconstruction(i) = 0;
                else
                   reconstruction(i) = 1;
                end
                 
            end
        
        