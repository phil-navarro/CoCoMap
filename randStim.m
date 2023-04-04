function A = randStim(N,T,prob_stim)
% This function produces a T x N binary random matrix with each entry
% having probability prob_stim of being '1'
%
% A - Matrix (for stimulation)
% N - Number of columns of matrix
%     (Total number of stimmed neurons)
% T - Number of rows of matrix
%     (Total number of trials)
% prob_stim - probability of neuron being stimmed per trial
%     (probability of neuron being stimmed per trial)
A = rand(T,N)<(prob_stim);