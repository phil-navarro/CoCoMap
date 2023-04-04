%% function B = blockComb(N,T,n_stim)
% This function produces a T x N block combinatorial matrix by sampling
% n_stim-sized subsets from a set of N elements.
%
% B - Matrix (for stimulation)
% N - Number of columns of matrix
%     (Total number of stimmed neurons)
% T - Number of rows of matrix
%     (Total number of trials)
% n_stim - Number of non-zero elements per row
%     (Number of neurons stimmed per trial)
function B = blockComb(N,T,n_stim)
%% Initilize variables
B = zeros(T,N);
% Determine number of trials needed to stim all cells
stims_per_round = ceil(N/n_stim);
% Calculate number of times all cells can be stimmed for a given number of
% trials
num_stim_rounds = floor(T/stims_per_round);
% Calculate leftover trials
partial_stim_rounds = rem(T,stims_per_round);
pad_size = stims_per_round*n_stim - N;
%% 
for i = 1:num_stim_rounds
    stim_round = zeros(stims_per_round,N);
    rand_stim = reshape([randperm(N) zeros(1,pad_size)],n_stim,[])';
    for j = 1:stims_per_round
        stim_round(j,nonzeros(rand_stim(j,:))) = 1;
    end
    B(1+(i-1)*stims_per_round:i*stims_per_round,:) = stim_round;
    
end

stim_round = zeros(stims_per_round,N);
rand_stim = reshape([randperm(N) zeros(1,pad_size)],n_stim,[])';
for k = 1:partial_stim_rounds
    stim_round(k,rand_stim(k,:)) = 1;
end
B(num_stim_rounds*stims_per_round+1:end,:) = stim_round(1:partial_stim_rounds,:);
end