% Created by Eugene M. Izhikevich, February 25, 2003
% Compressive Sensing + Group testing added by Phillip Navarro Oct, 2019
%
% This program simulates a partially observed neuronal network then tries
% to reconstruct the connections of that network using Compressive Sensing
%
%% Define experiment parameters
% Network Parameters
% TotalCells        : Total Number of cells in network
% N_obsCell         : Number of observed cells in network
% sparsity          : Probability of cell to cell connectivity
% synFailProb       : Probability current will not propagate to
%                     post-synaptic cells when pre-synaptic cell fires
% latency_mean      : mean of gaussian distributed latency of propagation
%                     of action potential to post-synaptic cell. Each
%                     connection has a fixed latency chosen randomly.
% latency_sd        : standard deviation of ...
% WS_beta           : Wattz-Strogatz Small world beta parameter. Sets local
%                     clustering of connectivity.
%                     0 = ring lattice. 1 = random graph
% frac_sin          : fraction of sinusoidal vs gaussian white input used
%                     to drive unobserved cells (unused)
%
% Stimulation Parameters
% test_samples      : Number of Trials
% fracStim          : Fraction of observed cells stimuated per trial
% offTargetProb     : Probability stimulation is off target. Extra cell is
%                     included in stimulation at this rate.
% phase             : phase of stim times based on frac_sin (unused)

% Reconstruction Parameters
% connect_thresh    : Weight Threshold for declaring a connection. Weights
%                     less than this value are not considered connections.
% lambda            : Sparsity-fit tradeoff used for basis pursuit. Weighs
%                     1-norm vs reconstruction error
%
%
% TotalCells = round(200./[1:-0.1:0.2 0.15:-0.05:0.05]);
TotalCells = 1000;
sparsity = [0.02:0.02:0.1] ; %cortico-cortico sparsity
%sparsity = 0.10;
%frac_sin = [0,0.5,1];
frac_sin = [0.0];
%synFailProb = [0:0.1:0.5];
synFailProb = 0;
%offTargetProb = [0:0.1:0.5];
offTargetProb = 0;
test_samples = [10:10:200];% 250:50:500 600:100:1000];
%test_samples = 100;

N_obsCell = [200];
lambda = 0.25;
%fractionStim = 0.1:0.1:0.5;
fractionStim = 0.1;

%latency_mean = 1.20;
latency_mean = 0.5;
%latency_sd = 0.60.*2.^(-2:1:2);
latency_sd = 0.0;
%phase = pi/4:pi/4:pi;
phase = pi/4;
WSbeta = 1;
%WSbeta = 0:0.2:1;
monteIter = 1;
scale_factor = 1.0;
% Whether stim is based on phase of signal or amplitude of network activity
phaseStim = 0; amplStim = 0;

% Results from simluation are stored in 'results' which is a cell array
% with dimensions   samples x parameter varied x metric
% % Metric indices are
% 1. False positive Rate. 2. False Negative Rate 3. True Correct Rate
% 4. Weighted Correct Rate 5. False Discovery Rate. 6. False Omission Rate
results = cell(length(test_samples),length(sparsity),length(1),2);
meanCorrect = zeros(length(test_samples),length(sparsity),length(1),5,2);
influences = cell(length(test_samples),length(sparsity),length(1),3);


%% Define Neural model
% Total Cells
for tot_c = 1:length(TotalCells)
    if length(TotalCells) > 1; tot_c, end
    for lat_m = 1:length(latency_mean)
        if length(latency_mean) > 1; lat_m, end
        for lat_sd = 1:length(latency_sd)
            if length(latency_sd) > 1; lat_sd, end
            %% Make cells
            % Excitatory neurons    Inhibitory neurons
            NeR = 0.8;              NiR = 1 - NeR;      % Ratio of excitatory:inhib
            Ne = round(NeR*TotalCells(tot_c));Ni = TotalCells(tot_c) - Ne;  % Number of excit/Inhib
            re=rand(Ne,1);          ri=rand(Ni,1);      % Randomize parameters
            a=[0.02*ones(Ne,1);     0.02+0.08*ri];      % Decay rate of conductance
            b=[0.2*ones(Ne,1);      0.25-0.05*ri];      % Sensitivity of conductance
            c=[-65+15*re.^2;        -65*ones(Ni,1)];    % Reset voltage of membrane
            d=[8-6*re.^2;           2*ones(Ni,1)];      % Reset conductance
            background_freq = 12;                       % Background Frequency of input to network
            
            % Calculate latency for each cell then convert into number of 0.5ms
            % time steps
            latency_raw = latency_mean(lat_m) + latency_sd(lat_sd)*randn(TotalCells(tot_c),TotalCells(tot_c));% Mean Latency 1.2 msec, 0.6 msec SD
            latency_steps = max(round(latency_raw * 2),1);
            
            for s=1:length(sparsity)
                if length(sparsity) > 1; s, end
                for beta = 1:length(WSbeta)
                    if length(WSbeta) > 1; beta, end
                    %% Connect cells to each other
                    IE_maxweights = scale_factor*[5, -10];
                    connect_thresh = 0.01*max(abs(IE_maxweights));%(1./sparsity(s));
                    %thalamic_sparsity = 0.2;    %thalamo-cortico sparsity, unused
                    
                    % Small-Worldness of Network
                    % Genrated unweighted digraph
                    hFull = WattsStrogatzDiG(Ne+Ni,round((Ne+Ni)*sparsity(s)),WSbeta(beta));
                    randWeights = [IE_maxweights(1)*rand(Ne+Ni,Ne),  IE_maxweights(2)*rand(Ne+Ni,Ni)];
                    shuffleTypes = randperm(TotalCells(tot_c));
                    % Weight digraph by multiplying with random weights.
                    S=randWeights(:,shuffleTypes).*hFull;
                    
                    % Label Excitatory cells for excitatory only case
                    E_idx = find(shuffleTypes<=Ne);
                    
                    % Assign latency to each connection
                    latency = latency_steps;% .* ( S ~= 0 );
                    max_delay = max(latency(:));
                    % Old S generation S=[5*rand(Ne+Ni,Ne).*(rand(Ne+Ni,Ne)<sparsity(s)),  -10*rand(Ne+Ni,Ni).*(rand(Ne+Ni,Ni)<sparsity(s))];
                    %tS_e = rand(Ne,1)<thalamic_sparsity;    % t-c excitory connections, Currently unused
                    %tS_i = rand(Ni,1)<thalamic_sparsity;    % t-c inhibitory connections, Currently unused
                    
                    % Iterate over experiment parameters
                    for a_r=1:length(fractionStim)
                        if length(fractionStim) > 1; a_r, end
                        for sf = 1:length(synFailProb)
                            if length(synFailProb) > 1; sf, end
                            for ot = 1:length(offTargetProb)
                                if length(offTargetProb)>1; ot, end
                                for n_o = 1:length(N_obsCell)
                                    if length(N_obsCell) > 1; n_o, end
                                    for fr_s = 1:length(frac_sin)
                                        if length(frac_sin)>1; fr_s, end
                                        %% Define current injection experiment parameters
                                        % Number of neurons with measureable signals
                                        N_obs = N_obsCell(n_o);
                                        % Number of measureable connections per neuron, average
                                        p = N_obs*sparsity(s);
                                        % NxT Binary matrix based on which N cells should be stimmed for which T
                                        % tests
                                        % Random iid test matrix
                                        
                                        cell_labels = randperm(TotalCells(tot_c));
                                        %cell_labels = E_idx(randperm(length(E_idx)));
                                        obs_labels = cell_labels(1:N_obs);
                                        %obs_labels = E_idx(cell_labels(1:N_obs));
                                        obs_S = S(obs_labels,obs_labels);
                                        
                                        
                                        % Stimulation Matrix
                                        T = max(test_samples);
                                        A = double(rand(T,N_obs)<(fractionStim(a_r)));
                                        % A = eye(T); Single cell stim test
                                        % Block combinatorial matrices are
                                        % also valid
                                        % A = blockComb(N_obs,T,10);
                                        
                                        % Voltage measurements @ post-synaptic
                                        % cells
                                        y = zeros(T,N_obs);
                                        test_every_x_msec = 50; % Perform a trial every X milliseconds
                                        
                                        % Current injection informed by literature
                                        % Using impulse responses for now
                                        for ph_s = 1:length(phase)
                                            if length(phase)>1; ph_s, end
                                            
                                            if phaseStim    % Phase-based stim
                                                stimPeriod = 1000/background_freq;
                                                phaseStimRads = (phase(ph_s))*1/(2*pi)*stimPeriod;
                                                Total_Timesteps = ceil(T*stimPeriod + 1);
                                                test_times = round(phaseStimRads:stimPeriod:Total_Timesteps);
                                            elseif amplStim % Membrane amplitude level based stim
                                                test_times = [];
                                            else            % Fixed interval stim
                                                Total_Timesteps = ceil(T*test_every_x_msec + 2);
                                                test_times = test_every_x_msec:test_every_x_msec:Total_Timesteps;
                                            end
                                            
                                            v = -65*ones(Ne+Ni,Total_Timesteps + max_delay);    % Initial values of v
                                            u = repmat(b,1,Total_Timesteps + max_delay).*v;     % Initial values of u
                                            I_stim = zeros(Ne+Ni,Total_Timesteps + max_delay);  % Contribution from stimulation experiment
                                            I_spont = zeros(Ne+Ni,Total_Timesteps + max_delay); % Spontaneous firing
                                            firings=[];             % spike timings
                                            tp_E = scale_factor*3.0; tp_I = scale_factor*1.2;   % thalamic input power, excitatory and inhibitory
                                            %tp_E = 5.0; tp_I = 2.0;
                                            %% Conduct current injection experiment
                                            test_number = 1;
                                            i = [];
                                            for t=1:Total_Timesteps            % simulation of 20*10000 ms
                                                I_sin = 1/sqrt(2/pi)*sin((t)*2*pi/(1000/background_freq));  % normalize area by half-gauss distribution expected value sqrt(2/pi) divided by sin expected value 2/pi
                                                I =[tp_E*I_sin*(frac_sin(fr_s)) + tp_E*(1-frac_sin(fr_s))*randn(Ne,1);...
                                                    tp_I*I_sin*(frac_sin(fr_s)) + tp_I*(1-frac_sin(fr_s))*randn(Ni,1)]; % thalamic input. Sinusoidal component + gaussian random noise
                                                I = I(shuffleTypes); % Make sure current input is given to appropriate type
                                                I(obs_labels) = 0;
                                                fired=find(v(:,t)>=30);    % indices of spikes
                                                firings=[firings; t+0*fired,fired];
                                                v(fired,t)=c(fired);
                                                u(fired,t)=u(fired,t)+d(fired);
                                                % Sum currents from synaptic projections
                                                % with probably of individual synapses failing
                                                % Used to generate figure
                                                %I_stim(:,t) = sum(S(:,intersect(fired,i)).*...                       % Sum of stimulated component
                                                %    (rand(size(S(:,intersect(fired,i))))>synFailProb(sf)),2);
                                                %I_spont(:,t) = I + sum(S(:,setdiff(fired, intersect(fired,i))).*...               % Sum of spontaeous firing
                                                %    (rand(size(S(:,setdiff(fired, intersect(fired,i)))))>synFailProb(sf)),2);
                                                %                                         I= I + sum(S(:,fired).*...
                                                %                                             (rand(size(S(:,fired)))>synFailProb(sf)),2);
                                                
                                                % Sum all incoming PSCs based on
                                                % arrival time and latency
                                                I_prop = accumarray([repmat((1:TotalCells(tot_c))',length(fired),1) ,reshape(latency(:,fired),[],1)] ,... % indices of elements to be added based on latency
                                                    reshape(S(:,fired), [],1).*... % values to be added based on which cells fired
                                                    (rand(size(reshape(S(:,fired), [],1))) > synFailProb(sf)),... % synaptic failure
                                                    [TotalCells(tot_c),max_delay]); % size of output array
                                                I_stim(:, t:t+max_delay-1) = I_stim(:, t:t+max_delay-1) + I_prop;
                                                I = I_stim(:,t) + I;
                                                %I = I_stim(:,t) + I_spont(:,t);
                                                % Update v based on A matrix
                                                v(:,t+1)=v(:,t)+0.5*(0.04*v(:,t).^2+5*v(:,t)+140-u(:,t)+I); % step 0.5 ms
                                                u(:,t+1)=u(:,t)+a.*(b.*v(:,t+1)-u(:,t));                    % for numerical stability
                                                % Force cells to fire.
                                                % Stims at preset intervals
                                                if(~amplStim && ismember(t,test_times) && test_number <= T)
                                                   i = A(test_number,:).*obs_labels; % Get indices of test cells
                                                    % Model off target stim by inclusing extra cell
                                                    if(rand < offTargetProb(ot))
                                                        i = [i obs_labels(randi(length(obs_labels)))];
                                                    end
                                                    i(i==0) = [];
                                                    v(i,t+1) = 30;                    % Set voltage to threshold
                                                    test_number = test_number + 1;
                                                elseif(amplStim && test_number <= T && norm(v(:,t)-c,2) < amplStimThresh)
                                                    %
                                                end
                                                
                                            end
                                            
                                            %plot(firings(:,1),firings(:,2),'.'); % plot raster
                                            pre_stim_avgNumSamples = 5;
                                            post_stim_avgNumSamples = max_delay;
                                            pre_stim_idx = coloncatrld((test_times - pre_stim_avgNumSamples + 1),(test_times));
                                            post_stim_idx = coloncatrld((test_times + 1),(test_times + post_stim_avgNumSamples + 1));
                                            
                                            pre_stim_amp_mat = mean(...
                                                reshape(v(obs_labels,pre_stim_idx),...
                                                length(obs_labels),pre_stim_avgNumSamples,[]),2);
                                            
                                            post_stim_amp_mat = reshape(v(obs_labels,post_stim_idx),...
                                                length(obs_labels),post_stim_avgNumSamples + 1,[]);
                                            
                                            %pre_stim_amp = squeeze(pre_stim_amp_mat);
                                            %post_stim_amp = squeeze(mean(post_stim_amp_mat));
                                            
                                            % Compute difference across pre-stim
                                            % and post stim period
                                            dv = diff([pre_stim_amp_mat, post_stim_amp_mat],1,2);
                                            
                                            dv_cs = cumsum(dv,2);
                                            dv_ratio = dv_cs(:,2:end,:)./dv_cs(:,1:end-1,:);
                                            % Calculate average decay
                                            decay_coeff = zeros(N_obs,1);
                                            dv_test2 = zeros(N_obs,1);
                                            for j = 1:N_obs
                                                [binCount, binEdges] = histcounts(dv_ratio(j,:,:),100,'BinLimits', [0.5 1]);
                                                [~,maxbinIdx] = max(binCount);
                                                decay_coeff(j) = mean(binEdges(maxbinIdx:maxbinIdx+1));
                                            end
                                            for j = 1:length(test_times)
                                                dv_test2(:,j) = sum(I_stim(obs_labels,test_times(j):test_times(j)+max_delay),2); % Ground truth for testing non varying latency
                                            end
                                            sum_ontoCell = sum(obs_S>0,2); % Ground truth for testing
                                            % Normalize change in amplitude by
                                            % expected decay using no prior
                                            % knowledge of system dynamics
                                            dv_normRise = [dv_cs(:,1,:) dv_cs(:,2:end,:) - dv_cs(:,1:end-1,:).*decay_coeff];
                                            
                                            % Only have excitatory cells
                                            % Sum rising edges only
                                            %dv_test = squeeze(sum(dv_normRise.*(dv_normRise>0.26),2));
                                            
                                            % For excit and inhib
                                            dv_test = squeeze(sum(dv,2));
                                            
                                            influence_BP = zeros(N_obs,N_obs);
                                            influence_BP2 = zeros(N_obs,N_obs);

                                            for ts = 1:length(test_samples)
                                                for lam = 1:length(lambda)
                                                    for mi = 1:monteIter
                                                        if monteIter == 1
                                                            rand_sample = 1:test_samples(ts);
                                                        else
                                                            rand_sample = randperm(T, test_samples(ts));
                                                        end
                                                        %% Decode connections from firing matrix
                                                        % Basis Pursuit with
                                                        % constraints
                                                        results_BP = zeros(N_obs, length(connect_thresh), 6);
                                                        for j = 1:N_obs
                                                            % Find trials where the measured cell was/was not stimulated
                                                            stim_trials = find(A(:,j)); no_stim_trials = setdiff(rand_sample,stim_trials);
                                                            %x_BP = netCS_Noise2(A(no_stim_trials,:),dv_test(j,no_stim_trials)', lambda(lam),setdiff(pos_label,j),setdiff(neg_label,j)); x_connect_BP = [x_BP obs_S(j,:)'];
                                                            %pos_label = find(obs_S(j,:)' > 0); neg_label = find(obs_S(j',:) < 0); % only connections tagged
                                                            pos_label = find(sum(obs_S)' > 0); neg_label = find(sum(obs_S)' < 0); % All cells tagged
                                                            %pos_label = []; neg_label = []; % No tags

                                                            x_BP = netCS_Noise2(A(rand_sample,:),dv_test2(j,rand_sample)', lambda(lam),setdiff(pos_label,j),setdiff(neg_label,j)); x_connect_BP = [x_BP obs_S(j,:)'];
                                                            x_connect_BP(j,:) = [];
                                                            influence_BP(j,:) = x_BP;
                                                            x_Sort2_BP = sortrows(x_connect_BP, 'descend', 'ComparisonMethod','abs');
                                                            x_ratio_BP = abs(x_Sort2_BP(1:end-1,1))./abs(x_Sort2_BP(2:end,1));
                                                            for ct = 1:length(connect_thresh)
                                                                cutoff_idx_BP = find(abs(x_Sort2_BP(:,1))<connect_thresh(ct),1);
                                                                xw_predict_BP = x_Sort2_BP(1:cutoff_idx_BP,2);
                                                                xw_total_BP = sum(abs(x_Sort2_BP(:,2)));
                                                                xw_percentCorrect_BP = sum(abs(xw_predict_BP))/xw_total_BP;
                                                                False_Positive_BP = sum(xw_predict_BP==0)/sum(x_Sort2_BP(:,2)==0);
                                                                False_Negative_BP = sum(x_Sort2_BP((1+cutoff_idx_BP):end,2)~=0)/sum(x_Sort2_BP(:,2)~=0);
                                                                True_Correct_BP =  (length(x_Sort2_BP) - sum(xw_predict_BP==0) - sum(x_Sort2_BP((1+cutoff_idx_BP):end,2)~=0)) ...
                                                                    ./length(x_Sort2_BP);
                                                                False_Negative_BP(isnan(False_Negative_BP)) = 0;
                                                                xw_percentCorrect_BP(isnan(xw_percentCorrect_BP)) = 1;
                                                                False_Discovery_BP = sum(xw_predict_BP==0)/length(xw_predict_BP);
                                                                False_Omit_BP = sum(x_Sort2_BP((1+cutoff_idx_BP):end,2)~=0)/length(x_Sort2_BP((1+cutoff_idx_BP):end,2));
                                                                results_BP(j,ct,1) = False_Positive_BP;% False Positive Rate
                                                                results_BP(j,ct,2) = False_Negative_BP;% False Negative Rate
                                                                results_BP(j,ct,3) = True_Correct_BP; % Fraction correct connections detected
                                                                results_BP(j,ct,4) = xw_percentCorrect_BP;% Weighted fraction correct connections detected
                                                                results_BP(j,ct,5) = False_Discovery_BP;
                                                                results_BP(j,ct,6) = False_Omit_BP;
                                                            end
                                                        end
                                                        meanCorrect(ts,s,1,1,1) = mean(squeeze(results_BP(:,:,1)));
                                                        meanCorrect(ts,s,1,2,1) = mean(squeeze(results_BP(:,:,2)));
                                                        meanCorrect(ts,s,1,3,1) = mean(squeeze(results_BP(:,:,3)));
                                                        meanCorrect(ts,s,1,4,1) = mean(squeeze(results_BP(:,:,4)));
                                                        meanCorrect(ts,s,1,5,1) = length(firings)./Total_Timesteps*1000;
                                                        results{ts,s,1,1} = results_BP;
                                                        
                                                        
                                                        % Basis Pursuit with no
                                                        % constraints
                                                        results_BP = zeros(N_obs, length(connect_thresh), 6);
                                                        for j = 1:N_obs
                                                            % Find trials where the measured cell was/was not stimulated
                                                            stim_trials = find(A(:,j)); no_stim_trials = setdiff(rand_sample,stim_trials);
                                                            %x_BP = netCS_Noise2(A(no_stim_trials,:),dv_test(j,no_stim_trials)', lambda(lam),setdiff(pos_label,j),setdiff(neg_label,j)); x_connect_BP = [x_BP obs_S(j,:)'];
                                                            %pos_label = find(obs_S(j,:)' > 0); neg_label = find(obs_S(j',:) < 0);
                                                            pos_label = []; neg_label = [];
                                                            x_BP = netCS_Noise2(A(rand_sample,:),dv_test2(j,rand_sample)', lambda(lam),setdiff(pos_label,j),setdiff(neg_label,j)); x_connect_BP = [x_BP obs_S(j,:)'];
                                                            x_connect_BP(j,:) = [];
                                                            influence_BP2(j,:) = x_BP;
                                                            x_Sort2_BP = sortrows(x_connect_BP, 'descend', 'ComparisonMethod','abs');
                                                            x_ratio_BP = abs(x_Sort2_BP(1:end-1,1))./abs(x_Sort2_BP(2:end,1));
                                                            for ct = 1:length(connect_thresh)
                                                                cutoff_idx_BP = find(abs(x_Sort2_BP(:,1))<connect_thresh(ct),1);
                                                                xw_predict_BP = x_Sort2_BP(1:cutoff_idx_BP,2);
                                                                xw_total_BP = sum(abs(x_Sort2_BP(:,2)));
                                                                xw_percentCorrect_BP = sum(abs(xw_predict_BP))/xw_total_BP;
                                                                False_Positive_BP = sum(xw_predict_BP==0)/sum(x_Sort2_BP(:,2)==0);
                                                                False_Negative_BP = sum(x_Sort2_BP((1+cutoff_idx_BP):end,2)~=0)/sum(x_Sort2_BP(:,2)~=0);
                                                                True_Correct_BP =  (length(x_Sort2_BP) - sum(xw_predict_BP==0) - sum(x_Sort2_BP((1+cutoff_idx_BP):end,2)~=0)) ...
                                                                    ./length(x_Sort2_BP);
                                                                False_Negative_BP(isnan(False_Negative_BP)) = 0;
                                                                xw_percentCorrect_BP(isnan(xw_percentCorrect_BP)) = 1;
                                                                False_Discovery_BP = sum(xw_predict_BP==0)/length(xw_predict_BP);
                                                                False_Omit_BP = sum(x_Sort2_BP((1+cutoff_idx_BP):end,2)~=0)/length(x_Sort2_BP((1+cutoff_idx_BP):end,2));
                                                                results_BP(j,ct,1) = False_Positive_BP;% False Positive Rate
                                                                results_BP(j,ct,2) = False_Negative_BP;% False Negative Rate
                                                                results_BP(j,ct,3) = True_Correct_BP; % Fraction correct connections detected
                                                                results_BP(j,ct,4) = xw_percentCorrect_BP;% Weighted fraction correct connections detected
                                                                results_BP(j,ct,5) = False_Discovery_BP;
                                                                results_BP(j,ct,6) = False_Omit_BP;
                                                            end
                                                        end
                                                        meanCorrect(ts,s,1,1,2) = mean(squeeze(results_BP(:,:,1)));
                                                        meanCorrect(ts,s,1,2,2) = mean(squeeze(results_BP(:,:,2)));
                                                        meanCorrect(ts,s,1,3,2) = mean(squeeze(results_BP(:,:,3)));
                                                        meanCorrect(ts,s,1,4,2) = mean(squeeze(results_BP(:,:,4)));
                                                        meanCorrect(ts,s,1,5,2) = length(firings)./Total_Timesteps*1000;
                                                        results{ts,s,1,2} = results_BP;
                                                        
                                                        % Store influence maps
                                                        influences{ts,s,1,1} = influence_BP;
                                                        influences{ts,s,1,2} = influence_BP2;
                                                        influences{ts,s,1,3} = obs_S;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
% Save output
sim_number = 1;
filename = strcat('SimLab_','IE_',num2str(sim_number),'_',date);
while(isfile(strcat(filename,'.mat')))
    sim_number = sim_number + 1;
    filename = strcat('SimLab_','IE_',num2str(sim_number),'_',date);
end
save(filename,'A','dv_test','TotalCells','sparsity','frac_sin','synFailProb','offTargetProb','latency_mean','latency_sd',...
    'test_samples', 'N_obsCell', 'lambda', 'fractionStim', 'phase', 'WSbeta','connect_thresh', ...
    'phaseStim', 'amplStim', 'meanCorrect','influences','results','IE_maxweights','tp_E','tp_I','background_freq')