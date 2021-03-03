function [H,i,D,Hpost,data]=ITDadultsmodel2(set, cue)
    % INPUT
    % firsttrial and lasttrial are the first and the last trial watched by the participant
    % nsequence is the sequence watched by the participant (from 1 to 16)
    
    % OUTPUT
    % H is Entropy of the sequence at each trial
    % i is how surprising each trial is
    % D is the Kullback-Leibler Divergence between previous and current
    %   trial, here used as measure of learning progress.
    % data returns the subset of trials that have been seen
    
    
    %% Code by F. Poli ( https://francescopoli.weebly.com/ ) and R. B. Mars ( http://www.rbmars.dds.nl/lab.html )
       
    %  Paper: Poli, F., Serino, G., Mars, R. B., & Hunnius, S. (2020). 
    %  Infants tailor their attention to maximize learning. Science Advances, 6(39), eabb5053.
    
    %% Define Input
    idx=find(set(:,1)==cue);
    data=[idx,set(idx,2)];
    
    elements=[1 2 3 4]; % the locations where target could appear
    
    % Prepare data
    [r, c] = size(data);
    if c>r, data = data'; end
    
    %% set things up for next for-loop
    trial_results.e = [];
    trial_results.P = []; % probability of target appeareance for each location, for each trial
    trial_results.p = []; % probabiilty of target appeareance for location where it actually appeared
    trial_results.I = []; % surprise of seeing the target in each location, for each trial
    trial_results.i = []; % surprise of seeing the target in location where it actually appeared
    trial_results.H = []; % Entropy of the sequence for each trial
    trial_results.D = []; % KL Divergence for each trial
    trial_results.Hpost = []; % Entropy of the sequence for each trial after observation
    
    %% Loop over trials
    
    for trial = 1:length(data)
        vector = data(1:trial,2); %all the targets that have been seen so far
        % if it's the first elements, our expectations are based only on the prior:
        if length(vector) == 1, p = ones(length(elements),1)/length(elements); end
        
        % at every trial, we compute:
        I = -log2(p); 
        i = I(vector(end,1),1); %surprise of the event that has been presented
        H = sum(p.*I); %entropy        
        
        % Update estimated probabilities
        p = [];  
        for element = 1:length(elements)
            p = [p; (length(find(vector==elements(element)))+1)/(length(vector)+length(elements))];
        end
        
        trial_results.e = [trial_results.e vector(end,1)];
        trial_results.P = [trial_results.P p];
        trial_results.p = [trial_results.p p(vector(end,1),1)];
        trial_results.I = [trial_results.I I];
        trial_results.i = [trial_results.i i];
        trial_results.H = [trial_results.H H];
        
        % once we have the updated probabilities, we can compute KL Divergence
        
        prevtrial = length(vector) - 1;
        if length(vector) == 1
            D = sum(p.*(log2(p./[0.25; 0.25; 0.25; 0.25])));
        else
            D = sum(p.*(log2(p./trial_results.P(:,(prevtrial))))); %KL divergence
        end
        % For saccadic latencies (i.e., before target appearance, we use H
        % as defined in line 86, but for look-away and looking time (after
        % target appearance) we compute H here (i.e., after observing the
        % current trial):
        Hpost = -sum(p.*log2(p)); % comment this line to keep H for saccadic latency
        
        trial_results.Hpost = [trial_results.Hpost Hpost];
        trial_results.D = [trial_results.D D];
     
    end
    
    i = [trial_results.i];
    H = [trial_results.H];
    D = [trial_results.D];
    Hpost = [trial_results.Hpost];
end

