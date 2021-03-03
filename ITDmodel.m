function [Hz,iz,Dz,Rz,data]=ITDmodel(firsttrial, lasttrial, nsequence, zedscore)
% firsttrial and lasttrial are the first and the last trial watched by the participant
% nsequence is the sequence watched by the participant
% zedscore should be set to 1 if trasformation of IT into zscores is required


% for each sequence, specify the target locations
seq1=ones(15,1)*4; %100 loc 4
seq2=[3 3 2 3 3 4 3 1 3 3 4 2 3 3 1]'; %60 loc 3
seq3=[1 1 1 2 1 1 1 1 4 1 1 1 3 1 1]'; %80 loc 1
seq4=ones(15,1)*2; %100 loc 2
seq5=[3 4 3 3 3 3 2 3 3 3 1 3 3 3 3]'; %80 loc 3
seq6=[4 1 4 4 3 4 4 1 2 4 4 3 4 3 4]'; %60 loc 4
seq7=[1 1 3 1 1 1 1 2 1 1 1 1 4 1 1]'; %80 loc 1
seq8=[2 2 1 3 2 2 3 2 4 2 2 3 2 2 4]'; %80 loc 2

seq9=ones(15,1)*3; %100 loc 3
seq10=[4 3 2 4 4 1 4 4 4 3 4 1 4 4 2]'; %60 loc 4
seq11=[3 2 2 2 2 1 2 2 2 2 2 3 2 2 2]'; %80 loc 2
seq12=ones(15,1)*1; %100 loc 1
seq13=[4 4 4 4 1 4 4 4 3 4 4 4 4 2 4]'; %80 loc 4
seq14=[2 1 1 4 1 3 1 1 1 3 1 1 2 4 1]'; %60 loc 1
seq15=[2 2 4 2 2 2 1 2 2 2 2 4 2 2 2]'; %80 loc 2
seq16=[3 4 3 2 3 3 4 1 3 3 4 3 3 4 3]'; %60 loc 3

% and specify the cue (and target) type (i.e. what shape it has)
seq1 = [ones(15,1)*101 seq1];
seq2 = [ones(15,1)*102 seq2];
seq3 = [ones(15,1)*103 seq3];
seq4 = [ones(15,1)*104 seq4];
seq5 = [ones(15,1)*105 seq5];
seq6 = [ones(15,1)*106 seq6];
seq7 = [ones(15,1)*107 seq7];
seq8 = [ones(15,1)*108 seq8];
seq9 = [ones(15,1)*101 seq9];
seq10= [ones(15,1)*102 seq10];
seq11= [ones(15,1)*103 seq11];
seq12= [ones(15,1)*104 seq12];
seq13= [ones(15,1)*105 seq13];
seq14= [ones(15,1)*106 seq14];
seq15= [ones(15,1)*107 seq15];
seq16= [ones(15,1)*108 seq16];

% create one single sequence of sub sequences
sequence = [seq1 seq2 seq3 seq4 seq5 seq6 seq7 seq8 seq9 seq10 seq11 seq12 seq13 seq14 seq15 seq16];
nsequence = nsequence*2;
sequence = sequence(firsttrial:lasttrial,nsequence-1:nsequence);

% set things up for later (the final output)
iall = [];
Hall = [];
Dall = [];
CHall = [];
LRall = [];
soltall = [];
solt2all = [];
Rall=[];


    data=sequence;
    elements=[1 2 3 4]; % the locations where target could appear
    
    % Prepare data
    [r, c] = size(data);
    if c>r, data = data'; end
    
    % set things up for next for-loop
    trial_results.e = [];
    trial_results.P = []; 
    trial_results.p = []; 
    trial_results.I = []; 
    trial_results.i = []; 
    trial_results.H = [];
    trial_results.CH = [];
    trial_results.D = [];
    trial_results.LR = [];
    trial_results.solt = [];
    trial_results.solt2 = [];
    trial_results.R = [];

for trial = 1:length(data)
        vector = data(1:trial,2); %all the targets that have been seen so far
        % if it's the first elements, our expectations are based only on the prior:
        if length(vector) == 1, p = ones(length(elements),1)/length(elements); end
        
        % at every trial, we compute:
        I = -log2(p); %complexity of every event (each location is a potential event)
        i = I(vector(end,1),1); %surprise of the event that has been presented
        H = sum(p.*I); %entropy        
        
        % Update estimated probabilities
        p = [];  
        for element = 1:length(elements)
        p = [p; (length(find(vector==elements(element)))+1)/(length(vector)+length(elements))];
        % +1 because in the prior there is one element of the same type; +4 because in the prior there are 4 elements
        % The influence of the prior should be sampled by a distribution or
        % set to a certain value based on Kidd et al. (2012, 2014)
        end
        
        trial_results.e = [trial_results.e vector(end,1)];
        trial_results.P = [trial_results.P p];
        trial_results.p = [trial_results.p p(vector(end,1),1)];
        trial_results.I = [trial_results.I I];
        trial_results.i = [trial_results.i i];
        
        % once we have the updated probabilities, we can compute KL
        % Divergence, Entropy and Cross-Entropy
        
        prevtrial = length(vector) - 1;
        if length(vector) == 1
            D = sum(p.*(log2(p./[0.25; 0.25; 0.25; 0.25])));
        else
            D = sum(p.*(log2(p./trial_results.P(:,(prevtrial))))); %KL divergence
        end
        H = sum(p.*log2(p));
        if trial<7
            subvector=vector;
        else
            subvector=vector(end-5:end);
        end
        inst=[nnz(subvector==1),nnz(subvector==2),nnz(subvector==3),nnz(subvector==4)];             
        R = sum(log(factorial(inst)));
        CH = H+D; %Cross-entropy
        
        trial_results.H = [trial_results.H H];
        trial_results.CH = [trial_results.CH CH];
        trial_results.D = [trial_results.D D];
        trial_results.R = [trial_results.R R];
%------------------------------------------------------------------------        
        % We also estimate the learning rate. 
        % At trial 1 there is no computable learning rate:
        if trial > 3
        
        % to estimate it, we fit a quadratic curve to the cross entropy
        fit = polyfit(1:length(vector),trial_results.CH,2); % 2 is the degree of the polynomial
        
        % estimation of points based on quadratic equation
        yest = polyval(fit,1:length(vector));
        
        % estimation of the learning rate from yest
        LR=yest(end-1)-yest(end);

        % plot the fitting for every additional trial of every sequence
%         figure
%         plot(1:length(vector),trial_results.CH,'o')
%         hold on 
%         plot(1:length(vector),yest)
        
        % Estimate LR based on change in CH instead of based on yest
        %LR = trial_results.CH(end-1)-trial_results.CH(end)
        
        %estimate of LR based on derivative of fit
        syms x % initialize equation
        f = fit(1)*(x^2)+fit(2)*(x)+fit(3); %create equation
        der=diff(f); % take the derivative
        solt=subs(der,trial); %compute the slope for the current trial
        solt2=subs(der,trial+1); %and for the next one
        solt = -double(vpa(solt)); % transform it into single numrical value
        solt2 = -double(vpa(solt2)); % transform it into single numerical value
        
        trial_results.LR = [trial_results.LR LR];
        trial_results.solt = [trial_results.solt solt];
        trial_results.solt2 = [trial_results.solt2 solt2];
        end
        
        %NOTE: now I estimated the learning rate based on the fitting of cross
        %entropy rather than directly on cross entropy because in this way all previous trials are
        %considered, and not just the difference between this and the
        %previous one.
        
        % NOTE2: the rate of change of the cross entropy corresponds to the
        % learning rate only for an ideal observer.
        
end
% Add the estimation of this sequence to previous sequences

%P = trial_results.P;
%p(seq) = trial_results.p;
%I(seq) = trial_results.I1;
iall = [iall; trial_results.i];
Hall = [Hall; trial_results.H];
CHall = [CHall; trial_results.CH];
Dall = [Dall; trial_results.D];
LRall = [LRall; trial_results.LR];
soltall = [soltall; trial_results.solt];
solt2all = [solt2all; trial_results.solt2];
Rall = [Rall; trial_results.R];

iz=iall';
Hz=Hall';
Dz=Dall';
Rz=Rall';
if zedscore
% compute zscores on single sequence
iz = zscore(iall)';
Hz = zscore(Hall)';
CHz = zscore(CHall)';
Dz = zscore(Dall)';
LRz = zscore(LRall)';
sloz = zscore(soltall)';
end

end

