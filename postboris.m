function firstlastbor = postboris(fix, stamps, lastseen, subjnumber)

% Import the data
xlsxname = strcat('C:/Users/u661121/Desktop/StatLearn\infantsdata\Infant',num2str(subjnumber),'\borised',num2str(subjnumber),'.xlsx');
hLA = xlsread(xlsxname);

% Clear temporary variables
clearvars raw;

% function that takes the boris output csv file of a participant (located
% in data) and his/her fixations matrix as computed by preRanalysis and
% returns the LA trials as estimated from handcoding in boris. It also
% require to know the subjnumber

%------------------------------------------------------------------------
% keep only cols that you need and convert them to array
% hLAfull = Boris11;
% hLA = table2array(hLAfull(:,1)); % time of the event
% behav=table2array(hLAfull(:,6)); %type of event
% behav2=table2array(hLAfull(:,9)); %if lookaway, start or end
% 
% for n=1:length(hLA)
%     if behav(n)=='look_away'
%         if behav2(n)=='START'
%             hLA(n,2)=1;
%         else
%             hLA(n,2)=2;
%         end
%     elseif behav(n)=='start'
%         hLA(n,2)=0;
%         tstart = hLA(n,1);
%         nstart = n;
%     else
%         hLA(n,2)=4;
%     end
% end
 
% set the start to the moment t = 0
idxstart= find(hLA(:,2)==0);
tstart = hLA(idxstart,1);
nstart = idxstart;
hLA(:,1)=hLA(:,1)-tstart; 
%remove the start line and what's before it (if anything)
%hLA = hLA(nstart+1:end,:);
%trasform time from sec to ms
%hLA(:,1)=hLA(:,1)./1000;

% KEEP LA only if longer than a certain minimum
minimumLA = 1000; %(ms)
% first remove the start, we don't need it anymore.
hLA(idxstart,:)=[];
%then find all the look-away ends (2)
ila=find(hLA(:,3)==2);
imin=[];
for n=ila(2):2:ila(end) %for every look-away end
    if hLA(n,1)-hLA(n-1)<minimumLA %check if LA time is less then minimum
        imin(n-1:n,1)=1; % create an idx
    end
end

% eliminate the look-aways shorter than 1 sec from hLA
hLA(find(imin==1),:)=[];

% now for each sequence, determine in what trials infants were not looking
% away. Take the time between sequence n beginning and sequence n+1
% beginning. 
% for every LA coded in boris, we can see the trials during which it
% happened. If trials happened during a LA, we tag them with a 0.
% Otherwise, they are a 1 (valid).
fix(:,8)=0;
for n=1:length(hLA)
    if hLA(n,3)==1 
        idxt = find(fix(:,2)<hLA(n,1));
        % if we don't find any fix before the current LA, it's because
        % infants were looking away before the first sequence even started.
        if isempty(idxt)
            % in that case, we have to check if they were still looking
            % away when the first trial started. So we find the moment in
            % which lookaway ends
            idxt2 = find(fix(:,2)<hLA(n+1,1));
            idxt2 = idxt2(end);
            % if this is empty as well, we do not worry about this LA
            % because it happened entirely before first seq started.
            % if it's not empty:
            if isempty(idxt2)==0
                % I set the first row of fix as the one the LA started at
                idxt=n;
                fix = [fix(1:idxt,:); [0 hLA(n,1) 0 NaN NaN 13 NaN 1]; fix(idxt+1:end,:)];
            end
        else %if we find a fix before LA, everything is normal and we check
            % when exactly that happens.
            idxt = idxt(end);
            fix = [fix(1:idxt,:); [0 hLA(n,1) 0 NaN NaN 13 NaN 1]; fix(idxt+1:end,:)];       
            % idxt is the moment the lookaway started. Now find the moment in
            % which lookaway ends
            idxt2 = find(fix(:,2)<hLA(n+1,1));
            idxt2 = idxt2(end);
            fix = [fix(1:idxt2,:); [0 hLA(n+1,1) 0 NaN NaN 13 NaN 1]; fix(idxt2+1:end,:)];
        end        
                %between LA start and end, fix8 is 1
        fix(idxt+1:idxt2+1,8)=1;
    end
end

k=[];
firstlastbor=[];
for k=1:lastseen
    if stamps(1,k)==stamps(1,k+16)
        startseq = 111;
        endseq = 111;
        lookaway = 111;
    else
    % I take only the current sequence, which starts at the moment in time
    % of its trial 1 end ends at the moment in time of next sequence trial 1
    thisfixstart = find(fix(:,2)==stamps(1,k));
    if k < lastseen
        thisfixend = find(fix(:,2)==stamps(1,k+1))-1;
    else % which means: k=lastseen, so we are in the last sequence
        thisfixend = length(fix);
    end
    thissfix=fix(thisfixstart:thisfixend,:);
    [nrows,~]=size(thissfix);
    [lengththissfix,~]=size(thissfix); % nrows and lengththisfix can be used interchangably
    % now I find the first moment within the sequence that has no Lookaway
    idxno = find(thissfix(:,8)==0);
    if isempty(idxno) %it means there are only LAs during the sequence
            startseq = 111;
            endseq = 111;
            lookaway = 111;
    % another possibility is that they looked but when no cues or targets were shown        
    elseif nnz(thissfix(idxno:end,6)==11)==0 || nnz(thissfix(idxno:end,6)==12)==0
            startseq = 111;
            endseq = 111;
            lookaway = 111;
    else
    idxno = idxno(1);
    n=idxno; %this is when the infant starts looking at the sequence
    % and I also trace back what trial this is
    % I look for the first cue or target presented after look back to seq
    while thissfix(n,6)~=11 && thissfix(n,6)~=12 && n<nrows
            n = n+1;
    end
    % I check how many cues or targets were presented before (i.e. an index
    % of which trial we are in when infant starts looking at sequence)
    if thissfix(n,6)==11
        startseq = nnz(thissfix(1:n,6)==11);
    elseif thissfix(n,6)==12
        startseq = nnz(thissfix(1:n,6)==12);
    end
    
    m=idxno;
    countcue=0;
    counttar=0;
    % now I check how many consecutive trials I have before a LA
    while thissfix(m,8)==0 && m<lengththissfix
        if thissfix(m,6)==11
            countcue=countcue+1; %how many cues
        elseif thissfix(m,6)==12
            counttar=counttar+1; %how many targets
        end        
        m=m+1;
    end
    nela = m; %(NExt LookAway)
    
 %If we reached the end of the sequence and still no LA was found, it means there was no LA.
 lookaway=1;
 if m==lengththissfix && thissfix(m,8)==0
    % in this case, lookaway is false
    lookaway=0;
 end
    
    % if countcue or counttar are bigger than 3, the infant has looked at
    % enough trials to analyze the sequence.
    % if countcue>counttar, 
 if counttar > 3 || countcue > 3 % THIS MUST CHANGE IF WE CONSIDER LESS THEN 3 TRIALS STILL VALID!
    % I look for the last cue or target presented after look away from seq
    while thissfix(m,6)~=11 && thissfix(m,6)~=12 
            m = m-1;
    end
    % remember that n indexes to the moment first cue or tar was looked.
    % if n and m are both cues, i.e.:
    if thissfix(n,6)==thissfix(m,6) && thissfix(m,6)==11
        % the participant looked until the cue of trial m, but not at the
        % target. This means that participant stopped looking
        % after the cue appeared, but before the target did. Since s/he has not
        % seen the target, his or her look away have been caused by what s/he
        % has seen in the last trial, meaning that we have to link lookaway to
        % the previous, and not this trial. thus,
        endseq = countcue-1;
    % if n is a target but m is a cue, i.e.:     
    elseif thissfix(n,6)~=thissfix(m,6) && thissfix(n,6)==12
        % something similar happens, but countcue in this case is wrong,
        % because it missed out the first cue. However, infant saw the
        % first target, meaning that s/he has info about its location.
        endseq = counttar-1;
    % if n is a cue or target but and m is a target
    elseif thissfix(m,6)==12
        % the participant looked away because of what s/he saw now
        endseq = counttar;
    end
 elseif counttar <= 3 && countcue <= 3 && nela <= lengththissfix
    % if I have less than 4 trials, there are two options: infant looked
    % away early during the sequence; infant looked at the sequence,
    % stopped looking e.g. whithin first trial, and then looked back at the
    % sequence. In this case we want to measure how many trials s/he looked at after
    % this look back.
    % I want to start from the moment infant looked away (nela) and see
    % whether when s/he looked back s/he looked for at least 4 trials.
    % So I crete a new thissfix that is a subsequence of thissfix which
    % starts from the moment infant looked away the first time.
    subbfix = thissfix(nela:end,:);
    [nsubrows,~]=size(subbfix); %this is the length of subfix
    [lengthsubbfix,~]=size(subbfix); %same thing
    idxsub = find(subbfix(:,8)==0);
    if isempty(idxsub) %it means there are only LAs during the sequence
            startseq = 111;
            endseq = 111;
            lookaway = 111;
    else
    idxsub = idxsub(1);
    n=idxsub;
    while subbfix(n,6)~=11 && subbfix(n,6)~=12 && n<nsubrows
            n = n+1;
    end
    % if we got to the end of the sequence and there is no cue or target
    % found, it means the sequence was over. In that case, we cannot use
    % this sequence. In such case, firstlastbor = 111 = invalid
    if n == nsubrows && subbfix(n,8)==0
        startseq =111;
        endseq=111;
        lookaway=111;
        % We continue only if that is false:
    else
        % we must do everything we did before but with subbfix instead of
        % thissfix.
        % I check how many cues or targets were presented before (i.e. an index
        % of which trial we are in when infant starts looking at sequence).
        % this still requires thissfix.
        if thissfix(n,6)==11
            startseq = nnz(thissfix(1:n,6)==11);
        elseif thissfix(n,6)==12
            startseq = nnz(thissfix(1:n,6)==12);
        end
    
        m=idxsub;
        countcue=0;
        counttar=0;
        % now I check how many consecutive trials I have before a LA
        while subbfix(m,8)==0 && m<nsubrows
            if subbfix(m,6)==11
                countcue=countcue+1;
            elseif subbfix(m,6)==12
                counttar=counttar+1;
            end        
        m=m+1;
        end
    nela = m;
    
    %If we reached the end of the sequence and still no LA was found, it means there was no LA.
    lookaway=1;
    if m==lengthsubbfix && subbfix(m,8)==0
    % in this case, noaway is true
    lookaway=0;
    end
    
    % if countcue or counttar are bigger than 3, the infant has looked at
    % enough trials to analyze the sequence.
    % if countcue>counttar, 
        if counttar > 3 || countcue > 3   
            % I look for the last cue or target presented after look away from seq
            while subbfix(m,6)~=11 && subbfix(m,6)~=12 
                m = m-1;
            end
            % remember that n indexes to the moment first cue or tar was looked.
            % if n and m are both cues, i.e.:
            if subbfix(n,6)==subbfix(m,6) && subbfix(m,6)==11
            % the participant looked until the cue of trial m, but not at the
            % target. This means that participant stopped looking
            % after the cue appeared, but before the target did. Since s/he has not
            % seen the target, his or her look away have been caused by what s/he
            % has seen in the last trial, meaning that we have to link lookaway to
            % the previous, and not this trial. thus,
            endseq = countcue-1+startseq-1;
            % if n is a target but m is a cue, i.e.:     
            elseif subbfix(n,6)~=subbfix(m,6) && subbfix(n,6)==12
            % something similar happens, but countcue in this case is wrong,
            % because it missed out the first cue. However, infant saw the
            % first target, meaning that s/he has info about its location.
            endseq = counttar-1+startseq-1;
            % if n is a cue or target but and m is a target
            elseif subbfix(m,6)==12
            % the participant looked away because of what s/he saw now
            endseq = counttar+startseq-1;
            end
        else
            % in the remaining cases, we do not have enough data STILL IF
            % WE ASSUME THAT TRIALS < 3 is a problem
            startseq = 111;
            endseq = 111;
            lookaway = 111;
        end
    end
    end
    end
    end
    end
 firstlastbor=[firstlastbor; startseq, endseq, lookaway];
end
%end the function
end
        
    
    
    
    