% -----------------------------------------------------------------
% Infant processing
% ----------------------------------------------------------------- 
%clearvars -except participant; close all; fclose('all'); clc;
%dbstop if error;
%commandwindow;

% switch to 0 if you don't want to save data at the end
saving = 0;
plotting=0;

% set the directory
cd 'C:/Users/u661121/Desktop/StatLearn'

% import the I2MC fixations files from each participant
folders.data = 'infantsoutput'; % the input of this algorithm is the output of I2MC
folders.output = 'Rinfantsinput'; % folder for output (will use structure in folders.data for saving output)

% create textfile and open for writing analyzed data (with saccades, dwell 
% times, look-away probability and anticipatory looking)
%fid = fopen(fullfile(folders.output,'allfix.txt'),'w');
%fprintf(fid,'Participant\Time\Duration\XPos\YPos\FixOrStamp\TypeAndPosition\AOIs\LookToTarget');

[fold,nfold] = FolderFromFolder(folders.data);

% Inizialize the output
allfix={};
allslat={};
allAL = [];
alldwell=[];
alllinger=[];
allboris=[];
allLA=[];
predictslat=[];
predictdwell=[];
predictlinger=[];
markasgood=[];

for e=1:nfold
%-------------------------------------------------------------------------
% Import data from text file.
% Initialize variables.
filename = char(strcat('C:/Users/u661121/Desktop/StatLearn\infantsoutput\',fold(e).name,'\allfixations.txt'));
delimiter = '\t';
startRow = 2;
% Read columns of data as text: For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% Close the text file.
fclose(fileID);
% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end
% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
% Create output variable
allfixations = cell2mat(raw); 
%-------------------------------------------------------------------------    
% Get only what you need from the I2MC output (i.e. participant number,
% beginning and end of saccades, their duration and location)
fixations = [allfixations(:,end-1) allfixations(:,1:5)];
stamps = load(char(strcat('C:/Users/u661121/Desktop/StatLearn\infantsdata\',fold(e).name,'\stampstimuli.mat')));
stamps = stamps.stampstimuli;

lasttrial = load(char(strcat('C:/Users/u661121/Desktop/StatLearn\infantsdata\',fold(e).name,'\lasttrial.mat')));
lasttrial = lasttrial.lasttrial;

% create a row for every fixation beginning (1) or end (2)
fix=zeros(length(fixations),7);
for n=1:2:2*length(fixations(:,1))
    m=round(n/2);
    fix(n,:)=[fixations(m,1:2) fixations(m,4:6) 1 0];
    fix(n+1,:)=[fixations(m,1) fixations(m,3:6) 2 0];
end
%-------------------------------------------------------------------------

% Add timestamps for cues (11) and targets (12)
% and also type of cue and location of target (col 7)

% for each sequence, specify the target locations
seq1=ones(15,1)*4; %100 loc 4
seq2=[3 3 2 3 3 4 3 1 3 3 4 2 3 3 1]'; %60 loc 3
seq3=[1 1 1 2 1 1 1 1 4 1 1 1 3 1 1]'; %80 loc 1
seq4=ones(15,1)*2; %100 loc 2
seq5=[3 4 3 3 3 3 2 3 3 3 1 3 3 3 3]'; %80 loc 3
seq6=[4 1 4 4 3 4 4 1 2 4 4 3 4 3 4]'; %60 loc 4
seq7=[1 1 3 1 1 1 1 2 1 1 1 1 4 1 1]'; %80 loc 1
seq8=[2 2 1 3 2 2 3 2 4 2 2 3 2 2 4]'; %60 loc 2
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
% Specify, for each sequence, if P is 100, 80 or 60.
ptype = [100, 60, 80, 100, 80, 60, 80, 60, 100, 60, 80, 100, 80, 60, 80, 60];

%take the sequence size
[nrows,ncols]=size(sequence);

allpred=[];
for k=2:2:ncols %n columns of sequence
    pred=mode(sequence(:,k));
    for n=1:nrows % of sequence
        if n == 1
            seqpred(n)=0;
        else
            if sequence(n,k)==pred
            seqpred(n)=1;
            else
            seqpred(n)=0;
            end
        end
    end
    allpred=[allpred;seqpred];
    prevpred=pred;
end
allpred=allpred';

% the target locations vary depending on participant number, which MAYBE is
% indicated by e. It is surely indicate by subjectnumber
subjnumber = load(char(strcat('C:/Users/u661121/Desktop/StatLearn\infantsdata\',fold(e).name,'\subjnum.mat')));
subjnumber = subjnumber.subjnum;
subjnumber = str2double(subjnumber);
for r=2:2:ncols %number of sequences
    for n = 1:nrows
        if ismember(subjnumber,2:4:1000)
            if sequence(n,r)==1
            sequence(n,r)=4;
            elseif  sequence(n,r)==2
            sequence(n,r)=3;
            elseif sequence(n,r)==3
            sequence(n,r)=1;
            elseif sequence(n,r)==4
            sequence(n,r)=2;
            end
        elseif ismember(subjnumber,3:4:1000)
            if sequence(n,r)==1
            sequence(n,r)=2;
            elseif  sequence(n,r)==2
            sequence(n,r)=1;
            elseif sequence(n,r)==3
            sequence(n,r)=4;
            elseif sequence(n,r)==4
            sequence(n,r)=3;
            end
        elseif ismember(subjnumber,4:4:1000)
            if sequence(n,r)==1
            sequence(n,r)=3;
            elseif  sequence(n,r)==2
            sequence(n,r)=4;
            elseif sequence(n,r)==3
            sequence(n,r)=2;
            elseif sequence(n,r)==4
            sequence(n,r)=1;
            end
        end
    end
end

% If I interrupted the exp before the last sequence, I reduce the sequence length

% lastseen = load(char(strcat('\\CNAS.RU.NL\U661121\Desktop\S_StatLearning\infantsdata\',fold(e).name,'\lastseen.mat')));
% lastseen = lastseen.lastseen(1);
% 
lastseen = ncols/2;
% sequence=sequence(:,1:lastseen*2);
% ptype=ptype(:,1:lastseen);
% %and also the lasttrial length
% lasttrial=lasttrial(1:lastseen,:);
% %and also the stamps length
% stamps=stamps(:,[1:lastseen ncols/2+1:ncols/2+lastseen]);
%---------------------------------------
% All the measures of information theory (ITs) - this must change if we
% implement also memory decay!
% I estimate ITs based on what they have seen. Moreover I can make two
% types of estimates: as if they whatched a single sequence or as if they
% watched many different sequences.

% what they have seen depends on postboris (the data analyzed after the
% handcoding of infants' LA).

% thisiall=iall(:,1:lastseen);
% thisDall=Dall(:,1:lastseen);
% thisLRall=LRall(:,1:lastseen);
% thisHall=Hall(:,1:lastseen);

%---------------------------------------

% I take the sequence size for thissequence
% [nrows,ncols]=size(sequence);

% check whether all the stamps are there. they must be 32 rows (16 for the
% cues followed by 16 for the targets)
if length(stamps)<32
    % if stamps are missing are usually the ones relative to targets for
    % sequences that have not been shown, while the ones for cues are still
    % there. so check how many rows of target stamps are missing:
    missingstamps = 32-length(stamps);
    % now add them. this is only to have the matrix of the correct size,
    % but firstlastbor will then take care of the fact that actually there
    % were no timestamps and the infant was actually not watching.
    correspcue=[16-missingstamps+1:16];
    stamps=[stamps, stamps(:,correspcue)];
end

% when an infant looks at every trial, lasttrial is not recorded, but it is
% nrows of the sequence
for n=1:ncols/2
     if lasttrial(n,1)==0 && stamps(1,n+16)>0 %*
         lasttrial(n,1)=nrows;
     end
end
%* there are actually two reasons why you can have a zero: the second is
%that the sequence changed before the first target is showed. In that case,
%the stamp for the first target is negative

% cut the sequences where participants stopped watching and make 2 2cols 
% vectors: one for type and location; one for cue and target onset.

thissequence=[];
thisstamps=[];
% thisi=[];
% thisD=[];
% thisLR=[];
% thisH=[];
for n=1:2:ncols
    if n==1
        m=1;
    else
        m=m+1;
    end
    thisseq=sequence(1:lasttrial(m,1),n:n+1);
    thissequence=[thissequence; thisseq];

    thisstamp=stamps(1:lasttrial(m,1),[m m+ncols/2]); % stamps are composed by:
    % all cue stamps first; all target stamps after (in correct order
    % 1:end). thus to match cue stamps with their target stamps I have to index
    % for a col number that is equal to half of the number of total cols in
    % the sequence
    thisstamps = [thisstamps; thisstamp];
    
%     thissur = thisiall(1:lasttrial(m,1),m);
%     thiskl = thisDall(1:lasttrial(m,1),m);
%     thisrate = thisLRall(1:lasttrial(m,1),m);
%     thish = thisHall(1:lasttrial(m,1),m);
%     thisi = [thisi; thissur];
%     thisD = [thisD; thiskl];
%     thisLR = [thisLR; thisrate];
%     thisH = [thisH; thish];
end 

% add timestamps
for j=1:2
    m=1;
    for n=1:length(thisstamps)
        while thisstamps(n,j) > fix(m,2) && m < length(fix)
            m=m+1;
        end
        if fix(m,6)==2
            stamp = [fix(1,1) thisstamps(n,j) NaN fix(m,4:5) j+10 thissequence(n,j)];
        else
            stamp = [fix(1,1) thisstamps(n,j) NaN -1 -1 j+10 thissequence(n,j)];
        end
        fix = [fix(1:m-1,:); stamp; fix(m:end,:)];
    end
end

% use postboris to fix the actual LAs found with handcoding
% it requires fix and lasttrial as input and returns a vector similar to 
% the first column of lasttrial, where it estimates, for each sequence, 
% at what trial infants looked away.
firstlastbor=[];
firstlastbor=postboris(fix, stamps, lastseen, subjnumber);

% now i create the LA matrix required for R
%borismatrix = boris2R(firstlastbor);

% Maybe: check from gazedata if stamp is NaN and notify that here if stamp
% is -1
%-------------------------------------------------------------------------

% create AOIs as squares for the 4 target locations and for cue location.
% cue location = 1
% 100%target location = 1
% 80%target most likely location = 2
% 60%target most likely location = 4
% remaining target location = 3
% none of the above = -1

xb=230; % half of the side of the square
xs=230;
ys=230;
yb=230;
leftXleft = 384 - xb;
leftXright = 384 + xs;
rightXleft = 1536 - xs;
rightXright = 1536 + xb;
upYup = 216 - ys;
upYdown = 216 + yb;
downYup = 864 - yb;
downYdown = 864 + ys;
cueXleft = 960 - xs;
cueXright = 960 + xs;
cueYup = 540 - xs;
cueYdown = 540 + xs;

AOI=-1.*ones(length(fix),1);
for n=1:length(fix)
    if leftXleft<fix(n,4) && fix(n,4)<leftXright
        if upYup<fix(n,5) && fix(n,5)<upYdown
            AOI(n)=2;
        elseif downYup<fix(n,5) && fix(n,5)<downYdown
            AOI(n)=3;
        end
    elseif rightXleft<fix(n,4) && fix(n,4)<rightXright
        if upYup<fix(n,5) && fix(n,5)<upYdown
            AOI(n)=4;
        elseif downYup<fix(n,5) && fix(n,5)<downYdown
            AOI(n)=1;
        end
    elseif cueXleft<fix(n,4) && fix(n,4)<cueXright
        if cueYup<fix(n,5) && fix(n,5)<cueYdown
            AOI(n)=0;
        end
    elseif isnan(fix(n,4)) || isnan(fix(n,5)) %% last added
        AOI(n)=NaN;                           %% last added 
    end
end
% add AOIs vector to fix as 8th column
fix = [fix AOI];

% ------------------------------------------------------------------------
% Extraction of Saccadic Latencies, Dwell times and Look-Away prob.
% ------------------------------------------------------------------------

% Remove the recording before the first cue onset
idx = find(fix(:,6)==11); %find cue onsets
idx = [1:idx(1)-1]; %take all the rows before first cue onset
fix(idx,:)=[]; %remove them

% with this loop i add fix(:,10) in which the trial number of each sequence
% is reported
ordinal=[];
lastn=[];
m=1;
for n=1:length(fix)
    if fix(n,7)>100
        ordinal=[ordinal;fix(n,7)];
        if length(ordinal)==1
           fix(n,10)=m;
        else
            if fix(n,7)==ordinal(end-1)
                m=m+1;
                fix(n,10)=m;
            else
                % I save in lastn the idx of last trial of every sequence
                lastn = [lastn; n];
                m=1;
                fix(n,10)=m;
            end
        end
    end
end

% adjust fix based on boris. Eliminate trials in which infants were not
% watching
% initialize for later 
ftocell = repmat(-111,16,2);
% Some sequences might missing. e.g., they go from 105 to 107 skipping 106.
% Create 2 vectors that specify the cue type and the sequence number for
% this particular subject.
tt=0;
tracktype=[];
for n=1:length(thissequence)-1
    if thissequence(n,1)~=thissequence(n+1,1)
        tt=tt+1;
        tracktype = [tracktype;thissequence(n,1),tt];
    end
end
    % this loop takes care of the last sequence (which cannot be compared
    % to the following one)
if thissequence(end,1)~=tracktype(end,1)
    tracktype = [tracktype;thissequence(end,1),tt+1];
end
for n=1:length(firstlastbor) %for every sequence
    % if the sequence is valid (where 111 equals invalid)
    if firstlastbor(n,1)<16
        if n<9
            % I get the sequence number for this specific subject
            idtrack=find(tracktype(:,1)==100+n);
            idtrack=tracktype(idtrack(1),2);
        else
            idtrack=find(tracktype(:,1)==100+n-8);
            idtrack=tracktype(idtrack(end),2);
        end
            % and I find the row associated to its starting point
        idseq=find(fix(:,10)==1);
        thisidseq=idseq(idtrack);

        % end when it ends (one row before beginning next one)
        if idtrack==tracktype(end,end)
            thisidseqend=length(fix);
        else
            thisidseqend=idseq(idtrack+1)-1;
        end
        % correct beginning based on boris
        idstart = firstlastbor(n,1);
        idstart = find(fix(thisidseq:thisidseqend,10)==idstart);
        idstart = thisidseq-1+idstart;
        % correct end based on boris.
        % First I find the beginning of the last boris trial
        idend = firstlastbor(n,2);
        idend = find(fix(thisidseq:thisidseqend,10)==idend);
        idend = thisidseq+idend(1)-1;
        % but now I have to find the end of the last boris trial
        % It corresponds to beginning +4 seconds, roughly. But it might
        % extend to beginning of the following trial, so I add 1.5 seconds
        % of the next cue to it. so it is 5.5 seconds, 5500 ms
        tend=fix(idend,2)+5500;
        idend = find(fix(:,2)<tend);
        idend = idend(end);
        % I sign fix11 as 1 if boris says they were looking
        fix(idstart:idend,11)=1;
        ftocell(n,1:2)=[idstart, idend];
    %else % if the sequence is invalid (111) then no 1 must be added.
    end
end

% So now I have two preliminary columns to look at before extracting every
% measurement, fix9 and fix11. Fix11 must be one because that's the part of
% fixations in which there is no look away and at the same time there a
% sequence that is being shown. Fix9 must be 1 because when it's zero the
% eye tracker did not detect the look to the target, and so estimation of
% eye behaviors is impossible (NaN)

% Let's take care of fix11 now.
% we can change the structure of fix, so that we have one cell for every
% sequence that have been seen. 
fixa=fix; % to save fix
fix={};
for n = 1:length(ftocell)
    if ftocell(n)>0
       fix{n}=fixa(ftocell(n,1):ftocell(n,2),:);
    end
end

% For every trial, I check whether the participant looked at the target at
% any moment between the cue offset and the target offset. If not, it means
% that s/he was not paying attention or that data are missing, so the trial
% cannot be analyzed
for w = 1:length(fix)
    tfix = fix(w);
    tfix = tfix{1,1};
    [lengthtfix,~]=size(tfix);
    if isempty(tfix)==0
    for k = 1:lengthtfix-1
    if tfix(k,6)==11
        % Given the current cue onset, I look for the next cue onset
        c = find(tfix(k+1:end,6)==11);
        if isempty(c) % the very last trial is not followed by another one, this fixes the problem
            [lengthtfix,~]=size(tfix);
            c = lengthtfix-k+1;
        end
        c = c(1)+k-1;
        a = tfix(k,2)+1000; %1000 is the cue length in ms. a is cue offset time
        m = k;
        % with the following loop I look for the first fixation that is after cue offset
        for n=k:c
            while a > tfix(m,2) && m < lengthtfix
                m=m+1;
            end
        end
        % now I check if in the range that goes from ((first fixation after cue 
        % offset (m))) target onset to time of next cue onset (excluded) (c), there is at
        % least one trial in which participant looked at the correct target
        % location
        b = find(tfix(m:c,6)==12); %find the moment of target onset within current range
        b = b+m-1; % row of the target location
        if length(b) > 1
            disp('something is wrong: there is more than 1 target marker in each trial')
            subjnumber
        end
        % create a new col and set 1 if participant looked at target
        if any(tfix(b:c,8)==tfix(b,7))
           tfix(k:c,9)=1;
        % if they have never looked at the target:
        % They cannot know what was going on if they were not looking, so
        % it's reasonable to discard the whole trial from the sequence. the
        % sequence will be tailored for the participant, and info theory
        % constructs will be computed on the sequence that has been seen
        % rather than on the standard sequence. This is implemented outside
        % the loop and can be commented out.
        end
    end
    end
    fix{w}=tfix;
    end
end

% Sometimes when eye-tracker looses the signal, it returns an impossible
% value for coordinates. This would result as -1 in column 8, i.e., the
% person was looking at a location different from the cue or target areas.
% However that's not correct. we thus add a new column (12) which specifies
% whether the -1 in row 8 should actually be considered as a NaN.

for w=1:length(fix)
    tfix = fix(w);
    tfix = tfix{1,1};
    if isempty(tfix)==0
        idxnanx = find(tfix(:,4)>1920);
        idxnany = find(tfix(:,5)>1080);
    end
    tfix(idxnanx,12)=1;
    tfix(idxnany,12)=1;
    fix{w}=tfix;
end


% Find saccadic latency to the target, where negative numbers (+ rt) mean
% anticipatory looking and positive (+ rt) mean saccadic reaction.
% rt is the time needed to initialize a saccade, so that if the saccade has
% been initialized before target appearence but it ends after target
% appearence, it is considered anticipatory anyway.
% Check whether participant is looking at correct target location when
% target appears. If so, go back until this is not true. That's the moment
% of the AL. If not, go forward until you find it. That's the moment of the
% saccadic reaction.

% measure saccadic latency of each trial
wslat={}; 
for w = 1:length(fix)
    prevslat = [];
    tfix = fix(w);
    tfix = tfix{1,1};
    [lengthtfix,~]=size(tfix);
    if isempty(tfix)==0
    for k = 1:lengthtfix-1
    if tfix(k,6)==11
        if tfix(k,9)==0 % If there is no look to target I cannot measure SL
            slat = NaN;
        elseif tfix(k,9)==1            
            % Given the current cue onset, I look for the next target onset
            t = find(tfix(k+1:end,6)==12);
            t = t(1)+k;
            % I check whether the person is already looking there at target
            % onset. If so, I go back, find when started looking there and take
            % the time. If not, I go forward and see when started to look there
            m=t;
            if tfix(t,7)==tfix(t,8)
                while tfix(t,7)==tfix(m,8)
                        m=m-1;
                end
                % if the row before the first look to target is a NaN, we
                % can't know whether the first look actually happened at
                % that point in time or before. If going even more back,
                % after a NaN, there is another look to the target, that is
                % the actual one.
                if isnan(tfix(m,8)) % if there is a NaN
                    mn=m;
                    % go to the beginning of NaNs
                    while isnan(tfix(mn,8)) && mn+1>k
                        mn=mn-1;
                    end
                    % check if the subject was already watching at the
                    % target
                    if tfix(t,7)==tfix(mn,8)
                        while tfix(t,7)==tfix(mn,8)
                            mn=mn-1;
                        end
                        m=mn;
                        slat=tfix(m+1,2)-tfix(t,2); %m+1 is the first trial in which subj looks at target
                    else
                        slat=NaN;
                    end
                else % if there is not a NaN
                    slat=tfix(m+1,2)-tfix(t,2); %m+1 is the first trial in which subj looks at target
                end
            else
                while tfix(t,7)~=tfix(m,8)
                        m=m+1;
                end   
            slat=tfix(m,2)-tfix(t,2); %in this case m is first trial in which subj looks at target
            end
        end
        prevslat = [prevslat; slat];
    end
    end
    end
    wslat{w}=prevslat;
end

%fix the length of sequences based on boris. This mismatch happens because
%some boris sequences end with a cue and so when looking for slats they are
%reported as NaNs while they are actually non existent.
%Also, arrange them in a single row.
prevslat=[];
for w=1:length(firstlastbor)
    if firstlastbor(w,1)<16 && firstlastbor(w,2)-firstlastbor(w,1)+1>3 %if it's not empty and if there are more than 3 trials
        % take the corresponding wslat and set the lastbor trial as last trial
        thisslat=wslat{w};
        lastbor=firstlastbor(w,2)-firstlastbor(w,1)+1;
        thisslat=thisslat(1:lastbor);
        prevslat=[prevslat;thisslat];
    end
end
        

% Extract Anticipatory Lookings from saccadic latencies.
wAL={};
tfix=[];
for w = 1:length(fix)
    prevalook = [];
    tfix = fix(w);
    tfix = tfix{1,1};
    [lengthtfix,~]=size(tfix);
    if isempty(tfix)==0
    for k = 1:lengthtfix-1
    if tfix(k,6)==11
        if tfix(k,9)==0 % If there is no look to target I do not want to measure Anticipatory look
            alook = [NaN, NaN]; %alook has two columns: first specify whether present absent, second whether correct or wrong
        elseif tfix(k,9)==1            
            % Given the current cue onset, I look for the next target onset
            t = find(tfix(k+1:end,6)==12);
            t = t(1)+k;
            % I check whether the person looked at any target position
            % after the cue offset. First I find an interval with all
            % saccades between cue onset + 375 and target onset+200 ms
            % first: cue onset +375 ms
            idx375 = tfix(k,2)+1000;
            n375 = find(tfix(:,2)>idx375);
            n375=n375(1);
            %second: target onset+200 ms
            idx200 = tfix(t,2)+200;
            n200 = find(tfix(:,2)<idx200);
            n200=n200(end);
            
            % now check which is the first target position that participant
            % looked at
            targettofind = 1;
            m=n375;
            while targettofind && m<n200
                if tfix(m,8)~= 0 && tfix(m,8)~= -1
                    targettofind=0;
                end
                m=m+1;
            end
            % if the target has been found, this is an anticipatory trial
            % and targettofind is now = 0. in that case the current alook
            % is 1 (present). otherwise it is 0 (absent)
            if targettofind == 0
                alook = 1;
            else
                alook = [0 NaN];
            end
            % and whether is matches the target position of this trial. In
            % that case the anticipation is correct. Otherwise it is wrong.
            % if alook = 1, check if correct or wrong
            if alook == 1
                if tfix(m,8)==mode(sequence(:,w*2))
                    alook = [tfix(m,8), 1]; %tfix(m,8)
                else
                    alook = [tfix(m,8), 0];
                end
            end
       
        end
     prevalook = [prevalook; alook];
    end
    end
    end
    wAL{w}=prevalook;
end

%fix the length of sequences based on boris. This mismatch happens because
%some boris sequences end with a cue and so when looking for slats they are
%reported as NaNs while they are actually non existent.
%Also, arrange them in a single row.
preAL=[];
for w=1:length(firstlastbor)
    if firstlastbor(w,1)<16 && firstlastbor(w,2)-firstlastbor(w,1)+1>3 %if it's not empty and if there are more than 3 trials
        % take the corresponding wslat and set the lastbor trial as last trial
        thisAL=wAL{w};
        lastbor=firstlastbor(w,2)-firstlastbor(w,1)+1;
        thisAL=thisAL(1:lastbor,:);
        preAL=[preAL;thisAL];
    end
end

% Measure dwell time on each trial
wdwell={};
tfix=[];
for w = 1:length(fix)
    predwell = [];
    tfix = fix(w);
    tfix = tfix{1,1};
    [lengthtfix,~]=size(tfix);
    if isempty(tfix)==0
    for k = 1:lengthtfix-1
    if tfix(k,6)==11
        if tfix(k,9)==0 % If there is no look to target I cannot measure dwell
            dwell = NaN;
        else
        % Given the current cue onset, I look for the next target onset
        t = find(tfix(k+1:end,6)==12);
        t = t(1)+k;
        % I check whether the person is already looking there at target
        % onset. If so, I go back, find when started looking there and take
        % the time. If not, I go forward and see when started to look there        
        m=t;
        n=t;
        if tfix(t,7)==tfix(t,8)
            while tfix(t,7)==tfix(m,8) && m > 0
                m=m-1;
            end
            while tfix(t,7)==tfix(n,8) && n<lengthtfix
                n=n+1;
            end
            dwell=tfix(n-1,2)-tfix(m+1,2);
         else
            while tfix(t,7)~=tfix(m,8) && m<lengthtfix
                m=m+1;
            end
            n = m;
            while tfix(t,7)==tfix(n,8) && n<lengthtfix
                n=n+1;
            end
            dwell=tfix(n-1,2)-tfix(m,2); % n is the first saccade that is not 
            %toward the target. n-1 is the end of the last sacc toward target
         end
       end
       predwell = [predwell; dwell];
    end
    end
    end
    wdwell{w}=predwell;
end  

%fix the length of sequences based on boris. This mismatch happens because
%some boris sequences end with a cue and so when looking for slats they are
%reported as NaNs while they are actually non existent.
%Also, arrange them in a single row.
predwell=[];
for w=1:length(firstlastbor)
    if firstlastbor(w,1)<16 && firstlastbor(w,2)-firstlastbor(w,1)+1>3 %if it's not empty and if there are more than 3 trials
        % take the corresponding wslat and set the lastbor trial as last trial
        thisdwell=wdwell{w};
        lastbor=firstlastbor(w,2)-firstlastbor(w,1)+1;
        thisdwell=thisdwell(1:lastbor);
        predwell=[predwell;thisdwell];
    end
end

% Measure lingering on each trial
wlinger={};
tfix=[];
for w = 1:length(fix)
    prelinger = [];
    tfix = fix(w);
    tfix = tfix{1,1};
    [lengthtfix,~]=size(tfix);
    if isempty(tfix)==0
    for k = 1:lengthtfix-1
    if tfix(k,6)==11
        if tfix(k,9)==0 % If there is no look to target I cannot measure dwell
           linger = NaN;
        else
        % Given the current cue onset, I look for the next target onset
        t = find(tfix(k+1:end,6)==12);
        t = t(1)+k;
        % And I determine the target offset
        idx = tfix(t,2)+1500; % 1500 is the length of target presentation
        % I also find the next cue onset
        n = find(tfix(k+1:end,6)==11);
        if isempty(n) % the very last trial is not followed by another one, this fixes the problem
            n = lengthtfix;
        else
            n = n(1)+k;
        end
        % And I add 375 ms and find the closest row
        idxn = tfix(n,2)+375;
        n = find(tfix(:,2)<idxn);
        n=n(end);
        % The last moment in which participant is looking
        % at target must be in a window between this target onset and 375ms
        % after next cue onset (to give enough time to go back to cue
        % location, but not enough to make an anticipatory saccade toward
        % the next target location)
        linger = find(tfix(t:n,8)==tfix(t,7));
        linger = tfix(linger(end)+t-1,2)-idx;
        end
        prelinger = [prelinger; linger];
    end
    end
    end
    wlinger{w}=prelinger;
end

%fix the length of sequences based on boris. This mismatch happens because
%some boris sequences end with a cue and so when looking for slats they are
%reported as NaNs while they are actually non existent.
%Also, arrange them in a single row.
prelinger=[];
for w=1:length(firstlastbor)
    if firstlastbor(w,1)<16 && firstlastbor(w,2)-firstlastbor(w,1)+1>3 %if it's not empty and if there are more than 3 trials
        % take the corresponding wslat and set the lastbor trial as last trial
        thisling=wlinger{w};
        lastbor=firstlastbor(w,2)-firstlastbor(w,1)+1;
        thisling=thisling(1:lastbor);
        prelinger=[prelinger;thisling];
    end
end
      
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%All that has been computed so far has been computed irrespective of how
%many trials in a sequence the infant has looked, but we want to discard
%sequences in which infants looked for less then 4 trials, because
% 1) it's not enough trials to have surprise- or kl-based expectations
% 2) computation of learning rate is imprecise for the first three trials

% We have thissequence, thisi, thisstamps, anticipatorytrials, predwell,
% prelinger and preslat which have to be adjusted.

% They have all the same length, so we find rows of sequences where
% ntrials<4 and we delete them.

% Those sequences are reported in firstlastbor

%------------------------------------------------------------------------    

% and we add them to the data of the other participants
% most of them has ONE row for every participant
 allfix{e}=fix; % FOR FIX ALSO SEQUENCES WITH NTRIAL<4 ARE KEPT
 allslat{e}=prevslat;
 allAL{e}=preAL;
 alldwell{e}=predwell;
 alllinger{e}=prelinger;
 allboris{e} = firstlastbor;
% create a sequence with only expected trials. I want to see whether
% saccadic latencies, dwell, lingering decrese with time, but this is true
% only when trials will be predictable. Conversely, AL can be higher even
% when the trial is unpredictable (because it is measured before t onset).
if length(prevslat)>19
    markasgood=[markasgood;subjnumber];
end

end

% For Look-Away survival analysis, data must be rearranged. This can be
% done with the boris2R function. The input is allboris, the output a
% rearranged version of it (see boris2R for details)

% first, keep only subjs that looked at least 20 trials
alln=[];
for n = 1:length(allslat)
    if length(allslat{n})>19
        alln=[alln,n];
    end
end
Roris=boris2R(allboris(alln),ptype,allpred);



if plotting
    Hd = quantileranks(Roris(:,7), 10); 
    Id = quantileranks(Roris(:,8), 10); 
    Dd = quantileranks(Roris(:,9), 10); 
    allpLA=[];
    for k=1:10
        idx=find(Hd==k);
        pLA= sum(Roris(idx,10))/length(idx);
        allpLA=[allpLA;pLA];
    end
    
figure    
plot(allpLA,'o')

[b,~,~,~,stats] = regress(Roris(:,7),[Roris(:,9),ones(length(Roris),1)]);
% plot(Roris(:,7),Roris(:,8),'o') %surprise and entropy
% plot(Roris(:,7),Roris(:,9),'o') %KL and entropy
% plot(Roris(:,9),Roris(:,8),'o') %KL and surprise


end

AL=[];
number=[];
for n=1:length(allAL)
    if length(allAL{n})>19
        thisAL=allAL{n};
        number=[number;n];
        
        
        AL=[AL;thisAL];
    end
end

notransf=1;
slat=[];
for n=1:length(allslat)
    if length(allslat{n})>19
        if notransf
            zslat=allslat{n};
        else
        thisslat=allslat{n};
        meanslat= nanmean(thisslat);% mean allslat taking into account NaN
        sdslat = nanstd(thisslat); % SD allslat taking into account NaN

        difslat = thisslat-meanslat;
        zslat=difslat./sdslat;
        end
        slat=[slat;zslat];
    end
end

linger=[];
for n=1:length(allslat)
    linger=[linger;alllinger{n}];
end

linger=[];
for n=1:length(alllinger)
    if length(alllinger{n})>19
        if notransf
            zlinger=alllinger{n};
        else
        thislinger=alllinger{n};
        meanslat= nanmean(thislinger);% mean allslat taking into account NaN
        sdslat = nanstd(thislinger); % SD allslat taking into account NaN

        difslat = thislinger-meanslat;
        zlinger=difslat./sdslat;
        end
        linger=[linger;zlinger];
    end
end

dwell=[];
for n=1:length(alldwell)
    if length(alllinger{n})>19
        if notransf
            zdwell=alldwell{n};
        else
            
        thisdwell=alldwell{n};
        meandwell= nanmean(thisdwell);% mean allslat taking into account NaN
        sddwell = nanstd(thisdwell); % SD allslat taking into account NaN

        difdwell = thisdwell-meandwell;
        zdwell=difdwell./sddwell;
        end
        dwell=[dwell;zdwell];
    end
end

% plot(Roris(:,9),dwell,'o')

% add slat, linger and dwell to Roris
Roris=[Roris, AL, slat, linger, dwell];

% remove outliers
for k=11:13
    Roris(:,k)=filloutliers(Roris(:,k),NaN,'mean');
end

% trasform IT measures in zscores, for every subject
ztrasform=0;
if ztrasform
    for n=1:max(unique(Roris(:,2)))
    idxs=find(Roris(:,2)==n);
    Roris(idxs,7)=zscore(Roris(idxs,7));
    Roris(idxs,8)=zscore(Roris(idxs,8));
    Roris(idxs,9)=zscore(Roris(idxs,9));
    Roris(idxs,11)=zscore(Roris(idxs,11));
    end
end

boris=[];
for n=1:length(allboris)
    %if length(alllinger{n})>19
    thisbor=allboris{n};
    for k=1:length(thisbor)
        if thisbor(k,1)~=111
            newborseq=[n,k,thisbor(k,1),thisbor(k,2),thisbor(k,3)];
            boris=[boris;newborseq];
        end
    end
end

if saving
    dlmwrite('Roris_nostd.csv',Roris);
    dlmwrite('LookAwayData.csv',boris);
end
