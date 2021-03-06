% -----------------------------------------------------------------
% adult processing
% ----------------------------------------------------------------- 
%clearvars -except participant; close all; fclose('all'); clc;
%dbstop if error;
%commandwindow;

%switch to 0 if you do not want to save data 

saving = 0;

cd 'C:/Users/u661121/Desktop/StatLearn'

% import the I2MC fixations files from each participant
folders.data = 'output'; % the input of this algorithm is the output of I2MC
folders.output = 'input'; % folder for output (will use structure in folders.data for saving output)

% create textfile and open for writing analyzed data (with saccades, dwell 
% times, look-away probability and anticipatory looking)
%fid = fopen(fullfile(folders.output,'allfix.txt'),'w');
%fprintf(fid,'participant\time\duration\xPos\yPos\fixOrStamp\typeAndPosition\aOIs\lookToTarget');

[fold,nfold] = FolderFromFolder(folders.data);

% Inizialize the output
allfix=[];
allslat=[];
allAL = [];
alldwell=[];
alllinger=[];
allLA=[];
predictslat=[];
predictdwell=[];
predictlinger=[];

for e=1:nfold

%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/giuliaserino/Desktop/adulti /input/Adult26/allfixations.txt
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2019/05/09 16:47:23

%% Initialize variables.
filename = char(strcat('C:/Users/u661121/Desktop/StatLearn/output/',fold(e).name,'/allfixations.txt'));
delimiter = '\t';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
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


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
allfixations = cell2mat(raw);
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
%-------------------------------------------------------------------------    
% Get only what you need from the I2MC output (i.e. participant number,
% beginning and end of saccades, their duration and location)
fixations = [allfixations(:,end-1) allfixations(:,1:5)];
%stamps = load(char(strcat('/Users/giuliaserino/Desktop/adulti /data/',fold(e).name,'/stampstimuli.mat')));
stamps = load(char(strcat('//CNAS.RU.NL/U661121/Desktop/S_StatLearning/data/',fold(e).name,'/stampstimuli.mat')));
stamps = stamps.stampstimuli;
stamps = [stamps(:,[1,4]); stamps(:,[2,5]); stamps(:,[3,6])];

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

% Load the sequence of stimuli that contains stimulus type (col 1) and
% stimulus location (col 2)
sequence1 = load('\\CNAS.RU.NL\U661121\Desktop\S_StatLearning/sequence1mod.mat');
sequence1=sequence1.sequence;
sequence2 = load('\\CNAS.RU.NL\U661121\Desktop\S_StatLearning/sequence2mod.mat');
sequence2=sequence2.sequence;
sequence3 = load('\\CNAS.RU.NL\U661121\Desktop\S_StatLearning/sequence3mod.mat');
sequence3=sequence3.sequence;

sequence2(:,1)=sequence2(:,1)+3;
sequence3(:,1)=sequence3(:,1)+6;
ind3=find(sequence3(:,1)==109);
sequence3(ind3,1)=103;
ind3=find(sequence3(:,1)==108);
sequence3(ind3,1)=104;
ind3=find(sequence2(:,1)==104);
sequence2(ind3,1)=108;
sequence = [sequence1, sequence2, sequence3];

%take the sequence size
[nrows,ncols]=size(sequence);

% the target locations vary depending on participant number, which MAYBE is
% indicated by e. It is surely indicate by subjectnumber
subjnumber = load(char(strcat('//CNAS.RU.NL/U661121/Desktop/S_StatLearning/data/',fold(e).name,'/subjnum.mat')));
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
sequence = [sequence(:,1:2); sequence(:,3:4); sequence(:,5:6)];

%INSERT HERE COMPUTATIONS BASED ON MODEL

% add timestamps     
for j=1:2
    m=1;
    for n=1:length(stamps)
        while stamps(n,j) > fix(m,2) && m < length(fix)
            m=m+1;
        end
        if fix(m,6)==2
            stamp = [fix(1,1) stamps(n,j) 0 fix(m,4:5) j+10 sequence(n,j)];
        else
            stamp = [fix(1,1) stamps(n,j) 0 -1 -1 j+10 sequence(n,j)];
        end
        fix = [fix(1:m-1,:); stamp; fix(m:end,:)];
    end
end

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
    end
end
% add AOIs vector to fix as 8th column
fix = [fix AOI];

% ------------------------------------------------------------------------
% Extraction of Saccadic Latencies, Dwell times and Look-Away prob.
% ------------------------------------------------------------------------

% For every trial, I check whether the participant looked at the target at
% any moment between the cue offset and the target offset. If not, it means
% that s/he was not paying attention or that data are missing, so the trial
% cannot be analyzed

todiscard = []; %initialize for later

for k = 1:length(fix)-1
    if fix(k,6)==11
        % Given the current cue onset, I look for the next cue onset
        c = find(fix(k+1:end,6)==11);
        if isempty(c) % the very last trial is not followed by another one, this fixes the problem
            c = length(fix)-k+1;
        end
        c = c(1)+k-1;
        a = fix(k,2)+750; %750 is the cue length in ms. a is cue offset time
        m = k;
        % with the following loop I look for the first fixation that is after cue offset
        for n=k:c
            while a > fix(m,2) && m < length(fix)
                m=m+1;
            end
        end
        % now I check if in the range that goes from ((first fixation after cue 
        % offset (m))) target onset to time of next cue onset (excluded) (c), there is at
        % least one trial in which participant looked at the correct target
        % location
        b = find(fix(m:c,6)==12); %find the moment of target onset within current range
        b = b+m-1; % row of the target location
        if length(b) > 1
            disp('something is wrong: there is more than 1 target marker in each trial')
        end
        % create a new col and set 1 if participant looked at target
        if any(fix(b:c,8)==fix(b,7))
           fix(k:c,9)=1;
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
% Remove the recording before the first cue onset
idx = find(fix(:,6)==11); %find cue onsets
idx = [1:idx(1)-1]; %take all the rows before first cue onset
fix(idx,:)=[]; %remove them

% Re-compute the sequence based on what has been seen
% % % % this is not necessary as it has already been done before based on lasttrial
% % % k=0;
% % % j=0;
% % % for n = 1:length(fix)
% % %     if fix(n,7)==101 || fix(n,7)==102 || fix(n,7)==103
% % %         k=k+1;
% % %         type(k,1)=fix(n,7);
% % %     elseif fix(n,7)==1 || fix(n,7)==2 || fix(n,7)==3 || fix(n,7)==4
% % %         j=j+1;
% % %         pos(j,1)=fix(n,7);
% % %         pos(j,2)=fix(n,9);
% % %     end
% % % end
% % % % now I have type of stimulus (col 1) position (col 2) and whether it was
% % % % looked at (col 3, where 1 is yes and 0 is no). 
% % % data = [type pos];
% % % % So I can:
% % % % extract the indexes of the trials that has not been missed, so that I can
% % % % keep only their surprise\KL when plotting the results
% % % survivors = data(:,3)==1;
% % % % Create a new sequence which contains only the stimuli that the
% % % % participant has seen. Attention: this implies that the absence of look to
% % % % the target is due to a real behavior of the participant rather than to
% % % % data loss.
% % % 
% % % % % %  data = data(data(:,3)==1); % % % COMMENTED OUT

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
prevslat = [];

for k = 1:length(fix)-1
    if fix(k,6)==11
        if fix(k,9)==0 % If there is no look to target I cannot measure SL
            slat = NaN;
        elseif fix(k,9)==1            
            % Given the current cue onset, I look for the next target onset
            t = find(fix(k+1:end,6)==12);
            t = t(1)+k;
            % I check whether the person is already looking there at target
            % onset. If so, I go back, find when started looking there and take
            % the time. If not, I go forward and see when started to look there
            m=t;
            if fix(t,7)==fix(t,8)
                while fix(t,7)==fix(m,8)
                        m=m-1;
                end
                % if the row before the first look to target is a NaN, we
                % can't know whether the first look actually happened at
                % that point in time or before. If going even more back,
                % after a NaN, there is another look to the target, that is
                % the actual one.
                % for adults there are no NaNs.
            slat=fix(m+1,2)-fix(t,2); %m+1 is the first trial in which subj looks at target
            else
                while fix(t,7)~=fix(m,8)
                        m=m+1;
                end
            slat=fix(m,2)-fix(t,2); %in this case m is first trial in which subj looks at target
            end
        end
        prevslat = [prevslat; slat];
    end
end

% Extract Anticipatory Lookings from saccadic latencies.
% CHANGE THIS! take it from newpreR
AnticipatoryTrials = prevslat<200; % if 1, the trial is anticipatory
prevalook = [];
check105 = 0;
for k = 1:length(fix)-1
    if fix(k,6)==11
        if fix(k,7)==105
            check105=1;
        end
        if fix(k,9)==0 % If there is no look to target I do not want to measure Anticipatory look
            alook = [NaN, NaN]; %alook has to column: first specify whether present absent, second whether correct or wrong
        elseif fix(k,9)==1            
            % Given the current cue onset, I look for the next target onset
            t = find(fix(k+1:end,6)==12);
            t = t(1)+k;
            % I check whether the person looked at any target position
            % after the cue offset. First I find an interval with all
            % saccades between cue onset + 375 and target onset+200 ms
            % first: cue onset +375 ms
            idx375 = fix(k,2)+375;
            n375 = find(fix(:,2)>idx375);
            n375=n375(1);
            %second: target onset+200 ms
            idx200 = fix(t,2)+200;
            n200 = find(fix(:,2)<idx200);
            n200=n200(end);
            
            % now check which is the first target position that participant
            % looked at
            targettofind = 1;
            m=n375;
            while targettofind && m<=n200
                if fix(m,8)~= 0 && fix(m,8)~= -1
                    targettofind=0;
                end
                m=m+1;
            end
            % if the target has been found, this is an anticipatory trial
            % and targettofind is now = 0. in that case the current alook
            % is 1 (present). otherwise it is 2 (absent)
            if targettofind == 0
                alook = 1;
            else
                alook = [0 NaN];
            end
        
            % and whether is matches the target position of this trial. In
            % that case the anticipation is correct. Otherwise it is wrong.
            % if alook = 1, check if correct or wrong
            if alook == 1
                if ismember(subjnumber,1:4:1000)
                        if fix(k,7)==101
                            corrtar = 3; %correct target location
                        elseif fix(k,7)==103 && check105==0
                            corrtar=1;
                        elseif fix(k,7)==102
                            corrtar=2;
                        elseif fix(k,7)==105
                            corrtar=3;
                        elseif fix(k,7)==106
                            corrtar=1;
                        elseif fix(k,7)==108
                            corrtar=4;
                        elseif fix(k,7)==103 && check105==1
                             corrtar=2;
                        elseif fix(k,7)==104
                            corrtar=1;
                        elseif fix(k,7)==107
                            corrtar=4;
                        end
                elseif ismember(subjnumber,2:4:1000)
                        if fix(k,7)==101
                            corrtar = 1; %correct target location
                        elseif fix(k,7)==103 && check105==0
                            corrtar=4;
                        elseif fix(k,7)==102
                            corrtar=3;
                        elseif fix(k,7)==105
                            corrtar=1;
                        elseif fix(k,7)==106
                            corrtar=4;
                        elseif fix(k,7)==108
                            corrtar=2;
                        elseif fix(k,7)==103 && check105==1
                             corrtar=3;
                        elseif fix(k,7)==104
                            corrtar=4;
                        elseif fix(k,7)==107
                            corrtar=2;
                        end
                elseif ismember(subjnumber,3:4:1000)
                        if fix(k,7)==101
                            corrtar = 4; %correct target location
                        elseif fix(k,7)==103 && check105==0
                            corrtar=2;
                        elseif fix(k,7)==102
                            corrtar=1;
                        elseif fix(k,7)==105
                            corrtar=4;
                        elseif fix(k,7)==106
                            corrtar=2;
                        elseif fix(k,7)==108
                            corrtar=3;
                        elseif fix(k,7)==103 && check105==1
                             corrtar=1;
                        elseif fix(k,7)==104
                            corrtar=2;
                        elseif fix(k,7)==107
                            corrtar=3;
                        end
                else 
                        if fix(k,7)==101
                            corrtar = 2; %correct target location
                        elseif fix(k,7)==103 && check105==0
                            corrtar=3;
                        elseif fix(k,7)==102
                            corrtar=4;
                        elseif fix(k,7)==105
                            corrtar=2;
                        elseif fix(k,7)==106
                            corrtar=3;
                        elseif fix(k,7)==108
                            corrtar=1; 
                        elseif fix(k,7)==103 && check105==1
                             corrtar=4;
                        elseif fix(k,7)==104
                            corrtar=3;
                        elseif fix(k,7)==107
                            corrtar=1;
                        end
                end
           if fix(m,8)==corrtar
           alook = [alook, 1];
           else
           alook = [alook, 0];
           end
        end
    end
    prevalook = [prevalook; alook];
end
end


% Measure dwell time on each trial
predwell = [];

for k = 1:length(fix)-1
    if fix(k,6)==11
        if fix(k,9)==0 % If there is no look to target I cannot measure dwell
            dwell = NaN;
        else
        % Given the current cue onset, I look for the next target onset
        t = find(fix(k+1:end,6)==12);
        t = t(1)+k;
        % I check whether the person is already looking there at target
        % onset. If so, I go back, find when started looking there and take
        % the time. If not, I go forward and see when started to look there        
        m=t;
       n=t;
        if fix(t,7)==fix(t,8)
            while fix(t,7)==fix(m,8) && m > 0
                m=m-1;
            end
            while fix(t,7)==fix(n,8) && n<length(fix)
                n=n+1;
            end
            dwell=fix(n-1,2)-fix(m+1,2);
         else
            while fix(t,7)~=fix(m,8) && m<length(fix)
                m=m+1;
            end
            n = m;
            while fix(t,7)==fix(n,8) && n<length(fix)
                n=n+1;
            end
            dwell=fix(n-1,2)-fix(m,2); % n is the first saccade that is not 
            %toward the target. n-1 is the end of the last sacc toward target 
         end
       end
       predwell = [predwell; dwell];
    end
end  

% Measure lingering on each trial
prelinger = [];
for k = 1:length(fix)-1
    if fix(k,6)==11
        if fix(k,9)==0 % If there is no look to target I cannot measure dwell
           linger = NaN;
        else
        % Given the current cue onset, I look for the next target onset
        t = find(fix(k+1:end,6)==12);
        t = t(1)+k;
        % And I determine the target offset
        idx = fix(t,2)+750; % 750 is the length of target presentation
        % I also find the next cue onset
        n = find(fix(k+1:end,6)==11);
        if isempty(n) % the very last trial is not followed by another one, this fixes the problem
            n = length(fix);
        else
            n = n(1)+k;
        end
        % And I add 375 ms and find the closest row
        idxn = fix(n,2)+375;
        n = find(fix(:,2)<idxn);
        n=n(end);
        % The last moment in which participant is looking
        % at target must be in a window between this target onset and 375ms
        % after next cue onset (to give enough time to go back to cue
        % location, but not enough to make an anticipatory saccade toward
        % the next target location)
        linger = find(fix(t:n,8)==fix(t,7));
        linger = fix(linger(end)+t-1,2)-idx;
        end
       
        prelinger = [prelinger; linger];
    end

end
%------------------------------------------------------------------------
% create a sequence with only expected trials. I want to see whether
% saccadic latencies, dwell, lingering decrese with time, but this is true
% only when trials will be predictable. Conversely, AL can be higher even
% when the trial is unpredictable (because it is measured before t onset).

%idx = [];
%or n = 1:length(sequence)
   % if sequence(n,1) == 101 && sequence(n,2) ~= 1
       % idx = [idx, n];
    %elseif sequence(n,1) == 102 && sequence(n,2) ~= 2
       % idx = [idx, n];
    %elseif sequence(n,1) == 103 && sequence(n,2) ~= 4
      %  idx = [idx, n];
   % end
%end
%predseq = sequence;
%predseq(idx,:)=[];

% Now take SL, dwell, linger and remove unexpected trials
%ppslat=prevslat;
%ppslat(idx)=[];

%ppdwell=predwell;
%ppdwell(idx)=[];

%pplinger=prelinger;
%pplinger(idx)=[];
%-------------------------------------------------------------------------
% estimate the cumulative number of anticipatory lookings and the learning
% rate from saccadic latencies. 

%-------------------------------------------------------------------------

 %allfix=[allfix;fix];
 allslat=[allslat; prevslat];
 %allAL = [allAL; prevalook]; % two rows for every participant
 alldwell=[alldwell; predwell];
 alllinger=[alllinger; prelinger];
 
% predictslat=[predictslat; ppslat'];
 %predictdwell=[predictdwell; ppdwell'];
 %predictlinger=[predictlinger; pplinger'];
end

%specify suubject number
allsubj=[];
for s = 1:e
    allsubj=[allsubj;ones(270,1)*s];
end

%run the ideal learner model for each sequence
[H1,i1,D1,Hpost1,data1]=ITDadultsmodel2(sequence1, 101);
[H2,i2,D2,Hpost2,data2]=ITDadultsmodel2(sequence1, 102);
[H3,i3,D3,Hpost3,data3]=ITDadultsmodel2(sequence1, 103);
% assemble them with the correct indexing
Hseq1=NaN(90,1);
Hseq1(data1(:,1))=H1;
Hseq1(data2(:,1))=H2;
Hseq1(data3(:,1))=H3;

iseq1=NaN(90,1);
iseq1(data1(:,1))=i1;
iseq1(data2(:,1))=i2;
iseq1(data3(:,1))=i3;

Dseq1=NaN(90,1);
Dseq1(data1(:,1))=D1;
Dseq1(data2(:,1))=D2;
Dseq1(data3(:,1))=D3;

Hpostseq1=NaN(90,1);
Hpostseq1(data1(:,1))=Hpost1;
Hpostseq1(data2(:,1))=Hpost2;
Hpostseq1(data3(:,1))=Hpost3;

%sequence 2
[H1,i1,D1,Hpost1,data1]=ITDadultsmodel2(sequence2, 105);
[H2,i2,D2,Hpost2,data2]=ITDadultsmodel2(sequence2, 106);
[H3,i3,D3,Hpost3,data3]=ITDadultsmodel2(sequence2, 108);
% assemble them with the correct indexing
Hseq2=NaN(90,1);
Hseq2(data1(:,1))=H1;
Hseq2(data2(:,1))=H2;
Hseq2(data3(:,1))=H3;

iseq2=NaN(90,1);
iseq2(data1(:,1))=i1;
iseq2(data2(:,1))=i2;
iseq2(data3(:,1))=i3;

Dseq2=NaN(90,1);
Dseq2(data1(:,1))=D1;
Dseq2(data2(:,1))=D2;
Dseq2(data3(:,1))=D3;

Hpostseq2=NaN(90,1);
Hpostseq2(data1(:,1))=Hpost1;
Hpostseq2(data2(:,1))=Hpost2;
Hpostseq2(data3(:,1))=Hpost3;

%sequence 3
[H1,i1,D1,Hpost1,data1]=ITDadultsmodel2(sequence3, 103);
[H2,i2,D2,Hpost2,data2]=ITDadultsmodel2(sequence3, 104);
[H3,i3,D3,Hpost3,data3]=ITDadultsmodel2(sequence3, 107);
% assemble them with the correct indexing
Hseq3=NaN(90,1);
Hseq3(data1(:,1))=H1;
Hseq3(data2(:,1))=H2;
Hseq3(data3(:,1))=H3;

iseq3=NaN(90,1);
iseq3(data1(:,1))=i1;
iseq3(data2(:,1))=i2;
iseq3(data3(:,1))=i3;

Dseq3=NaN(90,1);
Dseq3(data1(:,1))=D1;
Dseq3(data2(:,1))=D2;
Dseq3(data3(:,1))=D3;

Hpostseq3=NaN(90,1);
Hpostseq3(data1(:,1))=Hpost1;
Hpostseq3(data2(:,1))=Hpost2;
Hpostseq3(data3(:,1))=Hpost3;

% now put everything together
i = [iseq1;iseq2;iseq3];
H = [Hseq1;Hseq2;Hseq3];
Hpost = [Hpostseq1;Hpostseq2;Hpostseq3];
D = [Dseq1;Dseq2;Dseq3];

% now repeat all of them for each participant
i=repmat(i,e,1);
H=repmat(H,e,1);
Hpost=repmat(Hpost,e,1);
D=repmat(D,e,1);

% add time

Time90=repmat([1:90,1:90,1:90]',e,1);
Time270=repmat([1:270]',e,1);

% now put everything together:
data = [allsubj, allslat, alldwell, alllinger, i, H, Hpost, D, Time90, Time270];

dlmwrite('adult_statlearn.csv',data);
