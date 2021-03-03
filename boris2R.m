function borismatrix = boris2R(allboris,ptype,allpred) % postboris output
borismatrix = [];
allcountvalid=[]; %for later
for n=1:length(allboris)
    boris=allboris(n);
    boris=boris{:,:};
    countvalid=0; %for later
    newntrial=0; %for later
    for k=1:length(boris)
        % first of all, count how many trials the subject has already seen.
        % The number of seen trials is a good proxy of how much time the
        % infant has already spent looking at the video, and thus can be
        % added in the survival analysis as covariate to control for time
        if boris(k,1)~=111
            oldntrial=newntrial;
            newntrial=sum(boris(find(boris(1:k,2)~=111),2))-sum(boris(find(boris(1:k,1)~=111),1))+nnz(boris(1:k,1)~=111);
            thesetrials=[oldntrial+1:newntrial];
        end
        if boris(k,2)-boris(k,1)>2 %basically if I have at least 4 trials before LA
            countvalid=countvalid+1; %this counts how many valid sequences we have so far
            % which is also the sequence number
            newb=[boris(k,1):boris(k,2)';zeros(boris(k,2)-boris(k,1)+1,1)']';
            if boris(k,3)==1
                newb(boris(k,2)-boris(k,1)+1,2)=1; %Use this if you want a 1
                % only for the last trial
                %newb(:,2)=1; %otherwise use this
            end
        % Add the id number, sequence number and participant number n
        % THIS MIGHT CHANGE
        % to compute the id number, save the last countvalid of the
        % previous participant.
        newb(:,3)=sum(allcountvalid)+countvalid;
        newb(:,4)=n;
        newb(:,5)=countvalid;
        % rearrange newb:
        % col 1= id
        % col 2= participant number
        % col 3= sequence number
        % col 4= overall trial number
        % col 5= trial number in the sequence
        % col 6= LA in the sequence (yes=1)
        newb=[newb(:,3:5) thesetrials' newb(:,1:2)];
        %using first and last trial in allboris, and knowing the sequence
        %number, I want to estimate Information theory constructs.
        % I can do that with ITDmodel
        % ITD stands for Information Theory Dirichlet, bcecause it computes
        % information theory constructs starting from probabilities that
        % are updated following a dirichlet distribution
        % It requires sequence number and trial numbers (begin and end) as input
        [Hz,iz,Dz,Rz]=ITDmodel(boris(k,1),boris(k,2),k,0);
        newb=[newb(:,1:end-1),repmat(newb(end,5),length(newb(:,1)),1),Hz,iz,Dz,newb(:,end),Rz];
        % add also predictability (60, 80 or 100) and whether the trial is
        % predictable or unpredictable
        newb(:,end+1)=ptype(k);
        newb(:,end+1)=allpred(boris(k,1):boris(k,2),k);
        borismatrix=[borismatrix;newb];
        %countvalid.*ones(length(newb(:,1)),1)];
        end
    end
    
    allcountvalid=[allcountvalid, countvalid];
end

