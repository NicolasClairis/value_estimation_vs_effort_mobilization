%% designLearning75
% this script build the design matrix of the learning task with the
% following specifications:
%   - nt = 60 trials (24 gains / 24 loss / 12 neutrals)
%   - predictibility of cues: 75%
%   - valence & side of correct cue & predictibility orthogonalized 
%   - random permutation of the sequence of trials with the constraint of
%        decorrelated predictive sequences 

Valence = [ ones(1,24) repmat(3,1,24)]; % list of positive (1) and negative (3) pairs
% randomize order of each pair
indexValence = randperm(48);
npair_pn = Valence(indexValence); 

Side = [ repmat(-1,1,12) ones(1,12)]; % left/right side of the symbol
    % pseudo-randomization for each valence condition (prevent from random correlations)
    indexSide{1} = randperm(24);
    indexSide{2} = randperm(24);
    [r,p] = corr(indexSide{1}',indexSide{2}');
    while (r^2)>0.1
        indexSide{1} = randperm(24);
        indexSide{2} = randperm(24);
       [r,p] = corr(indexSide{1}',indexSide{2}');
    end
iPermut = 0;
for val = [1 3]
    iPermut = iPermut+1;
    side_pn(npair_pn == val) = Side(indexSide{iPermut}); 
end


Lottery = [ repmat(-1,1,3) repmat(1,1,9)];
     % pseudo-randomization for each valence condition (prevent from random correlations)
        indexLottery{1} = randperm(12);
        indexLottery{2} = randperm(12);
        indexLottery{3} = randperm(12);
        indexLottery{4} = randperm(12);
        [r,p] = corr([indexLottery{1}',indexLottery{2}',indexLottery{3}',indexLottery{4}']);
        r(1,1)=0; r(2,2)=0; r(3,3)=0; r(4,4)=0;
        while ~isempty(find(r.^2>0.2))
            indexLottery{1} = randperm(12);
            indexLottery{2} = randperm(12);
            indexLottery{3} = randperm(12);
            indexLottery{4} = randperm(12);
           [r,p] = corr([indexLottery{1}',indexLottery{2}',indexLottery{3}',indexLottery{4}']);
           r(1,1)=0; r(2,2)=0; r(3,3)=0; r(4,4)=0;
        end
iPermut = 0;
for val = [1 3]
    for side = [-1 1]
        iPermut = iPermut+1;
        lottery_pn(npair_pn==val & side_pn==side) = Lottery(indexLottery{iPermut});
    end
end

% Randomization for the neutral pair (12 trials)
% random & balanced side of cue (50%/50%)

npair_neu = 2*ones(1,12);  

Side = [ repmat(-1,1,6) ones(1,6)];
indexSide = randperm(12);
side_neu = Side(indexSide); 
lottery_neu = ones(1,12); 

% Create final trial vectors (2 consecutive blocks of 12 positive, 12 negative, and 6
% neutral pairs).
addneu = randperm(totaltrial,12);
side = nan(1,totaltrial);     side(addneu) = side_neu;         side(setdiff(1:totaltrial,addneu)) = side_pn;           % -1 = good on the left, 1 = good on the right
npair = nan(1,totaltrial);    npair(addneu) = npair_neu;       npair(setdiff(1:totaltrial,addneu)) = npair_pn;         % 1=gain 2=neutral 3=loss
lottery = nan(1,totaltrial);  lottery(addneu) = lottery_neu;   lottery(setdiff(1:totaltrial,addneu)) = lottery_pn;     % -1 = unlikely outcome (25%), % 1 = likely outcome (75%)

design = [side ; npair ; lottery ];
