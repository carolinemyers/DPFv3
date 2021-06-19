%% 00 Hello 

%This code was written by Caroline Myers, adapted from previous
%scripts by Mariel Roberts.

%%For CM collected data

%Purpose: to extract data from DPFv2 and put it in a structure for 
%further analysis. 

%Assumptions: Needs MGL to run, data (eye and mat) must be stored in 'data' by folder with observer
%name.  

%CORRECTED TO TAKE OUT NANS AND CORRECT NEGATIVES

%AMW2Comp = SubjectsStruct.s116PF_Exo_amw_S2  
%% 01 Init and Add to previous struct
%clear all
%close all

clearvars -except SubjectsStruct AMW2Comp
%SubjectsStruct = SubjectsStruct3
%% 02 Define vars and constants

subjects={'CL_Exo_Amb_Comb'};
SubNoString = '241'
Age = 24;
Heightin = 0;
Gender = 'm';
Session = 1;

Source = 'ExoAmbControl';
dataAvailability = 'eyebeh';
numLocations = 4;
Included = 1;
Notes = 'heightmissing';

ChildAdultCategory = 'adult';
ChildAdultNumeric = 3;

attentionCond = 'exo'

SubNo = str2num(SubNoString);
AAALastEditDate = '05.04.21'
%% 03 Pull out block info using MGL function

subj=1  

avg_num_testblocks = zeros(1,length(subjects));


Observer=subjects{subj};
Blocks=DPF_ExoAttnParseFiles(Observer);
num_blocks = size(Blocks,1);

for i=1:size(Blocks,1)
    if Blocks(i,2)<10
    eval(sprintf('load(''data/%s/%g_stim0%g.mat'')',Observer,Blocks(i,1),Blocks(i,2)));
    Data{i}=getTaskParameters(myscreen,task);
    Contrast(i)=stimulus.contrasts;
    NumFixBreak(i)=sum(stimulus.FixationBreak);
    fixBreak{i}=stimulus.FixationBreak;
    else
    eval(sprintf('load(''data/%s/%g_stim%g.mat'')',Observer,Blocks(i,1),Blocks(i,2)));
    Data{i}=getTaskParameters(myscreen,task);
    Contrast(i)=stimulus.contrasts;
    NumFixBreak(i)=sum(stimulus.FixationBreak);
    fixBreak{i}=stimulus.FixationBreak;
    end
end

%% 04 get subID
if Session == 1
    subjectID = strcat('s',SubNoString,Observer,'_S1');
    
elseif Session == 2
    subjectID = strcat('s',SubNoString,Observer,'_S2');
else
end

 %% 05 Get block info

BlockNums = Blocks(:,2)'; %get the block numbers
BlockTable = []; %preinit new table
for ii = 1:48 %number of trials in a block, this will have to be changed for other data
    BlockTable = [BlockTable;BlockNums];
end
%% 06 Pull out subject data
SubjectDataOriginal = [];
temp1=[];
for i=1:size(Blocks,1)
    fixIndex=logical(fixBreak{i}~=1);
    sum(fixIndex);
    temp1=[ temp1;Data{i}.randVars.targetOrientation(fixIndex)' Data{i}.response(fixIndex)' Data{i}.randVars.targetLocation(fixIndex)' Data{i}.randVars.ExoCueCondition(fixIndex)' Data{i}.reactionTime(fixIndex)'];

    SubjectDataOriginal=[ SubjectDataOriginal;Data{i}.randVars.targetOrientation(fixIndex)' Data{i}.response(fixIndex)' Data{i}.randVars.targetLocation(fixIndex)' Data{i}.randVars.contrast(fixIndex)' Data{i}.randVars.ExoCueCondition(fixIndex)' Data{i}.reactionTime(fixIndex)' Data{i}.blockNum(fixIndex)'];

end
temp1 = rmmissing(temp1,1);
SubjectDataOriginal = rmmissing(SubjectDataOriginal,1);

temp3=temp1(:,1)==temp1(:,2);

for ii = 1:height(temp1)
    if temp1(ii,5) < 0
        temp1(ii,5) = 0;
    else
    end
end
        

SubjectDataOriginal = horzcat(SubjectDataOriginal,temp3);



clearvars ii i 

%% 07 Get propcorrect and RT
% Location analysis: north=1, moves clockwise, west=4
for att=1:2 %1=valid, 2=neutral
for loc=1:4
    index=temp1(:,4)==att;
    TempData=temp1(index,:);
    index=TempData(:,3)==loc;
    TempData=TempData(index,:);

    TempResp=TempData(:,1)==TempData(:,2);
    TempRT=TempData(TempResp,5);            % only calculates RT for correct trials
    PropCorrect(att,loc,subj)=nanmean(TempResp);
   RT(att,loc,subj)=median(TempRT); % original - only since beginning of resp window
    if RT(att,loc,subj) <0 
        RT = 0 
    end
   %RT(att,loc,subj)=median(TempRT) + .820; % adds ms since stim onset
    clear index TempData TempResp;
end

end

%% 08 Define variables in table

SubjectDataOriginal = array2table(SubjectDataOriginal);
SubjectDataOriginal.Properties.VariableNames = {'targetOrientation', 'Response', 'TargetLoc', 'Contrast','ExoCueCondition', 'RT','Trial','rightwrong'};


%% IMPORTANT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_negO = find(SubjectDataOriginal.RT < 0); % creates an index of negative times
SubjectDataOriginal.RT(idx_negO) = 0; % transforms all negative RTs to zeros
SubjectDataOriginal = rmmissing(SubjectDataOriginal,'DataVariables','RT')

%% 09 add in intercardinal locations
%now add accuracy by loc
AccuracyValid = PropCorrect(1,:);
AccuracyNeutral = PropCorrect(2,:);

%add in Nans
AccuracyValid = [AccuracyValid(:,1),nan,AccuracyValid(:,2),nan,AccuracyValid(:,3),nan,AccuracyValid(:,4),nan];
AccuracyValid = array2table(AccuracyValid);
AccuracyValid.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

AccuracyNeutral = [AccuracyNeutral(:,1),nan,AccuracyNeutral(:,2),nan,AccuracyNeutral(:,3),nan,AccuracyNeutral(:,4),nan];
AccuracyNeutral = array2table(AccuracyNeutral);
AccuracyNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%% 10 Get MedianRTs

medianRTValid = RT(1,:);
medianRTNeutral = RT(2,:);

medianRTValid = [medianRTValid(:,1),nan,medianRTValid(:,2),nan,medianRTValid(:,3),nan,medianRTValid(:,4),nan];
medianRTNeutral = [medianRTNeutral(:,1),nan,medianRTNeutral(:,2),nan,medianRTNeutral(:,3),nan,medianRTNeutral(:,4),nan];

medianRTValid = array2table(medianRTValid);
medianRTValid.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

medianRTNeutral = array2table(medianRTNeutral);
medianRTNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

SubjectsStruct.(subjectID).medianRTValid = medianRTValid;
SubjectsStruct.(subjectID).medianRTNeutral = medianRTNeutral;
SubjectsStruct.(subjectID).Notes = Notes;

%% 11 now save outputs so far to struct 
SubjectsStruct.(subjectID).Source = Source;
SubjectsStruct.(subjectID).dataAvailability = dataAvailability;
SubjectsStruct.(subjectID).Session = Session;
SubjectsStruct.(subjectID).Included = Included;
SubjectsStruct.(subjectID).Gender = Gender;

%% 12 add observer initials
Observer = convertCharsToStrings(Observer);
ObserverTable = [];        
for ii = 1: height(SubjectDataOriginal) 
      ObserverTable = [ObserverTable;Observer];
end
ObserverTable = array2table(ObserverTable);
%SubjectDataOriginal.Observer = ObserverTable
SubjectDataOriginal = [SubjectDataOriginal  ObserverTable];
SubjectDataOriginal.Observer = SubjectDataOriginal.ObserverTable;

SubjectDataOriginal.ObserverTable=[];
%% 13 add height 
HeightinTable = [];        
for ii = 1: height(SubjectDataOriginal) 
      HeightinTable = [HeightinTable;Heightin];
end
HeightinTable = array2table(HeightinTable);
%SubjectDataOriginal.Heightin = HeightinTable
SubjectDataOriginal = [SubjectDataOriginal  HeightinTable];

SubjectDataOriginal.Heightin = SubjectDataOriginal.HeightinTable;
SubjectDataOriginal.HeightinTable=[];
%% 14 add SubNo
SubNoTable = [];        
for ii = 1: height(SubjectDataOriginal)
      SubNoTable = [SubNoTable;SubNo];
end
SubNoTable = array2table(SubNoTable);
SubjectDataOriginal = [SubjectDataOriginal  SubNoTable];

%SubjectDataOriginal.SubNo = SubNoTable
SubjectDataOriginal.SubNo = SubjectDataOriginal.SubNoTable;
SubjectDataOriginal.SubNoTable=[];
%% 15 Add age
AgeTable = [];        
for ii = 1: height(SubjectDataOriginal) 
      AgeTable = [AgeTable;Age];
end
AgeTable = array2table(AgeTable);
%SubjectDataOriginal.Age = AgeTable
SubjectDataOriginal = [SubjectDataOriginal  AgeTable];
SubjectDataOriginal.Age = SubjectDataOriginal.AgeTable;

SubjectDataOriginal.AgeTable=[];
%% 16 now make separate structs for valid and neutral

    indexValid =  SubjectDataOriginal.ExoCueCondition == 1;
    SubjectDataValid = SubjectDataOriginal(indexValid,:);
    clear indexValid

    indexNeutral =  SubjectDataOriginal.ExoCueCondition == 2;
    SubjectDataNeutral = SubjectDataOriginal(indexNeutral,:);
    clear indexNeutral
    
SubjectsStruct.(subjectID).SubjectDataValid = SubjectDataValid;    
SubjectsStruct.(subjectID).SubjectDataNeutral = SubjectDataNeutral;
SubjectsStruct.(subjectID).SubjectDataOriginal = SubjectDataOriginal;


%% 17 get weights

PropCorrectValid2 = [];
for ii = 1:numLocations
    index = SubjectDataValid.TargetLoc(:) == ii;
    PropCorrectValid2(ii)=sum(SubjectDataValid.rightwrong(index))/sum(index);
     NumTrialsPerLocValid(ii) = sum(index);
    WeightsValid(ii) = sum(index)/height(SubjectDataValid);
    %PropCorrect2.Location(ii)=sum(SubjectDataOriginal.TargetLoc(index))/sum(index);
    clear index
end

WeightsValid = [WeightsValid(1),nan,WeightsValid(2),nan,WeightsValid(3),nan,WeightsValid(4),nan];

PropCorrectNeutral2 = [];
for ii = 1:numLocations
    index = SubjectDataNeutral.TargetLoc(:) == ii;
    PropCorrectNeutral2(ii)=sum(SubjectDataNeutral.rightwrong(index))/sum(index);
     NumTrialsPerLocNeutral(ii) = sum(index);
    WeightsNeutral(ii) = sum(index)/height(SubjectDataNeutral);
    %PropCorrect2.Location(ii)=sum(SubjectDataOriginal.TargetLoc(index))/sum(index);
    clear index
end
WeightsNeutral = [WeightsNeutral(1),nan,WeightsNeutral(2),nan,WeightsNeutral(3),nan,WeightsNeutral(4),nan];


%% 18 Now correct accuracy - valid

%SubjectsStruct.(subjectID).PropCorrectNeutralCorrected = PropCorrectNeutral
%PropCorrectNeutralCorrected = []
clearvars ii 

AccuracyValid = table2array(AccuracyValid);
for ii = 1:width(AccuracyValid)
    if AccuracyValid(ii) <.50
        AccuracyValidCorrected(ii) = .50;
    else
        AccuracyValidCorrected(ii) = AccuracyValid(ii);
    end
    if ~isnan(AccuracyValidCorrected(ii))
     weightedValidCorrected(ii) = AccuracyValidCorrected(ii)*WeightsValid(ii);
    else
        weightedValidCorrected(ii) = nan;
    end
    end

AccuracyOverallValidWeighted = (nansum(weightedValidCorrected(:)))/(nansum(WeightsValid(:)));
%PropCorrectNeutralCorrected(ii)*
AccuracyOverallValidCorrectedWeighted = nanmean(AccuracyValidCorrected);


%% 19 Now correct accuracy - Neutral

%SubjectsStruct.(subjectID).PropCorrectNeutralCorrected = PropCorrectNeutral
%PropCorrectNeutralCorrected = []
clearvars ii 

AccuracyNeutral = table2array(AccuracyNeutral);
for ii = 1:width(AccuracyNeutral)
    if AccuracyNeutral(ii) <.50
        AccuracyNeutralCorrected(ii) = .50;
    else
        AccuracyNeutralCorrected(ii) = AccuracyNeutral(ii);
    end
    if ~isnan(AccuracyNeutralCorrected(ii))
     weightedNeutralCorrected(ii) = AccuracyNeutralCorrected(ii)*WeightsNeutral(ii);
    else
        weightedNeutralCorrected(ii) = nan;
    end
    end

AccuracyOverallNeutralWeighted = (nansum(weightedNeutralCorrected(:)))/(nansum(WeightsNeutral(:)));
%PropCorrectNeutralCorrected(ii)*
AccuracyOverallNeutralCorrectedWeighted = nanmean(AccuracyNeutralCorrected);



%% 20 Now define variables
AccuracyValid = array2table(AccuracyValid);
AccuracyValid.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

AccuracyValidCorrected = array2table(AccuracyValidCorrected);
AccuracyValidCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

SubjectsStruct.(subjectID).AccuracyOverallValidWeighted  = sum(SubjectDataValid.rightwrong)/length(SubjectDataValid.rightwrong);


%AccuracyNeutral = [AccuracyNeutral(:,1),nan,AccuracyNeutral(:,2),nan,AccuracyNeutral(:,3),nan,AccuracyNeutral(:,4),nan]
AccuracyNeutral = array2table(AccuracyNeutral);
AccuracyNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

AccuracyNeutralCorrected = array2table(AccuracyNeutralCorrected);
AccuracyNeutralCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


%define
SubjectsStruct.(subjectID).AccuracyValid = AccuracyValid;
SubjectsStruct.(subjectID).AccuracyValidCorrected = AccuracyValidCorrected;
SubjectsStruct.(subjectID).AccuracyOverallValidCorrectedWeighted = AccuracyOverallValidCorrectedWeighted;


SubjectsStruct.(subjectID).AccuracyNeutral = AccuracyNeutral;
SubjectsStruct.(subjectID).AccuracyNeutralCorrected = AccuracyNeutralCorrected;
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrectedWeighted = AccuracyOverallNeutralCorrectedWeighted;
SubjectsStruct.(subjectID).AccuracyOverallNeutralWeighted  = sum(SubjectDataNeutral.rightwrong)/length(SubjectDataNeutral.rightwrong);


%% 21 Now put everything else together
everythingElse = struct;

everythingElse.Age = Age;
everythingElse.BlockNums = BlockNums;
everythingElse.Blocks = Blocks;
everythingElse.Contrast = Contrast;


everythingElse.ContrastMean = mean(Contrast(:));
everythingElse.Data = Data;
everythingElse.Heightin = Heightin;
everythingElse.myscreen = myscreen;
everythingElse.Observer = Observer;
everythingElse.SubNo = SubNo;
everythingElse.task = task;
everythingElse.stimulus = stimulus;

SubjectsStruct.(subjectID).everythingElse = everythingElse;
SubjectsStruct.(subjectID).sf = stimulus.sf;
SubjectsStruct.(subjectID).ecc = stimulus.eccentricity;

SubjectsStruct.(subjectID).ChildAdultNumeric = ChildAdultNumeric;
SubjectsStruct.(subjectID).ChildAdultCategory = ChildAdultCategory;

%% 22 Get HVA VMA

HVANeutral(subj)=mean([PropCorrect(2,2,subj) PropCorrect(2,4,subj)])/mean([PropCorrect(2,1,subj) PropCorrect(2,3,subj)]);
HVAValid(subj)=mean([PropCorrect(1,2,subj) PropCorrect(1,4,subj)])/mean([PropCorrect(1,1,subj) PropCorrect(1,3,subj)]);
VMANeutral(subj)=PropCorrect(2,3,subj)/PropCorrect(2,1,subj);
VMAValid(subj)=PropCorrect(1,3,subj)/PropCorrect(1,1,subj);

SubjectsStruct.(subjectID).HVANeutral = HVANeutral;
SubjectsStruct.(subjectID).HVAValid = HVAValid;
SubjectsStruct.(subjectID).VMANeutral = VMANeutral;
SubjectsStruct.(subjectID).VMAValid = VMAValid;

% Contrasts(subj)=mean(Contrast);
% NumFixBreakTotal(subj) = sum(NumFixBreak);
% NumFixBreakMean(subj) = mean(NumFixBreak);

SubjectsStruct.(subjectID) = orderfields(SubjectsStruct.(subjectID));

SubjectsStruct.(subjectID) = orderfields(SubjectsStruct.(subjectID));
SubjectsStruct = orderfields(SubjectsStruct);
SubjectsStruct.(subjectID).NumTrialsLoc2 = table;
SubjectsStruct.(subjectID).NumTrialsLoc2.N = NumTrialsPerLocNeutral(1);
SubjectsStruct.(subjectID).NumTrialsLoc2.NE = NaN;
SubjectsStruct.(subjectID).NumTrialsLoc2.E = NumTrialsPerLocNeutral(2);
SubjectsStruct.(subjectID).NumTrialsLoc2.SE = NaN;
SubjectsStruct.(subjectID).NumTrialsLoc2.S = NumTrialsPerLocNeutral(3);
SubjectsStruct.(subjectID).NumTrialsLoc2.SW = NaN;
SubjectsStruct.(subjectID).NumTrialsLoc2.W = NumTrialsPerLocNeutral(4);
SubjectsStruct.(subjectID).NumTrialsLoc2.NW = NaN;


SubjectsStruct.(subjectID).everythingElse.NumFixBreak = NumFixBreak  
SubjectsStruct.AAALastEditDate = 'all corrections applied'
SubjectsStruct.(subjectID).attentionCond = attentionCond
%% 23 Save
cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/
save('SubjectsStruct.mat','SubjectsStruct');

cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/ExtractDataAdults/ExoAmb

disp(['Just finished with ',subjectID])