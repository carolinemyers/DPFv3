%% 00 Hello

%This script pulls out relevant trial info 
%for DPFv2 data collected by MG
%No attention here folks- sorry to disappoint

%% 01 Init
%SubjectsStruct2 = SubjectsStruct
clearvars -except SubjectsStruct
%% 02 Vars and constants

Observer= 'jjn';

SubNoString = '077'
Age = 5.54;
Heightin = 47;
Gender = 'm';
Session = 1;

Source = 'MGData';
dataAvailability = 'eyebeh';
numLocations = 4;
Included = 1;
Notes = '';

SubNo = str2num(SubNoString);
%% 03 Parse files
%Observer= {'bio'};
[Blocks,TestFiles] =DPFParseFilesCM(Observer);

for i=1:size(Blocks,1)
    if Blocks(i,2)<10
    eval(sprintf('load(''data/%s/%g_stim0%g.mat'')',Observer,Blocks(i,1),Blocks(i,2)));
    Data{i}=getTaskParameters(myscreen,task);
    Contrast(i)=stimulus.contrasts;
    else
    eval(sprintf('load(''data/%s/%g_stim%g.mat'')',Observer,Blocks(i,1),Blocks(i,2)));
    Data{i}=getTaskParameters(myscreen,task); %Data shoould be 1 x n cell where n = number of blocks
    Contrast(i)=stimulus.contrasts;
    end
end

%% 04 Pull out data from parsed files

% make loop so that block numbers can be concatenated 
BlockNums = Blocks(:,2)'; %get the block numbers
BlockTable = []; %preinit new table
for ii = 1:48 %number of trials in a block, this will have to be changed for other data
    BlockTable = [BlockTable;BlockNums];
end

%pull out relevant info (contrast, target ori, whatever you feel like)
SubjDataOriginal = [];
for i=1:size(Blocks,1)
    SubjDataOriginal=[ SubjDataOriginal;...
    Data{i}.randVars.targetOrientation'...
    Data{i}.response'...
    Data{i}.randVars.targetLocation'...
    Data{i}.randVars.contrast' ...
    Data{i}.reactionTime'...
    Data{i}.blockNum'...  #trial num
    BlockTable(:,i)];
end

SubjDataOriginal = array2table(SubjDataOriginal);
SubjDataOriginal.Properties.VariableNames = {'targetOrientation', 'Response', 'TargetLoc', 'Contrast', 'RT','Trial','Block'};

%% 05 add in whether they were correct
for ii = 1:height(SubjDataOriginal)
    if SubjDataOriginal.Response(ii) == SubjDataOriginal.targetOrientation(ii)
        SubjDataOriginal.rightwrong(ii) = 1;
    elseif SubjDataOriginal.Response ~= SubjDataOriginal.targetOrientation(ii)
        SubjDataOriginal.rightwrong(ii) = 0;
    end
end

%% 06 add observer initials
Observer = convertCharsToStrings(Observer);
ObserverTable = [];        
for ii = 1: height(SubjDataOriginal) 
      ObserverTable = [ObserverTable;Observer];
end
ObserverTable = array2table(ObserverTable);
SubjDataOriginal = [SubjDataOriginal  ObserverTable];

%% 07 add height 
HeightinTable = [];        
for ii = 1: height(SubjDataOriginal) 
      HeightinTable = [HeightinTable;Heightin];
end
HeightinTable = array2table(HeightinTable);
SubjDataOriginal = [SubjDataOriginal  HeightinTable];

%% 08 add SubNo
SubNoTable = [];        
for ii = 1: height(SubjDataOriginal)
      SubNoTable = [SubNoTable;SubNo];
end
SubNoTable = array2table(SubNoTable);
SubjDataOriginal = [SubjDataOriginal  SubNoTable];

%% 09 Add age
AgeTable = [];        
for ii = 1: height(SubjDataOriginal) 
      AgeTable = [AgeTable;Age];
end
AgeTable = array2table(AgeTable);
SubjDataOriginal = [SubjDataOriginal  AgeTable];

%% 10 Add variable names, Clear unneccessary variables and save variables to structure 

SubjDataOriginal.Properties.VariableNames = {'targetOrientation', 'Response', 'TargetLoc', 'Contrast', 'RT','Trial','Block','rightwrong', 'Observer','Heightin','SubNo','Age'};

clearvars AgeTable BlockTable HeightinTable i ii SubNoTable ObserverTable

%% 11 Now put everything else together
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
everythingElse.TestFiles = TestFiles;
%% 12 Now save data
%SubNoStr = num2str(SubNo);
if Session == 1
    subjectID = strcat('s',SubNoString,Observer,'_S1');
    
elseif Session == 2
    subjectID = strcat('s',SubNoString,Observer,'_S2');
else
end

%SubjectsStruct = struct;


SubjectsStruct.(subjectID).everythingElse = everythingElse;
SubjectsStruct.(subjectID).SubjectData = SubjDataOriginal;


%% 13 Location analysis
PropCorrect = [];
for ii = 1:numLocations
    index = SubjDataOriginal.TargetLoc(:) == ii;
    PropCorrect(ii)=sum(SubjDataOriginal.rightwrong(index))/sum(index);
     NumTrialsPerLoc(ii) = sum(index);
    Weights(ii) = sum(index)/height(SubjDataOriginal);
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end

totalNumTrials = sum(NumTrialsPerLoc(:));
SubjectsStruct.(subjectID).totalNumTrials = totalNumTrials;

Weights = [Weights(1),nan,Weights(2),nan,Weights(3),nan,Weights(4),nan];
NumTrialsPerLoc = [NumTrialsPerLoc(1),nan,NumTrialsPerLoc(2),nan,NumTrialsPerLoc(3),nan,NumTrialsPerLoc(4),nan];
%Overall accuracy
PropCorrectOverallWeighted = sum(SubjDataOriginal.rightwrong)/length(SubjDataOriginal.rightwrong);

PropCorrect = array2table(PropCorrect); %make table

% Location analysis: north=1, moves clockwise, west=4
PropCorrect.Properties.VariableNames = {'North','East','South','West'}; %name variables


%%%%Now add in intercardinals, and set them to NaNs
%these are all nans because this data is only cardinals 

%make table of intercardinals and fill with nans
intercardinalPropCorrect = table(nan,nan,nan,nan); 
intercardinalPropCorrect.Properties.VariableNames = {'Northeast','Southeast','Southwest','Northwest'};

PropCorrectNeutral = table;
PropCorrectNeutral = [PropCorrect.North ...
    intercardinalPropCorrect.Northeast...
    PropCorrect.East...
    intercardinalPropCorrect.Southeast...
    PropCorrect.South...
    intercardinalPropCorrect.Southwest...
    PropCorrect.West...
    intercardinalPropCorrect.Northwest];

%% 14 Correct to 50% 

%SubjectsStruct.(subjectID).PropCorrectNeutralCorrected = PropCorrectNeutral
%PropCorrectNeutralCorrected = []
for ii = 1:width(PropCorrectNeutral)
    if PropCorrectNeutral(ii) <.50
        PropCorrectNeutralCorrected(ii) = .50;
    else
        PropCorrectNeutralCorrected(ii) = PropCorrectNeutral(ii);
    end
    if ~isnan(PropCorrectNeutralCorrected(ii))
     weightedCorrected(ii) = PropCorrectNeutralCorrected(ii)*Weights(ii);
    else
        weightedCorrected(ii) = nan;
    end
    end

AccuracyOverallNeutralWeighted = (nansum(weightedCorrected(:)))/(nansum(Weights(:)));
%PropCorrectNeutralCorrected(ii)*
AccuracyOverallNeutralCorrected = nanmean(PropCorrectNeutralCorrected);

%AccuracyOverallNeutralCorrected = nanmean(PropCorrectNeutralCorrected);
%% 15 Label variables 

PropCorrectNeutral = array2table(PropCorrectNeutral);
PropCorrectNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


PropCorrectNeutralCorrected = array2table(PropCorrectNeutralCorrected);
PropCorrectNeutralCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


%%%%LABEL VARIABLES HERE%%%% 
SubjectsStruct.(subjectID).AccuracyNeutral = PropCorrectNeutral;
SubjectsStruct.(subjectID).AccuracyNeutralCorrected = PropCorrectNeutralCorrected;
SubjectsStruct.(subjectID).AccuracyOverallNeutral = PropCorrectOverallWeighted;
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrected = AccuracyOverallNeutralCorrected;
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrectedWeighted = AccuracyOverallNeutralWeighted;
SubjectsStruct.(subjectID).Notes = Notes;
SubjectsStruct.(subjectID).NumTrialsPerLoc = NumTrialsPerLoc;

%% 16 Now median RT 
clear index
clear ii 
%MedianRT = []
for ii = 1:4
    index = SubjDataOriginal.TargetLoc(:) == ii & SubjDataOriginal.rightwrong(:) == 1;
    MedianRT(ii) = median(SubjDataOriginal.RT(index));
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end

%%%%%%now add nans for intercardinal locations
MedianRT = [MedianRT(:,1),nan,MedianRT(:,2), nan,MedianRT(:,3),nan,MedianRT(:,4),nan];


MedianRT = array2table(MedianRT);
MedianRT.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


%%now add to sub structure
SubjectsStruct.(subjectID).medianRT = MedianRT;

%%%%%add source
SubjectsStruct.(subjectID).Source = Source;
SubjectsStruct.(subjectID).dataAvailability = dataAvailability;
SubjectsStruct.(subjectID).Session = Session;
SubjectsStruct.(subjectID).Included = Included;
SubjectsStruct.(subjectID).Gender = Gender;

everythingElse.stimulus = stimulus;

SubjectsStruct.(subjectID).everythingElse = everythingElse
SubjectsStruct.(subjectID).sf = stimulus.sf;
SubjectsStruct.(subjectID).ecc = stimulus.eccentricity;

%% 17 Now save
SubjectsStruct.(subjectID) = orderfields(SubjectsStruct.(subjectID))

SubjectsStruct = orderfields(SubjectsStruct)

cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts
save('SubjectsStruct.mat','SubjectsStruct');

%cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/DPFv2_Extract_MG_Data
%save('SubjectsStruct.mat','SubjectsStruct');