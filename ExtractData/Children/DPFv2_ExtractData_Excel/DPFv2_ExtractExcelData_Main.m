%% 00 Hello

%Caroline Myers

%This script pulls out relevant trial info 
%for DPFv2 data collected by Kate and Rishi
%No attention here folks- sorry to disappoint

%Assumes folder 'data' in same directory- contains data. 
%Saves output to structure- change path, probably. 
%% 01 Init

clearvars -except SubjectsStruct
%SubjectsStruct2 = SubjectsStruct
%SubjectsStruct = struct;
%% 02 Vars and constants
Observer= "BR";

SubNoString = '075';

SubNo = str2num(SubNoString);
Session = 1;
Source = '5811Grade';
dataAvailability = 'excel';
numLocations = 8;
Included = 1;
Notes = '';

%% other stuff
%import excel
eval(sprintf('Folder=dir(''data/%s'');',Observer));
FileNum=length(Folder);


for File=3:length(Folder)

if strcmp((Folder(File).name(end-4:end)),'.xlsx') == 1
    SubjDataOriginal = readtable(Folder(File).name);
else
end
end

SubjDataOriginal.Heightin = SubjDataOriginal.Height_in_;
SubjDataOriginal.Block = SubjDataOriginal.Run;


%% 08 add SubNo
SubNoTable = [];        
for ii = 1: height(SubjDataOriginal)
      SubNoTable = [SubNoTable;SubNo];
end
SubNoTable = array2table(SubNoTable);
SubjDataOriginal = [SubjDataOriginal  SubNoTable];


%% 10 Add variable names, Clear unneccessary variables and save variables to structure 

%SubjDataOriginal.Properties.VariableNames = {'targetOrientation', 'Response', 'TargetLoc', 'Contrast', 'RT','Trial','Block','rightwrong', 'Observer','Heightin','SubNo','Age'};

%clearvars AgeTable BlockTable HeightinTable i ii SubNoTable ObserverTable

%% 11 Now put everything else together
everythingElse = struct;

Age = SubjDataOriginal.Age(1,1);
BlockNums = unique(SubjDataOriginal.Run)';
Contrast = unique(SubjDataOriginal.Contrast)';

ContrastMean = mean(SubjDataOriginal.Contrast);
Heightin = SubjDataOriginal.Height_in_(1,1);


everythingElse.ContrastMean = ContrastMean;
everythingElse.Age = Age;
everythingElse.BlockNums = BlockNums;


everythingElse.Contrast = Contrast;
%everythingElse.Data = Data;
everythingElse.Heightin = Heightin;
%everythingElse.myscreen = myscreen;
everythingElse.Observer = Observer;
everythingElse.SubNo = SubNo;
%everythingElse.task = task;
%everythingElse.TestFiles = TestFiles;
%% 12 Now save data
%SubNoStr = num2str(SubNo);
if Session == 1
    subjectID = strcat('s',SubNoString,Observer,'_S1');
    
elseif Session == 2
    subjectID = strcat('s',SubNoString,Observer,'_S2');
else
end

%subjectID = strcat('s',SubNoStr,Observer);
%SubjectsStruct = struct;


SubjectsStruct.(subjectID).everythingElse = everythingElse;
SubjectsStruct.(subjectID).SubjectData = SubjDataOriginal;


%% 13 Location analysis
PropCorrect = [];
for ii = 1:numLocations
    index = SubjDataOriginal.Location(:) == ii;
    PropCorrect(ii)=sum(SubjDataOriginal.rightwrong(index))/sum(index);
    NumTrialsPerLoc(ii) = sum(index);
    Weights(ii) = sum(index)/height(SubjDataOriginal);
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end
totalNumTrials = sum(NumTrialsPerLoc(:));

%Overall accuracy
PropCorrectOverall = sum(SubjDataOriginal.rightwrong)/length(SubjDataOriginal.rightwrong);

PropCorrect = array2table(PropCorrect); %make table

% Location analysis: north=1, moves clockwise, west=4
PropCorrect.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'}; %name variables


%%%%Now add in intercardinals, and set them to NaNs
%these are all nans because this data is only cardinals 

%make table of intercardinals and fill with nans
 

PropCorrectNeutral = table2array(PropCorrect);


%% 14 Correct to 50% 

%SubjectsStruct.(subjectID).PropCorrectNeutralCorrected = PropCorrectNeutral
%PropCorrectNeutralCorrected = []
for ii = 1:width(PropCorrectNeutral)
    if PropCorrectNeutral(ii) <.50
        PropCorrectNeutralCorrected(ii) = .50;
    else
        PropCorrectNeutralCorrected(ii) = PropCorrectNeutral(ii);
        %Weights(ii) = PropCorrectNeutralCorrected(ii)*NumTrialsPerLoc(ii)
    end
    weightedCorrected(ii) = PropCorrectNeutralCorrected(ii)*Weights(ii);
%weightedCorrected(ii) = (PropCorrectNeutralCorrected(ii)*NumTrialsPerLoc(ii))/200
end
AccuracyOverallNeutralWeighted = (sum(weightedCorrected(:)))/(sum(Weights(:)));
%PropCorrectNeutralCorrected(ii)*
AccuracyOverallNeutralCorrected = nanmean(PropCorrectNeutralCorrected);
%% 15 Label variables 

PropCorrectNeutral = array2table(PropCorrectNeutral);
PropCorrectNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


PropCorrectNeutralCorrected = array2table(PropCorrectNeutralCorrected);
PropCorrectNeutralCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


%%%%LABEL VARIABLES HERE%%%% 
SubjectsStruct.(subjectID).AccuracyNeutral = PropCorrectNeutral;
SubjectsStruct.(subjectID).AccuracyNeutralCorrected = PropCorrectNeutralCorrected;
SubjectsStruct.(subjectID).AccuracyOverallNeutralWeighted = PropCorrectOverall;
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrected = AccuracyOverallNeutralCorrected;
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrectedWeighted = AccuracyOverallNeutralWeighted;
SubjectsStruct.(subjectID).totalNumTrials = totalNumTrials;
%% 16 Now median RT 
clear index
clear ii 
%MedianRT = []
for ii = 1:numLocations
    index = SubjDataOriginal.Location(:) == ii & SubjDataOriginal.rightwrong(:) == 1;
    MedianRT(ii) = median(SubjDataOriginal.RT(index));
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end

%%%%%%now add nans for intercardinal locations
%MedianRT = [MedianRT(:,1),nan,MedianRT(:,2), nan,MedianRT(:,3),nan,MedianRT(:,4),nan];


MedianRT = array2table(MedianRT);
MedianRT.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


%%now add to sub structure
SubjectsStruct.(subjectID).medianRT = MedianRT;
SubjectsStruct.(subjectID).SubjectDataNeutral = SubjectsStruct.(subjectID).SubjectData

%%% add gender
Gender = SubjDataOriginal.Gender(1,1);
% if isnan(Gender)
%     Gender = 'u';
% else
%     Gender = Gender;
% end


%%%%%add source
SubjectsStruct.(subjectID).Source = Source;
SubjectsStruct.(subjectID).dataAvailability = dataAvailability;
SubjectsStruct.(subjectID).Session = Session;
SubjectsStruct.(subjectID).Included = Included;
SubjectsStruct.(subjectID).Notes = Notes;
SubjectsStruct.(subjectID).NumTrialsPerLoc = NumTrialsPerLoc;
SubjectsStruct.(subjectID).Gender = Gender;

%% 17 Now save
cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/DPFv2_ExtractData_Excel
%SubjectsStruct2 = SubjectsStruct;
save('SubjectsStruct.mat','SubjectsStruct');