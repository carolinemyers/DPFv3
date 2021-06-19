%% 00 Hello

%Caroline Myers

%This script pulls out relevant trial info 
%for DPFv2 ADULT data collected by Katie and Rishi (KMAdults and RKAdults)
%No attention here folks- sorry to disappoint

%Assumes folder 'data' in same directory- contains data. 
%Saves output to structure- change path, probably. 

%last updated: 05.27.21
%% 01 Init

clearvars -except SubjectsStruct
%SubjectsStruct2 = SubjectsStruct
%SubjectsStruct = struct;
%% 02 Vars and constants
Observer= "dss2";

SubNoString = '357';

SubNo = str2num(SubNoString);
Session = 1;
Source = 'KMAdults';
dataAvailability = 'excel';
numLocations = 8;
Included = 1;
Notes = '';
ChildAdultCategory = 'adult';
ChildAdultNumeric = 3;

attentionCond = 'none'
sf = 4
ecc = 6
tilt = 30
%% other stuff
%import excel
%cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/DPFv2_ExtractData_Excel
%cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/Extract data/Adults/KMAdults
eval(sprintf('Folder=dir(''data/%s'');',Observer));
FileNum=length(Folder);


%for File=3:length(Folder)
for File = 10
if strcmp((Folder(File).name(end-4:end)),'.xlsx') == 1 || strcmp((Folder(File).name(end-3:end)),'.xls') == 1;
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
SubjectsStruct.(subjectID).SubjectDataNeutral = SubjDataOriginal;


%% 13 Location analysis
PropCorrect = [];
for ii = 1:numLocations
    index = SubjDataOriginal.Location(:) == ii;
    PropCorrect(ii)=sum(SubjDataOriginal.rightwrong(index))/sum(index);
    NumTrialsPerLoc(ii) = sum(index);
    WeightsNeutral(ii) = sum(index)/height(SubjDataOriginal);
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end

AccuracyNeutral = PropCorrect
clearvars ii
for ii = 1:width(WeightsNeutral)
weightedNeutralUncorrected(ii) = AccuracyNeutral(ii) * WeightsNeutral(ii)   
end

SubjectsStruct.(subjectID).AccuracyNeutralUncorrected = array2table(AccuracyNeutral);
SubjectsStruct.(subjectID).AccuracyNeutralUncorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

AccuracyOverallNeutralUncorrectedWeighted2 = nansum(weightedNeutralUncorrected)
SubjectsStruct.(subjectID).AccuracyOverallNeutralUncorrectedWeighted2 = AccuracyOverallNeutralUncorrectedWeighted2


totalNumTrials = sum(NumTrialsPerLoc(:));
NumTrialsNeutral = array2table(NumTrialsPerLoc)
NumTrialsNeutral.Properties.VariableNames = {'N','NE','E','SE','S','SW','W','NW'}
%Overall accuracy
PropCorrectOverall = sum(SubjDataOriginal.rightwrong)/length(SubjDataOriginal.rightwrong);

PropCorrect = array2table(PropCorrect); %make table

% Location analysis: north=1, moves clockwise, west=4
PropCorrect.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'}; %name variables

%%%%Now add in intercardinals, and set them to NaNs
%these are all nans because this data is only cardinals 

%make table of intercardinals and fill with nans
 
SubjectDataNeutral = SubjDataOriginal 

SubjectsStruct.(subjectID).SubjectDataNeutral = SubjectDataNeutral;
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
    weightedCorrected(ii) = PropCorrectNeutralCorrected(ii)*WeightsNeutral(ii);
%weightedCorrected(ii) = (PropCorrectNeutralCorrected(ii)*NumTrialsPerLoc(ii))/200
end
%% now weighted corrected

AccuracyNeutralCorrected = PropCorrectNeutralCorrected ;
SubjectsStruct.(subjectID).AccuracyNeutralCorrected = AccuracyNeutralCorrected;
clearvars ii
for ii = 1:width(WeightsNeutral)
weightedNeutralCorrected(ii) = AccuracyNeutralCorrected(ii) * WeightsNeutral(ii)   
end


AccuracyOverallNeutralCorrectedWeighted2 = nansum(weightedNeutralCorrected)
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrectedWeighted2 = AccuracyOverallNeutralCorrectedWeighted2

%%

AccuracyOverallNeutralWeighted = (sum(weightedCorrected(:)))/(sum(WeightsNeutral(:)));
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
SubjectsStruct.(subjectID).ChildAdultCategory = ChildAdultCategory;
SubjectsStruct.(subjectID).ChildAdultNumeric = ChildAdultNumeric
SubjectsStruct.(subjectID).Gender = SubjectDataNeutral.Gender{1}
%% 16 Now median RT 
clear index
clear ii 
%MedianRT = []
for ii = 1:numLocations
    index = SubjectDataNeutral.Location(:) == ii & SubjectDataNeutral.rightwrong(:) == 1;
    medianRTNeutral(ii) = median(SubjDataOriginal.RT(index));
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end

%%%%%%now add nans for intercardinal locations
%MedianRT = [MedianRT(:,1),nan,MedianRT(:,2), nan,MedianRT(:,3),nan,MedianRT(:,4),nan];


medianRTNeutral = array2table(medianRTNeutral);
medianRTNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


%%now add to sub structure
SubjectsStruct.(subjectID).medianRTNeutral = medianRTNeutral;
SubjectsStruct.(subjectID).SubjectDataNeutral = SubjectsStruct.(subjectID).SubjectDataNeutral

%%% add gender
Gender = SubjDataOriginal.Gender(1,1);


%%%%%add source
SubjectsStruct.(subjectID).Source = Source;
SubjectsStruct.(subjectID).dataAvailability = dataAvailability;
SubjectsStruct.(subjectID).Session = Session;
SubjectsStruct.(subjectID).Included = Included;
SubjectsStruct.(subjectID).Notes = Notes;
SubjectsStruct.(subjectID).NumTrialsPerLoc = NumTrialsPerLoc;
SubjectsStruct.(subjectID).NumTrialsNeutral = NumTrialsNeutral;
SubjectsStruct.(subjectID).Gender = string(Gender);

%% Get HVA and VMA

AccuracyNeutralCorrected = PropCorrectNeutralCorrected;
SubjectsStruct.(subjectID).AccuracyNeutralCorrected = AccuracyNeutralCorrected
%% HVA VMA uncorrected
HVANeutralUncorrected=mean([PropCorrectNeutral.East PropCorrectNeutral.West])/mean([PropCorrectNeutral.North PropCorrectNeutral.South]);
%HVAValid(subj)=mean([PropCorrectNeutral(1,2,subj) PropCorrectNeutral.South])/mean([PropCorrectNeutral(1,1,subj) PropCorrectNeutral(1,3,subj)]);
VMANeutralUncorrected=PropCorrectNeutral.South/PropCorrectNeutral.North;
%VMAValid(subj)=PropCorrectNeutral(1,3,subj)/PropCorrectNeutral(1,1,subj);

SubjectsStruct.(subjectID).HVANeutralUncorrected = HVANeutralUncorrected;
%SubjectsStruct.(subjectID).HVAValid = HVAValid;
SubjectsStruct.(subjectID).VMANeutralUncorrected = VMANeutralUncorrected;
%SubjectsStruct.(subjectID).VMAValid = VMAValid;

%% HVA VMA corrected
HVANeutralCorrected=mean([AccuracyNeutralCorrected.East AccuracyNeutralCorrected.West])/mean([AccuracyNeutralCorrected.North AccuracyNeutralCorrected.South]);
%HVAValid(subj)=mean([PropCorrectNeutral(1,2,subj) PropCorrectNeutral.South])/mean([PropCorrectNeutral(1,1,subj) PropCorrectNeutral(1,3,subj)]);
VMANeutralCorrected=AccuracyNeutralCorrected.South/AccuracyNeutralCorrected.North;
%VMAValid(subj)=PropCorrectNeutral(1,3,subj)/PropCorrectNeutral(1,1,subj);

SubjectsStruct.(subjectID).HVANeutralCorrected = HVANeutralCorrected;
%SubjectsStruct.(subjectID).HVAValid = HVAValid;
SubjectsStruct.(subjectID).VMANeutralCorrected = VMANeutralCorrected;
SubjectsStruct.(subjectID).sf = sf;
SubjectsStruct.(subjectID).ecc = ecc;
SubjectsStruct.(subjectID).tilt = tilt;

SubjectsStruct.(subjectID).HVANeutral = HVANeutralCorrected;
SubjectsStruct.(subjectID).VMANeutral = VMANeutralCorrected;

SubjectsStruct.(subjectID).numLocsTested = numLocations;
%% order fields
SubjectsStruct.(subjectID) = orderfields(SubjectsStruct.(subjectID));
SubjectsStruct = orderfields(SubjectsStruct);



%% 17 Now save
cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/
%SubjectsStruct2 = SubjectsStruct;
save('SubjectsStruct.mat','SubjectsStruct');

cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/ExtractData/Adults/KMAdults



