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
Observer= "HF";

SubNoString = '328';

SubNo = str2num(SubNoString);
Session = 2;
Source = 'CMv3';
dataAvailability = 'eyebeh';
numLocations = 4;
Included = 1;
Notes = '';
ChildAdultCategory = 'child';
ChildAdultNumeric = 1;
Gender = 'F';
Age = 9.443836;
sf = 4;
ecc = 6.4;
tilt = 20;
Heightin = 54;
attentionCond = 'exo';
%add attn cond
%% other stuff
%import excel
eval(sprintf('Folder=dir(''data/%s'');',Observer));
FileNum=length(Folder);


for File= 4   %3:length(Folder)

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

%% remove nans


SubjectDataOriginal = rmmissing(SubjDataOriginal,'DataVariables',{'rightwrong','RT'})

idx_negO = find(SubjectDataOriginal.RT < 0); % creates an index of negative times
SubjectDataOriginal.RT(idx_negO) = 0; % transforms all negative RTs to zeros
SubjectDataOriginal = rmmissing(SubjectDataOriginal,'DataVariables','RT')

clearvars idx_negO
SubjectsStruct.(subjectID).SubjectDataOriginal = SubjectDataOriginal
SubjectsStruct.(subjectID).everythingElse = everythingElse;
%SubjectsStruct.(subjectID).SubjectData = SubjDataOriginal;

%% Then, make 2 tables, one for valid, one for neutral
%%neutral
indexNeutral = SubjectsStruct.(subjectID).SubjectDataOriginal.CueCond == 2; %neutral
SubjectDataNeutral = SubjectsStruct.(subjectID).SubjectDataOriginal(indexNeutral,:);

%Valid
indexValid = SubjectsStruct.(subjectID).SubjectDataOriginal.CueCond == 1; %valid
SubjectDataValid = SubjectsStruct.(subjectID).SubjectDataOriginal(indexValid,:);

% indexB = strcmp(IncludedDataTableCleaned.Source,'MGData') == 1 & IncludedDataTableCleaned.sf == 6 ...
%     | strcmp(IncludedDataTableCleaned.Source,'MGDataAdults') == 1 & IncludedDataTableCleaned.sf == 6;
% B_IncludedDataTableCleaned6cpdNoDistractors = IncludedDataTableCleaned(indexB,:);
%% now do this for each 
%%note to CM- now you need to document each variable name you care about so
%%you can easily read them in your analysis script. Go to that and copy the
%%names tomorrow. 

%% Remember that the locations here are different. 1 = E, 2 = N, 3 = W, 4 = S 
PropCorrectNeutral = [];
for ii = 1:numLocations
    index = SubjectDataNeutral.Location(:) == ii;
    PropCorrectNeutral(ii)=sum(SubjectDataNeutral.rightwrong(index))/sum(index);
    NumTrialsPerLocNeutral(ii) = sum(index);
    WeightsNeutral(ii) = sum(index)/height(SubjectDataNeutral);
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end
WeightsNeutral = horzcat(WeightsNeutral(2),NaN, WeightsNeutral(1), NaN, WeightsNeutral(4), NaN,WeightsNeutral(3),NaN);
NumTrialsPerLocNeutral = horzcat(NumTrialsPerLocNeutral(2),NaN, NumTrialsPerLocNeutral(1), NaN, NumTrialsPerLocNeutral(4), NaN,NumTrialsPerLocNeutral(3),NaN);

clear index
totalNumTrialsNeutral = nansum(NumTrialsPerLocNeutral(:));

%Overall accuracy
PropCorrectOverallNeutral = sum(SubjectDataNeutral.rightwrong)/length(SubjectDataNeutral.rightwrong);

PropCorrectNeutral = horzcat(PropCorrectNeutral(2),NaN, PropCorrectNeutral(1), NaN, PropCorrectNeutral(4), NaN,PropCorrectNeutral(3),NaN);

PropCorrectNeutral = array2table(PropCorrectNeutral); %make table
PropCorrectNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%% Now do this for valid 
PropCorrectValid = [];
for ii = 1:numLocations
    index = SubjectDataValid.Location(:) == ii;
    PropCorrectValid(ii)=sum(SubjectDataValid.rightwrong(index))/sum(index);
    NumTrialsPerLocValid(ii) = sum(index);
    WeightsValid(ii) = sum(index)/height(SubjectDataValid);
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end
NumTrialsPerLocValid = horzcat(NumTrialsPerLocValid(2),NaN, NumTrialsPerLocValid(1), NaN, NumTrialsPerLocValid(4), NaN,NumTrialsPerLocValid(3),NaN);

WeightsValid = horzcat(WeightsValid(2),NaN, WeightsValid(1), NaN, WeightsValid(4), NaN,WeightsValid(3),NaN);

totalNumTrialsValid = nansum(NumTrialsPerLocValid(:));

%Overall accuracy
PropCorrectOverallValid = sum(SubjectDataValid.rightwrong)/length(SubjectDataValid.rightwrong);

PropCorrectValid = horzcat(PropCorrectValid(2),NaN, PropCorrectValid(1), NaN, PropCorrectValid(4), NaN,PropCorrectValid(3),NaN)

PropCorrectValid = array2table(PropCorrectValid); %make table
PropCorrectValid.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%% Now first get uncorrected weights- valid
AccuracyValid = table2array(PropCorrectValid);
for ii = 1:width(AccuracyValid)
    if ~isnan(AccuracyValid(ii))
     weightedValidUncorrected(ii) = AccuracyValid(ii)*WeightsValid(ii);
    else
        weightedValidUncorrected(ii) = nan;
    end
end
    
AccuracyOverallValidWeightedUncorrected = (nansum(weightedValidUncorrected(:)))/(nansum(WeightsValid(:)));

% and again for neutral
AccuracyNeutral = table2array(PropCorrectNeutral);
for ii = 1:width(AccuracyNeutral)
    if ~isnan(AccuracyNeutral(ii))
     weightedNeutralUncorrected(ii) = AccuracyNeutral(ii)*WeightsNeutral(ii);
    else
        weightedNeutralUncorrected(ii) = nan;
    end
end
    
AccuracyOverallNeutralWeightedUncorrected = (nansum(weightedNeutralUncorrected(:)))/(nansum(WeightsNeutral(:)));

%% Now get corrected weights- valid

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

AccuracyOverallValidCorrectedWeighted = (nansum(weightedValidCorrected(:)))/(nansum(WeightsValid(:)));
%PropCorrectNeutralCorrected(ii)*
%AccuracyOverallValidCorrectedWeighted = nanmean(AccuracyValidCorrected);


%% Now correct and weight neutral

clearvars ii 

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

AccuracyOverallNeutralCorrectedWeighted = (nansum(weightedNeutralCorrected(:)))/(nansum(WeightsNeutral(:)));

%% Now RT 

for ii = 1:numLocations
    index = SubjectDataNeutral.Location(:) == ii & SubjectDataNeutral.rightwrong(:) == 1;
    medianRTNeutral(ii) = median(SubjectDataNeutral.RT(index));
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end
medianRTNeutral = horzcat(medianRTNeutral(2),NaN, medianRTNeutral(1), NaN, medianRTNeutral(4), NaN,medianRTNeutral(3),NaN);


for ii = 1:numLocations
    index = SubjectDataValid.Location(:) == ii & SubjectDataValid.rightwrong(:) == 1;
    medianRTValid(ii) = median(SubjectDataValid.RT(index));
    %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
    clear index
end
medianRTValid = horzcat(medianRTValid(2),NaN, medianRTValid(1), NaN, medianRTValid(4), NaN,medianRTValid(3),NaN);



%% Now that we have everything, let's get it saved how we want for the struct. 
%accuracy neutral
SubjectsStruct.(subjectID).AccuracyNeutral = array2table(AccuracyNeutral);
SubjectsStruct.(subjectID).AccuracyNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%accuracy valid
SubjectsStruct.(subjectID).AccuracyValid = array2table(AccuracyValid);
SubjectsStruct.(subjectID).AccuracyValid.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%accuracy neutral corrected
SubjectsStruct.(subjectID).AccuracyNeutralCorrected = array2table(AccuracyNeutralCorrected);
SubjectsStruct.(subjectID).AccuracyNeutralCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%accuracy valid corrected
SubjectsStruct.(subjectID).AccuracyValidCorrected = array2table(AccuracyValidCorrected);
SubjectsStruct.(subjectID).AccuracyValidCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%Neutral Uncorrected overall
SubjectsStruct.(subjectID).AccuracyOverallNeutralWeightedUncorrected = AccuracyOverallNeutralWeightedUncorrected;
SubjectsStruct.(subjectID).AccuracyOverallNeutralWeighted = AccuracyOverallNeutralWeightedUncorrected;

%Valid Uncorrected overall
SubjectsStruct.(subjectID).AccuracyOverallValidWeightedUncorrected = AccuracyOverallValidWeightedUncorrected;
SubjectsStruct.(subjectID).AccuracyOverallValidWeighted = AccuracyOverallValidWeightedUncorrected;

%Neutral corrected overall
SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrectedWeighted = AccuracyOverallNeutralCorrectedWeighted;
%Valid corrected overall
SubjectsStruct.(subjectID).AccuracyOverallValidCorrectedWeighted = AccuracyOverallValidCorrectedWeighted;

%Num trials per loc neutral
SubjectsStruct.(subjectID).NumTrialsNeutral = array2table(NumTrialsPerLocNeutral);
SubjectsStruct.(subjectID).NumTrialsNeutral.Properties.VariableNames = {'N','NE','E','SE','S','SW','W','NW'};

%Num trials per loc valid
SubjectsStruct.(subjectID).NumTrialsValid = array2table(NumTrialsPerLocValid);
SubjectsStruct.(subjectID).NumTrialsValid.Properties.VariableNames = {'N','NE','E','SE','S','SW','W','NW'};


%MedianRTNeutral
SubjectsStruct.(subjectID).medianRTNeutral = array2table(medianRTNeutral);
SubjectsStruct.(subjectID).medianRTNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};


SubjectsStruct.(subjectID).medianRTValid = array2table(medianRTValid);
SubjectsStruct.(subjectID).medianRTValid.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};

%stim parameters

SubjectsStruct.(subjectID).tilt = tilt;
SubjectsStruct.(subjectID).sf = sf;
SubjectsStruct.(subjectID).ecc = ecc;
SubjectsStruct.(subjectID).ChildAdultCategory = ChildAdultCategory;
SubjectsStruct.(subjectID).ChildAdultNumeric = ChildAdultNumeric;
SubjectsStruct.(subjectID).Notes = Notes;
SubjectsStruct.(subjectID).Source = Source;
SubjectsStruct.(subjectID).Session = Session;
SubjectsStruct.(subjectID).everythingElse.BlockNums = BlockNums;
SubjectsStruct.(subjectID).everythingElse.Contrast = Contrast;
SubjectsStruct.(subjectID).everythingElse.ContrastMean = mean(Contrast);
SubjectsStruct.(subjectID).everythingElse.SubNo = SubNo;
SubjectsStruct.(subjectID).everythingElse.Heightin = Heightin;
SubjectsStruct.(subjectID).everythingElse.Observer = subjectID;
SubjectsStruct.(subjectID).dataAvailability = dataAvailability;
SubjectsStruct.(subjectID).attentionCond = attentionCond;
SubjectsStruct.(subjectID).Included = Included;
SubjectsStruct.(subjectID).everythingElse.Age = Age;
SubjectsStruct.(subjectID).numLocsTested = numLocations;
%% Last, add HVA and VMA

% HVANeutral(subj)=mean([PropCorrect(2,2,subj) PropCorrect(2,4,subj)])/mean([PropCorrect(2,1,subj) PropCorrect(2,3,subj)]);
% HVAValid(subj)=mean([PropCorrect(1,2,subj) PropCorrect(1,4,subj)])/mean([PropCorrect(1,1,subj) PropCorrect(1,3,subj)]);
% VMANeutral(subj)=PropCorrect(2,3,subj)/PropCorrect(2,1,subj);
% VMAValid(subj)=PropCorrect(1,3,subj)/PropCorrect(1,1,subj);
% 
% SubjectsStruct.(subjectID).HVANeutral = HVANeutral;
% SubjectsStruct.(subjectID).HVAValid = HVAValid;
% SubjectsStruct.(subjectID).VMANeutral = VMANeutral;
% SubjectsStruct.(subjectID).VMAValid = VMAValid;
%% Now save
SubjectsStruct.(subjectID) = orderfields(SubjectsStruct.(subjectID));
SubjectsStruct = orderfields(SubjectsStruct);

cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/
%SubjectsStruct2 = SubjectsStruct;
save('SubjectsStruct.mat','SubjectsStruct');

cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/ExtractData/Adults/CMv3/data
%% 14 Correct to 50% 
% 
% % %SubjectsStruct.(subjectID).PropCorrectNeutralCorrected = PropCorrectNeutral
% % %PropCorrectNeutralCorrected = []
% % for ii = 1:width(PropCorrectNeutral)
% %     if PropCorrectNeutral(ii) <.50
% %         PropCorrectNeutralCorrected(ii) = .50;
% %     else
% %         PropCorrectNeutralCorrected(ii) = PropCorrectNeutral(ii);
% %         %Weights(ii) = PropCorrectNeutralCorrected(ii)*NumTrialsPerLoc(ii)
% %     end
% %     weightedCorrected(ii) = PropCorrectNeutralCorrected(ii)*Weights(ii);
% % %weightedCorrected(ii) = (PropCorrectNeutralCorrected(ii)*NumTrialsPerLoc(ii))/200
% % end
% % AccuracyOverallNeutralWeighted = (sum(weightedCorrected(:)))/(sum(Weights(:)));
% % %PropCorrectNeutralCorrected(ii)*
% % AccuracyOverallNeutralCorrected = nanmean(PropCorrectNeutralCorrected);
% %% 15 Label variables 
% 
% PropCorrectNeutral = array2table(PropCorrectNeutral);
% PropCorrectNeutral.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};
% 
% 
% PropCorrectNeutralCorrected = array2table(PropCorrectNeutralCorrected);
% PropCorrectNeutralCorrected.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};
% 
% 
% %%%%LABEL VARIABLES HERE%%%% 
% SubjectsStruct.(subjectID).AccuracyNeutral = PropCorrectNeutral;
% SubjectsStruct.(subjectID).AccuracyNeutralCorrected = PropCorrectNeutralCorrected;
% SubjectsStruct.(subjectID).AccuracyOverallNeutralWeighted = PropCorrectOverallNeutral;
% SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrected = AccuracyOverallNeutralCorrected;
% SubjectsStruct.(subjectID).AccuracyOverallNeutralCorrectedWeighted = AccuracyOverallNeutralWeighted;
% SubjectsStruct.(subjectID).totalNumTrials = totalNumTrials;
% %% 16 Now median RT 
% clear index
% clear ii 
% %MedianRT = []
% for ii = 1:numLocations
%     index = SubjDataOriginal.Location(:) == ii & SubjDataOriginal.rightwrong(:) == 1;
%     MedianRT(ii) = median(SubjDataOriginal.RT(index));
%     %PropCorrect.Location(ii)=sum(SubjDataOriginal.TargetLoc(index))/sum(index);
%     clear index
% end
% 
% %%%%%%now add nans for intercardinal locations
% %MedianRT = [MedianRT(:,1),nan,MedianRT(:,2), nan,MedianRT(:,3),nan,MedianRT(:,4),nan];
% 
% 
% MedianRT = array2table(MedianRT);
% MedianRT.Properties.VariableNames = {'North','Northeast','East','Southeast','South','Southwest','West','Northwest'};
% 
% 
% %%now add to sub structure
% SubjectsStruct.(subjectID).medianRT = MedianRT;
% SubjectsStruct.(subjectID).SubjectDataNeutral = SubjectsStruct.(subjectID).SubjectData
% 
% %%% add gender
% Gender = SubjDataOriginal.Gender(1,1);
% % if isnan(Gender)
% %     Gender = 'u';
% % else
% %     Gender = Gender;
% % end
% 
% 
% %%%%%add source
% SubjectsStruct.(subjectID).Source = Source;
% SubjectsStruct.(subjectID).dataAvailability = dataAvailability;
% SubjectsStruct.(subjectID).Session = Session;
% SubjectsStruct.(subjectID).Included = Included;
% SubjectsStruct.(subjectID).Notes = Notes;
% SubjectsStruct.(subjectID).NumTrialsPerLoc = NumTrialsPerLocNeutral;
% SubjectsStruct.(subjectID).Gender = Gender;
% 
% %% 17 Now save
% cd /Volumes/purplab/EXPERIMENTS/1_Current_Experiments/Caroline/Caroline2/DPF_V2_(all)/DPFv2Scripts/DPFv2_ExtractData_Excel
% %SubjectsStruct2 = SubjectsStruct;
% save('SubjectsStruct.mat','SubjectsStruct');