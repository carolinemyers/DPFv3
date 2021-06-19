% need to edit so that requirement is not trial.end == 48, since now
% includes fix break trials so size of block is often > 48

function [TestBlocks]=DPF_ExoAttnParseFiles(Observer); 
eval(sprintf('Folder=dir(''data/%s'');',Observer));
if size(Folder,1)
FileNum=length(Folder); %for all items
for File=3:FileNum %excluding junk rows at top
    DataIndex(File)=strcmp(Folder(File).name(end-2:end),'mat'); %creates an array of every file that's .mat
end
Counter=0;
for File=3:FileNum
    if DataIndex(File) && ~strcmp(Folder(File).name(1),'.')
        Counter=Counter+1;
        eval(sprintf('load(''data/%s/%s'')',Observer,Folder(File).name));
     
        %now let's go through some checks. If the taskFileName (the task type)
     %is THRESHOLDING ('DPF_ExoAttnContrastThresh.m' (length = 27) don't
     %count it. If the taskFileName is 'Eyetrack_DPF_ExoAttnExp.m' (length
     %= 25) and the total number of trials = 48 (+ the trials where there
     %was a fixation break) count it! 
       
      if length(task{1, 1}.taskFilename) == 16 && stimulus.trialend==48 || length(task{1, 1}.taskFilename) == 31 && stimulus.trialend==48
          FullBlockIndex(Counter)=1;
          
      %If the taskFileName is 'Eyetrack_DPF_ExoAttnExp.m' (length
      %= 25) and the total number of trials = 48 (+ the trials where there
      %was a fixation break) count it!     
      elseif length(task{1, 1}.taskFilename) == 25 && stimulus.trialend==48+sum(stimulus.FixationBreak)
          FullBlockIndex(Counter)=1; 
      else
          FullBlockIndex(Counter)=0;
      end
        FileType{Counter,1}=task{1}.taskFilename; 
        FileType{Counter,2}=Folder(File).name(1:6); % places six digits of date (after file transformation to initials_XX for MS analyis
        FileType{Counter,3}=Folder(File).name(12:13); % places block number overall blocks of exp (including edf and non-main block)
        BlockPrep(Counter,1)=str2num(Folder(File).name(1:6));
        BlockPrep(Counter,2)=str2num(Folder(File).name(12:13));
    end
end
for i=1:length(FileType)
    if FullBlockIndex(i)
        checkFunction=strcmp(FileType{i,1},'DPF_ExoAttnExp.m')||strcmp(FileType{i,1},'DPF_ExoAttnExp_ExtendedTiming.m') ||strcmp(FileType{i,1},'Eyetrack_DPF_ExoAttnExp.m') ;
    TestIndex(i)=checkFunction;
    end
end
TestBlocks=BlockPrep(TestIndex,:);
else
    fprintf('Participant not found. Check Observer initials.\r')
TestBlocks=[];
end

