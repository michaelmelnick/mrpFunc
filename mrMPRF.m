function mrMPRF(subjectinitials,nSessions,workingdir,forceBugOff)


clear all; clc;



%mrMPRF
%The goal of this is to be able to combine multiple sessions and run the
%PRF on them, including our workaround for the 167 frame bug. Alignment
%will still need to be done manually, this is going to assume that you've
%downloaded and run mrSanity on both sessions, so that MotionComp and
%anatomical groups already exist


if ieNotDefined('forceBugOff') forceBugOff = 1;end

if ieNotDefined('subjectinitials') subjectinitials = '';end
if ieNotDefined('nSessions') nSessions = 2;end
if ieNotDefined('workingdir') workingdir = '/mri_work/';end

% setup params dialog
mprfParams = {};
mprfParams{end+1} = {'Subject Initials',subjectinitials,'type=string','Subjects initials, needs to match an existing set, all runs should be inside that set'};
mprfParams{end+1} = {'nSessions',nSessions,'minmax=[1 inf]','incdec=[-1 1]','Number of sessions to combine'};
mprfParams{end+1} = {'Working Directory',workingdir,'type=string','Base directory (usually /mri_work/)'};

% Display the initial dialogue to get some important information no
% matter what we're doing
initParams = mrParamsDialog(mprfParams,'Choose group parameters');

workingdir=strcat(initParams.Working_Directory,'/',initParams.Subject_Initials,'/');

%Get ready to present dialogue needed for each session

sessParams={};    

%Populate list with possible folders in initial folder, but accept paths to
%anywhere

cd(workingdir)


%List directories in current directory (subject folder) that start with 2
%(for dates starting with 20xx), sadly, this program won't work in 1000
%years
dList = dir('2*');

%Convert to cell

for pp = 1:length(dList)
   dCell{pp} = dList(pp).name; 
end

%sessParams{end+1} = {'Condition',condition,'type=popupmenu','Type of scan collected, only 3 options for now'};


for pp = 1:initParams.nSessions
    
        currSess = strcat('Session',num2str(pp));
        sessParams{end+1} = {currSess,dCell,'type=popupmenu',strcat('Full Date of Session ', num2str(pp))};
        
end

sessParams{end+1} = {'newDir','','type=string','New directory for this analysis to create'};

sessParams = mrParamsDialog(sessParams,'Choose session parameters');

%Check if the directory we want to make exists, if it doesn't create it, if
%it does exit with error

makeDir = exist(sessParams.newDir);

if makeDir ~= 7
%    mkdir(sessParams.newDir) 
%    mkdir('Raw');
%    mkdir('Raw/TSeries');
%    mkdir('Etc');
%    mkdir('Anatomy');
   %Make this a mlrDir

   makeEmptyMLRDir(sessParams.newDir,'description=PRF_Combined', strcat('subject=',initParams.Subject_Initials),'operator=PW','defaultParams=1');
      cd(sessParams.newDir)


else
   mrWarnDlg('New analysis directory already exists, assuming we want to continue in here')
   cd(sessParams.newDir);    
end

%CD to new directory and begin import of necessary files



%Copy Raw files then run mrInit

%List directories for raw data
rawFiles={};
for pp = 1:initParams.nSessions
    currSess = strcat('Session',num2str(pp));
    rawFiles{end + 1} = fullfile(workingdir, eval(sprintf('sessParams.%s', currSess)));
end

% %Copy raw files
% for pp = 1:initParams.nSessions
%     
%     evalstr = sprintf('cp -f %s/* %s', ...
%         fullfile(rawFiles{pp},'Raw/TSeries'), fullfile(workingdir, sessParams.newDir, 'Raw','TSeries'));
%     disp(sprintf('%s', evalstr));
%     [status, result] = system(evalstr);
% end

%Copy Anat files

for pp = 1:initParams.nSessions
    
    evalstr = sprintf('cp -f %s/* %s', ...
        fullfile(rawFiles{pp},'Anatomy'), fullfile(workingdir, sessParams.newDir, 'Anatomy'));
    disp(sprintf('%s', evalstr));
    [status, result] = system(evalstr);
    

    %Get rid of leftover .mat alignment files, etc.
    evalstr = sprintf('rm %s/*.mat', ...
        fullfile(workingdir, sessParams.newDir, 'Anatomy'));
    disp(sprintf('%s', evalstr));
    [status, result] = system(evalstr);
    
    dirOut = dir('Anatomy/run*');
    %We need to rename our anatomicals - this way there's no chance of
    %overwriting, and we know which scan they come from 
    
    for ss = 1:length(dirOut)
        currSess = strcat('Session',num2str(pp));
        fileExt = dirOut(ss).name(end-3:end);
        
        nameout = eval(sprintf('sessParams.%s', currSess));
        
        fullName = strcat(nameout, fileExt);
        
        evalstr = sprintf('mv Anatomy/%s Anatomy/%s', ...
            dirOut(ss).name, fullName);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr);
        
        

    end
    
    
end

v = newView;

%We're going to make two new local groups, then call MLR's built in
%function importGroupScans - this isn't automated, and will ask you to
%select the runs you want to grab from. Not worth the time to automate /
%steal code for this right now

for pp = 1:initParams.nSessions
    v = viewSet(v,'newGroup', strcat('MotionComp',num2str(pp)));
end

for pp = 1:initParams.nSessions

    importGroupScans
    
end


%Okay, now we have the entire motioncomp group from whatever number of
%sessions we wanted to combine. We MUST do an alignment here or further
%analysis is impossible

mrAlign

%Once alignment is finished, we can make averages, concatenation and
%finally PRF analysis, but we are going to work in the automated 167 frame
%bug workaround from mrSanity

keyboard



workingfile=pwd;

%We need to look at the motioncomped files for nFrames, etc., this will be
%done iteratively to treat each possible session separately


for pp = 1:initParams.nSessions
    
    %currSess = strcat('Session',num2str(pp));
    %nameout = eval(sprintf('sessParams.%s', currSess));
    nameout = strcat('MotionComp',num2str(pp));
    
    
    % get info about scans that live in Raw/TSeries
    scanDirName = fullfile(workingfile,nameout,'TSeries');
    scanFilenames = mlrImageGetAllFilenames(scanDirName,'mustNotHaveDotInFilename=1');
    nScans=0;scanNames = {};descriptions = {};totalFrames = {};nFrames = {};junkFrames = {};
    
    for i = 1:length(scanFilenames);
        % read the nifti header
        imageHeader = mlrImageHeaderLoad(fullfile(scanDirName,scanFilenames{i}));
        if imageHeader.nDim == 4
            nScans = nScans+1;
            scanNames{end+1} = scanFilenames{i};
            descriptions{end+1} = '';
            totalFrames{end+1} = imageHeader.dim(4);
            nFrames{end+1} = imageHeader.dim(4);
            junkFrames{end+1} = 0;
        else
            disp(sprintf('(mrInit) File %s is not a 4D Nifti file',scanFilenames{i}));
            
            
            
            
        end
    end
    
    [stimFileNames stimFileVols] = getStimFiles(workingfile);
    stimFileMatch = matchStimFiles(stimFileNames,stimFileVols,totalFrames);
    
    for i=1:nScans
        %Do the length of the stimfiles equal the length of the scans?
        equalvols(i,pp)=stimFileVols(i)==cell2mat(nFrames(i));
    end
    
end




mrWarnDlg('Looks like we may have run into the 167/168 frame bug in bar code - going to add a flag to make sure we start with a scan of 167, should work as a workaround for now')
pRFBarBug=1;
cd(workingfile)
svfile=strcat('Etc/pRFBarBug');
save(svfile,'pRFBarBug','stimFileVols','nFrames');

load('Etc/pRFBarBug.mat');

%Okay, let's grab the total number of conditions we need to check (i.e.
%averages we'll ultimately need to make). This means we need to know which
%scans were rings, wedges, bars, clockwise, etc. 

%Iterate for number of sessions
for pp = 1:initParams.nSessions
    
    nameout = strcat('MotionComp',num2str(pp));
    v=viewSet(v,'curgroup', nameout);
    nScans(pp) = viewGet(v,'nScans',nameout);
    
    groupNum = viewGet(v,'currentGroup');
    for iScan = 1:nScans(pp)
        
        %Get stimfile name for this particular scan in the motioncomp group
        stimfilename = viewGet(v,'stimFileName',iScan,groupNum);
        
        %Only try to grab the stim parameter if we have a stimfile there
        if ~isempty(stimfilename)
            barWidth(iScan,pp) = getParamfromStimfile(iScan,stimfilename,'stimulus.barWidth');
            stimulusType(iScan,pp) = getParamfromStimfile(iScan,stimfilename,'stimulus.stimulusType');
            stimDirection(iScan,pp) = getParamfromStimfile(iScan,stimfilename,'stimulus.direction');
        end
        
    end
    
end

stimFileVols_orig = stimFileVols;
totalBuggedSess = 0;
totalBuggedAffectedScans = 0;
totalSess = 0;
totalBuggedScans =0;
%Now we have the conditions that were run, how many averages will we need per group?

for pp = 1:initParams.nSessions
    
    stimWidths = unique(barWidth(:,pp));
    stimTypes = unique(stimulusType(:,pp));
    stimDirs  = unique(stimDirection(:,pp));
    
    %Get rid of zeros
    stimWidths(find(stimWidths==0))=[];
    stimTypes(find(stimTypes==0))=[];
    stimDirs(find(stimDirs==0))=[];

    
    cLengths = [length(stimWidths) length(stimTypes) length(stimDirs)];
    cVars = {'stimWidths', 'stimTypes', 'stimDirs','barWidth','stimulusType','stimDirection'};
    iConditions(pp) =  prod(cLengths);

    condWidths = repmat(stimWidths',1,iConditions(pp)/length(stimWidths));stimFileVols_orig
    condTypes = repmat(stimTypes', 1, iConditions(pp)/length(stimTypes));
    condDirs = repmat(stimDirs', 1, iConditions(pp)/length(stimDirs));
    
    condMat = vertcat (condWidths, condTypes, condDirs);
    condMat(3,:)=sort(condMat(3,:));
    condBugged{pp} = zeros(iConditions(pp),1);

    
    for iConds=1:iConditions(pp)
        
        if pp ==1
            stimFileVols = stimFileVols_orig(1:nScans(pp));
        else
            stimFileVols =  stimFileVols_orig(nScans(pp-1)+1:(nScans(pp)+nScans(pp-1)));
        end
        
        
        SLW = find(barWidth(:,pp)==condMat(1,iConds));
        SLT = find(stimulusType(:,pp)==condMat(2,iConds));
        SLD = find(stimDirection(:,pp)==condMat(3,iConds));
        
        %Now produce composite list
        
        scanList = intersect(SLW,SLT);
        scanList = intersect (scanList, SLD);
        
        
        %If we add a change here, and mess with the scanlist, we might be able
        %to fix the 167/168 problem. Let's see.
        
        
        %Okay, we're in a condition with a scanlist, let's end this loop
        %simply with binary indicators of whether each conditions in the
        %entire set HAS the bug
        

        
        if any(stimFileVols(scanList)==167)
            
            condBugged{pp}(iConds)=1;
        else
            condBugged{pp}(iConds)=0;
            
        end
        
        scanLists(iConds, pp) = length(scanList);
        allSL{pp}{iConds} = scanList;
        
        
        totalBuggedSess  = totalBuggedSess+condBugged{pp}(iConds);
        totalSess = totalSess+1;
        if condBugged{pp}(iConds)
            
            totalBuggedAffectedScans = totalBuggedAffectedScans + scanLists(iConds, pp);
            
        end
        
        %Now we need a specific instance of each bugged scan and where to find it
        %For this condition
        
        buggedScans{pp}{iConds}=find(stimFileVols==167);
        buggedScans{pp}{iConds}=intersect(buggedScans{pp}{iConds},scanList);
        
        totalBuggedScans = totalBuggedScans + length(buggedScans{pp}{iConds}); 
        
        
       
    end
end

%Okay, to this point we know whether each averaged condition will contain
%the bug or not. Now we need to determine the least cost in terms of data.
%I.e., if we ditch scans without the bug to include those with it, we lose
%data, but do we lose less than if we skip all runs with the bug? 


%Do we have a greater number of bugged or unbugged sessions?

totalScans = sum(sum(scanLists));
aveLost = totalSess-totalBuggedSess;



        
        
sprintf('Looks like of %i averages, %i *include* the bug, so we would lose %i average(s)', totalSess, totalBuggedSess, totalSess-totalBuggedSess)
sprintf('and of %i total scans, %i include the bug, so we would lose %i scans', sum(sum(scanList)), totalBuggedScans, totalBuggedScans)         


if aveLost>totalSess/2 || forceBugOff==1
    sprintf('We would lose a greater number of averages by accounting for the bug, so we will ditch bugged scans')
    
    
    %Here we'll just delete bugged scans from the possibilities, include an
    %option here to override and just ditch bugged scans
    for pp = 1:initParams.nSessions
        
        for iConds=1:iConditions(pp)
            aveList{pp}{iConds} = allSL{pp}{iConds};
            
            for q=1:length(buggedScans{pp}{iConds})
                
                if ~isempty(buggedScans{pp}{iConds})
                    whichScan = find(aveList{pp}{iConds}==buggedScans{pp}{iConds}(q));
                else
                    whichScan = [];
                end
                
                if isempty(whichScan)
                    continue
                else
                    
                    aveList{pp}{iConds}(whichScan)=[];
                    
                    
                end
                
                
            end
            
            
        end
        
        
    end
    
    
else
    sprintf('We lose less averages by accounting for the bug, starting bug workaround')
    
    %We need to build scanLists that start with the bug for each scanList
    %then build averages, easy
    
    for pp = 1:initParams.nSessions
%             allBugged = []
% 
%         
%         
%         for q = 1:length(buggedScans{pp})
%             allBugged = horzcat(allBugged, buggedScans{pp}{q}')
%         end
            for iConds=1:iConditions(pp)

                
                %Transfer scalist to aveList
                aveList{pp}{iConds} = allSL{pp}{iConds};
                
                %Rearrange
                if ~isempty(buggedScans{pp}{iConds})
                    whichScan = find(aveList{pp}{iConds}==buggedScans{pp}{iConds}(1));
                else
                    whichScan = [];
                end
                %If whichScan is already scan 1 move on, else move the scan
                
                if whichScan==1
                    continue
                elseif isempty(whichScan)
                    aveList{pp}{iConds}=[]
                else
                    
                    tempVal1 = aveList{pp}{iConds}(1);
                    tempVal2 = aveList{pp}{iConds}(whichScan);
                    
                    aveList{pp}{iConds}(1) = tempVal2;
                    aveList{pp}{iConds}(whichScan) = tempVal1;
                    
                    
                end
                
                
            end
        
        
    end
    
    
    
end



%Only thing left is to actually build the averages which are now contained
%in the aveList structure for either debugged or bugfixed versions

for pp = 1:initParams.nSessions
    
    
    
    
    for iConds=1:iConditions(pp)
        
        if aveLost>totalSess/2 || forceBugOff==1
            
            if ~isempty(aveList{pp}{iConds})
                if exist('scans','var')
                    clear scans
                end
                
                v=newView;
                nameout = strcat('MotionComp',num2str(pp));
                v=viewSet(v,'curgroup', nameout);
                
                scans.scanParams = viewGet(v,'scanparams',aveList{pp}{iConds});
                
                for i = 1:length(scans.scanParams)
                    
                    
                    scans.scanParams(i).nFrames=168;
                    scans.scanParams(i).junkFrames=0;
                    
                    v=viewSet(v,'updatescan',scans.scanParams(i),aveList{pp}{iConds}(i))
                end
            end
            
        else
            
            if ~isempty(aveList{pp}{iConds})
                
                if exist('scans','var')
                    clear scans
                end
                
                v=newView;
                nameout = strcat('MotionComp',num2str(pp));
                v=viewSet(v,'curgroup', nameout);
                
                scans.scanParams = viewGet(v,'scanparams',aveList{pp}{iConds});
                
                
                for i = 1:length(scans.scanParams)
                    
                    
                    scans.scanParams(i).nFrames=160;
                    scans.scanParams(i).junkFrames=7;
                    
                    v=viewSet(v,'updatescan',scans.scanParams(i),aveList{pp}{iConds}(i))
                    
                end
            end
            
        end
        
        if ~isempty(aveList{pp}{iConds})
            if exist('params','var')
                clear params
            end
            
            scanStr = strcat('scanList=', '[',num2str(aveList{pp}{iConds}'),']');
            
            [v params] = averageTSeries(v,[],'justGetParams=1','defaultParams=1',scanStr);
            v=averageTSeries(v,params);
            deleteView(v);
            
        end
        
    end
end






%%%%%%%%%%%%%%%%%%%%%%
%  get stimFileNums  %
%%%%%%%%%%%%%%%%%%%%%%
function [stimFileNames stimFileVols] = getStimFiles(workingfile)

stimfileDir = strcat(workingfile,'/','Etc');
dirList = dir(stimfileDir);
stimFileNums = [];stimFileNames = {};stimFileVols = [];
for i = 1:length(dirList);
    name = dirList(i).name;
    % name should be of form yymmdd_stimnn.mat (nn is sequence number)
    if (regexp(name,'\d\d\d\d\d\d_\w\w\w\w\d\d.mat'))
        % get the stimfile number
        stimFileNums(end+1) = str2num(name(12:13));
        % try to load and get some info
        load(sprintf(fullfile(stimfileDir,name)));
        if ieNotDefined('myscreen')
            disp(sprintf('(mrInit) %s does not look like an mgl stimfile (does not have a myscreen variable)',name));
            continue;
        end
        % make a string of some info myscreen
        stimfileStr = '';
        if isfield(myscreen,'volnum')
            stimfileStr = sprintf('%s[%i vols] ',stimfileStr,myscreen.volnum);
            stimFileVols(end+1) = myscreen.volnum;
        else
            stimFileVols(end+1) = nan;
        end
        if isfield(myscreen,'starttime')
            stimfileStr = sprintf('%s%s ',stimfileStr,myscreen.starttime);
        end
        if isfield(myscreen,'endtime')
            stimfileStr = sprintf('%s(End: %s) ',stimfileStr,myscreen.endtime);
        end
        stimFileNames{end+1} = sprintf('%s: %s',name,stimfileStr);
    end
end

function varval = getParamfromStimfile(scan,stimFileName,vartoget)

%Load up a stimfile for this particular scan

filetoload=fullfile(stimFileName);

%If we're trying to load an empty file, skip this and return empty
if ~isempty(filetoload)
    load(filetoload{1});
    
    str=sprintf(vartoget);
    %Try to get the value, if it's empty and we toss an error set this to empty
    try
        varval=eval(str);
    catch
        varval = [];
    end
else
    varval=[];
    
end

function stimFileMatch = matchStimFiles(stimFileNames,stimFileVols,totalFrames)

stimFileMatch = {};
% for now, just match in order
for i = 1:length(stimFileNames)
    % make sure there are enough vols (just check to see if there
    % are greater than 1/10 the total frames) since with tSense
    % etc there may be more totalFrames than vols, so this is
    % just to reject stimfiles with only a few volumes
    if (length(stimFileVols) < i) || (length(totalFrames)<i) || (stimFileVols(i) > totalFrames{i}/10)
        stimFileMatch{end+1} = putOnTopOfList(stimFileNames{i},{stimFileNames{:},'None'});
        % only make enough matches for each scan
        if length(stimFileMatch)==length(totalFrames)
            break;
        end
    end
end

% all unmatched scans
for i = length(stimFileMatch)+1:length(totalFrames)
    stimFileMatch{end+1} = putOnTopOfList('None',cellcat(stimFileNames,{'None'}));
end



