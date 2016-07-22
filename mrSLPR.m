%%Going to write a program that analyzes current PRF data for the response
%%of every voxel (say in an ROI if we want) to a single specified location
%%in space. That location will be determined by central coordinates and a
%%radius, and the response to the voxel at the time point at which there
%%was a stimulus present there. That will mean convolving with an HRF to
%%get the timing right. 

function [view d] = mrSLPR(view, params, varargin)

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end
if ieNotDefined('view'),view = newView;end
if ieNotDefined('xLoc'),xLoc = 0;end
if ieNotDefined('yLoc'),yLoc = 0;end
if ieNotDefined('radius'),radius = 2.5;end
if ieNotDefined('ROI'),ROI = [];end
if ieNotDefined('scanNum'),scanNum = 1;end
if ieNotDefined('doGLM'),doGLM = 1;end


 
% just return parameters
if justGetParams
  d = params;
  return
end

% get the group names
if ieNotDefined('groupName')
  groupName = viewGet(view,'groupNames');
else
  % if passed in name, put that on top of list to make it the default
  groupName = putOnTopOfList(groupName,viewGet(view,'groupNames'));
end


%Run  gui 

if ieNotDefined('stimParams')
    % setup params dialog
    slprParams = {};
    slprParams{end+1} = {'groupName',groupName,'Name of group from which to do eventRelated analysis'};
    slprParams{end+1} = {'scanNum',scanNum,'minmax=[1 inf]','incdec=[-1 1]','Which scan in the group to run on'};
    slprParams{end+1} = {'saveName','SLPR','File name to try to save as'};
    slprParams{end+1} = {'ROI',[],'Limit to ROI, leave empty to run on full scan'};
    slprParams{end+1} = {'xLoc',xLoc,'minmax=[-inf inf]','incdec=[-1 1]','Center of stimulus, X coordinate'};
    slprParams{end+1} = {'yLoc',yLoc,'minmax=[-inf inf]','incdec=[-1 1]','Center of stimulus, Y coordinate'};    
    slprParams{end+1} = {'radius',radius,'minmax=[0 inf]','incdec=[-1 1]','Stimulus radius'};
    slprParams{end+1} = {'cutoffPercent',0.1,'minmax=[0 inf]','incdec=[-.1 .1]','Cutoff percent of stimulus needed to activate RF'};
    slprParams{end+1} = {'doGLM',doGLM,'type=checkbox','Runs GLM in addition to deconvolution'};
    slprParams{end+1} = {'getStimParams',0,'type=checkbox','Whether to prompt for a stimfile to load, no need to specify stim params if clicked'};
end

% Display the initial dialogue to get some important information no
% matter what we're doing
initParams = mrParamsDialog(slprParams,'Choose analysis parameters');
params = initParams;


% if empty user hit cancel
if isempty(initParams)
  deleteView(view);
  return
end

workingdir=viewGet(view,'homedir');

if initParams.getStimParams==1
   stimPath = getPathStrDialog(workingdir,'Specify adaptation stimfile');
   load(stimPath, 'stimulus');
   fprintf(sprintf('\n Loading data from %s \n', stimPath))
   stimParams.x = stimulus.dots.xcenter;
   stimParams.y = stimulus.dots.ycenter;
   stimParams.r = stimulus.dots.rmax;
   clear stimulus
else
    
   stimParams.x = initParams.xLoc;
   stimParams.y = initParams.yLoc;
   stimParams.r = initParams.radius;  
end

params.stimParams=stimParams;
%Add stimulus location and radius to savename
initParams.saveName = strcat(initParams.saveName,'_',num2str(stimParams.x),'_',num2str(stimParams.y),'_',num2str(stimParams.r));


% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = mrParamsReconcile([],params);


if ~isfield(initParams,'groupName') || isempty(initParams.groupName)
    mrWarnDlg('No group set, defaulting to Concatenation');
    params.groupName = 'Concatenation';
end

view=viewSet(view, 'curGroup', initParams.groupName);
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',initParams.groupName);
if (groupNum ~= curGroup)
	mrWarnDlg(['Changing view to group: ',groupName]);
	view = viewSet(view,'currentGroup',groupNum);
end

% create the parameters for the overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = initParams.groupName;
r2.function = 'mrSLPR';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(view,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(view,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'eventRelatedPlot';
r2.mergeFunction = 'defaultMergeParams';


if initParams.doGLM
    % create the parameters for the overlay
    GLM.name = 'GLM';
    GLM.groupName = initParams.groupName;
    GLM.function = 'mrSLPR';
    GLM.reconcileFunction = 'defaultReconcileParams';
    GLM.data = cell(1,viewGet(view,'nScans'));
    GLM.date = dateString;
    GLM.params = cell(1,viewGet(view,'nScans'));
    GLM.range = [0 1];
    GLM.clip = [0 1];
    % colormap is made with a little bit less on the dark end
    GLM.colormap = hot(312);
    GLM.colormap = r2.colormap(end-255:end,:);
    GLM.alpha = 1;
    GLM.colormapType = 'setRangeToMax';
    GLM.interrogator = 'eventRelatedPlot';
    GLM.mergeFunction = 'defaultMergeParams';
end


%Extract raw time series
params.ROI = ROI;
if isempty(params.ROI)
   
    tSeries = loadTSeries(view);
    
else
    
    tSeries = loadROITSeries(view, params.ROI);
    
end

%Build specific stimvols

if ~isfield(initParams,'scanNum') || isempty(initParams.scanNum)
    mrWarnDlg('No scan number set, defaulting to scan 1');
    params.scanNum = 1;
    params.scanList=[1];
end

% get stimfiles
stimFile = viewGet(view,'stimfile',params.scanNum);

stimOut = pRFGetStimImageFromStimfile(stimFile);


for iScan = 1:length(stimOut)
    
    if length(stimOut{iScan})==1
        currStim = stimOut{iScan};
    else
        currStim = stimOut{iScan}{1};
    end
    
    %init vars
    stimOn=[];
    currIm=[];
    currCoords.x=[];
    currCoords.y=[];
    currCoords.xy=[];
    
    for iFrame = 1:length(currStim.t)
        
        %Extract image from the current frame, this is just a binary mask
        %telling us where the stimulus was, mapping onto actual degree
        %locations in .x and .y
        currIm = currStim.im(:,:,iFrame);
        
        %Make combined xy cell for isincircle
        currCoords.x = currStim.x(find(currIm==1));
        currCoords.y = currStim.y(find(currIm==1));
        currCoords.xy = {currCoords.x, currCoords.y};
        onPoints{iFrame}= pointsincircle(currCoords.xy, stimParams.r, [stimParams.x, stimParams.y]);
        
        %If there aren't any points in the circle move on
        % It may not be worth including a voxel if say, one point out of the
        % entire array is on. That is, if just one pixel out of the entire
        % screen enters the area, default will be .1% of voxels inside the stim
        % radius (which is a fair amount of stimulus, say about 1 degree)
        
        %If there are points in the circle just mark this frame as important
        %for our analysis
        
        if isempty(onPoints{iFrame}.in{1}) || length(onPoints{iFrame}.in{1}) < size(currStim.x,1)*size(currStim.x,2)*(params.cutoffPercent/100)
            stimOn(iFrame)=0;
        else
            stimOn(iFrame) = 1;
            
        end
        
        
        
    end
    
    %Store analysis 
    
    stimStore{iScan} = stimOn;
    stimTimes{iScan} = currStim.t(find(stimOn==1));
    
end

%Okay, now we have a timeseries for the analysis, and the stimtimes for
%every instance in which there was a stimulus crossing the area we're
%interested in

%Make two long stim traces
stimulusOn=[];
for iScans = 1:length(stimStore);
    stimulusOn=vertcat(stimulusOn, stimStore{iScans}');    
end
stimulusOff = (stimulusOn-1).^2;

%Next is basically a simple event related deconvolution at every voxel

frameperiod = viewGet(view, 'frameperiod');
%This analysis prefers its stimvolumes as framenumbers rather than just on
%or off

%Cell structure will contain just inverses of each other - stimulus on vs.
%stimulus off

stimvol{1} = find(stimulusOn==1);
stimvol{2} = find(stimulusOff==1);
stimNames{1} = 'On';
stimNames{2} = 'Off';


%This is our analysis loop - going to aim to do in in parallel for every
%voxel and that hopefully won't take too long. 

nVox = size(tSeries,1)*size(tSeries,2);

disp('Computing stimulus location population response...');
%initr2
allr2 = zeros(size(tSeries,1), size(tSeries, 2), size(tSeries,3));

waitHandle = mrWaitBar([0],['Progress'])

scanNum = params.scanNum;
numSlices = viewGet(view,'nSlices',scanNum);
  numVolumes = viewGet(view,'nFrames',scanNum);
  dims = viewGet(view,'dims',scanNum);
[numSlicesAtATime rawNumSlices numRowsAtATime precision] = getNumSlicesAtATime(numVolumes,dims);
  currentSlice = 1;
  ehdr = [];ehdrste = [];thisr2 = []; thisGLMr2=[];

    for i = 1:ceil(numSlices/numSlicesAtATime)
% calculate which slices we are working on
    thisSlices = [currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)];
    % set the row we are working on
    currentRow = 1;
    % clear variables that will hold the output for the slices we
    % are working on
    sliceEhdr = [];sliceEhdrste = [];sliceR2 = []; sliceGLMr2 = [];
    for j = 1:ceil(dims(1)/numRowsAtATime)
        % load the scan
        thisRows = [currentRow min(dims(1),currentRow+numRowsAtATime-1)];
        d = loadScan(view,scanNum,[],thisSlices,precision,thisRows);
        %Add stimvols to d structure
        d.stimvol = stimvol;
        d.stimNames = stimNames;
        % make a stimulation convolution matrix
        d = makescm(d,ceil(24/d.tr),1);
        % compute the estimated hemodynamic responses
        d = getr2(d);
        if initParams.doGLM
            %
            %Now do GLM if we're running it
            %GLM Specific settings
            d.designSupersampling = 1;
            for iRun=1:length(stimvol)
                d.stimDurations{iRun} = ones(size(stimvol{iRun}));
            end
            verbose = 1;
            %Make HRF
            [params, d.hrf] = hrfDoubleGamma([], frameperiod, [], 'defaultParams=1');
            params.spatialSmoothing=0;
            params.TFCE=0;
            %Design matrix
            glmd = makeDesignMatrix(d,[],verbose, scanNum);
            %Do analysis            
            [dGLM, out] = getGlmStatistics(glmd, params, verbose, precision, 1);%, computeTtests,computeBootstrap);
            
        end
        % update the current row we are working on
        currentRow = currentRow+numRowsAtATime;
        % if we are calculating full slice, then just pass that on
        if numRowsAtATime == dims(1)
            sliceEhdr = d.ehdr; 
            sliceEhdrste = d.ehdrste;
            sliceR2 = d.r2;
            
            if initParams.doGLM
               sliceGLMr2 =  out.r2;
            end
            
            % working on a subset of rows, cat together with what
            % has been computed for other rows
        else
            sliceEhdr = cat(1,sliceEhdr,d.ehdr);
            sliceEhdrste = cat(1,sliceEhdrste,d.ehdrste);
            sliceR2 = cat(1,sliceR2,d.r2);
            
            if initParams.doGLM
               sliceGLMr2 = cat(1, sliceGLMr2, out.r2); 
            end
            
        end
    end
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,sliceEhdr);
    ehdrste = cat(3,ehdrste,sliceEhdrste);
    thisr2 = cat(3,thisr2,sliceR2);
    if initParams.doGLM
        thisGLMr2 = cat(3, thisGLMr2, sliceGLMr2);
    end
    end
    
    % now put all the data from all the slices into the structure
    d.ehdr = single(ehdr);
    d.ehdrste = single(ehdrste);
    d.r2 = single(thisr2);
    if initParams.doGLM
        d.GLMr2 = single(thisGLMr2);
    end
  
  % get the actual size of the data (not just the size of the last
  % slice/set of rows we were working on).
  d.dim(1:3) = size(d.r2);
  
  
  % save the r2 overlay
  r2.data{scanNum} = d.r2;
  r2.params{scanNum} = initParams.scanNum;
 
  GLM.data{scanNum} = d.GLMr2;
  GLM.params{scanNum} = initParams.scanNum;

  % save other  parameters
  SLPR.d{scanNum}.ver = d.ver;
  SLPR.d{scanNum}.filename = d.filename;
  SLPR.d{scanNum}.filepath = d.filepath;
  SLPR.d{scanNum}.dim = d.dim;
  SLPR.d{scanNum}.ehdr = d.ehdr;
  SLPR.d{scanNum}.ehdrste = d.ehdrste;
  SLPR.d{scanNum}.nhdr = d.nhdr;
  SLPR.d{scanNum}.hdrlen = d.hdrlen;
  SLPR.d{scanNum}.tr = d.tr;
  SLPR.d{scanNum}.stimvol = d.stimvol;
  SLPR.d{scanNum}.stimNames = d.stimNames;
  SLPR.d{scanNum}.scm = d.scm;
  SLPR.d{scanNum}.expname = d.expname;
  SLPR.d{scanNum}.fullpath = d.fullpath;
  
  
  
mrCloseDlg(waitHandle);

% install analysis
SLPR.name = initParams.saveName;
SLPR.type = 'SLPR';
SLPR.range =[0 1];
SLPR.groupName =initParams.groupName;
SLPR.function = 'mrSLPR';
SLPR.reconcileFunction = 'defaultReconcileParams';
SLPR.mergeFunction = 'defaultMergeParams';
SLPR.guiFunction = 'None';
SLPR.params = params;
if initParams.doGLM
    SLPR.overlays = [r2 GLM];
else
    SLPR.overlays = r2;
end
SLPR.curOverlay = 1;
SLPR.date = dateString;
view = viewSet(view,'newAnalysis',SLPR);
if ~isempty(viewGet(view,'fignum'))
  refreshMLRDisplay(viewGet(view,'viewNum'));
end

% Save it
saveAnalysis(view,SLPR.name);

set(viewGet(view,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(d)
    SLPR.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(SLPR.d) == 1
    d = SLPR.d{1};
  else
    d = SLPR.d;
  end
end

function [points] = pointsincircle(xydata,radius,hk)
% POINTSINCIRCLE Identify points lying inside a circle

%   COMPLETELY VECTORIZED
%   Version 1.1 (Feb 29, 2012, 9:02 pm)
%   Version 1.0 (Feb 23, 2012) of this function was not uploaded

%   Find all points within a specified cut-off radius

%   Copyright, Sunil Anandatheertha Feb 23, 2012
%   http://soundcloud.com/sunilanandatheertha

%   INPUT:
%   xydata: a cell in the form xydata = {XvaluesOfCoordinates YvaluesOfCoordinates}.
%   radius: cut-off distance / cut-off radius
%   hk - h and k are the starting point of the radius vector.
%   OUTPUT:
%   points: this is a structure containing upto 2 values. 
%             points.in ---- cell of x- and y- coordinate data of all points inside the cut-off
%                                   radius
%             points.out ---- cell of x- and y- coordinate data of all points outside the cut-off
%                                   radius

dist = sqrt((hk(1) - xydata{1}(:)).^2   +   (hk(2) - xydata{2}(:)).^2); % distance calc.

% Find points inside circle
in = find(dist<radius);
points.in = {xydata{1}(in) xydata{2}(in)};

% Find points outside circle
% comment it out if its not needed
% its needed to run the testrun.m file though !!
out = find(dist>radius);  
points.out = {xydata{1}(out) xydata{2}(out)};

% ADDITIONAL --- USE TESTRUN.M TO SEE HOW TO USE THIS CODE



