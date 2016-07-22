% mrGLMdenoise.m
%
%      usage: mrGLMdenoise(v, varargin)
%         by: Michael Melnick
%       date: 06/05/16
%    purpose: Run GLMdenoise from Kendrick Kay's lab on current data
%
function v = mrGLMdenoise(v, stimVolFile, scanStart, scanStop, stimDur, varargin)

% % check arguments
% if ~any(nargin == [2:10])
%   help CB_adapat_anal
%   return
% end
% 
% get the input arguemnts

%%May add functionality to grab stimVols from MLR if it isn't specified in
%%a file...this may or may not be worth the time

addpath(genpath('~/GLMdenoise'))
getArgs(varargin, [], 'verbose=0');
if ieNotDefined('v'); v=newView;end
if ieNotDefined('scanNum'); scanNum = 1;end
if ieNotDefined('groupNum'); groupName = 'MotionComp';end
if ieNotDefined('stimVolFile'); stimVolFile = [];end

if ~isempty(stimVolFile)
   
    load(fullfile('Etc', stimVolFile));
    
else
    
    mrWarnDlg('No stimVolFile specified, try using mrDesignExport to grab it')
    
end

%load up fMRI data

for iScan = scanStart:scanStop
    
    v=viewSet(v, 'curGroup', groupName);
    groupNum = viewGet(v, 'curGroup');
    v = viewSet(v, 'curScan', iScan);
    frameperiod = viewGet(v, 'frameperiod');
    fname = viewGet(v,'tseriesPath',iScan,groupNum);
    funcData{iScan-scanStart+1} = cbiReadNifti(fullfile(fname),[],'single');
    nanreplace{iScan-scanStart+1} = isnan(funcData{1});
    funcData{iScan-scanStart+1}(nanreplace{iScan-scanStart+1})=0; 
    
    scanLength(iScan-scanStart+1) = size(funcData{iScan-scanStart+1},4);
    scansOut{iScan-scanStart+1} = funcData{iScan-scanStart+1}(:,:,:,1:end);
    
    condLength(iScan-scanStart+1) = length(Condition{iScan-scanStart+1});
end





[glmResults glmDenoisedata] = GLMdenoisedata(Condition, scansOut, stimDur, frameperiod);

glmDenoiseGroupName = 'DenoisedMotionComp';
% v = viewSet(v, 'newGroup', glmDenoiseGroupName);

% Open new view and set its group to the glmDenoise group name. Create the
% group if necessary.
viewGLM = newView;
glmDenoiseGroupNum = viewGet(viewGLM,'groupNum',glmDenoiseGroupName);
if isempty(glmDenoiseGroupNum)
  v = viewSet(v,'newgroup',glmDenoiseGroupName);
  glmDenoiseGroupNum = viewGet(viewGLM,'groupNum',glmDenoiseGroupName);
end
viewGLM = viewSet(viewGLM,'currentGroup',glmDenoiseGroupName);

for iScan=scanStart:scanStop
    iScan
    iFile = fullfile(viewGet(v, 'tseriesdir', 'MotionComp'), viewGet(v, 'tseriesfile', iScan, 'MotionComp'));
    % save tseries
    scanParams = viewGet(v,'scanParams',iScan, 'MotionComp');
    scanParams.description = ['glmDenoised ' scanParams.description];
    scanParams.fileName = [];
    scanParams.originalFileName{1} = viewGet(v,'tseriesfile',iScan);
    scanParams.originalGroupName{1} = 'MotionComp';    
    [viewGLM,tseriesFileName] = saveNewTSeries(viewGLM,glmDenoisedata{iScan-scanStart+1},scanParams,scanParams.niftiHdr);
    %Add a pause so we don't overwrite previous scan
%     disppercent(-inf,'Waiting for time counter');
%     while str2num(datestr(now,'SS'))~=0
%         disppercent((str2num(datestr(now,'SS'))/60), 'Waiting for time counter');
%     end
%     elapsedTime = disppercent(inf);
%     saveScanNum = viewGet(viewGLM,'nScans');   
        %Pause so we don't get identical names
        pause(1)
end

%We keep losing stimfiles in this and it's making me insane - here's a fix
%that pulls from the Raw group and sets it to the current
if ~isstruct('v')
    v=newView;
end

nScans = viewGet(v,'nScans',glmDenoiseGroupName);
v = viewSet(v, 'curGroup', glmDenoiseGroupName);

for iScan = scanStart:scanStop
    %Set ourselves to the Raw group and get its groupnumber
  v = viewSet(v, 'curGroup', 'Raw');
   groupNum = viewGet(v,'currentGroup');
%Get stimfile name for this particular scan in the raw group
  stimfilename = viewGet(v,'stimFileName',iScan,groupNum);
  %Change our view to the dewarpped group
  v = viewSet(v, 'curGroup', glmDenoiseGroupName);
  %Set stimfile to be the same as its Raw correlate
  %Need to split out the stimfilename output into fileparts - we apparently
  %automatically append directory
  if isstruct('stimfilename')
      if ~isempty(stimfilename{1})
          [a b c]=fileparts(stimfilename{1});
            viewSet(v,'stimfilename',strcat(b,c),iScan-scanStart+1,glmDenoiseGroupName);

      end
  end
end

