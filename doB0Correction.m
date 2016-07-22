% doB0Correction.m
%
%      usage: doB0Correction(v, magVol, phaseVol)
%         by: eli merriam
%       date: 03/13/14
%    purpose: applies field map correctin using the FreeSurfer script
%    epidewarp.fsl, which simply provides a nice wrapper for FSL's
%    PRELUDE and FUGE
%
% 1) convert magnitude volume from field map scans:  mri_convert ../20130619/7.gre_field_mapping/*0001.dcm -o Etc/run07.nii
% 2) convert phase difference volume:                mri_convert ../20130619/8.gre_field_mapping/*0001.dcm -o Etc/run08.nii
% 3) run doB0Correction on data in the Raw group:    v = newView; doB0Correction(v, 'Etc/run07.nii', 'Etc/run08.nii');
%
% for now, assume tediff = 2.46ms and esp = 0.77 ms (but we still need to confirm)
%

function retval = doB0Correction(v,magVol,phaseVol)

% check arguments
if ~any(nargin == [3])
  help doB0Correction
  return
end

tediff = 2.46;
esp = 0.77;
deWarppedGroupName = 'deWarped2';

% check that the inputs are legin
if ~isfile(magVol)
  disp(sprintf('The magnitude volume (%s) does not exist', magVol));
  return;
end

if ~isfile(magVol)
  disp(sprintf('The phase volume (%s) does not exist', phaseVol));
  return;
end

if system('which epidewarp.fsl')
  disp(sprintf('Can not find epidewarp.fsl in path, install FreeSurfer and check path'));
  return;
end

% create a new group
v = viewSet(v, 'newGroup', deWarppedGroupName);

v = viewSet(v, 'curGroup', deWarppedGroupName);

% number of scans in the Raw group
nScans = viewGet(v, 'nscans', 'Raw');

for iScan=1:nScans
  iFile = fullfile(viewGet(v, 'tseriesdir', 'Raw'), viewGet(v, 'tseriesfile', iScan, 'Raw'));
  evalstr = sprintf('epidewarp.fsl --mag %s --dph %s --epi %s --tediff %f --esp %f --vsm vsm.nii --epidw temp.nii --cleanup', ...
                magVol, phaseVol, iFile, tediff, esp);
  disp(sprintf('%s', evalstr));
  [status, result] = system(evalstr);
  if status
    disp(sprintf('uhoh, epidewarp.fsl failed...exiting'));
    return;
  end
  
  % save tseries
  scanParams = viewGet(v,'scanParams',iScan, 'Raw');
  scanParams.description = ['deWarpped ' scanParams.description];
  scanParams.fileName = [];
  scanParams.originalFileName{1} = viewGet(v,'tseriesfile',iScan);
  scanParams.originalGroupName{1} = 'Raw';
  [d,h] = cbiReadNifti('temp.nii');
  
  [v,tseriesFileName] = saveNewTSeries(v,d,scanParams,scanParams.niftiHdr);

  % % save evalstring for recomputing 
  % [pathstr, filename] = fileparts(tseriesFileName);
  % tseriesdir = viewGet(v, 'tseriesdir', deWarppedGroupName);
  % save(fullfile(tseriesdir,tseriesFileName),'evalstr','tseriesFileName');

  % clean up
  system(sprintf('rm -f vsm.nii.log.bak vsm.nii temp.nii vsm.nii.log'));
end


%We keep losing stimfiles in this and it's making me insane - here's a fix
%that pulls from the Raw group and sets it to the current
if ~isstruct('v')
    v=newView;
end
nScans = viewGet(v,'nScans',deWarppedGroupName);
v = viewSet(v, 'curGroup', deWarppedGroupName);

for iScan = 1:nScans
    %Set ourselves to the Raw group and get its groupnumber
  v = viewSet(v, 'curGroup', 'Raw');
   groupNum = viewGet(v,'currentGroup');
%Get stimfile name for this particular scan in the raw group
  stimfilename = viewGet(v,'stimFileName',iScan,groupNum);
  %Change our view to the dewarpped group
  v = viewSet(v, 'curGroup', deWarppedGroupName);
  %Set stimfile to be the same as its Raw correlate
  %Need to split out the stimfilename output into fileparts - we apparently
  %automatically append directory
  if isstruct('stimfilename')
      if ~isempty(stimfilename{1})
          [a b c]=fileparts(stimfilename{1});
            viewSet(v,'stimfilename',strcat(b,c),iScan,deWarppedGroupName);

      end
  end
end

