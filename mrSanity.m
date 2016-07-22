%MRI Pre-mrInit conversions and sanity checks, let's call him mrSanity

%Things this program can do:

%Download subject data from RCBI given subject initials and date of scan
%SI, SD

%Accept a subject name and condition
%SN, cond
%Build a new directory structure in mri_work for this subject if it doesn't
%already exist (Anatomy, TSeries, Etc)

%Convert data to useable NIFTI formats using system's mri_convert. Please
%make sure we're in tcsh, this won't work otherwise

%Run up a quick sanity check of type of stimulus and the parameters run.
%This should help when running mrInit and our initial stages. Should answer
%if we ran bars, changed widths, had a weird run that was too short, etc.

function mrSanity(subjectinitials,scandate,workingdir,rawdir,RCBIFetch,MRBuildDirs,MRIConvert,DoDewarp)

clear all; clc;

%             [sessionParams groupParams] = mrInit([],[],'justGetParams=1')
%             mrInit(sessionParams,groupParams);


condition={'Ring and Wedge Retinotopy','pRF','Ring and Wedge Retinotopy - Mask','Adaptation'};
if ieNotDefined('subjectinitials') subjectinitials = '';end
if ieNotDefined('scandate') scandate = '';end
if ieNotDefined('workingdir') workingdir = '/mri_work/';end
if ieNotDefined('rawdir') rawdir = '/mri_raw/';end
if ieNotDefined('RCBIFetch') RCBIFetch = 0;end
if ieNotDefined('MRBuildDirs') MRBuildDirs = 0;end
if ieNotDefined('MRIConvert') MRIConvert = 0;end
if ieNotDefined('DoDewarp') DoDewarp = 0;end
if ieNotDefined('CheckStimfileMatch') CheckStimfileMatch = 0;end
if ieNotDefined('RunmrInit') RunmrInit = 0;end
if ieNotDefined('HandleStimFiles') HandleStimFiles=0; end
if ieNotDefined('DoMotionComp') DoMotionComp=0; end
if ieNotDefined('DoAverage') DoAverage=0; end
if ieNotDefined('doAutorecon') doAutorecon=0; end
if ieNotDefined('doConcat') doConcat=0; end
if ieNotDefined('doAnalysis') doAnalysis=0; end

Output={'niftii hdr/img','niftii all in one'};


% setup params dialog
sanParams = {};
sanParams{end+1} = {'Subject Initials',subjectinitials,'type=string','Subjects initials, needs to match RCBIs record to fetch the data'};
sanParams{end+1} = {'Scan Date',scandate,'type=string','Format is: yyyymmdd'};
sanParams{end+1} = {'Condition',condition,'type=popupmenu','Type of scan collected, only 3 options for now'};
sanParams{end+1} = {'Output',Output,'type=popupmenu','Default is Nifti hdr/img, can change to nifti all in one'};
sanParams{end+1} = {'Working Directory',workingdir,'type=string','Base directory (usually /mri_work/)'};
sanParams{end+1} = {'Raw Directory',rawdir,'type=string','Raw directory (usually /mri_raw/)'};
sanParams{end+1} = {'RCBIFetch',RCBIFetch,'minmax=[0 1]','incdec=[-1 1]','Set to 1 will attempt to grab the data from the RCBI server'};

sanParams{end+1} = {'MRBuildDirs',MRBuildDirs,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to build the basic directory structure used for mrTools'};
sanParams{end+1} = {'MRIConvert',MRIConvert,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to convert all data and place in directory structure'};
sanParams{end+1} = {'DoDewarp',DoDewarp,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to convert Gradient Field Mapping images and run dewarpping using FSL'};
sanParams{end+1} = {'CheckStimfileMatch',CheckStimfileMatch,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to match scans to stimfiles'};
sanParams{end+1} = {'RunmrInit',RunmrInit,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will pass initials, etc. to mrInit and move on to stimfile matching'};
sanParams{end+1} = {'HandleStimFiles',HandleStimFiles,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will check for stim files in the working directory, and copy them if necessary'};
sanParams{end+1} = {'DoMotionComp',DoMotionComp,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to run motion compensation on dewarpped or raw datafiles'};
sanParams{end+1} = {'DoAverage',DoAverage,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to average the proper datafiles together'};
sanParams{end+1} = {'doConcat',doConcat,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to concatenate average runs together'};
sanParams{end+1} = {'doAnalysis',doAnalysis,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to perform relevant analysis (pRF for pRF, corAnal for Retinotopy'};
sanParams{end+1} = {'doAutorecon',doAutorecon,'minmax=[0 inf]','incdec=[-1 1]','Set to 1 will attempt to run autoreconall for the Freesurfer side'};


% Display the initial dialogue to get some important information no
% matter what we're doing
groupParams = mrParamsDialog(sanParams,'Choose group parameters');

workingdir=strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date);
rawdir=strcat(groupParams.Raw_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date);



if groupParams.RCBIFetch
    
    useSimpleSCP=0; %This flag will use an alternate, JAVA based library to get around some connection issues
    
    if ieNotDefined('RCBIuser') RCBIuser = 'mmelnick';end
    if ieNotDefined('RCBIpass') RCBIpass = '';end
    if ieNotDefined('RCBIdatadir') RCBIdatadir = '/rcbiData/Huxlin/Visual_Recovery/';end
    if ieNotDefined('RCBIserver') RCBIserver = 'hippocampus.rcbi.rochester.edu';end
    
    
    fetchParams = {};
    fetchParams{end+1} = {'RCBIuser',RCBIuser,'type=string','Username that will allow us to login to the RCBI, default is mmelnick for SSH access'};
    fetchParams{end+1} = {'RCBIpass',RCBIpass,'type=string','Password that will allow us to login to the RCBI, cleared after run'};
    fetchParams{end+1} = {'RCBIdatadir',RCBIdatadir,'type=string','defaults to /rcbiData/Huxlin/Visual_Recovery/'};
    fetchParams{end+1} = {'RCBIserver',RCBIserver,'type=string','defaults to hippocampus'};
    
    %Run dialogue for the fetch sequence
    fetchParams = mrParamsDialog(fetchParams,'Choose group parameters');
    
    %try fetch using scp
    
    %What if the raw directory already exists? Probably recommend
    %you don't bother with redownloading it
    if exist(rawdir,'dir')
        disp(sprintf('Woah, looks like you already have a directory with the raw data, going to check if you have all of the files you should have in the next step!'))
        rawdirExist=1;
    else
        rawdirExist=0;
    end
        
        remotedir=strcat(RCBIdatadir,groupParams.Subject_Initials,'/',groupParams.Scan_Date);
        if strcmp(groupParams.Condition,'pRF')
            stimdir='PRF';
        elseif strcmp(groupParams.Condition,'Ring and Wedge Retinotopy')
            stimdir='Retinotopy';
        elseif strcmp(groupParams.Condition,'Ring and Wedge Retinotopy - Mask')
            stimdir='RetinotopyMask';
        elseif strcmp(groupParams.Condition,'Adaptation')
            stimdir='Adaptation';
        end
        remotestimdir=strcat('/rcbiUsers/',fetchParams.RCBIuser,'/Desktop/',stimdir,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date(3:8));
        altremotestimdir=strcat('/rcbiUsers/',fetchParams.RCBIuser,'/Desktop/',stimdir,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date(1:8));
        
        %Make the raw directory unless it already exists - if it does we'll
        %check for its contents in a second
        if ~rawdirExist
            mkdir(rawdir);
        end
        
        if ~useSimpleSCP
        evalstr = sprintf('sshpass -p %s scp -r %s@%s:%s %s', ...
            strcat('''',fetchParams.RCBIpass,''''), fetchParams.RCBIuser, fetchParams.RCBIserver,remotedir, rawdir);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr)
        
        evalstr = sprintf('sshpass -p %s scp -r %s@%s:%s %s', ...
            strcat('''',fetchParams.RCBIpass,''''), fetchParams.RCBIuser, fetchParams.RCBIserver,remotestimdir, rawdir);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr)
        
        evalstr = sprintf('sshpass -p %s scp -r %s@%s:%s %s', ...
            strcat('''',fetchParams.RCBIpass,''''), fetchParams.RCBIuser, fetchParams.RCBIserver,altremotestimdir, rawdir);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr)
        
        else
           
           %Quick sanity check to make sure we're getting a proper SSH
           %connection using the toolbox
           try
               sshTest = ssh2_simple_command(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass','ls -la *ninjas*',1);
               
           catch
               mrErrorDlg('Uh oh, ssh connection seems to have failed. Make sure you are  connected to the U of R network, and check any JAVA errors you get. This could just be that you picked the wrong filename, or typed your password wrong');
           end
           
           %Okay, now let's get what we need from the server - we need to
           %start by listing out every file using SSH, as this version of
           %SCP doesn't do a good job with recursive transfers...
           
           if exist('sshTest','var')
               %make a file list structure
               evalstr = sprintf('ls %s',remotedir);
               sshDirs = ssh2_simple_command(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass',evalstr,1);
               
               %We have each subfile of the fMRI structure now, what we
               %need is to download recursively inside each folder.
               
               %But first! Did we have the raw folder before? Do we need to
               %redownload everything? Did we stop in the middle of an old
               %run and need to restart?
               keyboard
               dirListB=dir(rawdir);

               %Need to get rid of directory artifacts
               delIndex=[];
               for ii=1:length(dirListB);
                
                   if ii>length(dirListB)
                       break
                   end
                   if strcmp('.',dirList(ii).name) || strcmp('..',dirList(ii).name)
                       delIndex=[delIndex ii];
                       continue
                   end
               end
               if any(delIndex)
                   dirListB(delIndex)=[];
               end
               
               
               %if we have an empty folder (i.e. have never downloaded
               %anything) then skip ahead and start downloading. If we have
               %a partially populated folder we have more work to do - what
               %folders do we have, and are they populated properly?
               
               if length(dirListB)
                   
                   mrWarnDlg('Okay, I havent implemented this yet...it will eventually check the local directory for parity to the server. Sometime soon! Maybe just delete the raw folder for now and start over');
                   
               else
                   
                   for ff=1:length(sshDirs)
                       %List files inside the current folder
                       evalstr = sprintf('ls %s',fullfile(remotedir,sshDirs{ff}));
                       sshFiles = ssh2_simple_command(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass',evalstr,1);
                       
                       %Make the current internal folder
                       mkdir(fullfile(rawdir,sshDirs{ff}));
                       
                       %This has a tendency to stall - going to do a simple
                       %version that iterates through each file instead of
                       %feeding SCP a structure...
                       %Old version for posterity uses a structure:%                            %Get the files for that folder
%                            scp_simple_get(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass,sshFiles,fullfile(rawdir,sshDirs{ff}),fullfile(remotedir,sshDirs{ff}));

                       keyboard
                       for fff=1:length(sshFiles)
                           scp_simple_get(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass,sshFiles{fff},fullfile(rawdir,sshDirs{ff}),fullfile(remotedir,sshDirs{ff}));
                       end
                   end
               end
               
               %Okay, we should have all of the files we needed from the
               %RCBI, now try for stimfiles.
               
               %List files inside the current folder
               evalstr = sprintf('ls %s',fullfile(remotestimdir));
               sshFiles = ssh2_simple_command(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass',evalstr,1);
               
               %Make the current internal folder
               mkdir(fullfile(rawdir,groupParams.Scan_Date(3:8)));
               
               %Get the files for that folder
               scp_simple_get(fetchParams.RCBIserver,fetchParams.RCBIuser,fetchParams.RCBIpass,sshFiles,fullfile(rawdir,groupParams.Scan_Date(3:8)),fullfile(remotestimdir));
               
           else
               return
           end
           
           
        end
    
end

if groupParams.doAutorecon
   %Alright, we're hoping to build some code here that can pass data from
   %one computer to the supercomputer, run autorecon (from Freesurfer) on 3
   %anatomical datafiles, then potentially to ask for those files back at a
   %later time. Or we may sync them via Dropbox, let's see.
   
   %For now, let's just get a local version up and running in the
   %background. This can onl run after we've downloaded anatomical files.
   
   %Okay, first off, we need the raw directory and to find out where / how
   %many anatomical sets there are 
   
   rawdir=strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date);
   rawdirstruct=dir(rawdir);
   numDirs=length(rawdirstruct);
    keyboard
   whereAnat=[];
    
        for i=1:numDirs
            thisstring=strfind(rawdirstruct(i).name,'MPRAGE');
            if thisstring
                whereAnat=[whereAnat, i];
            end
        end
    %Make us a string for autorecon that has however many anatomical scans
    %we've got
    
    anatstr=[];
    for i=1:length(whereAnat)
        subjdir=strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date,'/',rawdirstruct(whereAnat(i)).name);
        subjdirstruct=dir(subjdir);
        
        
        %Need to get rid of directory artifacts
        deleted=0;
        for ii=1:length(subjdirstruct);
            if deleted
                ii=ii-1;
                deleted=0;
            end
            if ii>length(subjdirstruct)
                break
            end
            if strcmp('.',subjdirstruct(ii).name) || strcmp('..',subjdirstruct(ii).name)
                subjdirstruct(ii)=[];
                deleted=1;
                continue
            end
        end
        
        %Okay, sometimes parentheses are messing us up...this loop is to
        %fix that. 
        %Find any instances of parentheses...as far as I know these occur
        %in the MPRAGE scans as (192)
        try
            if ~isempty(strfind(rawdirstruct(whereAnat(1)).name,'('))
                for ii=1:length(whereAnat)
                    pInds(ii,1)=strfind(rawdirstruct(whereAnat(ii)).name,'(');
                    pInds(ii,2)=strfind(rawdirstruct(whereAnat(ii)).name,')');
                end
            end
        catch
        end
        %If there are any, we need to rename these directories and their
        %internal files, then redo everything we just did for valid
        %directories. 
        if exist('pInds','var')
            if any(pInds)
                for ii=1:length(whereAnat)
                    if any(pInds(ii,:))
                        newbasestr=rawdirstruct(whereAnat(ii)).name;
                        newbasestr(pInds(ii,pInds(ii,:)>0))=[];
                        %Rename the folder
                        movefile(fullfile(rawdir,rawdirstruct(whereAnat(ii)).name),fullfile(rawdir,newbasestr));
                        %Do the same thing for all files inside the folder
                        allDCMs=dir(fullfile(rawdir,newbasestr));
                        for iii=1:length(allDCMs)
                            clear newFilestr
                            pDInds{iii,1}=strfind(allDCMs(iii).name,'(');
                            pDInds{iii,2}=strfind(allDCMs(iii).name,')');
                            if ~isempty(pDInds{iii,1}) && ~isempty(pDInds{iii,2})
                                delInds=[pDInds{iii,1} pDInds{iii,2}];
                                newFilestr=allDCMs(iii).name;
                                newFilestr(delInds)=[];
                            elseif ~isempty(pDInds{iii,1})
                                newFilestr=allDCMs(iii).name;
                                newFilestr(pDInds{iii,1})=[];
                            elseif ~isempty(pDInds{iii,2})
                                newFilestr=allDCMs(iii).name;
                                newFilestr(pDInds{iii,2})=[];
                            end
                            if strcmp(allDCMs(iii).name,'.') || strcmp(allDCMs(iii).name,'..')
                                continue
                            else
                                movefile(fullfile(rawdir,newbasestr,allDCMs(iii).name),fullfile(rawdir,newbasestr,newFilestr));
                            end
                        end
                        
                        
                    end
                    
                    
                end
                
            end
        end
        
                
        anatstrs{i}=sprintf('-i %s',fullfile(rawdir,'/',rawdirstruct(whereAnat(i)).name,'/',subjdirstruct(1).name));
    end
    anatstr=[];
    for i=1:length(anatstrs)
        anatstr=sprintf('%s %s',anatstr,anatstrs{i});
    end
    
    nCores = feature('numCores');
    
    sidout=strcat(groupParams.Subject_Initials,groupParams.Scan_Date);
    if nCores>1
        paraStr = '-openmp';
        evalstr = sprintf('recon-all -sid %s %s -autorecon1 -autorecon2 %s %g', ...
            sidout, anatstr, paraStr, nCores);

    else
        evalstr = sprintf('recon-all -sid %s %s -autorecon1 -autorecon2', ...
            sidout, anatstr);

    end
    
   disp(sprintf('%s', evalstr));
   %Run in background so we can get on with the rest of this..., do some
   %conversions, etc.
   
   
   [status, result] = system(evalstr);
   
   
   
   
end




%MRI Conversion using mri_convert

if groupParams.MRBuildDirs
    
    %Build directory to use
    fulldir=strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date);
    rawdirstruct=dir(fulldir);
    
    %Does the main working directory exist yet?
    if ~exist(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date),'dir')
        mkdir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date));
    end
    %How about the Anatomy folder
    if ~exist(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/','Anatomy'),'dir')
        mkdir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date),'/Anatomy');
    end
    %Raw?
    if ~exist(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/','Raw'),'dir')
        mkdir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Raw'));
    end
    %Tseries
    if ~exist(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Raw/TSeries'),'dir')
        mkdir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Raw/TSeries'));
    end
    %Etc?
    if ~exist(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/','Etc'),'dir')
        mkdir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Etc'));
    end
    
end

if groupParams.MRIConvert
    
    %Build directory to use
    fulldir=strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date);
    rawdirstruct=dir(fulldir);
    
    numDirs=length(rawdirstruct);
    
    if strcmp(groupParams.Output,Output{1});
        outputtype='nifti1';
    elseif strcmp(groupParams.Output,Output{2});
        outputtype='nii';
    end
    whereGRE=[];
    
    if groupParams.DoDewarp
        
        
        
        
        for i=1:numDirs
            thisstring=strfind(rawdirstruct(i).name,'gre');
            if thisstring
                whereGRE=[whereGRE, i];
            end
        end
        
        deletedGRE=0;
        for i=1:length(whereGRE)
            thisstring=strfind(rawdirstruct(whereGRE(i-deletedGRE)).name,'(ac-pc)');
            if thisstring
                whereGRE(i-deletedGRE)=[];
                deletedGRE = deletedGRE+1;
            end
        end
    
    
        
        
        
        %Okay, sometimes parentheses are messing us up...this loop is to
        %fix that. 
        %Find any instances of parentheses...as far as I know these occur
        %in the MPRAGE scans as (192)
       
        %If there are any, we need to rename these directories and their
        %internal files, then redo everything we just did for valid
        %directories. 
        if exist('pInds','var')
            clear pInds
            
        end
        
        try
            for pp=1:length(whereGRE)
                if ~isempty(strfind(rawdirstruct(whereGRE(pp)).name,'('))
                    for ii=1:length(whereAnat)
                        pInds(ii,1)=strfind(rawdirstruct(whereAnat(ii)).name,'(');
                        pInds(ii,2)=strfind(rawdirstruct(whereAnat(ii)).name,')');
                    end
                end
            end
        catch
        end
        
        
        
        if exist('pInds','var')
            if any(pInds)
                for ii=1:length(whereGRE)
                    if any(pInds(ii,:))
                        newbasestr=rawdirstruct(whereGRE(ii)).name;
                        newbasestr(pInds(ii,pInds(ii,:)>0))=[];
                        %Rename the folder
                        movefile(fullfile(rawdir,rawdirstruct(whereGRE(ii)).name),fullfile(rawdir,newbasestr));
                        %Do the same thing for all files inside the folder
                        allDCMs=dir(fullfile(rawdir,newbasestr));
                        for iii=1:length(allDCMs)
                            clear newFilestr
                            pDInds{iii,1}=strfind(allDCMs(iii).name,'(');
                            pDInds{iii,2}=strfind(allDCMs(iii).name,')');
                            if ~isempty(pDInds{iii,1}) && ~isempty(pDInds{iii,2})
                                delInds=[pDInds{iii,1} pDInds{iii,2}];
                                newFilestr=allDCMs(iii).name;
                                newFilestr(delInds)=[];
                            elseif ~isempty(pDInds{iii,1})
                                newFilestr=allDCMs(iii).name;
                                newFilestr(pDInds{iii,1})=[];
                            elseif ~isempty(pDInds{iii,2})
                                newFilestr=allDCMs(iii).name;
                                newFilestr(pDInds{iii,2})=[];
                            end
                            if strcmp(allDCMs(iii).name,'.') || strcmp(allDCMs(iii).name,'..')
                                continue
                            else
                                movefile(fullfile(rawdir,newbasestr,allDCMs(iii).name),fullfile(rawdir,newbasestr,newFilestr));
                            end
                        end
                        
                        
                    end
                    
                    
                end
                
            end
        end
        
        
        %Do MRI Conversion for GRE mapping
        
        
        
        
        % 1) convert magnitude volume from field map scans:  mri_convert ../20130619/7.gre_field_mapping/*0001.dcm -o Etc/run07.nii
        iFile = strcat(fulldir,'/',rawdirstruct(whereGRE(1)).name,'/','*0001.dcm');
        outEtc=strcat(groupParams.Working_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Etc','/run',num2str(whereGRE(1)),'.nii');
        evalstr = sprintf('mri_convert %s -o %s', ...
            iFile, outEtc);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr);
        
        % 2) convert phase difference volume:                mri_convert ../20130619/8.gre_field_mapping/*0001.dcm -o Etc/run08.nii
        
        iFile = strcat(fulldir,'/',rawdirstruct(whereGRE(2)).name,'/','*0001.dcm');
        outEtc=strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Etc','/run',num2str(whereGRE(2)),'.nii');
        evalstr = sprintf('mri_convert %s -o %s', ...
            iFile, outEtc);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr);
        
    end
    
    %Do mri_convert for each 3D structural
    clear thisstring
    where3D=[];
    for i=1:numDirs
        thisstring=strfind(rawdirstruct(i).name,'MPRAGE');
        thisstring1=strfind(rawdirstruct(i).name,'t1');
        if ~isempty(thisstring) || ~isempty(thisstring1)
            where3D=[where3D, i];
        end
    end
    
    for i=1:length(where3D)
        iFile = strcat(fulldir,'/',rawdirstruct(where3D(i)).name,'/','*0001.dcm');
        out=strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Anatomy','/run',num2str(where3D(i)));
        outtype=outputtype;
        evalstr = sprintf('mri_convert %s %s -ot %s', ...
            iFile, out,outtype);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr);
        
    end
    
    %Same for 4D
    clear thisstring thisstring1
    where4D=[];
    for i=1:numDirs
        thisstring=strfind(rawdirstruct(i).name,'ep2d');
        if ~isempty(thisstring)
            where4D=[where4D, i];
        end
    end
    
    for i=1:length(where4D)
        if exist('dotinds','var');
            clear dotinds
        end
        iFile = strcat(fulldir,'/',rawdirstruct(where4D(i)).name,'/','*0001.dcm');
        allfiles=dir(strcat(fulldir,'/',rawdirstruct(where4D(i)).name,'/*.dcm'));
        fname=allfiles(1).name(1:(length(allfiles(1).name)-9));
        if ~isempty(find(fname=='.',1));
            dotinds=find(fname=='.');
            fname(dotinds)='_';
        end
        %         out=strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Raw/Tseries','/run',num2str(where4D(i)));
        out=strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Raw/Tseries','/',fname);
        
        outtype=outputtype;
        evalstr = sprintf('mri_convert %s %s -ot %s', ...
            iFile, out,outtype);
        disp(sprintf('%s', evalstr));
        [status, result] = system(evalstr);
        
    end
    
end
%Cool! Let's allow ourselves to pick the scans we're going to use now that
%we've converted everything, we're gonna get stim parameters too

%Need to index the files we've converted and added now, should be in
%working dir / Raw/Tseries

if strcmp(groupParams.Output,Output{1})
    extension='hdr';
    check_img_size=1;
elseif strcmp(groupParams.Output,Output{1})
    extension='nii'
    check_img_size=0;
end
scandirlist=dir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',(groupParams.Scan_Date),'/Raw/TSeries/*.',extension));

if check_img_size && exist('imgSize','var')
    imgList=dir(strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',(groupParams.Scan_Date),'/Raw/TSeries/*.','img'))
    for pp=1:length(imgList)
       imgSize(pp)=imgList(pp).bytes; 
    end
    incScans=imgSize>200000;
    scandirlist=scandirlist(incScans);
end

if groupParams.HandleStimFiles
    
    %Did we fetch stimfiles? Are they still in the Raw directory?
    
    if exist(strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date(3:8)),'dir')
        stimdir=strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date(3:8));
    else
        disp(sprintf('Oops, cannot find stimfiles in downloaded directory, can you tell me where to find the stimfiles?'))
        stimdir=input('Stimfile directory:');
    end
    
    %Okay, found stimfiles, copy to working directory
    evalstr = sprintf('cp -f %s/* %s', ...
        stimdir, strcat(groupParams.Working_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/','Etc'));
    disp(sprintf('%s', evalstr));
    [status, result] = system(evalstr);
    
end

if groupParams.CheckStimfileMatch
    
    %Let's go to the root of the working directory, where we need to be
    workingfile=strcat(groupParams.Working_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date);
    
    
    % get info about scans that live in Raw/TSeries
    scanDirName = fullfile(groupParams.Working_Directory,groupParams.Subject_Initials,groupParams.Scan_Date,'Raw','TSeries');
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
    % check to see if we got any good scans
    if nScans == 0
        disp(sprintf('(mrInit) Could not find any valid scans in Raw/TSeries'));
        sessionParams = [];groupParams = [];
        return
    end
    
    [stimFileNames stimFileVols] = getStimFiles(workingfile);
    stimFileMatch = matchStimFiles(stimFileNames,stimFileVols,totalFrames);
    
    %Okay, we've got information on both stimfiles and on scans, let's do some
    %quick comparisons
    
    %Do we have the same number of scans as stimfiles?
    scansequalstimfiles=nScans==length(stimFileVols);
    
    if scansequalstimfiles
        for i=1:length(stimFileVols)
            %Do the length of the stimfiles equal the length of the scans?
            equalvols(i)=stimFileVols(i)==cell2mat(nFrames(i));
        end
    else
        mrWarnDlg('Hey! Looks like you do not have the same number of stim files and scans, see if you can match the number - we will take the first set of matching stimfiles');
        mrWarnDlg('If you made changes to the directory type yes (in single quotes), if not type hit enter and this section will exit');
        changedstimfiles=input('Yes/no:');
        
        if strcmp(changedstimfiles,'yes') || strcmp(changedstimfiles,'Yes') || strcmp(changedstimfiles,'y')
            [stimFileNames stimFileVols] = getStimFiles(workingfile);
            stimFileMatch = matchStimFiles(stimFileNames,stimFileVols,totalFrames);
            
            for i=1:length(stimFileVols)
                %Do the length of the stimfiles equal the length of the scans?
                equalvols(i)=stimFileVols(i)==cell2mat(nFrames(i));
            end
        elseif isempty(changedstimfiles)
            mrWarnDlg('Looks like we have a bad match between stimfiles, quitting for now - do custom allocation via mrInit')
            return
            
        end
    end
    
    if any(equalvols==0)
        if strcmp(groupParams.Condition,'pRF')
            
            mrWarnDlg('Looks like we may have run into the 167/168 frame bug in bar code - going to add a flag to make sure we start with a scan of 167, should work as a workaround for now')
            pRFBarBug=1;
            cd(workingdir)
            svfile=strcat('Etc/pRFBarBug');
            save(svfile,'pRFBarBug','stimFileVols','nFrames');
        else
            mrWarnDlg('Hmm, some stimFiles do not have the same number of frames as your scans, going to suggest that these scans were not run properly')
            
            
            
            %     %Let's try and make everyone use 160 total frames
            %     for i=1:length(stimFileNames)
            %       nFrames{i} = 160;
            %       junkFrames{i} = 168-stimFileVols(i)+7;
            %     end
        end
        
        %
    end
    
    %One more sanity check, going to depend on what condition scan we ran
    
    %Retinotopy - did we scan in the usual order or not?
    
    if strcmp(groupParams.Condition,condition{1})
    end
    
    %pRF Which scans are 2 and which are 4 barWidth?
    if strcmp(groupParams.Condition,condition{1})
    end
    
    %Adaptation, what the hell did we do and where?
    if strcmp(groupParams.Condition,condition{1})
    end
    
end



if groupParams.RunmrInit
    if exist('workingdir','var')
        cd(workingdir)
    else
        workingdir=fullfile(groupParams.Working_Directory,groupParams.Subject_Initials,groupParams.Scan_Date);
        cd(workingdir)
    end
    
%      [sessionParams groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1')

    
    sessionParams.description = groupParams.Condition;
    sessionParams.subject = groupParams.Subject_Initials;
    sessionParams.operator = 'PW_probably';
    sessionParams.coil = 'Siemens Birdcage';
    sessionParams.magnet = 'Siemens Magnetom 3T';
    sessionParams.pulseSequence = 'RCBI_EPI';
    sessionParams.pulseSequenceText='';
    
    
    [sesParams gParams] = mrInit(sessionParams,[],'justGetParams=0')
    
    if isempty(gParams)
       mrWarnDlg('Hmm, we did not set any parameters, sometimes this happens if you are running matlab from tcsh, this like Bash better, try restarting matlab in bash') 
    end
    
end

if groupParams.DoDewarp
    
    %Move into the subject working directory
    if exist('workingdir','var')
        cd(workingdir)
    else
        workingdir=fullfile(groupParams.Working_Directory,groupParams.Subject_Initials,groupParams.Scan_Date);
        cd(workingdir)
    end
    
    %Okay, enumerate the Etc directory where we should have stored the
    %converted gradient echo images
    
    EtcNii=dir('Etc/*.nii');
    
    %What if we didn't get the files?
    if isempty(EtcNii)
        %Maybe let's try the conversion again...
        
        %Build directory to use
        fulldir=strcat(groupParams.Raw_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/',groupParams.Scan_Date);
        rawdirstruct=dir(fulldir);
        
        numDirs=length(rawdirstruct);
        
        if strcmp(groupParams.Output,Output{1});
            outputtype='nifti1';
        elseif strcmp(groupParams.Output,Output{2});
            outputtype='nii';
        end
        whereGRE=[];
        
        if groupParams.DoDewarp
            for i=1:numDirs
                thisstring=strfind(rawdirstruct(i).name,'gre');
                if thisstring
                    whereGRE=[whereGRE, i];
                end
            end
            
            %Do MRI Conversion for GRE mapping
            
            % 1) convert magnitude volume from field map scans:  mri_convert ../20130619/7.gre_field_mapping/*0001.dcm -o Etc/run07.nii
            
            iFile = strcat(fulldir,'/',rawdirstruct(whereGRE(1)).name,'/','*0001.dcm');
            outEtc=strcat(groupParams.Working_Directory,groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Etc','/run',num2str(whereGRE(1)),'.nii');
            evalstr = sprintf('mri_convert %s -o %s', ...
                iFile, outEtc);
            disp(sprintf('%s', evalstr));
            [status, result] = system(evalstr);
            
            % 2) convert phase difference volume:                mri_convert ../20130619/8.gre_field_mapping/*0001.dcm -o Etc/run08.nii
            
            iFile = strcat(fulldir,'/',rawdirstruct(whereGRE(2)).name,'/','*0001.dcm');
            outEtc=strcat(groupParams.Working_Directory,'/',groupParams.Subject_Initials,'/',groupParams.Scan_Date,'/Etc','/run',num2str(whereGRE(2)),'.nii');
            evalstr = sprintf('mri_convert %s -o %s', ...
                iFile, outEtc);
            disp(sprintf('%s', evalstr));
            [status, result] = system(evalstr);
            
        end
    end
    
    EtcNii=dir('Etc/*.nii');
    
    if isempty('EtcNii')
        mrWarnDlg('Yikes, looks like there is no gradient echo scans to work from, and we cannot convert it from Raw, double check that these scans exist, aborting for now');
        return
    end
    
    %Okay, if we're going to use the doB0Correction script we must have the
    %environmental variables from Freesurfer. Let's check for that
   
       v = newView; doB0Correction(v, fullfile('Etc',EtcNii(1).name), fullfile('Etc',EtcNii(2).name));
end

if groupParams.DoMotionComp
    v=newView;
    [v params] = motionComp(v,[], 'justGetParams=1','defaultParams=1');
    params.interpMethod='cubic';
    params.sliceTimeCorrection=0;
    v=motionComp(v,params);
    deleteView(v);
end

%If we're doing ring and wedge retinotopy, let's get things set:
%Note that we're now going to build averages for two cases in R+Ws. The
%first set of 3 averages will be Wedges, Rings, and MT Localizer. The 2nd
%set of 4 are going to be Wedges CW, Wedges CCW, Rings Contracting, Rings
%Expanding. The first set requires us to junk 8 frames. The 2nd set
%requires that we do not junk frames. The 1st set is reversed to combine
%averages, the 2nd set is not. This is in order to be able to run a PRF
%analysis on just Rings and Wedges. 



if strcmp(groupParams.Condition,condition{1}) && groupParams.DoAverage
    
    if ~isstruct('v')
        v = newView;
    end
    %Since we aren't assuming the nature of the scans (things can get
    %messed up during scanning) - we need to pull each stim type and its
    %direction from the stim files.
    gnames=viewGet(v,'groupnames');
    %Are we working with a motioncomped group?
    groupind=strfind(gnames,'MotionComp');
    if any(cell2mat(groupind))
       v=viewSet(v,'groupname','MotionComp'); 
    end
    nScans=viewGet(v,'nscans');
    
    
    for iScan=1:nScans
        %Get the stimfile for this scan
        retscan(iScan).stimFileName=viewGet(v,'stimfilename',iScan);
        if ~isempty(retscan(iScan).stimFileName)
            if ~isempty(retscan(iScan).stimFileName{1})
                retscan(iScan).stimtype = getParamfromStimfile(iScan,retscan(iScan).stimFileName,'stimulus.stimulusType');
                %load the stimfile to check for mtlocalizer just in case
                 currentfilename=getParamfromStimfile(iScan,retscan(iScan).stimFileName,'task{2}{1}.taskFilename');
                    if strcmp('mtlocalizer.m',currentfilename)
                        retscan(iScan).stimtype=4; %This will signal the MT Localizer as its own type of scan for averaging...
                    end
                if retscan(iScan).stimtype==3
                    mrWarnDlg('Oops - looks like these stimfiles are from a Bars only scan, try re-running mrSanity with pRF selected for type of analysis');
                    return
                elseif isempty(retscan(iScan).stimtype)
                    %Possible we have the MT localizer - let's look
                    currentfilename=getParamfromStimfile(iScan,retscan(iScan).stimFileName,'task{2}{1}.taskFilename');
                    if strcmp('mtlocalizer.m',currentfilename)
                        retscan(iScan).stimtype=4; %This will signal the MT Localizer as its own type of scan for averaging...
                    else
                    mrWarnDlg('No good - looks like this is not a valid stimfile or there was no retinotopy scan here, cannot find stimtype - quitting');
                    return
                    end
                end
                retscan(iScan).stimdir = getParamfromStimfile(iScan,retscan(iScan).stimFileName,'stimulus.direction');
            else
                retscan(iScan).noStim=1;
            end
        end
    end
    
    
      for ss=1:length(retscan)
          if ~isempty(retscan(ss).stimtype)
              allstimtypes(ss)= retscan(ss).stimtype;
          else
              allstimtypes(ss)=NaN;
          end
          if ~isempty(retscan(ss).stimdir)
              alldirections(ss)=retscan(ss).stimdir;
          else
              alldirections(ss)=NaN;              
          end
      end
      
      %Get number of unique stimtypes and number of unique directions
      allstimtypesX=unique(allstimtypes);
      alldirectionsX=unique(alldirections);
      allstimtypesXF=allstimtypesX(allstimtypesX>0);
      alldirectionsXF=alldirectionsX(alldirectionsX>-2);

      %We should have stimtypes 1 and 2 basically, and we're going to need
      %to loop through for directions -1 and 1, then stimtype 4 will be
      %mtlocalizer
     
      
      for ss=1:length(allstimtypesXF)
          for sss=1:length(alldirectionsXF)
              curstimlist=find(allstimtypes==allstimtypesXF(ss));
              curdirlist=find(alldirections==alldirectionsXF(sss));
              stimdirintersect=intersect(curstimlist,curdirlist);
              samestim{ss}{sss}=stimdirintersect;
          end
          
      end
      
      %Now we have stim lists for all rings and wedges runs, probably worth
      %noting the mtlocalizer ones as well.
            samestim{ss}{1}=find(allstimtypes==4);
            
            %Each cell has scans that we can average together (same stim
            %same direction)
       
            if ~isstruct('v')
                v=newView;
            end
            
            
    %Need to set the plan for reversing - the way we built things anything
    %in cell set 2 can be reversed
    %Make a set of all one direction, then a set of all the rest, we'll
    %reverse the first half
    
    %Number of different stimuli 
    for iAve=1:2
            scanlist=[];

    %Number of different directions possible 
    for ss=1:2
        scanlist=horzcat(scanlist,(samestim{iAve}{ss}));

    end
    
    %Now we'll plan to reverse half, shift all of it, and average it,
    %blammo.
    [v params] = averageTSeries(v,[],'justGetParams=1','defaultParams=1','scanList',scanlist);

    
    params.shiftList(:) = 2;
    params.reverseList(1:length(samestim{iAve}{1})) = 1;
    
    averageTSeries(v,params);

    
    end
    
    %Okay! Do the mtlocalizer

                scanlist=[];
scanlist=samestim{3}{1}
    [v params] = averageTSeries(v,[],'justGetParams=1','defaultParams=1','scanList',scanlist);

    
    params.shiftList(:) = 2;
    params.reverseList(1:length(samestim{3}{1})) = 0;
    
    averageTSeries(v,params);
      
end


%If we're doing pRF - let's do our averages
if strcmp(groupParams.Condition,condition{2}) && groupParams.DoAverage
    
     %Move into the subject working directory
    if exist('workingdir','var')
        cd(workingdir)
    else
        workingdir=fullfile(groupParams.Working_Directory,groupParams.Subject_Initials,groupParams.Scan_Date);
        cd(workingdir)
    end
    
    v=newView;
    v=viewSet(v,'curgroup', 'MotionComp');
    %for i=1:10; disp(viewGet(v,'description',i)); end
    
    %The way we typically run a pRF scan, we do 2 different conditions of
    %barWidth. Let's look for these in the stimfiles and thus determine which
    %scans to average together
    
%     v = viewSet(v, 'curGroup', 'MotionComp');
    nScans = viewGet(v,'nScans','MotionComp');
    
    for iScan = 1:nScans
        %Set ourselves to the MotionComp group and get its groupnumber
        v = viewSet(v, 'curGroup', 'MotionComp');
        groupNum = viewGet(v,'currentGroup');
        %Get stimfile name for this particular scan in the motioncomp group
        stimfilename = viewGet(v,'stimFileName',iScan,groupNum);
        
        %Only try to grab the stim parameter if we have a stimfile there
        if ~isempty(stimfilename)
            barWidth(iScan) = getParamfromStimfile(iScan,stimfilename,'stimulus.barWidth');
        else
            barWidth(iScan)=NaN;
        end
        
    end
    
    %How many types of widths are we dealing with?
    widthConditions=unique(barWidth);
    conds=widthConditions(widthConditions>-1);
    
    for iConds=1:length(conds)
        
        if ~isstruct(v)
            v=newView;
            v=viewSet(v,'curgroup', 'MotionComp');
        end
        
        scanList=find(barWidth==conds(iConds));
        
        %If we add a change here, and mess with the scanlist, we might be able
        %to fix the 167/168 problem. Let's see.
        
        if exist('Etc/pRFBarBug.mat','file')
            load('Etc/pRFBarBug.mat');
            
            if pRFBarBug
                %Okay, if we've gotten here, there was at least 1 file that
                %had the bug. If we do the workaround in one set of
                %averages and not another, we'll throw an error, so we need
                %the workaround to exist in every average or not at all.
                
                %Check for workaround in each list
                for listcheck=1:length(conds)
                    checklist(listcheck,:)=find(barWidth==conds(listcheck));
                end
                %Need to make sure we have the same number of stimfile
                %categories as scans, otherwise we'll throw an error
                if length(stimFileVols)<max(max(checklist))
                    checklist_orig=checklist;
                    checklist=checklist-(max(max(checklist))-length(stimFileVols));                   
                    
                end
                
                for listcheck=1:length(conds)
                    tmp=find(stimFileVols(checklist(listcheck,:))==167);
                    if ~isempty(tmp)
                        wherebugged(:,listcheck)=tmp;
                    else
                        continue
                    end
                end
                
                %Okay, if we've gotten this far we should be good,
                %unless one set of scans has no scans that need the
                %workaround, in that case, for now, we're just going to
                %pick all scans BUT the ones that need the workaround. 
                if any(wherebugged==0)
                    wherebugged(find(wherebugged>0))=checklist(wherebugged(find(wherebugged>0)));

                else
                    wherebugged=checklist(wherebugged);
                end
                %This is getting messy, let's make some structures to help
                %us clean things up
                
                for ss=1:length(conds)
                    scans(ss).group=ss;
                    scans(ss).scanlist=checklist(ss,:);
                    for sss=1:length(wherebugged)
                        if any(ismember(wherebugged(sss),scans(ss).scanlist))
                            scans(ss).buggedscans=wherebugged(sss);
                        end
                    end
                end
                
                %Okay, so we have the number of groups of scans, list of
                %scans in each set, and which scans are bugged. Now, is
                %there a scan in both that is bugged? If not we can't use
                %this.
                
                for ss=1:length(scans)
                    whichscansbugged(ss)=any(scans(ss).buggedscans);
                end
                if length(scans)==length(whichscansbugged) && ~any(whichscansbugged==0)
                    %Woo, okay, if we've done all this we can just
                    %rearrang the scans to start with a scan of 167
                    %frames, do it!
                    
                    for ss=1:length(scans)
                        scans(ss).scanlist_cat=horzcat(scans(ss).buggedscans(1),scans(ss).scanlist);
                        scans(ss).bugindex=find(scans(ss).scanlist==scans(ss).buggedscans(1));
                        scans(ss).scanlist_bugorder=scans(ss).scanlist_cat;
                        scans(ss).scanlist_bugorder(scans(ss).bugindex+1)=[];
                        
                    end
                    %If we dealt with a couple scans without stimfiles, we
                    %need to note this, our list will be shifted over
                    %relative to the actual scans in the set
                    if exist('checklist_orig','var')
                        %Add this file to the structure, then we can look
                        %for it, easy. 
                        for ss=1:length(scans)
                            scans(ss).scanlist_orig=checklist_orig(ss,:)
                            scans(ss).scanlist_orig_reordered=scans(ss).scanlist_bugorder+(max(max(checklist_orig))-length(stimFileVols));
                        end
                        
                    end
                    
                    %Okay, apparently we also need to change numframes to
                    %make this happen - use 160 numframes, junk 7 on all
                    %scans
                    
                    %Check number of frames on only the *actual* scan
                    %numbers we're using
                    for ss=1:length(scans(iConds).scanlist_orig_reordered)
                        scans(iConds).scanParams = viewGet(v,'scanparams',scans(iConds).scanlist_orig_reordered(ss)); 
                        scans(iConds).scanParams.nFrames=160;
                        scans(iConds).scanParams.junkFrames=7;

                        v=viewSet(v,'updatescan',scans(iConds).scanParams,scans(iConds).scanlist_orig_reordered(ss))
                    end
                else
                    mrWarnDlg('Yuh oh, my workaround for the frame bug is not going to work for this data set, we need at least one scan in each average to do it. Going to build an average that does not include the bugged datasets in hopes of maximizing available data ');
                    
%                     %Does the current condition we're in have a bugged
%                     %scan?
%                     if ~isempty(scans(iConds).buggedscans)
%                         %No bugged scans, build a normal scanlist of all
%                         %scans in the list
%                         
%                     else
%                         %Bugged scans listed, remove them from the list and
%                         %do the average that way. 
%                     end
%This should cause us to use the code below, which handles removal of any
%scan under 168 frames. We'll lose a scan or two but be able to run. 
                    clear scans
                    
                    
                end
            end
        end
        
        v=newView;
        v=viewSet(v,'curgroup', 'MotionComp');
        
        %If we ran into the pRF bug we'll have a scanlist to use, if not we
        %need to build a simple one
        
        if exist('scans','var')
            scanstouse=strcat('scanList=',mat2str(scans(iConds).scanlist_orig_reordered))
        else
            clear scanInds
            
            %Okay, there are X number of conditions, we're in the loop for
            %those, so pick the barWidth we're dealing with, pull it out
            %for each scan from each stimfile, then build the list from
            %there. 
            
            currWidth=widthConditions(iConds);
            for ii=1:nScans
                currSFN=viewGet(v,'stimfilename',ii);
                
                bwTmp{ii}=getParamfromStimfile(ii,currSFN,'stimulus.barWidth');
                if ~isempty(bwTmp{ii})
                    allBW(ii)=bwTmp{ii};
                end
                
            end
            scanInds=find(allBW==currWidth);
            
            %Verify that all scans have the right number of scans, ditch
            %scans that dont...
            
            for ii=1:length(scanInds)
                framesList(ii)=viewGet(v,'nframes',scanInds(ii));
            end
            delScans=find(framesList<168);
            scanInds(delScans)=[];
            scanstouse=strcat('scanList=',mat2str(scanInds));

            
        end;
        [v params] = averageTSeries(v,[],'justGetParams=1','defaultParams=1',scanstouse);
        v=averageTSeries(v,params);
        deleteView(v);
        
        
    end
    
end

if groupParams.doConcat && strcmp(groupParams.Condition,condition{2})
    %Verify that we have the average group'
    if ~isstruct('v')
        v = newView;
    end
    groupNames = viewGet(v,'groupNames');
    if any(strcmp('Averages',groupNames));
        v = viewSet(v, 'curGroup', 'Averages');
        groupNum = viewGet(v,'currentGroup');
        nScans = viewGet(v,'nScans');
        [v params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1','scanList',[1:nScans]);
        v=concatTSeries(v,params);
        
    end
end

if groupParams.doAnalysis
    if strcmp(groupParams.Condition,condition{2})
        %Okay, we're going to try to run a full pRF analysis - this won't
        %have an ROI restriction, so it might take a bit, but it should run
        %nicely. For now, we'll default to nelder-mead, and fitting the hdr
        %function
        if ~isstruct('v')
            v = newView;
        end
        
        %This is only true if we've built the ROI for the brain - this
        %own't always be true, commenting for now
%         v=loadROI(v,'Brain',1);
        
        gInfo = groupInfo;
        %Let's generally assume we have one concat scan that we're going to analyze
        v = viewSet(v,'curGroup','Concatenation');
        [v params] = pRF(v,[],'justGetParams=1','defaultParams=1','scanList=1')
        params.pRFFit.rfType='gaussian-hdr';
        groupNum = viewGet(v,'curGroup');
        params.groupNum=groupNum;
        [v d] = pRF(v,params,'justGetParams=0','scanList=1','defaultParams=1')%stGetParams=1','defaultParams=0','scanList=0');
    end
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




