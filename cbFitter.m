function cbFitter(dataFile, saveData, GOF, sessionSubsetStart, sessionSubsetEnd, subsetFromEnd, doPlot)

% New data fitter for Huxlin Lab
% MMelnick, December 2015
%
%cbFitter(dataFile, saveData, GOF, sessionSubsetStart, sessionSubsetEnd, subsetFromEnd)
%
%This function was written to analyze training data for CB subjects
%training at home on dot discrimination (fine or coarse). It uses psignifit
%to perform fitting and further analysis, so you need to install and
%configure psignifit via http://psignifit.sourceforge.net
%
%Required input for this function is the datafile you want to analyze. It
%analyzes .txt files that are delimited by semicolon characters (;).
%Easiest thing is to place datafiles in the same folder that you're running
%from, although it will also accept an entire path.
%
%saveData is by default set to 1, this will output the final table to a
%simple CSV file that can then be passed to a spreadsheet. Columns in the
%final table include Training Date, Fit Threshold, Percent Correct, Fit
%Slope, Stimulus X Location, Stimulus Y Location, and number of trials
%analyzed.
%
%GOF performs Goodness of Fit analysis and plot via psignifit. This will
%open a figure for every fittable (above 0.75 correct) threshold, so use
%with caution! A good sanity check to make sure your data is fitting
%properly.
%
%sessionSubsetStart & sessionSubsetEnd default to analyzing the entire log
%file input to the program. If you only want to analyze a subset of the
%logfile, you can these as the start and finish. I.e. if you have 100
%sessions in the log file but only want to analyze the last 50 sessions,
%sessionSubsetStart = 50, and sessionSubsetEnd = 100


% if ieNotDefined('dataFile') dataFile = '';end

if ieNotDefined('whichFitter') whichFitter = 2;end
if ieNotDefined('GOF') GOF = 0;end
if ieNotDefined('saveData') saveData = 1;end
if ieNotDefined('subsetFromEnd') subsetFromEnd = [];end
if ieNotDefined('doPlot') doPlot = 0;end



% dataIn = dlmread('ARFlogA.txt',';');
% dataIn = dlmread('JKMlogA-4.txt', ';');
% dataIn = dlmread('MAHlogA.txt', ';');



dataIn = dlmread(dataFile,';');

%Turn off warnings for now
warning off;

dates.n = sum(dataIn(:,1)>300); %How many training dates are compiled here?
dates.which = dataIn(dataIn(:,1)>300,1);
%Convert mat to cell

dates.index = find(dataIn(:,1)>300);

stim.locx = dataIn(dates.index,3);
stim.locy = dataIn(dates.index,4);

%How many total sessions have been run (in this log file)
all.n = length(dates.index);



warning('off')



%Check for psignifit installation
if exist('psignifit_version')
    psignifit_version
else
    warning('Psignifit does not appear to be installed, we need that for fitting, see psignifit.sourceforge.net')
    warning('Switching to internal Weibull fit')
    whcihFitter = 2;
end

%We need to loop through however many sessions we need to analyze

%columns: number, threshold, staircase (1-3), direction left/right, Reaction time, no 6, column 7: right
%or wrong

wBar = waitbar(0, 'Running fit on all available data');


if ieNotDefined('sessionSubsetStart') sessionSubsetStart = 1;end
if ieNotDefined('sessionSubsetEnd') sessionSubsetEnd = dates.n;end

if ~isempty(subsetFromEnd)
    sessionSubsetStart = dates.n-subsetFromEnd;
end

allData={};
global unfitN
unfitN = 0;

for i = sessionSubsetStart:sessionSubsetEnd-1
    wBar = waitbar(i/(sessionSubsetEnd-sessionSubsetStart),wBar, strcat('Running fit on all available data  ',num2str(i),' of ', num2str(sessionSubsetEnd-1)));
    %Current session label by date
    currSess.date = num2str(dates.which(i,:));
    
    %Which trials in the main file are we actually pointing to here
    currSess.trialIndex = dates.index(i)+1:dates.index(i+1)-1;
    
    %Current session number of trials should vary as necessary
    currSess.nTrials = max(dataIn(currSess.trialIndex, 1));
    currSess.pC = dataIn(dates.index(i), 6);
    currSess.stairMean = dataIn(dates.index(i), 7);
    
    %Difficulty levels tested
    currSess.diff = dataIn(currSess.trialIndex,2);
    %Correct / incorrect
    currSess.cc = dataIn(currSess.trialIndex, 7);
    %Palamedes implementation for fitting follows
    
    StimLevels = unique(currSess.diff)';
    
    %Number of positive responses (e.g., 'yes' or 'correct' at each of the
    %   entries of 'StimLevels'
    clear NumPos
    for pp = 1:length(StimLevels)
        NumPos(pp) = sum(currSess.cc(currSess.diff==StimLevels(pp)));
    end
    
    %Number of trials at each entry of 'StimLevels'
    clear OutOfNum
    for pp = 1:length(StimLevels)
        OutOfNum(pp) = sum(currSess.diff==StimLevels(pp));
    end
    
    %%%Add in a skip here to deal with sets of data that have chance
    %%%performance / no thresholds really
    
    if currSess.pC <75 || isempty(currSess.nTrials) || currSess.nTrials<60;
        currSess.fitThresh = 'unfit';
        thresh = 'unfit';
        currSess.CI = 'unfit';
        currSess.fitSlope = 'unfit';
        slope = 'unfit';
        currSess.lapse = 'unfit';
        
        allData{end+1} = currSess;
        unfitN = unfitN+1;
        
    else
        
        %%%If you're using psignifit set priors, etc., run psignifit bootstrap
        
        if whichFitter==1
            %         xyn = [flip(StimLevels') NumPos' OutOfNum'];
            xyn = [StimLevels' NumPos' OutOfNum'];
            
            if xyn(1,1) == 0
                xyn(:,1) = xyn(:,1)+1;
            end
            
            xyn(:,1) = xyn(:,1)./361;
            xyn(:,1) = 1./xyn(:,1);
            
            %Log transform
            xyn(:,1) = log(xyn(:,1));
            priors.m_or_a = 'None';
            
            %If this is fine use a uniform prior
            if max(currSess.diff)<91
                priors.w_or_b = 'Uniform(0.0001,10';
                %Otherwise improper prior
            else
                %         priors.w_or_b = 'Uniform(.000001,0.1)';
                priors.w_or_b = 'None';
                
                %         priors.m_or_a = 'Gauss(150,100)';     
            end
            
            %     priors.lambda = 'Beta(1.5,12)';
            priors.lambda = 'Uniform(0,0.1)';
            
            %     priors.gamma  = 'Uniform(0,.1)';
            
            clear results
            %     results = BootstrapInference ( xyn, priors, 'sigmoid', 'logistic', 'core', 'ab','nafc', 2);
%             results = BootstrapInference ( xyn, priors, 'sigmoid', 'gumbel_l', 'core', 'weibull','nafc', 2);
%             
%             results = BayesInference ( xyn, priors,'sigmoid', 'logistic','core', 'ab','nafc', 2);
%             
            pilot = BootstrapInference ( xyn, priors, 'sigmoid', 'gauss', 'core', 'mw0.1','nafc', 2, 'cuts', [0.75 0.8, 0.82], 'samples', 500);
            results = BayesInference ( xyn, priors, 'sigmoid', 'gauss', 'core', 'mw0.1','nafc', 2, 'cuts', [0.75 0.8, 0.82], 'pilot', pilot.mcestimates);
            
            
                    
            thresh = getThres(results, 3);
            
            thresh = 10^(360*(1/thresh));
            slope  =  getSlope(results, 3);
            
            CI = getCI(results, 3,0.95);
            
            
            
            threshRange = min(currSess.diff):max(currSess.diff);
            %     threshFlipped = flip(threshRange);
            
            %     thresh = threshFlipped(find(threshRange==round(thresh)));
            
            %If we're using our own Weibull fitting program, do so here, may add
            %bootstrapping
            
            
        elseif whichFitter==2
            
            
            
            %         xyn = [StimLevels' NumPos' OutOfNum'];
            
            %Do SL transforms
            StimLevels2 = StimLevels+1; %Move away from 0
            
            
            if max(StimLevels2)<92
                fineDir = 1;
                SLINL = StimLevels2;
                
            else
                fineDir = 0;
                SLIN = 1./(StimLevels2./360); %Invert and normalize
                SLINL = log(SLIN); %log transform data
                
            end
            
            PC = NumPos./OutOfNum; %Percent correct
            
            [fitresult, gof] = createFit(SLINL, PC, OutOfNum, 1, i, sessionSubsetStart, sessionSubsetEnd, doPlot);
            
            fitOut = coeffvalues(fitresult);
            thresh75 = fitOut(2);
            slope = fitOut(3);
            lapse = fitOut(1);
            
            syms x
            assume(x, 'real')
            warning('ON');
            numericSol = vpasolve(lapse-exp(-exp(2*slope*thresh75/log(2)*(log(x)-log(thresh75))+log(log(2)))) == 0.75, x, 0.7);
            
%             if isempty(numericSol)
%                 keyboard
%                 numericSol = exp(0.240227/(slope*thresh75))*thresh75;
%                 
% %                 numericSol = real(2^((1.5708i)/(thresh75*slope))*exp(0.127024/(thresh75*slope))*thresh75*log(0.346574/(thresh75*slope))*(lapse-0.75));
%                 
%                 if fineDir
%                    thresh = numericSol;
%                    
%                 end
%             end
            
            if isempty(numericSol)
                currSess.fitThresh = 'unfit';
                thresh = 'unfit';
                currSess.CI = 'unfit';
                currSess.fitSlope = 'unfit';
                slope = 'unfit';
                currSess.lapse = 'unfit';
                allData{end+1} = currSess;
                continue
            end
            if ~fineDir
                thresh75 = 10^numericSol; %Remove log transform
                thresh = double((1/thresh75)*360);
            else
                thresh = sym2poly(numericSol);
            end
            allCI = confint(fitresult);
            CI = allCI(:,2);
            
            if ~fineDir
                %CI's need to be transformed
                CI = 10.^CI;
                CI = (1./CI)*360;
            end
            lapse = 1-lapse;
            
            %Need to shift to CI's to numerical solution for threshold
            if ~fineDir
                
                oldThresh =  numericSol;
                oldThresh = 10^oldThresh; %Remove log transform
                oldThresh = double((1/oldThresh)*360);
            else
                                oldThresh =  numericSol;
            end
            CIint(1) = abs(CI(1) - oldThresh);
            CIint(2) = abs(CI(2) - oldThresh);
            
            if ~isstr(thresh)
                CI(2) = thresh+CIint(2);
                CI(1) = thresh-CIint(1);
            else
                CI = 'unfit';
            end
            
            
        elseif whichFitter==2
            
            
            
        end
        
        
        
        disp('done:')
        if ~isstr(thresh)
            message = sprintf('Threshold estimate: %6.4f',thresh);
            disp(message);
            message = sprintf('Slope estimate: %6.4f',slope);
            disp(message);
            message = sprintf('Lapse estimate: %6.4f\r',lapse);
            disp(message);
        else
            message = sprintf('Threshold not fittable');
            disp(message);
        end
        
        %Store our fit
        
        
        currSess.fitThresh = thresh;
        currSess.fitSlope = slope;
        currSess.CI = CI;
        currSess.lapse = lapse;
        
        allData{end+1} = currSess;
        
        
        try
            if GOF
                figH = figure;
                GoodnessOfFit(results)
                waitfor(figH);
            end
        catch
            disp('Goodness of Fit not functional for this fit, probably too messy to even try this one')
        end
        
        clear currSess NumPos OutOfNum StimLevels;
    end
    
end

%Init vars

Date={};
Thresh={};
Slope={};
CI = {};
CIlow={};
CIhigh={};
Lapse ={};
PC=[];
nTrials = [];
clear threshRange
clear threshFlipped
for pp = 1:length(allData)
    
    
    if pp>1
        if strcmp(allData{pp}.date, allData{pp-1}.date)
            continue
        end
    end
    
    if length(allData{pp}.date)<8
        Date{end+1} = allData{pp}.date(1:end);
    else
        Date{end+1} = allData{pp}.date(1:8);
    end
    
    
    PC(end+1) = allData{pp}.pC;
    nTrials{end+1} = allData{pp}.nTrials;
    CI{end+1} = allData{pp}.CI;
    
    if isstr(allData{pp}.fitThresh)
        Thresh{end+1} = allData{pp}.fitThresh;
        CIlow{end+1} = 'unfit';
        CIhigh{end+1} = 'unfit';
        Lapse{end+1} = allData{pp}.lapse;
        Slope{end+1} = allData{pp}.fitSlope;
        
    else
        Thresh{end+1} = sprintf('%4.2f',allData{pp}.fitThresh);
        Slope{end+1} = sprintf('%4.2f',allData{pp}.fitSlope);
        if allData{pp}.CI(1) <0
            CIlow{end+1} = 0;
        else
            CIlow{end+1} = sprintf('%4.2f',allData{pp}.CI(1));
        end
        CIhigh{end+1} = sprintf('%4.2f',allData{pp}.CI(2));
        Lapse{end+1} = sprintf('%1.2f',allData{pp}.lapse);
        
       
        
    end
    
end

close(wBar)
Lapse = Lapse';
Thresh = Thresh';
CIlow = CIlow';
CIhigh = CIhigh';
PC = PC';
Slope = Slope';
X = stim.locx(1:length(Thresh));
Y = stim.locy(1:length(Thresh));
nTrials=nTrials';

%Make a quick loop through the Date variable and look for duplicates,
%modify as necessary

if length(Date) ~= length(unique(Date))
    
    for ff = 1:length(Date)
        currDate = Date(ff);
        dupes = sum(strcmp(currDate,Date));
        
        if dupes >1 
            dupeInds = strcmp(currDate,Date);
            dupeDates = Date(dupeInds);
            dupeFind = find(dupeInds==1);
            for qq = 2:length(dupeDates)
                Date(dupeFind(qq)) = strcat(Date(dupeFind(qq)),'-',num2str(qq));
            end
            
            
        end
        
    end
    
end

CBTrainingTable = table(Thresh, CIlow, CIhigh, PC,Slope,Lapse,X, Y,nTrials,'RowNames',Date)
keyboard
if saveData
    
    writetable(CBTrainingTable,strcat(dataFile,'_FitOn_',date,'.csv'),'FileType','text','WriteRowNames', 1)
    
    
end
warning('on')



function [fitresult, gof] = createFit(logSL, PC, OON, plotData, subplotN, subplotStart, subplotEnd, doPlot)
%CREATEFIT(TESTL,PC)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : testl
%      Y Output: PC
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 12-Jan-2016 11:40:57


%% Fit: 'Weibull Fit'.
[xData, yData, weights] = prepareCurveData( logSL, PC, OON );

% Set up fittype and options.
ft = fittype( 'lapse-exp(-exp(2*s*m/log(2)*(log(x)-log(m))+log(log(2))))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0.7 0 eps];
opts.Upper = [1 Inf 200];
% opts.StartPoint = [0.392227019534168 0.655477890177557];
opts.Weights = OON./max(OON);

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
if plotData
    global unfitN
    % Plot fit with data.
    plotSubFig=-1;
    if subplotN
        
        totalPlots = subplotEnd-subplotStart;
        %         if totalPlots>20
        %
        %             if (subplotN-subplotStart)+1>20
        %
        %                 plotSubFig=plotSubFig+1;
        %                 figure
        %             end
        %             subplot(5, 4,  ((subplotN-subplotStart)+1)-unfitN-(plotSubFig*20))
        %
        %         else
        %
        %             subplot(ceil(totalPlots/4), 4,  (subplotN-subplotStart)+1)
        %
        %         end
        
        if totalPlots>20 && (subplotEnd-subplotN)<21 && doPlot
            
            
            subplot(5, 4,  (subplotEnd-subplotN)+1)
            h = plot( fitresult, xData, yData);
  
            %     set(h(1),'MarkerSize',20)
            mSize = 10; %base marker size
            %     legend( h, 'Percent Correct vs Stimulus Levels', 'Weibull fit', 'Location', 'NorthEast' );
            % Label axes
            xlabel('Stimulus Intensity')
            ylabel PC
            grid on
            hold on
            %Redo plot with markersize proportional to nTrials weight
            for i=1:length(logSL);
                sz = mSize*(OON(i)/max(OON));
                plot(logSL(i),PC(i),'ko','MarkerFaceColor','b','MarkerSize',sz);
            end
            
            if max(logSL)==91
                
                %Make custom X marks
                %         xV = 1:10:91;
                %         xAx = 1./(xV./91); %Invert and normalize
                %
                %
                %         xAxi = log(xAx); %log transform data
                %         xAxi(end) = 0.05;
                %         set(gca, 'XTick', xAxi);
                %         axis([0 95 0.5 1]);
                %
                %         set(gca, 'XTickLabel', {'0','10', '20','30', '40','50','60','70', '80', '90'});
                set(gca, 'FontSize', 12);
                %         set(gca, 'XScale', 'log');
                legend HIDE
                
            else
                
                
                %Make custom X marks
                xV = 1:40:361;
                xAx = 1./(xV./361); %Invert and normalize
                
                
                xAxi = log(xAx); %log transform data
                xAxi(end) = 0.05;
                set(gca, 'XTick', flip(xAxi));
                set(gca, 'XTickLabel', {'0','40', '80','120', '160','200','240','280', '320', '360'});
                set(gca, 'FontSize', 12);
                set(gca, 'XScale', 'log');
                legend HIDE
                axis([0.05 5.886 0.5 1]);
                
                
            end
            %         else
            %             subplot(ceil(totalPlots/4), 4,  (subplotN-subplotStart)+1)
            %
        end
        
    else
        
        figure( 'Name', 'Weibull Fit' );
        
    end
    
    h = plot( fitresult, xData, yData);
              hold on
            xrange = min(xData):max(xData);
            plot(xrange,repmat(0.75, 1,length(xrange)), 'b', 'LineWidth', 3);
    %     set(h(1),'MarkerSize',20)
    mSize = 10; %base marker size
    %     legend( h, 'Percent Correct vs Stimulus Levels', 'Weibull fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel('Stimulus Intensity')
    ylabel PC
    grid on
    hold on
    %Redo plot with markersize proportional to nTrials weight
    for i=1:length(logSL);
        sz = mSize*(OON(i)/max(OON));
        plot(logSL(i),PC(i),'ko','MarkerFaceColor','b','MarkerSize',sz);
    end
    if max(logSL)==91
        
        %Make custom X marks
        %         xV = 1:10:91;
        %         xAx = 1./(xV./91); %Invert and normalize
        %
        %
        %         xAxi = log(xAx); %log transform data
        %         xAxi(end) = 0.05;
        %         set(gca, 'XTick', xAxi);
        %         axis([0 95 0.5 1]);
        %
        %         set(gca, 'XTickLabel', {'0','10', '20','30', '40','50','60','70', '80', '90'});
        set(gca, 'FontSize', 12);
        %         set(gca, 'XScale', 'log');
        legend HIDE
        
    else
        
        
        %Make custom X marks
        xV = 1:40:361;
        xAx = 1./(xV./361); %Invert and normalize
        
        
        xAxi = log(xAx); %log transform data
        xAxi(end) = 0.05;
        set(gca, 'XTick', flip(xAxi));
        set(gca, 'XTickLabel', {'0','40', '80','120', '160','200','240','280', '320', '360'});
        set(gca, 'FontSize', 12);
        set(gca, 'XScale', 'log');
        legend HIDE
        axis([0.05 5.886 0.5 1]);
        
        
    end
    
    
end



