%% Code by Jace Curran
% Last published 03/06/2025 on jacecurran.com

%% What's all this then?
% This code uses PicoHarp .dat files ONLY to analyze fluorescence time
% dynamics, quantum beat frequency and amplitude, and boost behavior.

% Why does fit plotting in R2024b look so chunky? The fit object itself
% seems to be the same, but when plotted the fit looks hilariously low-res.
% Until I understand why, I'm using R2023.

%% High Level Static Variables
% Hands off!
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
warning('off','curvefit:fit:noStartPoint');
TIDamp = 0;

%% High Level Tuneable Variables
% TODO: The smoothing window should not be used combined with the linNorm
% function. Smoothing introduces a high variance in peak height, which
% destroys the data when normalization occurs.
smoothingWindow = 0;    % Set window size for averaging / smoothing
beatDisp = 30;          % Beat data will be graphed from 0.1 to beatDisp ns
analysisMode = 0;       % Switch for boost (0) or beat (1) analysis
fftMode = 0;            % Switch for cos (0) or FFT (1) fitting
beatDecay = 0;          % Switch for standard cos (0) or TID-type fitting (1)
fitDivide = 0;          % For beat fitting, divide by fit (0) or another trial (1)?
linNorm = 0;            % For boost only: should the active trial be linearly shifted
                        % so that the peak aligns with the control trial?
intMode = 0;            % Switch for int mode; display sum of photons per unit time
                                % (only works in boost mode)
% NOTE: intMode = 0 has been depreciated, may not work.
fitOverride = 0;        % Skip fitting to avoid issues with wide bin trials.
                        % Fine for boost analysis, may cause issues with beats

%% File IO and Fluorescence Time-Dynamics
% Open the file via explorer window, parse enough just to get the time step
% (recorded by the PicoHarp), then read all data as a table.
[fileName,filePath] = uigetfile('.dat');
trash = readlines(join([filePath fileName]));
timeStep = double(extractBefore(trash(9),char(9)));
allData = table2array(readtable(join([filePath,fileName]),'NumHeaderLines',10));

timeArray = 0:timeStep:((length(allData)-1)*timeStep);

% Optional Smoothing
if smoothingWindow
    allData = smoothdata(allData,'movmean',smoothingWindow);
end

% Fitting with logarithmic error is... complicated (to me). Instead,
% we'll simply have a parallel data structure that contains the log of
% each photon bin. We also find the peak of each data set; we'll set t=0 to
% be the lowest peak value (negative before, positive after). Looking
% ahead, the domain of the fit function is (-1,inf]; we'll actually
% truncate the log time array to begin 1 step after -1.
% TODO: Could implement a shifting algo to make all peaks the same
[~,peakPos] = max(allData);
minPeak = min(peakPos);
logData = log(allData);
logTimes = -(timeStep*minPeak):timeStep:timeStep*(length(allData)-1-minPeak);
n1Index = find(logTimes>-1,1);
logTShort = logTimes(n1Index:end);

% Cool, now we define a cell array, with each cell containing a fit object.
% We have n such fits / cells where n is the number of trials.

if ~fitOverride
    fitArray = cell(1,width(allData));
    for i = 1:width(logData)
        timeToFit = logTimes(minPeak:end);
        dataToFit = logData(minPeak:end,i);
        fitEquation = 'log(a1*exp(-1*x/b1)+c1*exp(-1*x/d1)+o1*(1+x)^p1+q1*(1+x)^r1)';
        fitTypeOpt = fittype(fitEquation,'independent','x');
        lowerBounds = [0 0 0 0 0 -Inf 0 -Inf];
        startPointArray = [35000 1.2 11000 4.0 800000 -5.1 7000 -0.35];
        [fitArray{i},~] = fit(timeToFit',dataToFit,fitTypeOpt,'Lower',...
            lowerBounds,'StartPoint',startPointArray);
    end
end
    
    % First Plot is lin-log, second plot is log-log
    figure(1);
    plot(logTimes,logData);
    legend();
    xlabel("Time (ns)");
    ylabel('10^n Photon Events');
    title('Semi-log Fluorescence Data');
    grid on;
    
    figure(2);
    plot(logTimes,logData);
    legend();
    xlabel("Time (ns)");
    ylabel('10^n Photon Events');
    title('Log-Log Fluorescence Data');
    set(gca,'Xscale','log');
    grid on;

% Luke Code: Truncate all data at first zero
% [row, col] = size(allData);
% zeroInd = [];
% for i = 1:col
%     zeroInd = [zeroInd, find(allData(:,i) == 0, 1)];
% end
% allData = allData(1:min(zeroInd - 1),:);

if analysisMode
    % First, all data trials need to be divided by their fit. We feed each
    % trials fit function the log time array, and then divide (log
    % properties make it look like subtraction) to find the beat
    % information without the background exp/power law. 
    divData = zeros(length(allData)-n1Index+1,width(allData));
    if fitDivide
        divTrial = str2num(input('Which trial should be used as the control? ','s'));
    end
    for i = 1:width(allData)
        fitFunk = fitArray{i};
        if fitDivide
            divData(:,i) = exp(logData(n1Index:end,i)-logData(n1Index:end,divTrial));
        else
            divData(:,i) = exp(logData(n1Index:end,i)-fitFunk(logTShort));
        % DELETE. For Xallan trials only
        % divData(:,i) = exp(logData(n1Index:end,i)-logData(n1Index:end,9));
        end
    end

    actTrial = str2num(input('Which trial would you like to analyze?: ','s'));

    while actTrial == 0 || actTrial > width(allData)
        actTrial = str2num(input(['Please enter a trial number between 1 and ', ...
            num2str(width(allData)),': '],'s'));
    end
    while actTrial > 0

        % First step, just plot the divided data
        figure(3);
        plot(logTShort,divData(:,actTrial));
        % title(['Trial ',num2str(actTrial),' Data / Fit']);
        grid on;
        xlim([0.1 beatDisp]);
        xlabel('Time after peak (ns)');
        ylabel('Photon Flux Ratio');

        if fftMode
            % We want the user to specify the range of FFT analysis
            fftStart = str2num(input('Where to begin FFT analysis? (ns): ','s'));
            fftEnd = str2num(input('Where to end FFT analysis? (ns): ','s'));
            fftStartInd = find(logTShort > fftStart,1);
            fftEndInd = find(logTShort > fftEnd,1);

            % The following code is... inspired by MATLAB documentation
            fftData = divData(fftStartInd:fftEndInd,actTrial);
            fftTimes = fftStart:1:fftEnd;
            fftArray = fft(fftData);
            fftFreqs = (0:length(fftArray)-1)/length(fftArray);
            fftFreqs = fftFreqs(2:floor(end/2))./timeStep;
            fftArray = fftArray(2:floor(end/2));
            figure(3);
            plot(fftFreqs,abs(fftArray));
            grid on;
            set(gca,'XMinorGrid','on','YMinorGrid','on');
            title('Frequency Spectra of Beat Data');
            xlabel('Frequency (GHz)');
            ylabel('Relative Strength');
            xlim([0 20]);

            % Last thing we'll do is print the top 3 frequencies and
            % their relative amplitudes
            [mFreqs,mInds] = maxk(abs(fftArray),3);
            fftTable = table(round([fftFreqs(mInds(1)) fftFreqs(mInds(2)) fftFreqs(mInds(3))],4)', ...
                round([mFreqs(1) mFreqs(2) mFreqs(3)],4)', ...
                'VariableNames',{'Frequency (GHz)','Relative Amp'},'RowNames',{'1st','2nd','3rd'});
            disp(fftTable);

        else
            % If we're not running an FT, we'll fit with a cos function

            % We'll ask the user for the region to fit, and an approximation for 
            % frequency, then apply a cosine fit.
            bfStart = str2num(input('Where to begin beat analysis? (ns): ','s'));
            bfEnd = str2num(input('Where to end beat analysis? (ns): ','s'));
            bfOsc = str2num(input('How many oscillations in this span? ','s'));
            bfStartInd = find(logTShort > bfStart, 1);
            bfEndInd = find(logTShort > bfEnd, 1);
            fitTimes = logTShort(bfStartInd:bfEndInd);
            fitData = divData(bfStartInd:bfEndInd,actTrial);
            bfFreqLin = bfOsc / (bfEnd - bfStart);
            bfEquation = 'a2*cos(2*pi*b2*x+c2) + d2';
            bfOptions = fittype(bfEquation,'independent','x');
            spArray = [0.1 bfFreqLin 0 0];
            lowerBF = [0 0 -3.14159 -Inf];
            upperBF = [Inf Inf 3.14159 Inf];
            [beatFit,~] = fit(fitTimes',fitData,bfOptions,'StartPoint',spArray ...
                ,'Lower',lowerBF,'Upper',upperBF);
    
            if beatDecay
                % If we have an exponentially decaying cos fit, we'll lock all
                % previous parameters, and apply an exponential decay. All we
                % need is an assumed 'noTID' amplitude once per run (get 
                % that from a known symmetry direction).
    
                % TODO: For now, the TID analysis limits are inherited, one
                % block up. We could have different limits
                if TIDamp == 0
                    TIDamp = str2num(input('Enter the assumed no-dephasing amplitude: ','s'));
                end
                bfTIDEq = 'a3*exp(-x*f3)*cos(2*pi*b3*x+c3) + d3';
                bfTIDOp = fittype(bfTIDEq,'independent','x');
                lbTIDArray = [TIDamp beatFit.b2 beatFit.c2 beatFit.d2 0];
                ubTIDArray = [TIDamp beatFit.b2 beatFit.c2 beatFit.d2 Inf];
                [TIDfit,~] = fit(fitTimes',fitData,bfTIDOp, ...
                    'Lower',lbTIDArray,'Upper',ubTIDArray);
    
                % Pretty much done, we just plot and display data
                figure(4);
                clf;
                hold on;
                plot(TIDfit,logTShort,divData(:,actTrial),'-');
                plot(TIDfit);
                xline(bfStart,'--');
                xline(bfEnd,'--');
                hold off;
                legend('off');
                grid on;
                % title('Beat (TID) Fit and Data');
                xlabel('Time after peak (ns)');
                ylabel('Photo Flux Ratio');
                xlim([0.1 beatDisp]);
    
                % Parameter and error display
                conf = confint(beatFit);
                confTID = confint(TIDfit);
                % One quick note; freezing in the previous cosine function 
                % means we need to use the old confidence intervals for error
                % for amplitude, freq, and phase.
                dataTable = table(round([abs(TIDfit.a3) TIDfit.b3 TIDfit.c3 TIDfit.f3],4)', ...
                    round([0.5*(conf(2,1)-conf(1,1)) 0.5*(conf(2,2)-conf(1,2)) 0.5*(conf(2,3)-conf(1,3)) 0.5*(confTID(2,5)-confTID(1,5))],4)', ...
                    'VariableNames',{'Value','+/-'},'RowNames',{'Amplitude','Frequency (GHz)','Phase','k (TID)'});
                disp(dataTable);
    
            else
    
                % Pretty much done, we just plot and display data
                figure(4);
                clf;
                hold on;
                plot(beatFit,logTShort,divData(:,actTrial),'-');
                xline(bfStart,'--');
                xline(bfEnd,'--');
                hold off;
                legend('off');
                grid on;
                % title('Beat Fit and Data');
                xlabel('Time after peak (ns)');
                ylabel('Photon Flux Ratio');
                xlim([0.1 beatDisp]);
        
                % We used to save fit data and display on program close, but
                % the way I process data, there's really no need. Instead, we'll
                % print freq amp and phase here, along with error.
                conf = confint(beatFit);
                dataTable = table(round([abs(beatFit.a2) beatFit.b2 beatFit.c2],4)', ...
                    round([0.5*(conf(2,1)-conf(1,1)) 0.5*(conf(2,2)-conf(1,2)) 0.5*(conf(2,3)-conf(1,3))],4)', ...
                    'VariableNames',{'Value','+/-'},'RowNames',{'Amplitude','Frequency (GHz)','Phase'});
                disp(dataTable);
    
            end

        end
        
        % All we have to do now is repeat the loop until the user quits
        % (anything that isn't a positive integer will end the program)
        actTrial = str2num(input('Which trial would you like to analyze? ("n" to quit): ','s'));
        if actTrial >= 0
            while actTrial == 0 || actTrial > width(allData)
                actTrial = str2num(input(['Please enter a trial number between 1 and ', ...
                    num2str(width(allData)),': '],'s'));
            end
        end
    end

    disp('Thank you for using BBA2024!');

else
    % Boost Analysis
    % For the boost, we don't divide by each trial's fit, but by a 
    % control trial (usually either B=0 or Xallan angle)
    conTrial = str2num(input('Please identify the control trial: ','s'));
    while conTrial == 0 || conTrial > width(allData)
        conTrial = str2num(input(['Please enter a trial number between 1 and ', ...
            num2str(width(allData)),': '],'s'));
    end

    actTrial = str2num(input('Which trial would you like to analyze?: ','s'));
    while actTrial == 0 || actTrial > width(allData)
        actTrial = str2num(input(['Please enter a trial number between 1 and ', ...
            num2str(width(allData)),': '],'s'));
    end

    while actTrial > 0

    actData = allData(:,actTrial);
    conData = allData(:,conTrial);
    % Normalization: we'll linearly boost the peak of the act trial
    % to match the con trial. This likely does almost nothing.
    if linNorm
        [actMax,actInd] = max(actData);
        [conMax,conInd] = max(conData);
        actData = actData ./ (actMax / conMax);
    end

    % Truncate the data from 0 ns onwards
    zeroInd = find(logTimes>0,1);
    boostTimes = logTimes(zeroInd:end);
    actData = actData(zeroInd:end);
    conData = conData(zeroInd:end);

    % For a quick debug, we'll display the control and active trials,
    % shifted such that the peak occurs at t = 0
    figure(3);
    clf;
    hold on;
    plot(boostTimes,conData);
    plot(boostTimes,actData);
    hold off;
    legend('Control Trial','Active Trial');
    xlabel("Time (ns)");
    ylabel('10^n Photon Events');
    title('Boost Preview');
    grid on;
    set(gca,'YScale','log');
    xlim([0 30]);

    % If intMode, create a parallel array for both con and act trials, and
    % make that the sum of all photons received by that point. Then plot.
    % Or better yet, have MATLAB do the algo for you.
    if intMode
        figure(4);
        clf;
        hold on;
        plot(boostTimes,cumsum(conData));
        plot(boostTimes,cumsum(actData));
        xlabel("Time (ns)");
        ylabel('10^n Cumulative Sum of Photon Events');
        title('Boost Preview - Cumulative Sum');
        grid on;
        set(gca,'YScale','log','XScale','log');
        % set(gca,'XScale','log');
        % xlim([0 60]);

        % We'll pull out two different boost occurances:
        % 0-1 ns for historic reasons, 1-3 ns and 1-5 ns to avoid IRF
        ind1ns = find(boostTimes>=1,1);
        ind3ns = find(boostTimes>=3,1);
        ind5ns = find(boostTimes>=5,1);
        oneBoost = sum(actData(1:ind1ns)) / sum(conData(1:ind1ns));
        threeBoost = sum(actData(ind1ns:ind3ns)) / sum(conData(ind1ns:ind3ns));
        fiveBoost = sum(actData(ind1ns:ind5ns)) / sum(conData(ind1ns:ind5ns));
        disp(['0-1 ns boost ratio: ', num2str(oneBoost,'%.4f')]);
        disp(['1-3 ns boost ratio: ', num2str(threeBoost,'%.4f')]);
        disp(['1-5 ns boost ratio: ', num2str(fiveBoost,'%.4f')]);
        xline(boostTimes(ind1ns),'--');
        xline(boostTimes(ind3ns),'--');
        xline(boostTimes(ind5ns),'--');
        legend('Control Trial','Active Trial','');
        hold off;
    else

        % Now for analysis; we'll actually have the user identify the crossing
        % point. I know, but there's enough variation in how these trials look
        % that an automated algorithm would be a little tricky
        cPoint = str2num(input('Identify the crossing point ("n" for none): ','s'));
        % TODO: Catch some bad input here
    
        % Then, we'll step through and calculate the boost ratio (sum of all
        % photons received in region in active trial/control trial) and max
        % boost (largest difference in single t photo count between active and
        % control trials). This will be done for two different regions if a
        % crossing point has been identified.
    
        % TODO: Right now, the 'max boost' quant is useless, especially at t>>0
        if cPoint > 0
            cPointInd = find(boostTimes>cPoint,1);
            actSum1 = sum(actDataS(1:cPointInd));
            actSum2 = sum(actDataS(cPointInd:end));
            conSum1 = sum(conDataS(1:cPointInd));
            conSum2 = sum(conDataS(cPointInd:end));
            boostRat1 = actSum1/conSum1;
            boostRat2 = actSum2/conSum2;
            boostRatO = (actSum1+actSum2)/(conSum1+conSum2);
            maxBoost1 = 0;
            maxBoost2 = 0;
            maxBoostO = 0;
            for i = 1:length(actDataS(1:cPointInd))
                if abs(actDataS(i)/conDataS(i)) > abs(maxBoost1)
                    maxBoost1 = actDataS(i)/conDataS(i);
                end
            end
            for i = cPointInd:length(actDataS(cPointInd:end))
                if abs(actDataS(i)/conDataS(i)) > abs(maxBoost2)
                    maxBoost2 = actDataS(i)/conDataS(i);
                end
            end
            if abs(maxBoost1) > abs(maxBoost2)
                maxBoostO = maxBoost1;
            else
                maxBoostO = maxBoost2;
            end
            boostTable = table([cPoint cPoint cPoint]',round([boostRat1 boostRat2 boostRatO],4)' ...
                ,[maxBoost1 maxBoost2 maxBoostO]', 'VariableNames', ...
                {'Crossing Point' 'Boost Ratio' 'Max Boost'},'RowNames', ...
                {'1st Region' '2nd Region' 'Overall'});
            disp(boostTable);
        else
            actDataS = actData;
            conDataS = conData;
            actSum = sum(actDataS);
            conSum = sum(conDataS);
            boostRat = actSum/conSum;
            maxBoost = 0;
            for i = 1:length(actDataS)
                if abs(actDataS(i)/conDataS(i)) > abs(maxBoost)
                    maxBoost = actDataS(i)/conDataS(i);
                end
            end
            boostTable = table("None",round(boostRat,4),maxBoost,'VariableNames', ...
                {'Crossing Point' 'Boost Ratio' 'Max Boost'},'RowNames',"Value");
            disp(boostTable);
        end
    end

    % All we have to do now is repeat the loop until the user quits
    % (anything that isn't a positive integer will end the program)
    actTrial = str2num(input('Which trial would you like to analyze? ("n" to quit): ','s'));
        if actTrial > 0
            while actTrial == 0 || actTrial > width(allData)
                actTrial = str2num(input(['Please enter a trial number between 1 and ', ...
                    num2str(width(allData)),': '],'s'));
            end
        end
    end

    disp('Thank you for using BBA2024!');
end



