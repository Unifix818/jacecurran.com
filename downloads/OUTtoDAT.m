%% Code by Jace Curran
% Last published 03/06/2025 on jacecurran.com

% This program will accept a .out file from PTUtoOUT.m 
% (either T2 or T3 mode) and parse the resulting data into a .dat file.
% This is most useful if you're using my TDPLAnalyzer.m code, otherwise
% you could use a different, more useful final format.

% Note that the only things required are:
% Channels and TimeTags from T2 or
% DTime from T3
[fName,fPath] = uigetfile('*.out')

% First step, save all lines in an array. If 'CHN' is detected, we're 
% dealing with a T2 file; else, T3.
text = readlines(join([fPath,fName]));
t3Res = 0;
t2Width = 4;
if contains(text(4),'CHN') || contains(text(5),'CHN')
    mode = 'T2'
    t2Width = str2num(input('Enter bin size in ps (min 4, please use powers of 2): ','s'));
    % Doesn't actually check for pow of 2, just checks if even. You'd need
    % a log base 2 or some junk
    if t2Width < 4 || mod(t2Width,2)~=0
        t2Width = 4;
    end
else
    mode = 'T3'
    t3Res = str2num(input('Enter the bin size in ps used for T3 (hopefully included in info.dat): ','s'));
end

keepGoing = 'y';

% The processing in T2/T3 is very different, so we split here
if strcmp(mode,'T3')
    % Bins is a 2D array, allocated to be 4096 x 16 for speed, but we'll cut
    % down at the end
    bins = zeros(4096,16);
    trialInd = 1;

    while strcmpi(keepGoing,'y')
        % First, we have a list of DTimes. Kill the first line, I don't want it.
        % Take DTime * t3Res to get actual time of arrival in ps.
        % NOTE: This is where the limiting factor of the chosen resolution
        % comes into play. Data limit = 4096*t3Res. So, for the highest
        % resolution bin size of 4 ps, the data will only extend to 16.380 ns.
        text(1) = [];
        text = str2double(text);
        text = text .*(t3Res);
        text = rmmissing(text); % Removes errors
        
        % Now we have arrival time in ps, but we need a histogram to be
        % compatible with PicoHarp analysis. Two cools things: numBin is
        % already set and static (4096), and the arrival times have t3Res steps
        % between.
        for i=1:length(text)
            binNum = text(i)/t3Res + 1;
            if ~binNum
                binNum = 1;
            end
            bins(binNum,trialInd) = bins(binNum,trialInd) + 1;
        end

        load gong.mat;
        y=y/5;
        sound(y);
        % Ask user if there's an additional OUT to tack onto this DAT
        keepGoing = input('Add another OUT as an additional trial? (y or n): ','s');
        if strcmpi(keepGoing,'y')
            trialInd = trialInd + 1;
            [fName,fPath] = uigetfile('*.out')
            text = readlines(join([fPath,fName]));
        end
    end
    bins(:,trialInd+1:end) = [];
else
    % Notes on why the bins array is sized this way are further down
    bins = zeros(1.7e7,8);
    trialInd = 1;

    while strcmpi(keepGoing,'y')
        % In T2, we have a list of TimeTags for Ch0 (laser) and Ch1 (sample).
        % First, kill the first line, I don't want it.
        % Go through each line, constantly updating the latest Ch0 time tag
        % When you see a Ch1 time tag, take the relative difference to the last
        % Ch0 time tag, that's our DTime array.
        dTime = zeros(2,1);
        text(1) = [];
        lastZero = 0;
        for i=1:length(text)
            if contains(text(i),'CHN 0')
                num = erase(text(i),['CHN 0',' ']);
                lastZero = str2num(num); % numnum
            elseif contains(text(i),'CHN 1')
                num = str2num(erase(text(i),['CHN 1',' ']));
                dTime(end+1) = num - lastZero;
            end
        end
        dTime(1:2) = [];
        % Now, take dTime * 4 to get true dTime in ps.
        dTime = dTime.*4;
    
        % Okay, now create the histogram. The largest timeTag received
        % will depend on the rep rate of the laser, but fundamentally
        % the data is stored in 24 bits. 24 bits = 68 us, the largest
        % conceivable dTime. If we use the highest possible resolution (you
        % could go higher but I'm implementing a limit here of 4 ps, that's
        % 1.67e7 bins. So we should be covered no matter what! with 1.7e7 bins.
        % If I've done some math wrong here, the program will throw an error,
        % so don't worry about it.
        % Then, we'll take a look at the end for when the empty bins begin, and
        % delete all bins until the last filled bin. That's our arrival time
        % histogram!
        for i = 1:length(dTime)
            binNum = floor(dTime(i)/t2Width);
            if ~binNum
                binNum = 1;
            end
            bins(binNum,trialInd) = bins(binNum,trialInd) + 1;
        end
        if length(bins) > 1.7e7
            disp('ERROR: Histogram overflow. Run this code again with a larger bin size.');
        end
        lastVal = find(bins(:,trialInd),1,'last');
        bins(lastVal+1:end) = [];

        load gong.mat;
        y=y/5;
        sound(y);
        % Ask user if there's an additional OUT to tack onto this DAT
        keepGoing = input('Add another OUT as an additional trial? (y or n): ','s');
        if strcmpi(keepGoing,'y')
            trialInd = trialInd + 1;
            [fName,fPath] = uigetfile('*.out')
            text = readlines(join([fPath,fName]));
        end
    end
    bins(:,trialInd+1:end) = [];
end

outText="";

endFile = replace(fName,'.out','.dat')
fid = fopen(append(fPath,endFile),'wt');

% Header Text
outText = outText + "#PicoHarp 300 Histogram Data from T2/T3 Mode \n";
outText = outText + "#channels per curve\n";
outText = outText + "65536\n";
outText = outText + "#display curve no.\n";
outText = outText + "0\t1\n";
outText = outText + "#memory block no.\n";
outText = outText + "0\t1\n";
outText = outText + "#ns/channel\n";
if strcmp(mode,'T3')
    outText = outText + append(num2str(t3Res/1000),"\t");
    outText = outText + append(num2str(t3Res/1000),"\n");
else
    outText = outText + "0.0040\t0.0040\n";
end
outText = outText + "#counts\n";

fprintf(fid,outText);
outText = '';

% Convert to file writing rather than building one big string
for row = 1:length(bins)
    for col = 1:width(bins)
        outText = outText + string(bins(row,col)) + "\t";
    end
    outText = outText + "\n";
    fprintf(fid,outText);
    outText = "";
end

fclose(fid);