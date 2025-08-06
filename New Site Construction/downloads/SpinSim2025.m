%% Code by Jace Curran
% Last published 03/06/2025 on jacecurran.com

%% What's all this then?
% This code simulates the triplet exciton energies, singlet (S=0)
% overlaps, quantum beat frequencies and amplitudes, and expected
% fusion / fluorescence boost using the following assumptions:
% 1) The unit cell contains two differently-oriented molecules
% 2) The material is capable of undergoing singlet fission
% 3) The Hamiltonian for the system is dependent on zero-field spin
% energies, Zeeman shift due to an external magentic field, and
% dipole-dipole interaction between two triplets.

% NOTE: As of 1/6/25, the AB fix is now sometimes required again

%% High Level Static Variables
configRGB = {[1 0 0] [0 0 1] [148 51 255]/256 [0 1 0]};
configTitles = ["AA" "BB" "AB" "Averaged"];

% Bohr Magneton includes a pre-baked gfactor of 2
bohrMagneton = 2.8e-3;                  % GHz / Gauss
% Luke bM: bohrMagneton = 1.78e-2;
expBoostPoints = false;                 % Do Not Change
simBoostPoints = false;                 % Do Not Change
set(groot,'DefaultLineLineWidth',1);

%% High Level Tuneable Variables
% Magnetic Field Angles (W.R.T. Crystal axis (axis between rotated
% molecules). Defined s.t.      <0,0> or <0,90> is c
%                               <90,0> is a
%                               <90,90> is b
% Unused in magMode 1 (angle)
magTheta =  90;
magPhi =    90;
% Xallan Angle approx at (39.75,0)

% Magnetic Field Strength (max value if in mag mode, only value if in angle
% mode)
magStrMax = 5000;

% Loop through ||B|| (0) or angle (1)?
magMode = 0;

% Beat Prediction Bool (speeds up processing if left off)
predictMode = 0;

% How much of the predicted boost is the AB case?
ABBoost = 0.50;

% Should we try to clean up degenerate state projections?
degenFix = 0;

% Resolution Parameters (larger num = higher res) 
% (too low may lead to eigen-flipping)
nSteps = 1000;

% EPR (ZFS) Parameters in GHz (See notes)
% Rubrene Molecular
dZFS = 1.6564;
eZFS = -0.11872;

% % Tetracene Set 1
% dZFS = -0.1858;
% eZFS = 0.7435;
% 
% % Tetracene Set 2
% dZFS = 1.5590;
% eZFS = -0.1559;

% % Tetracene Alt E
% eZFS = -0.1199;

% Tait-Bryan Angles (ZYX), rotates between Crystal axes to MolA axes (See notes)
% Rubrene Angles
tbAlpha = 31;
tbBeta  = 0;
tbGamma = 0;

% % Tetracene Angles
% tbAlpha = 0;
% tbBeta = 26;
% tbGamma = 0;

% Nearest neighbor vector
nnVec = [1 0 0];
nnVec = nnVec/norm(nnVec);
nnX = nnVec(1);
nnY = nnVec(2);
nnZ = nnVec(3);

% Interaction Constant Strength (GHz)
% intStr = 0.008;
% intStr = 8e-8;
intStr = 8e-10;
% intStr = 0;
% Wang 2015 recommends 0.008, historically we've used 8e-8, lowest 
% value that avoids eigenflipping seems to be ~10^-12

%% ZFS and DP-DP Terms

% Rotation Matrix for mol A to crystal
yawRotA = [cosd(tbAlpha) -sind(tbAlpha) 0; sind(tbAlpha) cosd(tbAlpha) 0; 0 0 1];
pitRotA = [cosd(tbBeta) 0 sind(tbBeta); 0 1 0; -sind(tbBeta) 0 cosd(tbBeta)];
rolRotA = [1 0 0; 0 cosd(tbGamma) -sind(tbGamma); 0 -sind(tbGamma) cosd(tbGamma)];
rotA = yawRotA*pitRotA*rolRotA;

% Rotation Matrix for mol B to crystal
yawRotB = [cosd(tbAlpha) sind(tbAlpha) 0; -sind(tbAlpha) cosd(tbAlpha) 0; 0 0 1];
pitRotB = [cosd(tbBeta) 0 -sind(tbBeta); 0 1 0; sind(tbBeta) 0 cosd(tbBeta)];
rolRotB = [1 0 0; 0 cosd(tbGamma) sind(tbGamma); 0 sind(tbGamma) cosd(tbGamma)];
rotB = yawRotB*pitRotB*rolRotB;

% Construction of the ZFS term for AA, AB, BB, and DL (delocalized)
% The ZFS term is 3 dimensional (9 x 9 x 4), where each page is
% configuration (AA, AB, BB, DL)
% The ZFS matrix has been cyclically rotated such that abc = xyz
% zeroMini = [dZFS/3-eZFS 0 0; 0 dZFS/3+eZFS 0; 0 0 -2*dZFS/3];
% zeroMini = 2*pi*[-2*dZFS/3 0 0; 0 dZFS/3-eZFS 0; 0 0 dZFS/3+eZFS];
zeroMini = [-2*dZFS/3 0 0; 0 dZFS/3-eZFS 0; 0 0 dZFS/3+eZFS];
zeroMiniAA = rotA'*zeroMini*rotA;
zeroMiniBB = rotB'*zeroMini*rotB;
zeroMiniAvg = 0.5*(zeroMiniAA+zeroMiniBB);
ZFSTerm = zeros(9,9,4);
ZFSTerm(:,:,1) = kron(zeroMiniAA,eye(3)) + kron(eye(3),zeroMiniAA);
ZFSTerm(:,:,2) = kron(zeroMiniBB,eye(3)) + kron(eye(3),zeroMiniBB);
ZFSTerm(:,:,3) = kron(zeroMiniAA,eye(3)) + kron(eye(3),zeroMiniBB);
ZFSAVG = (1/4)*(ZFSTerm(:,:,1)+ZFSTerm(:,:,2)+2*ZFSTerm(:,:,3));
ZFSTerm(:,:,4) = ZFSAVG;

% % Compare to the 'bad' averaging technique
% zeroMiniBadAvg = 1/2*zeroMiniAA+1/2*zeroMiniBB;
% ZFSBadAvg = kron(zeroMiniBadAvg,eye(3)) + kron(eye(3),zeroMiniBadAvg);
% ZFSTerm(:,:,4) = ZFSBadAvg;

% Interaction Term 
intZMini = [0 1-3*nnZ^2 3*nnY*nnZ; -1+3*nnZ^2 0 -3*nnX*nnZ; -3*nnY*nnZ 3*nnX*nnZ 0];
intYMini = [0 3*nnY*nnZ 1-3*nnY^2; -3*nnY*nnZ 0 3*nnX*nnY; -1+3*nnY^2 -3*nnX*nnY 0];
intXMini = [0 -3*nnX*nnZ 3*nnX*nnY; 3*nnX*nnZ 0 1-3*nnX^2; -3*nnX*nnY -1+3*nnX^2 0];
intMat = intStr*[zeros(3,3) intZMini intYMini; -intZMini zeros(3,3) intXMini;...
    -intYMini -intXMini zeros(3,3)];

%% Zeeman Shift + Full Hamiltonian Assembly
% IFF angle mode is selected, we'll need to to build an array of theta and
% phi values.

% This will not work in the BC version without converting to the new
% 'cab' order

if magMode
    sPR = floor(nSteps/3);      % Steps per range
    phiBC = 90*ones(1,sPR);
    thetaBC = linspace(90,0,sPR);
    thetaCA = linspace(0,90,sPR);
    phiCA = zeros(1,sPR);
    thetaAB = 90*ones(1,sPR);
    phiAB = linspace(0,90,sPR);
    thetaArray = [thetaBC thetaCA thetaAB];
    phiArray = [phiBC phiCA phiAB];
    nSteps = length(thetaArray);
end

% The Zeeman Shift term is the only non-constant part of the hamiltonian.
% We'll loop and calculate Zeeman for each magStr / angle. Each Hamiltonian
% matrix will be a different page in a 3D array (9 x 9 x nStep) per each
% configuration (AA to DL), s.t. the final data structure is (9 x 9 x nStep
% x 4).
if magMode
    hamFull = zeros(9,9,nSteps,4);
    loopArray = linspace(-90,180,nSteps);
    for i = 1:nSteps
        zeeman = findHZee(thetaArray(i),phiArray(i),magStrMax,bohrMagneton);
        for config = 1:4
            hamFull(:,:,i,config) = ZFSTerm(:,:,config) + zeeman + intMat;
        end
    end
else
    hamFull = zeros(9,9,nSteps,4);
    loopArray = linspace(0,magStrMax,nSteps);
    for i = 1:nSteps
        zeeman = findHZee(magTheta,magPhi,loopArray(i),bohrMagneton);
        for config = 1:4
            hamFull(:,:,i,config) = ZFSTerm(:,:,config) + zeeman + intMat;
        end
    end
end

% The crux of the code is John D'Errico's eigenshuffle function, which uses
% knowledge of the previous and next page to maintain a consistent order
% of eigenvalues and eigenvectors
eigVals = zeros(9,nSteps,4);
eigVecs = zeros(9,9,nSteps,4);
for config = 1:4
    [eigVecs(:,:,:,config),eigVals(:,:,config)] = eigenshuffle(hamFull(:,:,:,config));

    % Non-eigenshuffle alternative
    % for i = 1:nSteps
    %     [eigVecs(:,:,i,config),temp] = eig(hamFull(:,:,i,config));
    %     eigVals(:,i,config) = diag(temp);
    % end
end

%% Analysis and Plotting
% First stop is the eigenvalues = energy levels
figure(1);
clf;
t1 = tiledlayout('flow');
title(t1, 'Eigenergies of 9 Triplet Pair States');
for config = 1:4
    nexttile;
    plot(loopArray,real(eigVals(:,:,config)),'Color',configRGB{config});
    title(configTitles(config));
    if magMode
        xlabel('Angle');
        xline(-90,'-','b');
        xline(0,'-','c');
        xline(90,'-','a');
        xline(180,'-','b');
    else
        xlabel('Magnetic Field Strength (Gauss)');
    end
    ylabel('Energy (GHz)');
    set(gca,'XMinorGrid','on');
end

% Next, we look at the probability to measure each of the nine abitrary
% spin states as spin 0: S = 1/sqrt(3)[|xx>+|yy>+|zz>]. First, we find
% innerProd by summing the 1st, 5th, and 9th eiegenvectors. Then we find
% the measurable probability by squaring and adjusting by 1/3.
s0Vecs = zeros(9,nSteps,4);
s0Prob = zeros(9,nSteps,4);
figure(2);
clf;
t2 = tiledlayout('flow');
title(t2, 'Probability of S=0');
for config = 1:4
    nexttile;
    for i = 1:nSteps
        s0Vecs(:,i,config) = eigVecs(1,:,i,config) + eigVecs(5,:,i,config) + eigVecs(9,:,i,config);
        % If degenFix is flipped, we look for any degenerate states (by
        % looking at energy) for each val of magnetic field. If found, the
        % singlet probabilities are overwritten with the sum and difference
        % of probabilities instead.
        if degenFix
            for j = 1:9
                for k = j+1:9
                    if ismembertol(eigVals(j,i,config),eigVals(k,i,config),1e-4)
                        temp1 = s0Vecs(j,i,config);
                        temp2 = s0Vecs(k,i,config);
                        s0Vecs(j,i,config) = (1/sqrt(2))*(temp1 + temp2);
                        s0Vecs(k,i,config) = (1/sqrt(2))*(temp1 - temp2);
                    end
                end
            end
        end
        s0Prob(:,i,config) = (1/3).*abs(s0Vecs(:,i,config)).^2;
    end
    
        % AB Fix: There's something weird happening with the AB singlet states
        % at extremely low field (under xGauss). Until this is fixed/understood,
        % the s0Vecs and s0Prob for the AB config will be overwritten s.t.
        % 0 to xG = xG.
        xGauss = 1;
        loopInd = find(loopArray>xGauss,1);
        s0VecsInsert = zeros(9,loopInd) + s0Vecs(:,loopInd,3);
        s0ProbInsert = zeros(9,loopInd) + s0Prob(:,loopInd,3);
        s0Vecs(:,1:loopInd,3) = s0VecsInsert;
        s0Prob(:,1:loopInd,3) = s0ProbInsert;
    
        plot(loopArray,s0Prob(:,:,config),'Color',configRGB{config});
        title(configTitles(config));
        if magMode
            xlabel('Angle');
            xline(-90,'-','b');
            xline(0,'-','c');
            xline(90,'-','a');
            xline(180,'-','b');
        else
            xlabel('Magnetic Field Strength (Gauss)');
        end
        ylabel('P(S=0))');
        set(gca,'XMinorGrid','on','YMinorGrid','on');
end

% Figure plotting for ZFTID Paper (1/6/25)
figure(10);
t10 = tiledlayout('vertical');
nexttile;
plot(loopArray,s0Prob(:,:,1),'Color',configRGB{1});
set(gca,'XMinorGrid','on','YMinorGrid','on');
xlabel('Magnetic Field Strength (Gauss)');
ylabel('Singlet Probability');
ylim([0 0.7]);
nexttile;
plot(loopArray,s0Prob(:,:,3),'Color',configRGB{2});
set(gca,'XMinorGrid','on','YMinorGrid','on');
xlabel('Magnetic Field Strength (Gauss)');
ylabel('Singlet Probability');
ylim([0 0.7]);

% figure(11);
% plot(loopArray,s0Prob(:,:,3),'Color',configRGB{3});
% set(gca,'XMinorGrid','on','YMinorGrid','on');
% xlabel('Magnetic Field Strength (Gauss)');
% ylabel('Singlet Probability');

% Next Step: The spin-wavefunction in a non-specific basis is the LC
% C1|1> + C2|2> + ... C9|9>. These coefficients can be found by
% constraining the wavefunction at t=0 to be the S=0 wavefunction, which
% basically makes this a 9x9 system of equations.
% INEFF: Could this be done with a page-wise division, removing the i loop?
psiCoeff = zeros(9,nSteps,4);
s0Vector = [1/sqrt(3) 0 0 0 1/sqrt(3) 0 0 0 1/sqrt(3)];
% Equipartition vector
ePVector = 1/sqrt(9)*eye(1,9);
for config = 1:4
    for i = 1:nSteps
        psiCoeff(:,i,config) = eigVecs(:,:,i,config)\s0Vector';
        % psiCoeff(:,i,config) = eigVecs(:,:,i,config)\ePVector';
        psiCoeff(:,i,config) = psiCoeff(:,i,config)/norm(psiCoeff(:,i,config));
    end
end

% Now that we have the coefficients solved for, we carry out the entire
% time evolution and projection onto the singlet state in the generalized
% basis. There are two 'kinds' of terms to come out of the algebra: 9 'like
% like' terms, and 36 'like unlike' cross terms. The 9 LL terms will have a
% cancelled exponential, and their coeff*vecs will give a 'boost' in
% fluorescence, averaged over t >> beat period.
boostPred = zeros(nSteps,4);
figure(3);
clf;
t3 = tiledlayout('flow');
title(t3, 'Predicted Boost per Config');
for config = 1:4
    for i = 1:nSteps
        for j = 1:9
            % boostPred(i,config) = boostPred(i,config) +...
            % (1/3)*abs(psiCoeff(j,i,config))^2*abs(s0Vecs(j,i,config))^2;
            boostPred(i,config) = boostPred(i,config) + ...
            abs(psiCoeff(j,i,config))^4;
        end
    end
    boppo = boostPred(1,config);
    boostPred(:,config) = boostPred(:,config) / boppo;
    nexttile;
    plot(loopArray,boostPred(:,config),'Color',configRGB{config});
    title(configTitles(config));
    if magMode
        xlabel('Angle');
        xline(-90,'-','b');
        xline(0,'-','c');
        xline(90,'-','a');
        xline(180,'-','b');
    else
        xlabel('Magnetic Field Strength (Gauss)');
    end
    ylabel('Predicted Boost');
    set(gca,'XMinorGrid','on','YMinorGrid','on');
end

% We're also interested in comparing the observed boost effect to the
% effect of heterofusion, or trapping-driven triplet fusion in Tetracene
% by Merrifield 71. Here, the coefficients still matter, but in a
% completely different algebraic form. Note that this assumes that the
% rates of fission, fusion (outside the spin character), and ballistic
% collision are independent of field.
kM = 2.8e9;
kS = 1.1e9;
kT = 1.7e9;

heteroPred = zeros(nSteps,4);
figure(31);
clf;
t3 = tiledlayout('flow');
title(t3, 'Merrifield Heterofusion Rate per Config');
for config = 1:4
    for i = 1:nSteps
        for j = 1:9
            heteroPred(i,config) = heteroPred(i,config) + ...
            kS*abs(psiCoeff(j,i,config))^2/((kM+(1/2)*kT)+(kS-(3/2)*kT)*abs(psiCoeff(j,i,config))^2);
        end
        heteroPred(i,config) = (1/9)*heteroPred(i,config);
    end
    bippo = heteroPred(1,config);
    heteroPred(:,config) = heteroPred(:,config) / bippo;
    nexttile;
    plot(loopArray,heteroPred(:,config),'Color',configRGB{config});
    title(configTitles(config));
    if magMode
        xlabel('Angle');
        xline(-90,'-','b');
        xline(0,'-','c');
        xline(90,'-','a');
        xline(180,'-','b');
    else
        xlabel('Magnetic Field Strength (Gauss)');
    end
    ylabel('Predicted Heterofusion Rate');
    set(gca,'XMinorGrid','on','YMinorGrid','on');
end


% Point Storage for Pred Boost Comparison (optional)
% 10 deg B to A
% xExpBoost = [200 400 600 800 1000 1200 1400 1600 1800];
% yExpBoost = [0.9212 0.9322 0.9453 0.9597 0.9721 1.1003 1.0994 1.0889 1.0955];
% 10 deg B to C
% xExpBoost = [200 400 600 800 1000 1200 1400 1600 1800];
% yExpBoost = [0.9050 0.9161 0.9434 0.9597 0.9868 1.0244 1.0385 1.0621 1.0714];
% expBoostPoints = true;
% Xallan Angle, 0-40ns
% xExpBoost = [200 400 600 800 1000 1250 1500 1750 2000 2500 3000];
% yExpBoost = [0.9119 0.9261 0.9262 0.9875 0.9728 1.0623 1.0812 1.1667 1.1826 1.2591 1.4213];
% expBoostPoints = true;
% Xallan Angle, 0-Crossing Point
% xExpBoost = [800 1000 1250 1500 1750 2000 2500 3000];
% yExpBoost = [1.0492 1.0318 1.1263 1.1520 1.2327 1.2509 1.3286 1.5147];
% expBoostPoints = true;
% Near Zero, mag||b
% xExpBoost = [68 95 115 138 158 177 197 218 240 258 276 300 350 400 450 500 550 600];
% yExpBoost = [0.9213 0.8879 0.8763 0.9010 0.9200 0.8930 0.9047 0.9107 0.9069... 
%     0.9036 0.9169 0.8377 0.8473 0.8231 0.9031 0.9277 0.9409 0.9698];
% expBoostPoints = true;

% Point Storage for Pred Boost Comparison (optional) (BBA2024 - Normalized)
% 15 deg off B, BC Plane
% xExpBoost = [140 180 220 240 260 300 340];
% yExpBoost = [0.9087 0.9004 0.9152 0.9089 0.9120 0.9175 0.9063];
% expBoostPoints = true;
% 30 deg off B, BA Plane
% xExpBoost = [140 180 220 240 260 300 340];
% yExpBoost = [0.9087 0.9004 0.9001 0.8996 0.8961 0.8961 0.9169];
% expBoostPoints = true;
% 10 deg off B, BA Plane
% xExpBoost = [200 400 600 800 1000 1200 1400 1600 1800];
% yExpBoost = [0.9458 0.9549 0.9630 0.9850 0.9869 1.0152 1.0176 1.0581 1.0828];
% expBoostPoints = true;
% 10 deg off B, BC Plane
% xExpBoost = [200 400 600 800 1000 1200 1400 1600 1800];
% yExpBoost = [0.9335 0.9440 0.9703 0.9852 1.0108 0.9986 1.0281 1.0351 1.0304];
% expBoostPoints = true;
% 34.5 deg off C, CA Plane (Xallan Angle)
% xExpBoost = [200 400 600 800 1000 1250 1500 1750 2000 2500 3000];
% yExpBoost = [0.8971 0.9097 0.9396 0.9876 1.0155 1.0811 1.0988 1.1803 1.2117 1.2943 1.2931];
% expBoostPoints = true;
% Mag||b
% xExpBoost = [68 95 115 138 158 177 197 218 240 258 276];
% yExpBoost = [0.9234 0.9281 0.9187 0.8858 0.8998 0.8952 0.8714 0.8958 0.8936 0.8856 0.8744];
% expBoostPoints = true;

% Point Storage for Pred Boost Comparison Sept 2024
% Solenoid, Mag||c
% xExpBoost = [25 50 75 100 125 150 175 200 225 250 275 300 325 350 375];
% yExpBoost = [0.986 0.9822 0.9654 0.9841 0.9790 0.9759 0.9697 0.9442 0.9677 0.9682 0.9533 0.9551 0.9667 0.9718 0.9720];
% expBoostPoints = true;

% % Mag||b, Max Boost Ratio Comparison
% xExpBoost = [1000 2000 3000 4000];
% yExpBoost = [1.1665 1.3604 1.5162 1.5899];
% expBoostPoints = true;
% % Mag||b, Sim Boost Comparison (Full Sim - No Dephasing)
% xSimBoost = [50 100 150 200 250 300 350 400 450 500 600 700 800 900 1000 ...
%     1200 1400 1600 1800 2000 2500];
% ySimBoost = [0.9574 0.8678 0.7981 0.7767 0.7954 0.8645 0.9412 1.0191 1.0919 1.1569 ...
%     1.2661 1.3476 1.4090 1.4552 1.4907 1.5402 1.5719 1.5933 1.6083 1.6192 1.6364];
% simBoostPoints = true;

% Mag||b, Max Boost Ratio Comparison - 1 ns (Saimonth Data)
% xExpBoost = [10 20 30 40 50 60 70 80 90 100 110 120 130 140 200 300 400 500 ... 
%     600 700 800 900 1000 1200 1400 1600 1800 2000 2250 2500 2750 3000 3250 3500 3750 4000];
% yExpBoost = [0.9938 0.9872 0.9888 0.9691 0.9533 0.9502 0.9226 0.9179 0.9096 0.8877 ... 
%     0.8877 0.8824 0.8666 0.8728 0.7705 0.7805 0.8118 0.8375 0.8820 0.9103 0.9535 0.9881 0.9915 ...
%     1.0602 1.1240 1.1688 1.1987 1.2314 1.2557 1.2861 1.2983 1.3139 1.3354  1.3359 1.3925 1.381];
% expBoostPoints = true;
% % Mag||b, Max Boost Ratio Comparison - 1 ns to MAX
% xSimBoost = [100 200 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000 ...
%     2250 2500 2750 3000 3250 3500 3750 4000];
% ySimBoost = [0.8187 0.7705 0.7805 0.8118 0.8375 0.8820 0.9103 0.9535 1.0384 1.033 ...
%     1.0926 1.1554 1.2002 1.2227 1.2532 1.2785 1.3183 1.3419 1.3631 1.3976 ...
%     1.4027 1.4687 1.4641];
% simBoostPoints = true;

% Mag||c, Max Boost Ratio Comparison
% xExpBoost = [500 1500 2000 2500 3000];
% yExpBoost = [0.9944 1.3579 1.4644 1.5232 1.5430];
% expBoostPoints = true;
% Mag along Xallan Angle (35degfromC,CA), Max Boost Ratio
% xExpBoost = [200 400 600 800 1000 1250 1500 1750 2000 2500 3000];
% yExpBoost = [0.9897 1.0065 1.0147 1.0605 1.0978 1.1699 1.1965 1.2794 1.3172 1.4115 1.4182];
% expBoostPoints = true;
% WTF is this
% xExpBoost = [500 1000 2300 3600];
% yExpBoost = [0.9896 1.2641 1.6629 2.7685];
% expBoostPoints = true;
% Mag||a, Max Boost Ratio Comparison
% xExpBoost = [250 500 1000 1500 2000 3000 4100];
% yExpBoost = [1.1766 1.3832 1.2082 1.0345 0.98122 1.1596 2.0132];
% expBoostPoints = true;
% Mag||a, Max Boost Ratio Comparison (Take two!)
% xExpBoost = [250 500 1000 2000 4000];
% yExpBoost = [0.9487 0.9313 0.9575 1.1272 1.3192];
% expBoostPoints = true;

% Zach Points, Mag || C
% xExpBoost = [100 200 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000 2250 2500 2750 3000];
% yExpBoost = [0.9689 1.0424 1.0367 0.9616 0.9366 1.0091 1.0573 1.0353 1.064 1.089 1.0982 1.1276 1.1299 1.1613 1.1225 1.1517 1.1731 1.2164 1.2306];
% expBoostPoints = true;

% % 4-13-25, Mag||b, Boost at 1ns
% xExpBoost = [200 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000];
% yExpBoost = [0.8100 0.8044 0.8384 0.8849 0.9038 0.9367 0.9652 1.0207 1.0379 1.1117 1.167 1.2093 ...
%     1.2298 1.2603 1.3437 1.3933 1.4091 1.4607];
% expBoostPoints = true;
% % 4-13-25, Mag||b, Boost at 5ns
% xExpBoost = [200 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000];
% yExpBoost = [0.7461 0.7418 0.7715 0.7976 0.8097 0.8423 0.8675 0.9120 0.9294 0.9965 1.0585 1.1127 ...
%     1.1452 1.1940 1.3135 1.3939 1.4434 1.5112];
% expBoostPoints = true;
% 4-13-25, Mag||b, MAX Boost
% xExpBoost = [200 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000];
% yExpBoost = [0.7234 0.7177 0.7480 0.7603 0.7708 0.7922 1.0077 1.0762 1.0751 1.1436 1.1839 1.2164 ...
%     1.2385 1.2609 1.3706 1.4325 1.4635 1.5243];
% expBoostPoints = true;

% % Data Source: 2025-05-09 for Mag||c (Unicorn Crystal)
% xExpBoost = [20 40 60 80 100 120 140 160 180 200 250 300 350 400 450 500 600 700 800 ...
%     900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500 5000];
% yExpBoost = [0.9986 0.9882 0.9991 0.9702 0.9780 0.9852 0.9545 0.9429 0.9507 ...
%     0.9206 0.9123 0.8938 0.8954 0.9022 0.9021 0.9178 0.9437 0.9654 1.0174 1.0447 ... 
%     1.0784 1.1098 1.1355 1.1911 1.2024 1.2014 1.2329 1.2609 1.2717 1.2766 1.2789 1.3008];
% expBoostPoints = true;

% If we assume the TID assumption, that the AB config is twice as likely as
% the AA and BB configs, we can throw together a naive predicted boost
% for a real rubrene crystal.
expBoost = (ABBoost)*boostPred(:,3)+((1-ABBoost)/2)*boostPred(:,1)+((1-ABBoost)/2)*boostPred(:,2);
% expBoost = (11/40)*boostPred(:,1) + (11/40)*boostPred(:,2) + (22/40)*boostPred(:,3);
% expBoost = boostPred(:,4);
figure(4);
clf;
hold on;
plot(loopArray,expBoost(:)/expBoost(1),'Color','black');
if expBoostPoints
    scatter(xExpBoost,yExpBoost,'filled');
end
if simBoostPoints
    scatter(xSimBoost,ySimBoost,'filled');
end
yline(1,"--");
hold off;
title('Predicted Boost - Crystal Average');
if magMode
    xlabel('Angle');
    xline(-90,'-','b');
    xline(0,'-','c');
    xline(90,'-','a');
    xline(180,'-','b');
else
    xlabel('Magnetic Field Strength (Gauss)');
end
ylabel('Predicted Boost Ratio');
set(gca,'XMinorGrid','on','YMinorGrid','on');

% One more observable! The 36 cross terms have an exponential that doesn't
% cancel, e^(E1-E2). Those differences in energy will become quantum beats.
% The amplitude of those terms, (2/3)*C1*C2*v1*v2 will be the relative
% amplitude.
if predictMode
    freqArray = zeros(36,nSteps,4);
    ampArray = zeros(36,nSteps,4);
    for config = 1:4
        counter = 0;
        for j = 1:9 
            for k = (j+1):9
                counter = counter + 1;
                for i = 1:nSteps
                    freqArray(counter,i,config) = eigVals(j,i,config)-eigVals(k,i,config);
                    % ampArray(counter,i,config) = (2/3)*abs(psiCoeff(j,i,config)* ...
                    %     psiCoeff(k,i,config)*s0Vecs(j,i,config)*s0Vecs(k,i,config));
                    ampArray(counter,i,config) = abs(2*abs(psiCoeff(j,i,config))^2*abs(psiCoeff(k,i,config))^2);
                end
            end
        end
    end
    
    % This cutoff function kinda solves the issue in the next comment
    % without having to sort the array; simply don't plot points with
    % amplitudes below a certain threshold.
    ampArray(ampArray < 0.00001) = NaN;
    
    % TODO: QB Prediction is misleading with that whole 'marker size'
    % thing. Could be an issue with sorting.
    figure(5);
    clf;
    t5 = tiledlayout('flow');
    title(t5, 'Quantum Beat Predictions')
    for config = 1:4
        nexttile;
        hold on;
        for i = 1:counter
            scatter(loopArray,freqArray(i,:,config),(ampArray(i,:,config)+0.0001).^2*50,configRGB{config},'filled');
            % scatter(loopArray,real(freqArray(i,:,config)),ampArray(i,:,config).^2*500,configRGB{config},'filled');
        end
        hold off;
        title(configTitles(config));
        if magMode
            xlabel('Angle');
            xline(-90,'-','b');
            xline(0,'-','c');
            xline(90,'-','a');
            xline(180,'-','b');
        else
            xlabel('Magnetic Field Strength (Gauss)');
        end
        ylabel('Predicted Beat Frequency (GHz)');
        set(gca,'XMinorGrid','on','YMinorGrid','on');
    end
end

%% Functions
% Calculates Zeeman Shift given angle, str, and the Bohr magneton
function zeemanFull = findHZee(magTheta,magPhi,magStr,bohrMag)
    % Calculate Magnetic Field Components
    magX = magStr*sind(magTheta)*cosd(magPhi);
    magY = magStr*sind(magTheta)*sind(magPhi);
    magZ = magStr*cosd(magTheta);

    % Construct 3x3 mini matrices for each component
    magXMini = magX*eye(3);
    magYMini = magY*eye(3);
    magZMini = magZ*eye(3);

    % Construct 3x3 mini matrix for main diagonal
    magDiagMini = [0 -magZ magY; magZ 0 -magX; -magY magX 0];

    % Now construct full hZeeman term according to Tapping layout
    zeemanFull = 1i*bohrMag*[magDiagMini -magZMini magYMini;... 
        magZMini magDiagMini -magXMini; -magYMini magXMini magDiagMini];
end

%% Notes and References

% Common EPR Parameters
% Rubrene (Biaggio Fit)
% dZFS = 1.6564;
% eZFS = -0.11872;
% Tetracene
% dZFS = 1.5589;
% eZFS = -0.15589;
% Anthracene
% dZFS = 2.0805;
% eZFS = -0.250624;
% Napthalene
% dZFS = 2.89897;
% eZFS = -0.47667;

% Tait-Bryan Angles
% Rubrene       (0,31,0)
% Tetracene     (26,0,90)   Check?
% Anthracene    (0,125,0)   Check?
