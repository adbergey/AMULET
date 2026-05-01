% ADB 2/20/2026

% This script will generate all of the experimental result plots from 
% AMULET: Acoustic Metastructure for Direction-of-Arrival Estimation 
% Underwater using a Single Hydrophone

% Ensure that AMULET_Data.mat is in the same folder as this script

clc
clear
close all
tic
%% Load in the data

load("AMULET_Data.mat")

% Setup 
chirp_duration = 0.1;
fs =  192000;
f_start = 1000; %[Hz]
f_stop = 88000; %[Hz]
t = (0:1/fs:chirp_duration-1/fs)';
TX_chirp = chirp(t,f_start,chirp_duration,f_stop);
angles = 0:359;
numAngles = length(angles);

%% Figure 2(a-c)

%Initialize the matching score matrices
match_scores_a = zeros(numAngles);
match_scores_b = zeros(numAngles);
match_scores_c = zeros(numAngles);

%Compute how similar each unknown signature is to each reference signature
for u = 1 : numAngles
    for c = 1 : numAngles
        match_scores_a(u,c) = max(abs(xcorr(amulet_Data.noMS1(:,c),amulet_Data.noMS2(:,u))));
        match_scores_b(u,c) = max(abs(xcorr(amulet_Data.owlet1(:,c),amulet_Data.owlet2(:,u))));
        match_scores_c(u,c) = max(abs(xcorr(amulet_Data.amuletPressFit1_1m(:,c),amulet_Data.amuletPressFit2_1m(:,u))));
    end
end

%Plot the result
figure
subplot(1,3,1)
imagesc(match_scores_a)
title('No Metastructure')
xlabel('Reference Angle (Degrees)')
ylabel('Unknown Angle (Degrees)')
clim([0 1])

subplot(1,3,2)
imagesc(match_scores_b)
title('Owlet-type Structure')
xlabel('Reference Angle (Degrees)')
ylabel('Unknown Angle (Degrees)')
clim([0 1])

subplot(1,3,3)
imagesc(match_scores_c)
title('AMULET-type')
xlabel('Reference Angle (Degrees)')
ylabel('Unknown Angle (Degrees)')
colorbar
clim([0 1])

%% Figure 8

figure
plot(abs(amulet_Data.pier_55cm_pressFit))
hold on
plot(abs(amulet_Data.pier_1m_pressFit))
legend('55 cm Separation', '1 m Separation')
title('Signatures are Distance Agnostic')


%% Figure 9

comp_mat_noMS = zeros(numAngles);
comp_mat_owlet = zeros(numAngles);
comp_mat_amulet = zeros(numAngles);


%Compute the similarity score of each signature to its neighbors
for i = 1:numAngles
    for j = i:numAngles

        temp_cor = abs(xcorr(amulet_Data.noMS1(:,i),amulet_Data.noMS1(:,j)));
        comp_mat_noMS(i,j) = max(temp_cor);
        comp_mat_noMS(j,i) = comp_mat_noMS(i,j);

        temp_cor = abs(xcorr(amulet_Data.owlet1(:,i),amulet_Data.owlet1(:,j)));
        comp_mat_owlet(i,j) = max(temp_cor);
        comp_mat_owlet(j,i) = comp_mat_owlet(i,j);

        temp_cor = abs(xcorr(amulet_Data.amuletPressFit_55cm(:,i),amulet_Data.amuletPressFit_55cm(:,j)));
        comp_mat_amulet(i,j) = max(temp_cor);
        comp_mat_amulet(j,i) = comp_mat_amulet(i,j);

    end
end

%Variable setup
diag_width = 10; %How many elements each direction to extend beyond the MD

goodness_score_noMS = zeros(numAngles,1);
goodness_score_owlet = zeros(numAngles,1);
goodness_score_amulet = zeros(numAngles,1);

%Excluding the immediate neighbors up to diag_width, compute 
for i = 1:numAngles
    excludeInds = mod((i-diag_width-1:i+diag_width-1)+360,360)+1;
    inds = angles+1;
    inds(excludeInds) = [];
    goodness_score_noMS(i) = comp_mat_noMS(i,i)/mean(comp_mat_noMS(i,inds));
    goodness_score_owlet(i) = comp_mat_owlet(i,i)/mean(comp_mat_owlet(i,inds));
    goodness_score_amulet(i) = comp_mat_amulet(i,i)/mean(comp_mat_amulet(i,inds));
end

%Convert to dB
db_score_noMS = 20*log10(goodness_score_noMS);
db_score_owlet = 20*log10(goodness_score_owlet);
db_score_amulet = 20*log10(goodness_score_amulet);

%Plot the result
figure
plot(db_score_noMS)
hold on
plot(db_score_owlet)
plot(db_score_amulet)
title('Goodness Metric')
ylabel('Goodness (dB)')
xlabel('Angle (degrees)')
legend('No MS','Owlet','AMULET')

%Print the average goodness score
disp("No Metastructure Average Goodness Score of " + mean(db_score_noMS) + "dB")
disp("Owlet-type Average Goodness Score of " + mean(db_score_owlet) + "dB")
disp("AMULET-type Average Goodness Score of " + mean(db_score_amulet) + "dB")


%% Figure 11(a) Design Variation

% NOTE on Figure 11. To save time and line of code, the final numbers are
% precomputed. However the raw data is all included and an example of how
% the numbers were computed is included.







% Load in the precomputed results
results = amulet_Data.design_variation;
avg = results(1,:);
upper95 = results(2,:);
lower95 = results(3,:);
goodness = results(4,:);
x = ["No MS", "Owlet", "Water Filled Slit-Swirl","Water Filled Swirl", "Air Filled Slit Swirl", "Amulet"]; 


figure
subplot(2,1,1)
bar(x,goodness)
ylabel("Goodness (dB)")
set(gca, 'FontName', 'Arial', 'FontSize',20)


subplot(2,1,2)
upperOffset = upper95-avg;
lowerOffset = avg-lower95;
yMax = round(max(upper95)*1.3,1);
bar(x,avg)
hold on
er = errorbar(1:length(x),avg,lowerOffset,upperOffset,'LineWidth',2,'CapSize',10);
er.LineStyle = 'none';
er.Color = '[0 0 0]';
ylabel("Error (Degrees)")
ylim([0,yMax])
ytix = get(gca, 'YTick');
set(gca, 'FontName', 'Arial', 'FontSize',20)


%% Figure 11(b) Diameter Variation

% Load in the precomputed results
results = amulet_Data.diameter_variation;
avg = results(1,:);
upper95 = results(2,:);
lower95 = results(3,:);
goodness = results(4,:);
x = ["4.5","6.2","7.7"]; 


figure
subplot(2,1,1)
bar(x,goodness)
ylabel("Goodness (dB)")
set(gca, 'FontName', 'Arial', 'FontSize',20)


subplot(2,1,2)
upperOffset = upper95-avg;
lowerOffset = avg-lower95;
yMax = round(max(upper95)*1.3,1);
bar(x,avg)
hold on
er = errorbar(1:length(x),avg,lowerOffset,upperOffset,'LineWidth',2,'CapSize',10);
er.LineStyle = 'none';
er.Color = '[0 0 0]';
ylabel("Error (Degrees)")
ylim([0,yMax])
ytix = get(gca, 'YTick');
set(gca, 'FontName', 'Arial', 'FontSize',20)
xlabel("Structure Diameter (cm)")

%% Figure 11(c) Material Variation

% Load in the precomputed results
results = amulet_Data.material_variation;
avg = results(1,:);
upper95 = results(2,:);
lower95 = results(3,:);
goodness = results(4,:);
x = ["Resin","PLA","Seal-coated PLA"];


figure
subplot(2,1,1)
bar(x,goodness)
ylabel("Goodness (dB)")
set(gca, 'FontName', 'Arial', 'FontSize',20)


subplot(2,1,2)
upperOffset = upper95-avg;
lowerOffset = avg-lower95;
yMax = round(max(upper95)*1.3,1);
bar(x,avg)
hold on
er = errorbar(1:length(x),avg,lowerOffset,upperOffset,'LineWidth',2,'CapSize',10);
er.LineStyle = 'none';
er.Color = '[0 0 0]';
ylabel("Error (Degrees)")
ylim([0,yMax])
ytix = get(gca, 'YTick');
set(gca, 'FontName', 'Arial', 'FontSize',20)


%% Figure 12

%Initialize the matching score matrices
match_scores_a = zeros(numAngles);
match_scores_b = zeros(numAngles);
match_scores_c = zeros(numAngles);
match_scores_d = zeros(numAngles);
match_scores_e = zeros(numAngles);
match_scores_f = zeros(numAngles);

%Compute how similar each unknown signature is to each reference signature
for u = 1 : numAngles
    for c = 1 : numAngles
        match_scores_a(u,c) = max(abs(xcorr(amulet_Data.amulet_desktop1(:,c),amulet_Data.amulet_desktop2(:,u))));
        match_scores_b(u,c) = max(abs(xcorr(amulet_Data.amulet_saltwater1(:,c),amulet_Data.amulet_saltwater2(:,u))));
        match_scores_c(u,c) = max(abs(xcorr(amulet_Data.coated_amulet_pier1(:,c),amulet_Data.coated_amulet_pier2(:,u))));
        match_scores_d(u,c) = max(abs(xcorr(amulet_Data.amuletPressFit1_1m(:,c),amulet_Data.amuletPressFit_55cm(:,u))));
        match_scores_e(u,c) = max(abs(xcorr(amulet_Data.amulet_pier1(:,c),amulet_Data.coated_amulet_pier2(:,u))));
        match_scores_f(u,c) = max(abs(xcorr(amulet_Data.amulet_saltwater1(:,c),amulet_Data.coated_amulet_pier2(:,u))));
    end
end

%Plot the result
figure
subplot(2,3,1)
imagesc(match_scores_a)
title('Desktop Tank')
xlabel('Reference Angle (Degrees)')
ylabel('Unknown Angle (Degrees)')
clim([0 1])

subplot(2,3,2)
imagesc(match_scores_b)
title('Saltwater Tank')
xlabel('Reference Angle (Degrees)')
ylabel('Unknown Angle (Degrees)')
clim([0 1])

subplot(2,3,3)
imagesc(match_scores_c)
title('Lake Pier')
xlabel('Reference Angle (Degrees)')
ylabel('Unknown Angle (Degrees)')
colorbar
clim([0 1])

subplot(2,3,4)
imagesc(match_scores_d)
title('Across Distance')
xlabel('0.5m Separation (Degrees)')
ylabel('1m Separation (Degrees)')
colorbar
clim([0 1])

subplot(2,3,5)
imagesc(match_scores_e)
title('Across Days')
xlabel('Day 1 (Degrees)')
ylabel('Day 2 (Degrees)')
colorbar
clim([0 1])

subplot(2,3,6)
imagesc(match_scores_f)
title('Across Testbeds')
xlabel('Saltwater Tank (Degrees)')
ylabel('Lake Pier (Degrees)')
colorbar
clim([0 1])


%% Figure 13(a)
% NOTE: For simplicity all of the errors have already been computed for
% this example. To see how they were computed, refer to the next section. 

%Setup
numPaths = 10;
errorBins = zeros(181,1);
allErrors = [];

figure
subplot(1,4,1:3)
hold on

for i = 1 : numPaths
    errors = amulet_Data.("saltwater_tracking_errors"+num2str(i));

    %Remove outliers
    p = prctile(abs(errors),[5,95],"all");
    errors(abs(errors) > p(2)) = [];
    
    %Create a violin plot to show the error distribution
    violinplot(i,abs(errors))
    medError = median(abs(errors));

    e = errorbar(i,min(abs(errors)),0,medError,'LineWidth',3,'CapSize',30);
    e.Color = [0,0,0];
    e = errorbar(i,min(abs(errors)),0,max(abs(errors)),'LineWidth',3,'CapSize',30);
    e.Color = [0,0,0];

    counts = histc(abs(errors),0:180);
    errorBins = errorBins+counts;

    %Append the errors to the master list for the cdf
    allErrors = [allErrors; abs(errors)];

end

xlabel('Path #')
ylabel('Error (degrees)')

yMax = 10;
ylim([0,yMax])


%Compute CDF
totalPts = sum(errorBins);
cdf = cumsum(errorBins)/totalPts;
maxInd = yMax;
[F,x] = ecdf(allErrors);

%Plot CDF
subplot(1,4,4)
plot(F,x)
xlabel('CDF')
xlim([0,1])


%% Figure 13(b-c)

%Define trajectory
traj_num = 1; % Change this to any number 1-10 

% (Note that there is a small timing offset between the start of the motion
% and the start of the recording for some trajectories. This can be
% compensated for via the initial grndTrthT value below)

trackSig = amulet_Data.("saltwater_tracking"+num2str(traj_num));
numChirps = amulet_Data.tracking_numChirps(traj_num);
pathTraveled = amulet_Data.("saltwater_tracking_path"+num2str(traj_num));

%Define the calibration signatures
ref_sigs = amulet_Data.amulet_saltwater1;

%Setup
txNumSamps = 115200;
chirpsPerSec = fs/txNumSamps; 
degPerS = 1.48;

keepLowerLimit = 25;
keepUpperLimit = 125;
corWidth = keepUpperLimit+keepLowerLimit+1;


% Identify the first correlation peak in the recording (used to align all
% subsequent chirps)
[cor, lags] = xcorr(trackSig(1:fs),TX_chirp);
[pks, locs] = findpeaks(abs(cor(floor(length(cor)/2):end)), 'MinPeakDistance', chirp_duration*3*fs*.8);

corrOffset = 500;  %To ensure the entire peak is included)

%Isolate the first received chirp with some buffer
firstCor = abs(cor(locs(2)-corrOffset+floor(length(cor)/2):locs(2)+2*corrOffset+floor(length(cor)/2)));

%Narrow in on the received chirp and find its starting index
[~,locs1] = findpeaks(abs(firstCor),'MinPeakDistance', 60,'MinPeakHeight',6);
firstCor = firstCor(locs1(1)-keepLowerLimit:locs1(1)+keepUpperLimit);
startInd = lags(locs(2)-corrOffset+locs1(1)-keepLowerLimit+floor(length(cor)/2));

%Define the indexes of the regularly spaced received impusle responses
testChunk = trackSig(startInd:startInd+chirp_duration*fs-1);
temp = xcorr(testChunk,TX_chirp);
[pks1,locs1] = findpeaks(abs(temp),'MinPeakDistance', 60, 'MinPeakHeight',6);
corKeepRange = locs1(1)-keepLowerLimit:locs1(1)+keepUpperLimit;

%Setup variables
matchScores = zeros(numChirps,length(angles));
predictions = zeros(numChirps,1);

%Compare with each point
for u = 1 : numChirps

    %Use the known chirp to get the next test signature
    testChirp = trackSig(startInd+(u-1)*txNumSamps:startInd+(u-1)*txNumSamps+fs*chirp_duration-1);
    temp = xcorr(testChirp,TX_chirp); 
    testIR = temp(corKeepRange);
    normTestIR = testIR/sqrt(sum(testIR.^2));

    %Match the test signature to the reference signatures
    for c = 1 : length(angles)
        matchScores(u,c) = max(abs(xcorr(ref_sigs(:,c),normTestIR)));
    end

    %Predict the position based on the maximum similarity
    [~, closeInd] = max(matchScores(u,:));
    prediction = closeInd-1;
    predictions(u) = prediction;

end

%Shift early negative predictions to positive
if predictions(1) > 300
    predictions(1) = predictions(1)-360;
end

% OPTIONAL: Ignore discontinuities by duplicating current point
jumpthresh = 75;
for i = 1:numChirps-1
    if abs(wrapTo180(predictions(i)-predictions(i+1))) > jumpthresh
        %predictions(i+1) = predictions(i);
    end
end


%Account for wrap around errors
predictions = rad2deg(unwrap(deg2rad(predictions)));

%Plot the experimental results
trackingT = 0:1/chirpsPerSec:numChirps/chirpsPerSec-1/chirpsPerSec;

figure
imagesc(flipud(matchScores'))
title('Matching scores at each point along the trajectory')
colorbar
clim([0,1])

figure
plot(trackingT,predictions)
title('Angular Trajectory')
hold on


% Compute the ground truth trajectory
trackingT = 0:1/chirpsPerSec:numChirps/chirpsPerSec - 1/chirpsPerSec;

grndTrthT = 5; %This initial point is equal to the time between starting the recording script and the movement script
grndTrth = pathTraveled(1);
prevTime = grndTrthT;

for i = 2 : length(pathTraveled)
    grndTrth = [grndTrth, pathTraveled(i)];
    grndTrthT = [grndTrthT, prevTime + abs(pathTraveled(i)-pathTraveled(i-1))/degPerS];
    prevTime = prevTime + abs(pathTraveled(i)-pathTraveled(i-1))/degPerS;

end

plot(grndTrthT,grndTrth)
legend('AMULET','Trajectory')
gndTrthPts = interp1(grndTrthT,grndTrth,trackingT);

%Compute the error for each prediction
errors = predictions-gndTrthPts';

%Account for initial offset errors
errors(isnan(errors)) = [];

disp("The average angular error is " + num2str(round(mean(abs(errors)),2)) + " degrees")



%% Figure 14(b)

%Define trajectory
trackSig = amulet_Data.pier_8mtracking;
numChirps = amulet_Data.tracking_numChirps(11);
pathTraveled = amulet_Data.saltwater_tracking_path5;

%Define the calibration signatures
ref_sigs = amulet_Data.coated_amulet_pier1;

%Setup
chirpsPerSec = fs/txNumSamps; %Ensure these numbers are accurate  %4.7OG
degPerS = 1.41;

keepLowerLimit = 25;
keepUpperLimit = 125;
corWidth = keepUpperLimit+keepLowerLimit+1;


% Identify the first correlation peak in the recording (used to align all
% subsequent chirps)
[cor, lags] = xcorr(trackSig(1:fs),TX_chirp);
[pks, locs] = findpeaks(abs(cor(floor(length(cor)/2):end)), 'MinPeakDistance', chirp_duration*3*fs*.8);


corrOffset = 500;  %To ensure the entire peak is included


%Manually set it, at least for initial testing
locs = 267165;
firstCor = abs(cor(locs-corrOffset:locs+2*corrOffset));


[pks1,locs1] = findpeaks(abs(firstCor),'MinPeakDistance', 60,'MinPeakHeight',0.5);


firstCor = firstCor(locs1(1)-keepLowerLimit:locs1(1)+keepUpperLimit);

startInd = lags(locs(1)-corrOffset+locs1(1)-keepLowerLimit);

testChunk = trackSig(startInd:startInd+chirp_duration*fs-1);
temp = xcorr(testChunk,TX_chirp);
[pks1,locs1] = findpeaks(abs(temp),'MinPeakDistance', 60, 'MinPeakHeight',0.5);

%corKeepRange = locs1(1)-keepLowerLimit:locs1(1)+keepUpperLimit;
corKeepRange = 19200:19350;

%Setup variables
matchScores = zeros(numChirps,length(angles));
predictions = zeros(numChirps,1);

%Compare with each point
for u = 1 : numChirps

    %Use the known chirp to get the next test signature
    testChirp = trackSig(startInd+(u-1)*txNumSamps:startInd+(u-1)*txNumSamps+fs*chirp_duration-1);
    temp = xcorr(testChirp,TX_chirp); 
    testIR = temp(corKeepRange);
    normTestIR = testIR/sqrt(sum(testIR.^2));

    %Match the test signature to the reference signatures 
    % TODO switch to xcorr2
    for c = 1 : length(angles)
        matchScores(u,c) = max(abs(xcorr(ref_sigs(:,c),normTestIR)));
    end

    %Predict the position based on the maximum similarity
    [~, closeInd] = max(matchScores(u,:));
    prediction = closeInd-1;
    predictions(u) = prediction;

end

%Shift early negative predictions to positive
if predictions(1) > 300
    predictions(1) = predictions(1)-360;
end

% OPTIONAL: Ignore discontinuities by duplicating current point
jumpthresh = 75;
for i = 1:numChirps-1
    if abs(wrapTo180(predictions(i)-predictions(i+1))) > jumpthresh
        %predictions(i+1) = predictions(i);
    end
end


%Account for wrap around errors
predictions = rad2deg(unwrap(deg2rad(predictions)));

%Plot the experimental results
trackingT = 0:1/chirpsPerSec:numChirps/chirpsPerSec-1/chirpsPerSec;

figure
imagesc(flipud(matchScores'))
title('Matching scores at each point along the trajectory')
colorbar
clim([0,1])

figure
plot(trackingT,predictions)
title('Angular Trajectory')
hold on


% Compute the ground truth trajectory
trackingT = 0:1/chirpsPerSec:numChirps/chirpsPerSec - 1/chirpsPerSec;

grndTrthT = 0; %This initial point is equal to the time between starting the recording script and the movement script
grndTrth = pathTraveled(1);
prevTime = grndTrthT;

for i = 2 : length(pathTraveled)
    grndTrth = [grndTrth, pathTraveled(i)];
    grndTrthT = [grndTrthT, prevTime + abs(pathTraveled(i)-pathTraveled(i-1))/degPerS];
    prevTime = prevTime + abs(pathTraveled(i)-pathTraveled(i-1))/degPerS;

end

plot(grndTrthT,grndTrth)
legend('AMULET','Trajectory')
gndTrthPts = interp1(grndTrthT,grndTrth,trackingT);

%Compute the error for each prediction
errors = predictions-gndTrthPts';

%Account for initial offset errors
errors(isnan(errors)) = [];

disp("The average angular error is " + num2str(round(mean(abs(errors)),2)) + " degrees")


%% Figure 14(c)

%Load in the data
normUpperIRs = amulet_Data.simul_upperFreqs_IRs;
normLowerIRs = amulet_Data.simul_lowerFreqs_IRs;
trackSig = amulet_Data.pier_simulTracking;
%trackSig = rxSig;

%Specify the chirp parameters
fStart = 1000;
fStop = 88000;
fMid = 45000;
fs = 192000;
duration = 0.2;

numChirps = 400;


singleChirpT = (0:1/fs:duration-1/fs)';
singleChirp=chirp(singleChirpT,fStart,duration,fStop);

%Parameters of the recording
%txNumSamps = length(singleTX); %How many samples between start of consequtive TX chirps
txNumSamps = 153600; 


%Specify the density of calibration points
angleInc = 1; %Degrees between each sample
angles = 0:angleInc:360-angleInc;


%Specify the ground truth path information
pathTraveled = [0,180,160,180,160,180,0];


waitSpots = [0];
waitTime = 0;

chirpsPerSec = fs/txNumSamps; %Ensure these numbers are accurate  %4.7OG
%degPerS = 1.5*2.72/2/1.1;   %math says 1.59
degPerS = 1.45; %1,
%degPerS = 1.52; %3,2,4,8,9

keepLowerLimit = 25;
keepUpperLimit = 125;
corWidth = keepUpperLimit+keepLowerLimit+1;

t = singleChirpT;

%Define the normalized chirp
p = [4.383e-09,-4.126e-05,0.1138]; %Obtained from real data
normFunc = polyval(p,1:length(t));
[minv,mini] = min(normFunc);
normFunc(1:mini) = minv; %Set the minimum

normChirp = singleChirp./normFunc';
scaleFact = max(abs(normChirp));
normChirp = normChirp./scaleFact;


%Define the smooth start chirp signal
startSig = cos(2*pi*t*fStart);
stopSig = cos(2*pi*t*fStop);
midSig = cos(2*pi*t*fMid);

smoothStart = [startSig; singleChirp; stopSig];
smoothStart1 = [startSig; singleChirp(1:length(singleChirp)/2); midSig];
smoothStart2 = [midSig; singleChirp(length(singleChirp)/2:end); stopSig];

%Define the normalized smooth start chirp signal
normStartSig = startSig/minv/scaleFact;
normStopSig = stopSig/normFunc(end)/scaleFact;
normMidSig = midSig/scaleFact/normFunc(floor(length(normFunc)/2));
normSmooth = [normStartSig; normChirp; normStopSig];

normSmooth1 = [normStartSig; normChirp(1:length(singleChirp)/2); normMidSig];
normSmooth2 = [normMidSig; normChirp(length(singleChirp)/2:end); normStopSig];


% Identify the first correlation peak in the recording (used to align all
% subsequent chirps)
%temp = trackSig(85000:123400-1);
temp = trackSig(1:48000);

freqClearance = 6000;

numSamps = length(temp);
fshift = (-numSamps/2:numSamps/2-1)*(fs/numSamps);
mask1 = (abs(fshift) < fMid-freqClearance) & (abs(fshift) > fStart); %testing this
mask2 = (abs(fshift) < fStop) & (abs(fshift) > fMid+freqClearance);

tempFFT = fftshift(fft(temp));
tempFilt = ifft(ifftshift(tempFFT.*mask2'));

chirpLower = singleChirp(1:length(singleChirp)/2);
chirpUpper = singleChirp(length(singleChirp)/2:end);




corKeepRange = 64607:64768;
startInd = 1;



% Match the Chirps to Predict Angular Trajectory

%Setup variables
matchScores1 = zeros(numChirps,length(angles));
errors1 = zeros(numChirps,1);
predictions1 = zeros(numChirps,1);
matchScores2 = matchScores1;
errors2 = errors1;
predictions2=predictions1;



%Compare with each point
for u = 15 : numChirps-10 

    testChirp = trackSig(startInd+(u-1)*txNumSamps:startInd+(u-1)*txNumSamps+48000-1);
    testFFT = fftshift(fft(testChirp));
    upperFilt = ifft(ifftshift(testFFT.*mask2'));
    lowerFilt = ifft(ifftshift(testFFT.*mask1'));

    temp1 = xcorr(lowerFilt,chirpLower); %This expects to be honed in on IR
    %[~,maxi1] = max(abs(temp1));
    %testIR1 = temp1(maxi1-keepLowerLimit:maxi1+keepUpperLimit);
    testIR1 = temp1(corKeepRange);
    normTestIR1 = testIR1/sqrt(sum(testIR1.^2));

    temp2 = xcorr(upperFilt,chirpUpper); %This expects to be honed in on IR
    %testIR2 = temp2(corKeepRange);
    [~,maxi2] = max(abs(temp2));
    testIR2 = temp2(maxi2-keepLowerLimit:maxi2+keepUpperLimit);
    normTestIR2 = testIR2/sqrt(sum(testIR2.^2));

    for c = 1 : length(angles)
        matchScores1(u,c) = max(abs(xcorr(normLowerIRs(:,c),normTestIR1)));
        matchScores2(u,c) = max(abs(xcorr(normUpperIRs(:,c),normTestIR2)));
    end

    [~, closeInd] = max(matchScores1(u,:));
    prediction = (closeInd-1)*angleInc;
    predictions1(u) = prediction;

    [~, closeInd] = max(matchScores2(u,:));
    prediction = (closeInd-1)*angleInc;
    predictions2(u) = prediction;

end

%Shift early negative predictions
if predictions1(1) > 300
    predictions1(1) = predictions1(1)-360;
end

% Ignore discontinuities by duplicating current point
jumpthresh = 75;
for i = 1:numChirps-1
    if abs(wrapTo180(predictions1(i)-predictions1(i+1))) > jumpthresh
        %predictions(i+1) = predictions(i);
    end
end

%Account for wrap around errors
predictions1 = rad2deg(unwrap(deg2rad(predictions1)));


%Shift early negative predictions
if predictions2(1) > 300
    predictions2(1) = predictions2(1)-360;
end

% Ignore discontinuities by duplicating current point
jumpthresh = 75;
for i = 1:numChirps-1
    if abs(wrapTo180(predictions2(i)-predictions2(i+1))) > jumpthresh
        %predictions(i+1) = predictions(i);
    end
end

%Account for wrap around errors
predictions2 = rad2deg(unwrap(deg2rad(predictions2)));



%Plot the results
trackingT = 0:1/chirpsPerSec:numChirps/chirpsPerSec-1/chirpsPerSec;


figure
plot(trackingT,predictions1)
title('Angular Trajectory')
hold on
plot(trackingT,predictions2)


% Compute the ground truth trajectory
trackingT = 0:1/chirpsPerSec:numChirps/chirpsPerSec - 1/chirpsPerSec;

grndTrthT = 10; %This initial point is equal to the time between starting the recording script and the movement script
grndTrth = pathTraveled(1);
prevTime = grndTrthT;

for i = 2 : length(pathTraveled)
    grndTrth = [grndTrth, pathTraveled(i)];
    grndTrthT = [grndTrthT, prevTime + abs(pathTraveled(i)-pathTraveled(i-1))/degPerS];
    prevTime = prevTime + abs(pathTraveled(i)-pathTraveled(i-1))/degPerS;

    if ismember(i,waitSpots)
        grndTrth = [grndTrth, pathTraveled(i)];
        grndTrthT = [grndTrthT, prevTime+waitTime];
    end
end

plot(grndTrthT,grndTrth)
plot(grndTrthT,grndTrth+37.3)

gndTrthPts = interp1(grndTrthT,grndTrth,trackingT);

errors1 = predictions1-gndTrthPts'-37.3;
errors2 = predictions2-gndTrthPts';

%Account for initial offset errors
errors1(isnan(errors1)) = [];
errors2(isnan(errors2)) = [];

lowerFreqAvgEr = mean(abs(errors1));
upperFreqAvgEr = mean(abs(errors2));

disp("The average error of the lower frequencies is " + num2str(round(lowerFreqAvgEr,1)) + " degrees")
disp("The average error of the upper frequencies is " + num2str(round(upperFreqAvgEr,1)) + " degrees")




%% Figure 15

% NOTE: for simplicity, the final errors are provided after the analysis.
% To re-run the analysis you can uncomment the following block of code and
% comment out the load final results from the datafile line


%{
meanErrors = zeros(size(upperLims));
upper95 = meanErrors;
lower95 = meanErrors;

filtSigs1 = zeros(size(rxSigs1));
filtSigs2 = zeros(size(rxSigs1));
IRs1 = zeros(keepUpperLimit+keepLowerLimit+1,length(angles));
IRs2 = IRs1;
normIRs1 = IRs1;
normIRs2 = IRs1;

maskLen = length(rxSigs);

for i = 1 : length(lowerLims)

    %Define the mask
    numSamps = maskLen;
    fshift = (-numSamps/2:numSamps/2-1)*(fs/numSamps);
    mask = (abs(fshift) < upperLims(i)) & (abs(fshift) > lowerLims(i));

    %FFT all data for filtering
    rxFFTs1 = fft(rxSigs1);
    rxFFTs2 = fft(rxSigs2); 

    %Apply the filter to each signal
    for a = 1 : length(angles)
        %Filter out non-chirp freqs
        rxFFT = fftshift(rxFFTs1(:,a));
        filtSigs1(:,a) = ifft(ifftshift(rxFFT.*mask'));
        
        rxFFT = fftshift(rxFFTs2(:,a));
        filtSigs2(:,a) = ifft(ifftshift(rxFFT.*mask'));
    end


    %Use the max bandwidth to set the IR range
    if i == 1
        [cor1,lags1] = xcorr(filtSigs1(:,1),singleChirp);
    
        corOffset = floor(length(cor1)/2);
        [pks,locs] = findpeaks(abs(cor1),'MinPeakHeight',5.5,'MinPeakDistance',10);
        corKeepRange = locs(1)-keepLowerLimit : locs(1)+keepUpperLimit;
    end


    %Grab all of the IRs
    cors = xcorr2(filtSigs1,singleChirp);
    IRs1 = cors(corKeepRange-(length(cor1)-length(cors)),:);

    cors = xcorr2(filtSigs2,singleChirp);
    IRs2 = cors(corKeepRange-(length(cor1)-length(cors)),:);


    %Normalize the IR energies
    for a = 1:length(angles)
        normIRs1(:,a) = IRs1(:,a)/sqrt(sum(IRs1(:,a).^2));
        normIRs2(:,a) = IRs2(:,a)/sqrt(sum(IRs2(:,a).^2));
    end


    %Matching Scores
    matchScores = zeros(length(angles));

    errors = zeros(length(angles),1);
    predictions = zeros(length(angles),1);

    %Compute the score for each pair of angles
    for u = 1 : length(angles)
        for c = 1 : length(angles)
             matchScores(u,c) = max(abs(xcorr(normIRs1(:,c),normIRs2(:,u))));
        end
    
        %Compare the prediciton to the theoretical 
        target = angles(u);
        [~, closeInd] = max(matchScores(u,:));
        prediction = (closeInd-1)*angleInc;
        
        predictions(u) = prediction;
        errors(u) = abs(wrapTo180(target-prediction));
    end

end
%}



% Load in the final results instead of computing them
meanErrors = amulet_Data.ablation_errors; %Comment out to recalculate

%Analysis setup parameters
BWStep = 5000;
fMax = 88000;
fMin = 1000;
minBW = 1000;
maxBW = fMax-fMin;

numBWs = length(fMin+minBW:BWStep:fMax) + 1;
upperLims = [fMax (fMin+minBW):BWStep:fMax ones(1,numBWs)*fMax];
lowerLims = [fMin ones(1,numBWs)*fMin (fMax-minBW):-BWStep:fMin];
tests = [lowerLims;upperLims];

BWs = upperLims-lowerLims;

upperBWInds = [1, length(upperLims):-1:numBWs+4];
lowerBWInds = [4:numBWs, 1];



%Generate output plots
figure
plot(BWs(lowerBWInds),meanErrors(lowerBWInds),'o-')
hold on
plot(BWs(upperBWInds),meanErrors(upperBWInds),'+-')
legend('Lower Frequencies','Upper Frequencies')
xlabel('Bandwidth (kHz)')
ylabel('Mean Error (degrees)')
title('Ablation Study')


%% Figure 16(a)

% Load in the precomputed results
results = amulet_Data.depth_variation;
x = ["0.8","1.6","2.4"]; 


figure
bar(x,results)
ylabel("Error (Degrees)")
ylim([0 round(max(max(results))*1.2,1)])
set(gca, 'FontName', 'Arial', 'FontSize',20)
xlabel("Evalutation Depth (m)")
lgd = legend('0.8m', '1.6m', '2.4m');
title(lgd, 'Calibration Depth')

%% Figure 16(b)

% Load in the precomputed results
results = amulet_Data.modulation_variation;
avg = results(1,:);
upper95 = results(2,:);
lower95 = results(3,:);
x = ["FSK","OFDM","Chirp"]; 

figure
upperOffset = upper95-avg;
lowerOffset = avg-lower95;
yMax = round(max(upper95)*1.3,1);
bar(x,avg)
hold on
er = errorbar(1:length(x),avg,lowerOffset,upperOffset,'LineWidth',2,'CapSize',10);
er.LineStyle = 'none';
er.Color = '[0 0 0]';
ylabel("Error (Degrees)")
xlabel("Modulation Scheme")
ylim([0,yMax])
ytix = get(gca, 'YTick');
set(gca, 'FontName', 'Arial', 'FontSize',20)


%% Figure 16(c)

% NOTE: The SNR simulation adds in AWGN, so the results may have
% minor varitions compared with what is shown in the paper

%Load in the data
raw_signals = amulet_Data.amulet_saltwater_rawSig1;
cali_signatures = amulet_Data.amulet_saltwater2;

[len,~] = size(raw_signals);

noisy_signals = raw_signals;
SNRs = -4:-2:-20;
meanErrors = zeros(length(SNRs),1);
matchScores = zeros(360);
predictions = zeros(360,1);
errors = zeros(360,1);
IRs = zeros(keepLowerLimit+keepUpperLimit+1,360);
normIRs = IRs;


for i = 1 : length(SNRs)

    for a = 1:360
        noisy_signals(:,a) = raw_signals(:,a) + awgn(ones(len,1)*sqrt(sum(abs(raw_signals(:,a)).^2)/len), SNRs(i),'measured')-sqrt(sum(abs(raw_signals(:,a)).^2)/len);
    end

    cors2 = xcorr2(noisy_signals,TX_chirp);

    for a = 1:360
        [~,maxi] = max(abs(cors2(:,a)));
        IRs(:,a) = cors2(maxi-keepLowerLimit:maxi+keepUpperLimit,a);
        normIRs(:,a) = IRs(:,a)/sqrt(sum(IRs(:,a).^2));
    end

    
    for u = 1 : 360
        for c = 1 : 360

            matchScores(u,c) = max(abs(xcorr(cali_signatures(:,c),normIRs(:,u))));

        end

        %Compare the prediciton to the theoretical 
        target = u-1;
        [~, closeInd] = max(matchScores(u,:));
        prediction = (closeInd-1)*angleInc;
        
        predictions(u) = prediction;
        errors(u) = abs(wrapTo180(target-prediction));
    end

    meanErrors(i) = mean(errors);

end


figure
plot(SNRs,meanErrors)
ylabel('Mean Errors (Degrees)')
xlabel('SNR')


%% Figure 16(d)


%Specify the interpolation spacings to plot
interp_spacings = [1,5,10,20,30,50,90];

%Specify the data to use
cali_signatures = amulet_Data.amulet_saltwater1; %will be interpolated
test_signatures = amulet_Data.amulet_saltwater2; %will be estimated

interpIRs = zeros(size(test_signatures));
overall_errors = zeros(length(interp_spacings),1);

for itsp = 1:length(interp_spacings)
    for i = 1 : interp_spacings(itsp) : numAngles
       
        %Set the anchor IRs to interp between
        bottomIR = cali_signatures(:,i);
    
        %Handle wrap around for the last case
        if i+interp_spacings(itsp) < numAngles
            topIR = cali_signatures(:,i+interp_spacings(itsp));
        else
            topIR = interpIRs(:,1);
        end
    
        %Store the lower anchor in the output
        interpIRs(:,i) = bottomIR;
    
        %Linear interpolation between the anchors
        for j = 1:interp_spacings(itsp)-1
            interpIRs(:,i+j) = (topIR-bottomIR)*j/(interp_spacings(itsp))+bottomIR;
        end
    end

    %Reset the storage variables
    matchScores = zeros(numAngles);
    interpMatchScores = matchScores;
    errors = zeros(numAngles,1);
    interpErrors = errors;
    predictions = zeros(numAngles,1);
    interpPredictions = predictions;
    
    
    %Normalize the interpolated calibration signal energy
    normInterpCaliIRs = zeros(size(test_signatures));
    for i = 1:numAngles
        normInterpCaliIRs(:,i) = interpIRs(:,i)/sqrt(sum(interpIRs(:,i).^2));
    end
    
    
    
    for u = 1 : numAngles
        for c = 1 : numAngles
    
            interpMatchScores(u,c) = max(abs(xcorr(normInterpCaliIRs(:,c),test_signatures(:,u))));

        end
    
        %Compare the prediciton to the theoretical 
        target = u-1;
    
        [~, closeInd] = max(interpMatchScores(u,:));
        prediction = closeInd-1;
        
        interpPredictions(u) = prediction;
        interpErrors(u) = abs(wrapTo180(target-prediction));
    end


    overall_errors(itsp) = mean(abs(interpErrors));

end

figure
plot(interp_spacings,overall_errors)
title('Average Error vs Interpolation Spacing')
xlabel('Interpolation Spacing (Degrees)')
ylabel('Average Error (Degrees)')

toc