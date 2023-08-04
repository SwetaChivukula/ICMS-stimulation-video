% clear all, close all, clc
% find data snippet in samples based on electrode 10's neural data
% Create data matrix numElectrodex x numSamplesInSnippet
% Get the video heatmap code section working
%% Load the Neural Data
% ADD CODE HERE - update file names and locations
% Open the neural data file
fileLocation = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk\NSP Data';
fileName = 'RP01_ICMS_2-Contact_Characterization_2023-01-06-15-06-32-152_NSP1006.ns5';
fileStr = strcat(fileLocation,filesep,fileName);
openNSx(fileStr,'read'); %make sure NPMK folder is added to MATLAB's path
% Open the sensory data file
fileLocationSensory = 'C:\Users\chivu\OneDrive\Documents\MATLAB\research_project_graczyk';
fileNameSensory = 'RP01_ICMS_2-Contact_Characterization_2023-01-06-15-06-32-152.mat';
fileStrSensory = strcat(fileLocationSensory,filesep,fileNameSensory);
load(fileStrSensory);
%% Load neural data
neuralDataUnmapped = (double(NS5.Data)./4).*10^(-6);
% May need to break this up into smaller steps to avoid MATLAB/computer
% overload
%% Organize the time
fs = 30000;
dt = 1/fs;
%TendSamples = length(NS5.Data(1,:));
%TendSec = TendSamples/fs;
TendSec = NS5.MetaTags.DataDurationSec;
TendSamples = NS5.MetaTags.DataPoints;
TstartSec = 0;
TstartSamples = 1;
timeAll = TstartSec:dt:TendSec-dt;
%% Map the neural data from channel-indexing to electrode-indexing
electrodeNumber = zeros(128,1);
zeroChar = double('0');
for i = 1:128
   elecChar =  NS5.ElectrodesInfo(i).Label(7:9);
   elecDouble = double(elecChar);
   elecDouble2 = elecDouble(elecDouble ~= 0);
   numDigits = size(elecDouble2,2);
   if(numDigits == 1)
       elecNumber = elecDouble2-zeroChar;
   elseif(numDigits == 2)
       elecNumber = sum([10 1] .* (elecDouble2-zeroChar));
   elseif(numDigits == 3)
       elecNumber = sum([100 10 1] .* (elecDouble2-zeroChar));
   else
       disp('Error')
   end
   electrodeNumber(i) = elecNumber;
end
elecChannel = [electrodeNumber,[1:128]'];
elecChannelSorted = sortrows(elecChannel,1);
neuralData = neuralDataUnmapped([elecChannelSorted(:,2);129;130],:);
%% Plot the Channel-specific information for visualization
%Stimulated electrode = 10
%Stimulated electrode = 123
%Synch pulse - high during stim, started a little early - channel 130
%Monitor port - stimulation waveform, precise timing - channel 129
for electrode = 129 %1:128
    figure
    plot(timeAll,neuralData(electrode,:));
    xlabel('Time (secs)')
    ylabel('Voltage (uV)')
    titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
    title(titleString);
    hold on;    
end
synchOutput = neuralData(129,:);
monitorPortOutput = neuralData(130,:);

%% Determine duration of stim artifact
thresh = 0.0000133; % determined from plot of electrode 129
for electrode = 129
    stimDuration = timeAll(neuralData(electrode,:) >= thresh);
end
% use for more precise threshold without outliers
% stimDuration(1,length(stimDuration));
% stimDuration(1,1);
% stimDuration = stimDuration(1,length(stimDuration)) - stimDuration(1,1);

%% Use sample and hold to extract artifact from monitor port
counter = 1;
sample = 0;
holdData = neuralData(130,:);
for counter = 1:length(timeAll)
    if timeAll(1,counter) == stimDuration(1,1)
       if sample == 0
          sample = holdData(1,(counter-1)); % determine sample value
          index = counter; % set the index for first thresh cross
       end
    else
       counter = counter + 1;
    end
end

for counting = 1:length(stimDuration)
    holdData(1,index) = sample;
    index = index + 1;
end

% plot with the sample and hold
for electrode = 130 %1:128
    figure
    plot(timeAll,holdData(1,:));
    xlabel('Time (secs)')
    ylabel('Voltage (uV)')
    titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
    title(titleString);
    hold on;    
end

%sample = 1.4550e-04
%index = 663886
%% Extract the stimulation onset times
figure
plot(timeAll,monitorPortOutput);
thresholdMonitorPort = 6400; % Need to adjust for each file
MonitorPortBinary = monitorPortOutput<thresholdMonitorPort; % Create a logical array of above and below threshold samples
yline(thresholdMonitorPort);
%% Extract data snippets just before, during, and just after stimulation
% Try to keep track of which electrode was stimulated when
dataTime = 22:27; % CHANGE TO ACTUAL DATA (22 seconds to 27 seconds)
neuralData_new = neuralData;
for c = 660066:810081
    if timeAll(c) >= 22 && timeAll(c) <= 27
        neuralData_new(electrode,c) = 0;
    end
end
plot(timeAll, neuralData(129,:))
s1 = find(timeAll == 21.5);
s2 = find(timeAll == 23.5);
plotting_new = neuralData(:, s1:s2);
hold off
figure
plot(timeAll, neuralData_new(electrode, :))
 xlabel('Time without snippet (secs)')
 ylabel('Voltage (uV)')
 titleString = strcat("Voltage Plot for Electrode ",num2str(electrode));
 title(titleString);
% Set data = number of electrodes x number of samples for snipped
%% Video of stimulation spread across time
% ADD CODE HERE: update the time vector to be in samples and make it work
binSize = 1; %samples
dt = binSize/fs;
timeVec = -1+dt:dt:3; % change to 0:numSamples; eventually redo with seconds
figure(1)
subplot(2,1,2);
hold on;
%plot(timeVec,monitorPortOutput);
stimElectrode = neuralData(10, s1:s2);
timeVec = (1:length(stimElectrode))/fs;
plot(timeVec,stimElectrode); %neural data at electrode 10, 1 x samples vector
%xline(0,'--k','LineWidth',2); % Add a line at time 0 to indicate the stimulus presentation
%xline(2,'--k','LineWidth',2); % Add a line at time 1 to indicate the stimulus end
c = xline(0,'r','LineWidth',1);
titleStr = strcat("Electrode 10");
title(titleStr);
xlabel('Time (sec)');
% xlabel('Samples');
ylabel('Voltage (V)');
set(gca,'FontSize',16)
data = plotting_new(1:128,:);
%% stim video
starting = 1;
for timePoint = 1:50:length(timeVec)
subplot(2,1,1);
plotData = data(:,timePoint);
% Map data spatially
dataImage = reshape(plotData,8,16); % reshape the data into an array shape
dataImage(1:8,1:8)=flip(dataImage(1:8,1:8)',2); % reshape 1/2 the data into the right spatial layout
dataImage(1:8,9:16)=flip(dataImage(1:8,9:16)',2); % reshape the other 1/2
cmin = -0.01; % need to tune
cmax = 0.01; % need to tune
clims = [cmin cmax]; % set the colorbar limits
[x1,y1] = meshgrid(1:size(dataImage,2), 1:size(dataImage,1)); % create first set of sample points
[finerx,finery] = meshgrid(1:0.01:size(dataImage,2), 1:0.01:size(dataImage,1)); % create second set of sample points
newData = interp2(x1,y1,dataImage,finerx,finery); % interpolate the data points
imagesc(newData,clims)
% imagesc(dataImage,clims); % plot the image scaled by the data with colorbar limits
% imagesc(dataImage) % plot the image without scaling
%optional - apply smoothing
colorbar % show the colorbar
titleStr = strcat("Heatmap of ICMS Stimulation Spread");
title(titleStr);
set(gca,'FontSize',16)
xlabel('Lateral Array   -   Medial Array');
set(gcf,'position',[100,100,1000,1000])
% smoothing
% cmin = -2;
% cmax = 4;
% clims = [cmin cmax]; % set the colorbar limits
% [x1,y1] = meshgrid(1:size(dataImage,2), 1:size(dataImage,1)); % create first set of sample points
% [finerx,finery] = meshgrid(1:0.01:size(dataImage,2), 1:0.01:size(dataImage,1)); % create second set of sample points
% newData = interp2(x1,y1,dataImage,finerx,finery); % interpolate the data points
% imagesc(newData,clims); % plot the image scaled by the data with colorbar limits
% colorbar % show the colorbar
% titleStr = strcat("Heatmap Rate = 60"," Depth = 5");
% title(titleStr);
% set(gca,'FontSize',16)
% xlabel('Lateral Array   -   Medial Array');
% set(gcf,'position',[100,100,1000,1000])
% % reset the axis labels to original size
% set(gca, 'XTick', linspace(1,size(finerx,2),size(x1,2))); 
% set(gca, 'YTick', linspace(1,size(finerx,1),size(x1,1)));
% set(gca, 'XTickLabel', 1:size(x1,2));
% set(gca, 'YTickLabel', 1:size(x1,1));

subplot(2,1,2); % Can specify subplot size here [l,b,w,h]
delete(c);
plot(timeVec,stimElectrode);
c = xline(timeVec(timePoint),'r','LineWidth',1);
pause(0.01)
% Can use the timer function to be more efficient...
end
% if you have bonus time, could you show a video of the average spread of
% stimulation artifact for each electrode across time?
%% Clean up
close all;