%% Automated Seizure Detection Using Short-Time Energy on All Channels
clear; close all; clc;

% Sampling frequency (adjust as needed)
Fs = 480;

%% 1) Load EEG Data
% Read the CSV file (assumed to have a header row) and extract columns:
% Column 1: Elapsed Time, Columns 2: FP1, 3: FP2, 4: O1, 5: O2.
rawData = readmatrix('EEGrelaxed.csv', 'Range','2:1048576');
timeData = rawData(:,1);          % Elapsed Time (assumed in seconds)
EEGdata = rawData(:,2:5);         % All available channels: FP1, FP2, O1, O2

% In case "timeData" is not valid, create a time vector based on Fs
if isempty(timeData) || max(timeData) < 1
    N = size(EEGdata,1);
    t = (0:N-1)/Fs;
else
    t = timeData;
end

% Ensure data is numeric and remove any non-finite values
EEGdata(~isfinite(EEGdata)) = 0;

%% 2) Short-Time Energy Computation (for Seizure Detection)
% Parameters for windowing
windowSec = 1;                    % Window length in seconds
windowLen = round(windowSec * Fs); % Window length in samples
overlap = round(0.5 * windowLen);    % 50% overlap
step = windowLen - overlap;
numWindows = floor((size(EEGdata,1) - windowLen) / step) + 1;

% Pre-allocate an array to store energy for each window and channel
energy = zeros(numWindows, size(EEGdata,2));

% Calculate energy for each window and each channel
for ch = 1:size(EEGdata,2)
    for win = 1:numWindows
        idx_start = (win-1)*step + 1;
        idx_end = idx_start + windowLen - 1;
        segment = EEGdata(idx_start:idx_end, ch);
        energy(win, ch) = sum(segment.^2);
    end
end

% Compute a detection threshold for each channel based on baseline statistics.
% Here, we use mean + 2*std as a simple threshold.
threshold = mean(energy) + 2*std(energy);

% Detect seizure windows: if any channel's energy exceeds its threshold, flag that window.
seizureWindows = any(energy > threshold, 2);

% Create a time vector for the center of each window
winTime = (((0:numWindows-1)*step) + windowLen/2) / Fs;

%% 3) Plot Short-Time Energy for Each Channel with Thresholds
figure('Name','Short-Time Energy per Channel');
for ch = 1:size(EEGdata,2)
    subplot(4,1,ch);
    plot(winTime, energy(:, ch), '-o','LineWidth',1.5);
    hold on;
    yline(threshold(ch), 'r--','LineWidth',1.2);
    xlabel('Time (s)');
    ylabel(['Energy, Ch ' num2str(ch)]);
    title(['Short-Time Energy for Channel ' num2str(ch)]);
    grid on;
end

%% 4) Plot Seizure Detection Results
% Plot a binary indicator showing which windows were flagged
figure('Name','Seizure Detection');
plot(winTime, seizureWindows, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Detection (1 = Seizure, 0 = Normal)');
title('Automated Seizure Detection (Any channel exceeding threshold)');
grid on;
