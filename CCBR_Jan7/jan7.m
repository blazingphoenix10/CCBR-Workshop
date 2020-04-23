%% CCBR 7th January 2020 - Assignment 1

% Part 1: Spike train of a neuron
% Modeling neuron spikes as a Memoryless Poisson process

%{
clear all

FR = 0;                             % Desired Firing rate - number of spikes per second
Dt = 1e-3;                          % 1 ms
time = 200;                         % Time in seconds for which the expt is carried out for

i = 1;                              % Setting for loop index
MISI = zeros(100,1);
DISI = zeros(100,1);

for FR=1:100
    timePoints = time/Dt;               % Total number of datapoints under study
    R = rand(timePoints, 1);            % Creating a random vector
    Lambda = FR*Dt;                     % Threshold for a spike to occur
    Spikes = zeros(timePoints, 1);      % Creating a zero vector of length of number of data points
    Spikes(R<Lambda) = 1;               % Setting values below threshold to fire
    ObservedFR = sum(Spikes==1)/time;   % Observed Firing rate
    SpikeIndices = find(Spikes==1);     % Finding indices where the neuron fires
    SpikeTimes = SpikeIndices*Dt;       % Identifying times associated with the firing of the neuron
    ISI = diff(SpikeTimes);             % Inter Spike Interval
    MISI(i) = mean(ISI);
    DISI(i) = std(ISI);
    i = i + 1;
end

figure; loglog(MISI, DISI, '.')
xlabel("Mean of Inter Space Interval");
ylabel("Std Dev of Inter Space Interval");

binSize = 20;                           % Entire data is divided into 20 equal parts
figure; hist(ISI,binSize);
title("Plot at 100 Hz firing rate")
xlabel("Time in ms")
ylabel("Counts")
%}

%%

% Part 2: Autocorrelation

%{
[correlation, lags] = xcorr(Spikes,'coeff');            % coeff normalizes the Spikes array autocorrelation values
plot(correlation, lags,'.')                                 % Symmetric about lag = 0
%}

% Why are we even doing autocorrelation?
% ?

%%

clear all

% Part 3: *Periodic* spike train

T = 20;
Dt = 1e-3;
timePoints2 = T/Dt;
Spikes2 = zeros(timePoints2, 1);
Spikes2(1:20:end) = 1;

% Q1(5 points) How many spikes would you expect in this spike train of 20seconds?
disp("Answer for Question 1: We can expect 20000/20 spikes i.e 1000 spikes, with the last spike at 19,981") 

% Q2(10 points) How would you check that?(function/code in MATLAB)
disp("Answer for Question 2")
sum(Spikes2 == 1)

% Q3(15 points) What is the mean firing rate for this spike train? Store this value in a
% variable called MFR.

disp("Answer for Question 3")
MFR = sum(Spikes2==1)/T

[correlation2,lags2] = xcorr(Spikes2,"coeff");

figure; plot(correlation2,lags2,'.');
xlabel("Normalized autocorrelation values");
ylabel("Lag values");

%% HW problems

%{
Q4.(30 points) Given the actual data from a hippocampal neuron, compute the
ISI histogram and autocorrelation. 
Load this data using the command load('HW1Data_SpikeTrain.mat'); make sure that you "set path" to the 
folder containing this file. Alternatively, use function cd('FolderLocationForTheFile')
Start with the Spiketimes (ST) and the time bins(T_bc, in seconds). 
Compute ISI histogram of the Spiketimes. 
Next use, "hist" function to compute Spike counts in each time bin.
Then compute autocorrelation. 

What is the time bin(Spike counts are provided as vector SC to cross check )? 

How long was this recording made? If the time bins are in seconds, 
at what resolution(in seconds) is the experiment being recorded?(write code to obtain this)
%}

clear all

% cd("CCBR_Jan7")
load("HW1Data_SpikeTrain.mat");
ISI2 = diff(ST);
% ISI Histogram of Spike Times
subplot(1,3,1);
hist(ISI2,T_bc);
title("Histogram of ISI")
xlabel("Time in ms")
ylabel("Counts")

%[SC1, edges] = hist(ISI2,T_bc);
% Using hist function to find spike counts
[SC1, edges] = hist(ST,T_bc);

[correlation3,lags3] = xcorr(SC,"coeff");
subplot(1,3,2);
plot(correlation3,lags3,'.');
title("For given SC");
xlabel("Normalized autocorrelation values");
ylabel("Lag values");

[correlation4,lags4] = xcorr(SC1,"coeff");
subplot(1,3,3);
plot(correlation4,lags4,'.');
title("For calculated SC");
xlabel("Normalized autocorrelation values");
ylabel("Lag values");

% What is the time bin(Spike counts are provided as vector SC to cross check )? 
disp("Bin Size is")
binSize = T_bc(2) - T_bc(1)

% How long was this recording made? If the time bins are in seconds,
% at what resolution(in seconds) is the experiment being recorded?(write code to obtain this)
disp("Time taken for the recording")
T = T_bc(end) + 0.005 % Since hist function gives the mid point of the bin

%%

clear all

%{
Q5. (40 points) Create 2 spike trains, at a resolution of 1ms for total
duration of 5second. Let the first spike train fire every 500ms, and the second one
fire every 700ms(both have their first spike at t=1ms). 

What would you expect to see in the autocorrelation of each of the spike trains?

Use "xcorr" function and verify your expectations. Next use the xcorr function 
to compute the cross-correlation between these spike trains. Comment in 5 lines about what you see, and how
& why it is expected or unexpected.
%}

T = 5;
dT = 1e-3;
TimePoints = T/dT;

% What would you expect to see in the autocorrelation of each of the spike trains?
% Symmetric behaviour around lag 0, with lag 0 correlation value being the
% highest, and the correlation value decreasing as we move away from lag 0
% The difference between the 500ms autocorrelation and the 700ms
% autocorrelation plot will be the rate of decrease in correlation values.

Spikes1 = zeros(TimePoints, 1);
Spikes1(1:500:end) = 1;
[correlation1,lags1] = xcorr(Spikes1,"coeff");
subplot(1,4,1);
plot(correlation1,lags1,'.');
title("Firing every 500ms");
xlabel("Normalized autocorrelation values");
ylabel("Lag values");

Spikes2 = zeros(TimePoints, 1);
Spikes2(1:700:end) = 1;
[correlation2,lags2] = xcorr(Spikes2,"coeff");
subplot(1,4,2);
plot(correlation2,lags2,'.');
title("Firing every 700ms");
xlabel("Normalized autocorrelation values");
ylabel("Lag values")

[correlation3,lags3] = xcorr(Spikes1,Spikes2);
subplot(1,4,3);
plot(correlation3,lags3,'.');
title("Cross correlation plot");
xlabel("Autocorrelation values");
ylabel("Lag values")

[correlation4,lags4] = xcorr(Spikes1,Spikes2,"coeff");
subplot(1,4,4);
plot(correlation3,lags3,'.');
title("Normalized cross correlation plot");
xlabel("Normalized Autocorrelation values");
ylabel("Lag values")

% Why do you see such a behviour for the cross correlation plot?
% It should meet every 3500 ms, the auto correlation doesn't seem to be
% symmetric either.