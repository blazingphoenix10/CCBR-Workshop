% RUPAK K V
% BE15B023 IIT MADRAS

%% 8th Jan Tutorial
%{
%Create PSTH
clear all
clc
T=5; % Time for 1 trial in secs
N=100; % This is the number of trials
Dt=1e-3; %Dt is the discretization we are using, 10^-3 means 1ms intervals
Timepoints=T/Dt; % Total time points 
TimeVector=Dt:Dt:5;
FR=5; %We set the desired firing rate of the neuron as 5Hz, ie 5 spikes per second
%ie 5 spikes per 1000 data points
Lambda=FR*Dt; %Lambda, which is a threshold
Rand=rand([Timepoints,N]); %Create a matrix of random numbers
Spikes=zeros([Timepoints,N]); % Initiate a matrix with all zeros
%Note now we have 5000 datapoints in 1 trial and 100 such trials total, so
%the dimension of Spikes is 5000 x 100
Spikes(Rand<Lambda)=1; % We are going to consider spikes to be binary, 
%and when the random number exceeds the threshold, 
%we switch "Spikes" at that time to 1
figure;
subplot(3,1,1)%break down a single figure into 3 subplots, arranged into 3 rows, and 1
%columns, then plot out the first of these subplots as below.
[Index_time,Index_trial]=find(Spikes==1);%We want the index corresponding to time and 
%trial of when the spike happened, ie when Spikes was =1. So we use the
%find function to get those indices
plot(Index_time,Index_trial,'.','markersize',1)%Set size of points to 1, so that they 
%do not overlap with each other

% Compute its Tuning curve 
TuningCurve=sum(Spikes');% We take the transpose to ensure that we are 
%summing across trials
subplot(3,1,2);%Still using the same figure, with 3 subplots and
%now we are plotting the second one out of these.
plot(TimeVector,TuningCurve);
xlabel('Time in seconds')
ylabel('Total Spikes in 100 trials')
%Note that we can bin the TuningCurve above, to attain a lower sampling
%rate and a smoother, less noisy tuning curve
%Hence instead of computing the Tuningcurve at 1ms resolution(5000
%datapoints) let us compute at 100ms resolution, by binning every 100 data
%points together
NewTimeVector=TimeVector(100:100:end);%Take only every 100th datapoint
for i=1:Timepoints/100%This is a loop which will compute the sum of 
    %spikes in jumps of 100 time-points
NewTuningCurve(i)=sum(TuningCurve((100*(i-1)+1):100*i));
end
subplot(3,1,3)
plot(NewTimeVector,NewTuningCurve);
xlabel('Time in seconds')

%% Define sparsity and depth of modulation
temp=makeSparsity(NewTuningCurve);
fprintf(['Untuned spike matrix sparsity= ' num2str(temp) '\n'])
%}

%%
%Q1. Create 4 sinusoids(over the range of 0 to pi with a resolution of 0.1)
%with different scaling factors, and with different offsets,
%and compute their sparsities.
%eg. a1_eg=sin(0:0.1:pi)*3 +40;%where 40 is the offset, and 3 is the
%scaling factor. Comment in 4 lines about how sparsity is affected by both
%these parameters. (15 points)

a1=sin(0:0.1:pi)*3 + 40;
sa1=makeSparsity(a1);
fprintf(['a1 sparsity= ' num2str(sa1) '\n'])

a2=sin(0:0.1:pi)*6 + 40;
sa2=makeSparsity(a2);
fprintf(['a2 sparsity= ' num2str(sa2) '\n'])

a3=sin(0:0.1:pi)*3;
sa3=makeSparsity(a3) + 80;
fprintf(['a3 sparsity= ' num2str(sa3) '\n'])

a4=sin(0:0.1:pi)*6 + 80;
sa4=makeSparsity(a4);
fprintf(['a4 sparsity= ' num2str(sa4) '\n'])

% let s=scaling factor and o=offset
% When s,o are changed by the same factor, the sparsity does not change.
% Sparsity is directly proportional to scaling factor.
% Sparsity is inversely proportional to Offset.

%% PSTH of a tuned neuron
%{
% Create tuned neuron, with varying lambda, find its sparsity
%Create PSTH
clear all
T=5; % Time for 1 trial in secs
N=100; % This is the number of trials
Dt=1e-3; %Dt is the discretization we are using, 10^-3 means 1ms intervals
Timepoints=T/Dt; % Total time points 
TimeVector=Dt:Dt:5;

%We set the desired firing rate of the neuron in each trials with the
%following criteria
%1. Fires 0Hz away from Time t=2.5 seconds
%2. Has a peak firing rate of 20Hz at t=2.5 seconds
%3. Firing rate varies sinusoidally 
Equispaced_Angles=linspace(0,pi,Timepoints);% We create 5000 datapoints from 0 to pi
%which we can use as the angles whose "sin" is taken to compute FR below
FR=0 + 20*sin(Equispaced_Angles); % Check that the minimum value is 0 and maximum is 20
Lambda=FR*Dt; %Lambda, which is a threshold
%Since we want the same threshold for all trials we repeat this vector 100
%times
Lambda_AllTrials=repmat(Lambda,100,1);
Rand=rand([Timepoints,N]); %Create a matrix of random numbers
Spikes=zeros([Timepoints,N]); % Initiate a matrix with all zeros
%Note now we have 5000 datapoints in 1 trial and 100 such trials total, so
%the dimension of Spikes is 5000 x 100
%Note that Lambda on the other hand is a matrix of 100 x 5000 so we will
%have to take its transpose below 
Spikes(Rand<transpose(Lambda_AllTrials))=1; % We are going to consider spikes to be binary, 
%and when the random number exceeds the threshold,we switch "Spikes" at that time to 1
figure;
subplot(3,1,1)%break down a single figure into 2 subplots, arranged into 2 rows, and 1
%columns, then plot out the first of these subplots as below.
[Index_time,Index_trial]=find(Spikes==1);%We want the index corresponding to time and 
%trial of when the spike happened, ie when Spikes was =1. So we use the
%find function to get those indices
plot(Index_time,Index_trial,'.','markersize',1)%markersize determines the size of the dots
% Compute its Tuning curve 
TuningCurve=sum(Spikes');% We take the transpose to ensure that we are summing across trials
subplot(3,1,2);%Still using the same figure, with 2 subplots in 2 row, 1 columns and
%now we are plotting the second one out of these.
plot(TimeVector,TuningCurve);
xlabel('Time in seconds')
ylabel('Total Spikes in 100 trials')
%Note that we can bin the TuningCurve above, to attain a lower sampling
%rate and a smoother, less noisy tuning curve
%Hence instead of computing the Tuningcurve at 1ms resolution(5000
%datapoints) let us compute at 100ms resolution, by binning every 100 data
%points together
NewTimeVector=TimeVector(100:100:end);%Take only every 100th datapoint
for i=1:Timepoints/100
NewTuningCurve(i)=sum(TuningCurve((100*(i-1)+1):100*i));
end
subplot(3,1,3)
plot(NewTimeVector,NewTuningCurve);
xlabel('Time in seconds')

%Next, compute the sparsity of this tuning curve
temp=makeSparsity(NewTuningCurve);
fprintf(['Tuned spike matrix sparsity= ' num2str(temp) '\n'])
%}
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Next we are going to go through the code related to shuffling a spike matrix to
%obtain a null-distribution and eventually quantify significance of tuning
ShuffledSpikeMatrix=ShuffleTrials(Spikes);

%Next we compute the tuning curve of this Shuffled spike matrix and compute
%its sparsity
 TuningCurve_shuffled=sum(ShuffledSpikeMatrix');% We take the transpose to ensure that we are summing across trials
figure
subplot(2,1,1)
plot(TimeVector,TuningCurve_shuffled);

NewTimeVector=TimeVector(100:100:end);%Take only every 100th datapoint
for i=1:Timepoints/100
NewTuningCurve_shuffled(i)=sum(TuningCurve_shuffled((100*(i-1)+1):100*i));
end
subplot(2,1,2)
plot(NewTimeVector,NewTuningCurve_shuffled);
xlabel('Time in seconds')
temp=makeSparsity(NewTuningCurve_shuffled);

fprintf(['Shuffled spike matrix sparsity= ' num2str(temp) '\n'])
%}


%% Q2
%{
Q2. Create an artificial spike train, for 100 trials with the following 
 expected firing rate profile over 5 seconds of recording per trial at 1ms
 resolution.
a. Baseline firing rate of 2Hz. 
b. Sinusoidal profile from time points of 2 to 3sec, with baseline firing rate of
 2Hz but a peak firing rate of 10Hz
Plot this expected firing rate profile. (10 points)
Compute the sparsity of the Spikes matrix obtained for this neuron, after
finding the average tuning curve across 100 trials and performing summation
to 100ms resolution as done in class (NewTuningCurve) (10 points)
Next perform 100 iterations of shuffling analysis on this spikematrix 
(hint=use a for loop with the makeSparsity function inside, and store all
outputs of the shuffled spikematrices in a vector) (15 points)
Compare the actual spikematrix sparsity with the 100 shuffled spikematrix
sparsities. Create a new vector called  
Sparsities=["100 shuffled sparsities";"Actual sparsity"]. Hence it would
have 101 entries with the last entry being the actual sparsity. 
use "zscore" function in MATLAB and obtain a zscore or a measure of
signficance for the actual sparsity(ie last entry in the vector 
"Sparsities" above). (40 points)
%}

clear all
T=5; % Time for 1 trial in secs
N=100; % This is the number of trials
Dt=1e-3; %Dt is the discretization we are using, 10^-3 means 1ms intervals
Timepoints=T/Dt; % Total time points 
TimeVector=Dt:Dt:5;

FR=ones(1,5000)*2;
Equispaced_Angles=linspace(0,pi,1000);
FR(2001:3000)=2 + 10*sin(Equispaced_Angles);

figure;
plot(TimeVector, FR);
xlabel('Time in seconds')
ylabel('FR')

Lambda=FR*Dt;
Lambda_AllTrials=repmat(Lambda,100,1);
Rand=rand([Timepoints,N]);
Spikes=zeros([Timepoints,N]);
Spikes(Rand<transpose(Lambda_AllTrials))=1;
figure;
subplot(3,1,1)
[Index_time,Index_trial]=find(Spikes==1);
plot(Index_time,Index_trial,'.','markersize',1)
TuningCurve=sum(Spikes');
subplot(3,1,2);
plot(TimeVector,TuningCurve);
xlabel('Time in seconds')
ylabel('Total Spikes in 100 trials')

NewTimeVector=TimeVector(100:100:end);%Take only every 100th datapoint
for i=1:Timepoints/100
NewTuningCurve(i)=sum(TuningCurve((100*(i-1)+1):100*i));
end
subplot(3,1,3)
plot(NewTimeVector,NewTuningCurve);
xlabel('Time in seconds')

temp1=makeSparsity(NewTuningCurve);
fprintf(['Tuned spike matrix sparsity= ' num2str(temp1) '\n'])

for i=1:100
ShuffledSpikeMatrix=ShuffleTrials(Spikes);

TuningCurve_shuffled=sum(ShuffledSpikeMatrix'); 

NewTimeVector=TimeVector(100:100:end);
for i=1:Timepoints/100
NewTuningCurve_shuffled(i)=sum(TuningCurve_shuffled((100*(i-1)+1):100*i));
end

temp=makeSparsity(NewTuningCurve_shuffled);
Sparsities(i)=temp;
end
Sparsities(101) = temp1;

Z_Actual_sparsity = zscore(Sparsities(101));
fprintf(['zscore of actual sparsity= ' num2str(Z_Actual_sparsity) '\n'])

%% lastly, let us quickly look at bandpassing the LFP signal in theta range
%{
clear all
load('SampleLFP_Speed.mat')
Freq=1e6/median(diff(Time));%In Hz
%The time vector "Time" is provided in microseconds, hence we scale by 1e6 to
%obtain frequency in Hz above
theta_min=4;%frequency range of theta band 
theta_max=12;%upper and lower bounds, in Hz

%Next we design a butterworth filter
[q,w] = butter(2,2*[theta_min theta_max]/Freq,'bandpass');

filtered_LFP = filtfilt(q,w,(LFP_raw));
Theta_power=abs(hilbert(filtered_LFP));

%Next we are going to bin the running speed on a logarithmic scale from
%1cm/sec to 100cm/sec
%Plot the histogram of Running_speed, what we see is that there is a lot of
%occurences of low speed, and fewer occurences of high speed.
SpeedBinsEdges=linspace(0,100,50);%Function linspace will now create 50 linearly sampled
%spaced bins starting from 0 to 100 

for i=1:(length(SpeedBinsEdges)-1)%We have 50 speedbinEdges as above, so we will compute 
    %median theta power in each of these ranges, leading to 49 values,
    %since 50 binsedges correspond to 49 bins
Indices=RunningSpeed_cm_sec<SpeedBinsEdges(i+1) & RunningSpeed_cm_sec>SpeedBinsEdges(i);
%We first find the indices corresponding to datapoints when the running
%speed is in the given range
MedThetaPower(i)=median(Theta_power(Indices));%We find the median theta power for all
%datapoints when the running speed was in this particular range
BinCenter(i)=(SpeedBinsEdges(i+1)+SpeedBinsEdges(i))/2;% For each of plotting later, we also 
%compute the bin center for each of these speed bins
end

figure;plot(BinCenter,MedThetaPower);
%}