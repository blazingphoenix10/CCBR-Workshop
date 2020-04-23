%% 
%  1/7/2019 CCBR MATLAB SESSION 
%%

% part 0 MATLAB as a calculator
% Using the + - exp(x) sin cos functions
% Using plot, figure and subplot functions
% Create a vector, create a matrix and access its entries
% Read up the "rand" function to understand how functions work
% Create a new script and use %% to create and run different subsections.

clear all

%Part 1 Creating vector of 0 and 1 using random numbers
R=rand(10000,1);

figure; 
hist(R,50)
Lambda=0.1;
Binary=R<Lambda;
Count1=sum(Binary==1);
Count0=sum(Binary==0);

Lambda=0.9;
Binary=R<Lambda;
Count1=sum(Binary==1);
Count0=sum(Binary==0);

%%
%Part 2 Code to create a spike train
T=200; % Total Time for simulation in secs
Dt=1e-3; %Dt is the discretization we are using, 10^-3 means 1ms intervals
Timepoints=T/Dt; % Total time points created would be 200 thousand, since each second is 1000 data points
%Since each time point is binary- 0 or 1 the maximum firing rate we can
%have is 1000Hz, which is atleast an order of magnitude more than typical
%neurons, so we are fine
FR=5; %We set the desired firing rate of the neuron as 5Hz, ie 5 spikes per second ie 5 spikes per 1000 data points
Lambda=FR*Dt; %Lambda, which is a threshold or spike probability %at any datapoint.
Rand=rand([Timepoints,1]); %Create a matrix of random numbers
Spikes=zeros([Timepoints,1]); % Initiate a matrix with all zeros, which are spikes
Spikes(Rand<Lambda)=1; % We are going to consider spikes to be binary, 
%and when the random number exceeds the threshold, 
%we switch "Spikes" at that time to 1
Observed_FR=sum(Spikes==1)/T;%This will have the units of spikes per second, or events per
%Second, hence in Hz

x_1 = zeros(1,100)
y_1 = zeros(1,100)
i=1
for FR = 1:100
    Lambda=FR*Dt; %Lambda, which is a threshold or spike probability %at any datapoint.
    Rand=rand([Timepoints,1]); %Create a matrix of random numbers
    Spikes=zeros([Timepoints,1]); % Initiate a matrix with all zeros, which are spikes
    Spikes(Rand<Lambda)=1; % We are going to consider spikes to be binary, 
    %and when the random number exceeds the threshold, 
    %we switch "Spikes" at that time to 1
    Observed_FR=sum(Spikes==1)/T;%This will have the units of spikes per second, or events per
    %Second, hence in Hz
    SpikeIndices=find(Spikes==1); % Find function will give us the indices of the vector Spikes which are 1, ie when spike happened
    SpikeTimes=SpikeIndices*Dt; % Since 1000 indices correspond to 1 second, we scale the indices by Dt, to get time of spikes in sec
    ISI=diff(SpikeTimes); % the diff function will give us the difference between consencutive elements, hence the vector ISI holds the
    x_1(i) = mean(ISI)
    y_1(i) = std(ISI)
    i= i+1
end

figure;loglog(x_1,y_1,'.')
xlabel('Mean(ISI)')
ylabel('STD(ISI)')
%Comments and future work
%Self-Study. Change lambda, to see how firing rates change.

%% ISI histogram
SpikeIndices=find(Spikes==1); % Find function will give us the indices of the vector Spikes which are 1, ie when spike happened
SpikeTimes=SpikeIndices*Dt; % Since 1000 indices correspond to 1 second, we scale the indices by Dt, to get time of spikes in sec

ISI=diff(SpikeTimes); % the diff function will give us the difference between consencutive elements, hence the vector ISI holds the
% difference in seconds, between 2 consecutive spikes. Hence we get the
% Inter spike intervals.
figure; [y,x]=hist(ISI,20);plot(x,y)
xlabel('Time in ms')
ylabel('Counts')
    
%% Autocorrelation
% We will use the inbuilt xcorr function in MATLAB
% Use CTRL+D on the function name to read up more details, or Gooogle 
[correlation,lags]=xcorr(Spikes,'coeff');
plot(correlation,lags)
%A function in MATLAB takes inputs, in this case the Spikes vector and the 
%String 'coeff' which determines the normalization. The quantities inside
%the square brackets, to the left of "=" are the outputs of this function. Correlation is the
%value of cross correlation, normalized such that the correlation at lag zero is
%set to unity. lags is a vector of the amount of shift used to compute the
%corresponding "correlation". Hence verify that their dimensions are the
%same and value of correlation corresponding to lag=0 is 1.
% If using the latest trial version of MATLAB xcorr function might be
% unavailable to you, in that case use the xcorrCCBR function provided 
%% Let us now create another spike train which is periodic
% let us create a spike train for 20 seconds, which fires every 20 ms at a
% resolution of 1ms.
T=20;
Dt=1e-3; % Let us stick to the units of seconds, so 1ms=10^-3 seconds
Spikes2=zeros(T/Dt,1);
%Next we would not be using the Poisson properties but artificially making a
%spike every 20ms ie after 20 data points each, with the first spike at
%t=1ms
Spikes2(1:20:end)=1;
%We are asking MATLAB to start at index =1 go ahead in jumps of index=20
%and keep doing this till it reaches the end of the Spikes2 vector. For all the values
%corresponding to these indices we want the value stored in Spikes to
%become 1 instead of the zero we started with.

%Q1(5 points) How many spikes would you expect in this spike train of 20seconds?
%Q2(10 points) How would you check that?(function/code in MATLAB)
%Q3(15 points) What the mean firing rate for this spike train? Store this value in a
%variable called MFR.

%Now we will compute the autocorrelation of this spike train, what do you
%expect to see?
[correlation2,lags2]=xcorr(Spikes2,'coeff');
figure;plot(lags2,correlation2)
xlabel('Lag')
ylabel('normalized autocorrelation')
%xcorr function with a single vector as 
%input computes the autocorrelation, i.e. the xcorr of a vector with itself
%By default, the lags2 would hold values of indices corresponding to the
%shift, but we would like to convert it to time. Hence we should scale it
%by Dt to get lags_Time in the units of seconds

lags_Time=lags2*Dt; %Each timepoint corresponds to Dt amount of time
figure;plot(lags_Time,correlation2)
xlabel('Time in seconds')
ylabel('normalized autocorrelation')

%Why do we see a traingle?---------------- Finite data effect

%Why does the x-axis range from -20 to +20?

%Let us zoom in around +-100ms, how many peaks would you expect therein?
xlim([-0.1 0.1])

% Now lets quickly create another spike train, with the same properties
% but whose firing starts at the 11th ms.
T=20;
Spikes3=zeros(T/Dt,1);
Spikes3(11:20:end)=1;

%We are going to compute the cross correlation between these spike trains,
%what would you expect to see?
[cross_corr,lags_cc]=xcorr(Spikes2,Spikes3,'coeff');
figure;plot(lags_cc,cross_corr)
xlabel('lag')
ylabel('normalized autocorrelation')

%% HW problems
%{
Q4.(30 points) Given the actual data from a hippocampal neuron, compute the ISI histogram and
autocorrelation. Load this data using the command
load('HW1Data_SpikeTrain.mat'); make sure that you "set path" to the folder containing 
 %this file. Alternatively, use function cd('FolderLocationForTheFile')
Start with the Spiketimes (ST) and the time bins(T_bc, in seconds). Compute ISI
histogram of the Spiketimes. Next use, "hist" function to compute Spike counts in each time bin. Then compute autocorrelation. What is the
time bin(Spike counts are provided as vector SC to cross check )? 
%How long was this recording made? If the time bins are in seconds, 
%at what resolution(in seconds) is the experiment being recorded?(write code to obtain this)
Q5. (40 points) Create 2 spike trains, at a resolution of 1ms for total
duration of 5second. Let the first spike train fire every 500ms, and the second one
fire every 700ms(both have their first spike at t=1ms). 
% What would you expect to see in the autocorrelation of each of the spike trains?
% Use "xcorr" function and verify your expectations. Next use the xcorr function 
%to compute the cross-correlation between these spike trains. Comment in 5 lines about what you see, and how
& why it is expected or unexpected.
%}
    