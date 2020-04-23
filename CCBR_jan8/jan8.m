% 8th Jan CCBR tutorial
% Histogram 10ms bin

FR = 5; % Firing rate - 5Hz
Dt=1e-3; %Dt is the discretization we are using, 10^-3 means 1ms intervals
T=5; % Total Time for simulation in secs
Timepoints=T/Dt; % Total time points created would be 200 thousand, since each second is 1000 data points
N = 100 % Number of trials
i = 1
spikeMatrix = zeros(Timepoints,N) % Dimension of this matrix is now 5000 X 100

for trials = 1:N
    Lambda=FR*Dt; %Lambda, which is a threshold or spike probability at any datapoint.
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
    spikeMatrix(:,i) = Spikes;
    
    % Fig 1: x axis: Time, y axis: Trial number, A dot if neuron fires in that trial
    subplot(3,1,1)
    plot(SpikeIndices, i, '.'); hold on
    xlabel('Time')
    ylabel('Trial Number')  
    i = i+1;
    
end

%%

