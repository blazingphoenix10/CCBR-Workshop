function [ SpikeMatrix_Shuffled ] = ShuffleTrials( SpikeMatrix )
%Shuffle trials with random shifts in each trial to obtain a null
%distribution.
%Let us say that the response of the neuron is "tuned" i.e. it responds in
%non-uniform fashion with higher firing rate at a certain time and
%lower,baseline firing at other times, in all trials. If this tuning is
%stable, i.e. occurs repeated and reliably in each trial, then in all
%trials, the firing rate would be highest at the same location(or time
%points) and lower at others. To quantify this alignment of high firing at
%the same location, we will shuffle each trial by some random amounts, so
%that the overall firing rate is the same, the small-timescale statistics
%of spike train is unaffected, but we are "simulating the null hypothesis"
%viz across trials, high firing occurences are NOT aligned.

[r,c]=size(SpikeMatrix);%rows correspond to trials, and columns correspond to 
% time points in a given trial. We would expect c to be same as "N" which
% was the number of trials we started with

% One can use debugging points inserted using F12 to pause the execution of
% code at certain lines, to examine outputs.(refer to the red/grey dot on the left)
SpikeMatrix_Shuffled=SpikeMatrix;

for i=1:c%Repeat the same logic for each trial.
    RandNumber=randi(r);%Create an random integer which is the amount of shift 
    %we would be introducing
SpikeMatrix_Shuffled(:,i)=circshift(SpikeMatrix(:,i),RandNumber);
%Circshift ensures that the dimension of the matrix does not change, and
%the new, shifted vector is wrapped around
end

end

