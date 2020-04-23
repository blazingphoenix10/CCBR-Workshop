%% Load actual brain data, and compute a place cell tuning
clear all
load('Data1D.mat');
%{
This has 2 sessions worth of data, and each session has the following
information-

beh- holds the time vector(in microseconds) and behavior of the subject like 
position( in cm, on a linear track, 0 to +225) on the
1D field and running speed(which could be derived by using
diff(position) and time). Also is provided a boolean called "Forward" which
is 1 when the rat is running from 0 to +225 ie running forward. 0
otherwise

LFP- Time and raw LFP signal around 1 recording wire aka electrode.

UIF- single unit information. ST= timings of individual spikes from a
particular neuron. in microseconds
(ignore the other sections for brevity)
%}
% ses=1,ID=1 is directional place field, ses=2,ID=1 is bidirectional
ses=1;% We will only work on the first session for now.
ID=1; %Lets find the spatial tuning of first unit recorded

SpikeTimes=Data1D(ses).UIF(ID).ST;%Note how data inside a structure is accessed,
%with A.B indicating the B subfield inside A, and then referencing with
%numbers in parenthesis for multiple subfields
TimeVector=Data1D(ses).beh.Time;
Dt=median(diff(TimeVector))/1e6;%since time vector is in units of microseconds
Position=Data1D(ses).beh.Position;%is in cm

figure; 
SpikeIndices=nearestpoint(SpikeTimes,TimeVector);%Nearestpoint works as its name suggests
% It finds the indices in TimeVector corresponding to SpikeTimes which are
% closest to it.
x=Position(SpikeIndices);
y=TimeVector(SpikeIndices);%Note that vector y and vector SpikeTimes should have very
%small error between them, and this error should be <= the
%resolution of recording TimeVector (why?)
plot(x,y,'.')
% So what we see is a neuron which fires heavily around 100cm, ie near the
% center of the track, and also again at 1 end but not the other. Do these
% spikes come when he is running towards +220cm or when he is returning
% from there?
%%
%Now let us consider if the spikes happened when the rat was going left to
%right or was coming back, ie classify based on the boolean Forward
Forward=Data1D(ses).beh.Forward;
figure; 
subplot(1,2,1)

forwardRunningSpikes=Forward(SpikeIndices);%We would need the information of whether
%a given spike occured when the rat is running forward or coming back, so
%we store that as another boolean forwardRunningSpikes.

plot(x(forwardRunningSpikes),y(forwardRunningSpikes),'.')

subplot(1,2,2)
plot(x(~forwardRunningSpikes),y(~forwardRunningSpikes),'.')
%Self-study, change timevector to be in seconds, rather than the
%microseconds it is right now, and label the axis
%% Next, we want to compute a tuning curve as we had done for the PSTH earlier
%Note now that for a rat running on a linear track, the occupancy of time
%at the ends of the track is not the same as that at the center of the
%track (why?- rat runs quickly through the center, but sits down and recieves rewards
% at the ends of the track)
% So our tuning curve would be such that 
%Fr_x=(spikes at location x)/(amount of time spent in location x);%Hence
%will be in Hz as expected

SpatialBins=linspace(min(Position),max(Position),26);

%Spatial corresponds to "space" ie spatial bins are bins for space, in this
%case space is just 1D, but for your Homework it would be 2D
Occupancy=hist(Position,SpatialBins)*Dt;
SpikesBinned=hist(Position(SpikeIndices),SpatialBins);
TuningCurve=SpikesBinned ./ Occupancy; % ./ works to compute bin-by-bin division
%So we are finding number of spikes divided by occupancy in each bin
figure;plot(SpatialBins,TuningCurve)
xlabel('Position (cm)')
ylabel('Firing rate (Hz)')

%% 
%HW:
%Q1. Similar to this analysis, load actual brain data, and compute theta
%phase modulation in 30 bins of 18 degrees each for theta phase
%Theta phases are stored in Data1D(ses).LFP.ThetaPhase. Comment in 2 lines
%about what you see. Seperate Forward and return trials in different
%subplots. (25 points)

%Q2 Compute sparsity for the 2 tuning curves of place tuning obtained in
% the section above(for forward and backward direction seperately)
%Also compute the sparsity for theta phase modulation.(10 points)

%Q3. Use the data given, Data2D.mat,to load a session with 7 units recorded
%during a rat running on a 2D table. Plot out 2 place cells, as a scatter
%plot with spikes as red stars (r*) and occupancy/behavior as black dots (k.)
%(25 points)

% Q4 Use the histcn and surf function to plot a colored 2D place field
%Note that the occupancy is inside a circle, but you will be creating a 
%rectangular grid for binning though "histcn" function. Hence, make sure 
%that for bins with zero occupancy(like the red regions in the image
%ZeroOccupancyillustration.png), you reset the firing rate to NaN.(30 points)

% Q5. Read Nature Neuro paper provided in the zip file,
% "NatureNeuro2014_MehtaLab.pdf" and describe in your own words(7 lines) 
%the idea of a "Motif" and its use to compute phase precession.  (10 points)

%% Use spike position and spike phase, to check precession
% Instead of the yaxis being time, we will plot yaxis as theta phase,
% obtained from the LFP by bandpassing in theta range, computing phase with
% angle(hilbert(filtered_LFP)) instead of abs(..) we used in the last
% tutorial to compute magnitude or power of signal in theta band
%Let us also only use the data from running epochs

ThetaPhase=Data1D(ses).LFP.ThetaPhase;
Forward=Data1D(ses).beh.Forward;
SpeedBoolean=Data1D(ses).beh.Speed>10;
figure; 
subplot(1,2,1)
SpikeIndices=nearestpoint(SpikeTimes,TimeVector);%Nearestpoint works as its name suggests
% It finds the indices in TimeVector corresponding to SpikeTimes which are
% closest to it.
x=Position(SpikeIndices);
y=ThetaPhase(SpikeIndices);%Note that vector y and vector SpikeTimes should have very
%small error between them, and this error should be <= the
%resolution of recording TimeVector (why?)
forwardRunningSpikes=Forward(SpikeIndices);
RunningSpikes=SpeedBoolean(SpikeIndices);

plot(x(forwardRunningSpikes & RunningSpikes),y(forwardRunningSpikes & RunningSpikes),'.')

subplot(1,2,2)
plot(x(~forwardRunningSpikes & RunningSpikes),y(~forwardRunningSpikes & RunningSpikes),'.')

%% Convert this into a color plot(if time permits)
PositionBins=linspace(0,225,30);
ThetaPhaseBins=linspace(-pi,pi,50);
%Read the description of the file histcn, it requires the input data as an
%(M x N) array, so we will concatenate the ThetaPhase and Position vectors

%Also note and convince yourself that for precession to work it is CRITICAL 
%that forward and backwards running trials are seperately treated, 
%otherwise the precessions in those 2 directions would negate each other.
OccupancyBoolean=Forward & SpeedBoolean;
Occupancy_2D=histcn([Position(OccupancyBoolean) ThetaPhase(OccupancyBoolean);],...
    PositionBins,ThetaPhaseBins);
%next we will bin the spikes in 2D using the same x and y from above which
%are positions at which spikes happened(x) and thetaphase at which spikes
%happened(y)
SpikesBoolean=forwardRunningSpikes & RunningSpikes;%use running data only, in forward
Spikes2D=histcn([x(SpikesBoolean) y(SpikesBoolean);],...
    PositionBins,ThetaPhaseBins);

ColoredPrecession=Spikes2D ./ Occupancy_2D;
figure;
subplot(1,2,1)
imagesc(PositionBins,ThetaPhaseBins,ColoredPrecession)


%Now literally copy paste the code and use backward direction data
OccupancyBoolean=~Forward & SpeedBoolean;%Use ~Forward ie backward running data
Occupancy_2D=histcn([Position(OccupancyBoolean) ThetaPhase(OccupancyBoolean);],...
    PositionBins,ThetaPhaseBins);
%next we will bin the spikes in 2D using the same x and y from above which
%are positions at which spikes happened(x) and thetaphase at which spikes
%happened(y)
SpikesBoolean=~forwardRunningSpikes & RunningSpikes;%use running data only, in forward
Spikes2D=histcn([x(SpikesBoolean) y(SpikesBoolean);],...
    PositionBins,ThetaPhaseBins);

ColoredPrecession=Spikes2D ./ Occupancy_2D;

subplot(1,2,2)%Plot it in the second subplot
imagesc(ThetaPhaseBins,PositionBins,ColoredPrecession)
imagesc(PositionBins,ThetaPhaseBins,ColoredPrecession)