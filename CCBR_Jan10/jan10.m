% B P Kailash - BE15B007


clear all

load("Data1D.mat");

ses = 1;
ID = 1;

%%
% Self study: Change time vector to seconds and graph it.

SpikeTimes = (Data1D(ses).UIF(ID).ST);
RatsPosition = Data1D(ses).beh(ID).Position;
Direction = Data1D(ses).beh(ID).Forward;
TimeVector = Data1D(ses).beh(ID).Time;
Dt = median(diff(SpikeTimes))/1e6; % Since time vector is in ms we divide by 1e6

SpikeIndices = nearestpoint(SpikeTimes, TimeVector);
% Nearest point finds the indices of TimeVector closest to SpikeTimes

x = RatsPosition(SpikeIndices);
y = TimeVector(SpikeIndices);
y2 = TimeVector(SpikeIndices)/1e6;
% y and SpikeTimes should have very less difference between them why?
% Difference should be less than resolution of time vector if the nearest
% point algorithm works as expected.

% Now we take into consideration the direction, it's boolean vector forward
% 1 if forward, 0 if backward
forwardRunningSpikes = Direction(SpikeIndices);

figure;
subplot(1,3,1);
plot(x,y2,'.');
title("Both");
xlabel("Position in cm")
ylabel("Time in seconds");
subplot(1,3,2);
plot(x(forwardRunningSpikes), y2(forwardRunningSpikes),'.');
title("Forward");
xlabel("Position in cm");
ylabel("Time in seconds");
subplot(1,3,3);
plot(x(~forwardRunningSpikes), y2(~forwardRunningSpikes),'.');
title("Backward");
xlabel("Position in cm");
ylabel("Time in seconds");

% Here we note that most of the spikes at around 100 cm are when the mouse
% is moving towards the 220 cm end i.e moving in the forward direction

%%
% Next step is getting the tuning curve
% Firing rate is the number of spikes / Time taken

ls = linspace(min(RatsPosition),max(RatsPosition),26);
% Why 26? 25 bins of 9 cm because the experiment is for 225cms
Occupancy = hist(RatsPosition, ls)*Dt;
% plot(Occupancy)
% Scaling above by amount of time taken by each neuron
SpikesBinned = hist(RatsPosition(SpikeIndices),ls);
% plot(SpikesBinned)
TuningCurve= ls./Occupancy; % ./ works to compute bin-by-bin division
%So we are finding number of spikes divided by occupancy in each bin
figure;plot(ls,TuningCurve)
xlabel('Position (cm)')
ylabel('Firing rate (Hz)')

%% 

%Q1. Similar to this analysis, load actual brain data, and compute theta
%phase modulation in 30 bins of 18 degrees each for theta phase
%Theta phases are stored in Data1D(ses).LFP.ThetaPhase. Comment in 2 lines
%about what you see. Seperate Forward and return trials in different
%subplots. (25 points)
ses = 1;
ID = 1;

load("Data1D.mat");
Theta = Data1D(ses).LFP.ThetaPhase;
SpikeTimes = (Data1D(ses).UIF(ID).ST);
RatsPosition = Data1D(ses).beh(ID).Position;
Direction = Data1D(ses).beh(ID).Forward;
TimeVector = Data1D(ses).beh(ID).Time;
Dt = median(diff(SpikeTimes))/1e6; % Since time vector is in ms we divide by 1e6

SpikeIndices = nearestpoint(SpikeTimes, TimeVector);
x = RatsPosition(SpikeIndices);
y = Theta(SpikeIndices);

forwardRunningSpikes = Direction(SpikeIndices);

figure;
subplot(1,3,1);
plot(x,y,'.');
title("Both");
xlabel("Postion")
ylabel("Theta Phase (degrees)");
subplot(1,3,2);
plot(x(forwardRunningSpikes), y(forwardRunningSpikes),'.');
title("Forward");
xlabel("Postion")
ylabel("Theta Phase (degrees)");
subplot(1,3,3);
plot(x(~forwardRunningSpikes), y(~forwardRunningSpikes),'.');
title("Backward");
xlabel("Postion")
ylabel("Theta Phase (degrees)");

ls2 = linspace(-(3/2)*pi,(3/2)*pi,30);
Dt2 = median(diff(Theta));
Occupancy = hist(Theta, ls2)*Dt2;
SpikesBinned = hist(Theta(SpikeIndices),ls2);
TuningCurve= SpikesBinned./Occupancy;
figure;plot(ls2,TuningCurve)
xlabel('Theta Phase')
ylabel('Firing rate (Hz)')

%%

%Q2 Compute sparsity for the 2 tuning curves of place tuning obtained in
% the section above(for forward and backward direction seperately)
%Also compute the sparsity for theta phase modulation.(10 points)

disp("Sparsity of forward running spikes")
makeSparsity(forwardRunningSpikes)
disp("Sparsity of backward running spikes")
makeSparsity(~forwardRunningSpikes)
disp("Sparsity of Theta phase")
makeSparsity(Theta)

%%
%Q3. Use the data given, Data2D.mat,to load a session with 7 units recorded
%during a rat running on a 2D table. Plot out 2 place cells, as a scatter
%plot with spikes as red stars (r*) and occupancy/behavior as black dots (k.)
%(25 points)

load("Data2D.mat");
SpikeTimes1 = Data2D.UIF(1).ST;
SpikeTimes2 = Data2D.UIF(2).ST;
TimeVector = Data2D.beh.Time;
SpikeIndices1 = nearestpoint(SpikeTimes1, TimeVector);
SpikeIndices2 = nearestpoint(SpikeTimes2, TimeVector);
Position_x = Data2D.beh.Position(:,1);
Position_y = Data2D.beh.Position(:,2);

x1 = Position_x(SpikeIndices1);
y1 = Position_y(SpikeIndices1);

x2 = Position_x(SpikeIndices2);
y2 = Position_y(SpikeIndices2);

figure;
plot(Position_x, Position_y, "k.")
hold on;
plot(x1, y1, "r*")
hold on;
plot(x2, y2, "r*")
xlabel("x position")
ylabel("y position")
% zlabel("Time in seconds")
title("Plotting behavior(position), red star - spike, behaviour - black dot")

%%
% Q4 Use the histcn and surf function to plot a colored 2D place field
%Note that the occupancy is inside a circle, but you will be creating a 
%rectangular grid for binning though "histcn" function. Hence, make sure 
%that for bins with zero occupancy(like the red regions in the image
%ZeroOccupancyillustration.png), you reset the firing rate to NaN.(30 points)

xbins = linspace(-100, +100, 25);
ybins = linspace(-100, +100, 25);
mat = histcn([Position_x Position_y;], xbins, ybins);
mat1 = histcn([x1 y1;], xbins, ybins);
mat2 = histcn([x2 y2;], xbins, ybins);
mat(mat == 0) = NaN;
% surf(mat1, mat2)
precession1 = mat1 ./ mat;
precession2 = mat2 ./ mat;

figure;
subplot(1,2,1)
surf(precession1)
title("Precession for Place cell 1")
xlabel("xbins")
ylabel("ybins")
zlabel("Firing rate")
subplot(1,2,2)
surf(precession2)
xlabel("xbins")
ylabel("ybins")
zlabel("Firing rate")
title("Precession for Place cell 2")
%%
% Q5. Read Nature Neuro paper provided in the zip file,
% "NatureNeuro2014_MehtaLab.pdf" and describe in your own words(7 lines) 
%the idea of a "Motif" and its use to compute phase precession.  (10 points)

%{
A Motif is a neuronal firing pattern whose firing rate is sustained beyond
a threshold based on the peak firing rate for a pre-defined period of time.
 
Theta phase of spikes in a motif field is indicative of a animal's
location, this is known as phase precession.
%}