clear all
close all
clc

addpath('src/');

% this script runs an example of the moving window stretching method

% u0 = reference trace at 6 km/s
% u1 = perturbed trace at 5.94 km/s (a -1.0% perturbation)

%% load the data and plot section

load('exampleData/traces.mat');

npts = numel( u0 ); % number of points in traces

tArray = ( 0 : npts - 1 ) .* dt; % time vector

figure;
plot(tArray,u0); hold on;
plot(tArray,u1);
xlabel('Time [s]'); ylabel('Amplitude [a.u.]');
legend('u0','u1'); legend boxoff;
xlim([1 3]);
ylim([-5e-2 5e-2]);

%% do moving window stretching

winLength = 0.5;     % [s] length of moving window
tStep     = dt * 10; % [s] make a measurement every 'tStep'

dVmax     = 0.02;    % set maximum for stretch parameter search (0.5 = 50%)
dV        = dVmax/6; % sample interval from -dVmax:dV:dVmax for epsilon values

[ ccArray, dtot, tSamp ] = movingWinStretch( u0, u1, dt, winLength, tStep, dVmax, dV );

% tSamp is the array that gives the window center and changes with tStep

figure;
subplot( 2, 1, 1 )
plot( tSamp, ccArray ); ylabel('Corr. Coeff.'); ylim([0.9 1]);
subplot( 2, 1, 2)
plot( tSamp, dtot ); ylabel('\epsilon'); ylim([-dVmax dVmax]);
xlabel('Time [s]');