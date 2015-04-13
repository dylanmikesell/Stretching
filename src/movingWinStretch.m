function [C,epsArray,tSamp] = movingWinStretch(u0,u1,dt,twin,tstep,dVmax,dV)
%
% USAGE: [C,epsArray,tSamp] = movingWinStretch(u0,u1,dt,twin,dVmax,dV)
%
% INPUT:
%   u0    = reference trace
%   u1    = perturbed trace
%   dt    = sample interval [s]
%   twin  = length of analysis window [s]
%   tstep = time step of windows [s]
%   dVmax = maximum epsilon value to test (i.e. 0.5 = 50%)
%   dV    = sampel interval for epsilon vector (i.e. epsilon = -dVmax:dV:dVmax)
% OUTPUT:
%   C        = maximum correlation coefficient in each window
%   epsArray = epsilon value for each window corresponding to maximum correlation coefficient in C
%   tSamp    = time vector containing the center time of each analyzed window
%
% This function computes the moving window stretching method to estimate
% velocity changes.
%
% see also corrcoef.m, interp1.m
% 
% Written by Dylan Mikesell (mikesell@mit.edu)
% Last modified 13 April 2015
%
% Guts of algorithm taken from an example by Berenice Froment

if nargin < 4
    twin = 1;
    fprintf('Window length not set (default set to %0.2f (s))\n',twin);
elseif nargin == 4
    fprintf('Window length set to %0.2f (s).\n',twin);
end

if numel(u1) ~= numel(u0)
    error('MATLAB:accumMovingWinStretch','traces are not the same length');
end

npts  = numel( u1 );         % number of data points in the traces
nstep = round( tstep / dt);  % number of samples to step through windows
nwin  = round( twin / dt );  % number of samples in the window
nwin2 = floor( nwin / 2 );   % half the sample window
eps   = -dVmax : dV : dVmax; % dV array to test
neps  = numel( eps );        % number of dV values to test

nMeas = numel( nwin2 + 1 : nstep : npts - nwin2 - 1 ); % number of measurements

ccArray = zeros( neps, nMeas ); % allocate correlation coefficient vector
tSamp   = zeros( 1, nMeas ); % allocate window center time vector

for ii = 1 : neps % test each epsilon value
    
    fprintf( 'Testing epsilon = %0.4f\n', eps(ii) );
    
    % stretch the reference trace for the given eps_vec
    tPrime =  ( 1 : npts ) .* ( 1 + eps(ii) ); % the stretched time axis
    u0int  = interp1( 1 : npts, u0, tPrime, 'pchip', 0 ); % stretched trace
    
    cnt = 0; % counter for windows
    
    for jj = nwin2 + 1 : nstep : npts - nwin2 - 1 % loop through windows 

        cnt = cnt + 1; % next window -> update counter
        tSamp(cnt) = (jj-1) * dt; % [s] center of window  
        
        startIdx = jj - nwin2; % window starts half width before center of window
        stopIdx  = jj + nwin2 - 1; % window ends half width after center of window
               
        u0tmp = u0int( startIdx : stopIdx ); % locally windowed trace
        u1tmp = u1(    startIdx : stopIdx ); % locally windowed trace
        
        % zero lag correlation of this window
        toto = corrcoef( u0tmp, u1tmp ); % equation 6 in Hadziioannou (2012).
        ccArray( ii, cnt ) = toto( 2 ); % take cross correlation coefficient
        
    end
    
end

% get the max dV within in each window
[C, CCidx] = max( ccArray, [], 1 );
epsArray = eps( CCidx );

return