clear all
close all
clc

indir  = '/Users/dmikesell/Dropbox/Documents/research/CodaWaves/SUFDMOD_dt_Test/single_traces';

a      = dir([indir '/*.mat']);
ntrc   = numel(a);

for ii = 1:ntrc
    
    fprintf( 'Loading %s\n', a(ii).name );
    b = load([indir '/' a(ii).name],'-mat');
    
    dt         = b.dt;
    tmax       = b.tmax;
    dv(ii)     = b.dv;
    npts       = b.npts;
    data(:,ii) = b.trc;
    tvec       = b.tvec;

end
clear b;

u0 = data(:,2); % 6 km/s
u1 = data(:,1); % 5.94 km/s

save('./traces.mat','u0','u1','dt');