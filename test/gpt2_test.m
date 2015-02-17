dmjd = 56141.d0;
% Latitude [-90,90]
dlat = [48.20e0; 89.20e0; -89.20e0];
% Longtitude [-pi:pi] or [0:2pi]
dlon = [16.37e0; 16.37e0; 16.37e0];
% Height
hell = [156.0e0;  156.0e0;  156.0e0 ];
% NR of stations
nstat = 3;

it = 0;

dlat.*pi/180.d0;
dlon.*pi/180.d0;

[p,T,dT,e,ah,aw,undu] = gpt2 (dmjd,dlat,dlon,hell,nstat,it)