function Define_Global_Vars
% This defines variables used in all functions and saves them to a file for
% easier access later

% Global variable structure
glbl = struct;

% Init
glbl.dir_top = 'C:\Research\Projects\PhaseMatch\999_FinalV2_0721\';
addpath([glbl.dir_top 'Functions']);

% Strings
glbl.dir_data     = [glbl.dir_top 'Data\'];
glbl.dir_info     = 'F:\Datasets\PODAAC\';
glbl.file_bath    = 'F:\Datasets\HYCOM\depth_GLBc0.04_27.a';
glbl.file_dshore  = 'F:\Datasets\GEBCO\GEBCO_gridone_distance_to_coastline_v2.0m.nc';
glbl.file_mdlgrid = 'F:\Datasets\HYCOM/HYCOM25_Grid.mat';
glbl.file_basins  = 'F:\Datasets\NOAA\range_area.msk';
glbl.dir_mdl      = 'F:\Datasets\HYCOM\PhaseMatch\';
glbl.dir_altm     = 'F:\Datasets\Altimetry\Jason_StandardCorrections_fromEd\';
glbl.str_mdl      = 'DA_Modern';                             % HYCOM output labels
glbl.str_altm     = 'j3';                                    % Altimeter mission labels
glbl.str_tide     = {'M2','S2','N2','K1','O1'};              % Tidal constituent labels
glbl.str_astroarg = 'ray';                                   % Astronomical argument source
glbl.str_raytide  = 'F:\Datasets\Altimetry\RayOld\trimmed.mat';  % Richard Ray's stationary tide model
glbl.str_ztide    = 'F:\Datasets\Altimetry\Zaron2019\HRET_v8.1_compressed.nc';

% Numerical
glbl.rE            = 6371;           % Radius of Earth (km) for distance-along-track calculations
glbl.dist_lin_step = 6.6;            % Step size (km) of linearly-spaced sample points 
glbl.altm_tref     = datenum('1985-01-01');      % Altimeter data is referenced to this time
glbl.tide_tref     = datenum('2016-01-01');      % Tidal phases will be referenced to this time
glbl.itide_zaron   = [1 2 4 5];                  % Indices of tidal constituents in Zaron 2019 baroclinic tide model
glbl.ibasin_land   = [1 3 5 7 9 11 13 15 17 19]; % These indices in file_basins are coastal/land points
glbl.mdl_fac       = 100/9.8;   % Convert HYCOM ssh values to cm
glbl.time_lim      = [736695.999 737791];   % Time range for altimeter missions

% "Valid" altimeter point parameters
glbl.bath_min      = 1500;   % Minimum bathymetric depth
glbl.dshore_min    = 12;     % Minimum distance from shore
glbl.seglength_min = 30; %1200;   % Segments must be at least this long
glbl.gap_max       = 30;     % Minimum size of gaps within segments
glbl.basin_incl    = [2 4 6 8 10 12 14 16 18 20 21]; % Valid basin indices (excluding seas, etc.)
% Jay does butterworth filter on available points (no nanfill or something to make data evenly-spaced)
% Also Jay fits with x=index not x=along-track-distance
% Brian uses 2 butterworth filters, (2,0.02) and (2,0.32)

% Options for boxes (longitude range, latitude range, basin indices, name)
glbl.box_info = {[140  240],     [ 40  60],  [8 10 12],  'N  PACIFIC';   ...
            [120  170],     [ 05  40],  [8 10 12],  'LUZON';        ...
            [170  210],     [ 05  40],	[8 10 12],  'HAWAII';       ...
            [210  280],     [ 05  40],  [8 10 12],  'NW PACIFIC';   ...
            [120  280],     [-05  05],  [8 10 12],  'EQ PACIFIC';  ...
            [150  180],     [-30 -05],  [8 10 12],  'TAHITI';       ...
            [180  230],     [-30 -05],	[8 10 12],  'SC PACIFIC';	...
            [230  290],     [-30 -05],  [8 10 12],  'SW PACIFIC';   ...
            [150  290],     [-50 -30],  [8 10 12 14 18 2 6],  'S PACIFIC';    ...
            [ 40  100],     [-05  30],  [14 16 18], 'N INDIAN';     ...
            [ 30   70],     [-30 -05],  [14 16 18], 'MADAGASCAR';   ...
            [ 70  130],     [-30 -05],  [14 16 18], 'W INDIAN';     ...
            [ 20  150],     [-50 -30],  [14 16 18 2 6 8 12], 'S INDIAN';     ...
            [300  360],     [ 50  66],  [2 4 6],    'N ATLANTIC';   ...
            [280  320],     [ 20  50],  [2 4 6],    'GULF STREAM';  ...
            [320  360],     [ 20  50],  [2 4 6],    'NW ATLANTIC';  ...
            [280  380],     [-20  20],  [2 4 6],    'EQ ATLANTIC';  ...
            [290  380],     [-50 -20],  [2 4 6 14 18 8 12],    'S ATLANTIC';   ...
            [ 20  380],     [-66 -50],  [20 21],    'ANTARCTIC';    ...
           };   
       
% Tidal frequencies IN RADIANS PER SECOND  
glbl.tide_omega = [1.932274 2.000000 1.895982 1.002738 0.929536].*(2*pi);    % UNITS OF RADIANS PER DAY!

% Tidal arguments
[h0,s0,p0,N0,~,~]=ray_arguments(year(glbl.tide_tref),day(glbl.tide_tref));
glbl.tide_chi = [2*h0-2*s0,        0, 2*h0-3*s0+p0,     h0+90,                h0-2*s0-90          ];	% Ephemerides 'chi' (in degrees)
glbl.tide_f   = [1-0.037*cosd(N0), 1, 1-0.037*cosd(N0), 1.006+0.115*cosd(N0), 1.009+0.187*cosd(N0)]; % Amplitude nodal factors 'f' (unitless)
glbl.tide_nu  = [-2.1*sind(N0),    0, -2.1*sind(N0),    -8.9*sind(N0),        10.8*sind(N0)       ]; % Phase nodal factors 'nu' (degrees?)

% Save
save([glbl.dir_top 'global.mat'],'glbl');

end