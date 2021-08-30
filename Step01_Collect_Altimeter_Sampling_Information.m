% Collect altimeter standard tracks
% Arin Nelson
% Rewritten 12/12/2020
% V3 07/25/2021
%==========================================================================
clear; clc; close all; addpath('Functions');

% Switches
Switch     = zeros(9,1);
Switch(01) = 0;         % Gather altimeter track & sampling point info
Switch(02) = 0;         % Interpolate bathymetry from HYCOM bathymetry file & generate land-sea mask
Switch(03) = 0;         % Interpolate distance-from-shore from GEBCO
Switch(04) = 0;         % Interpolate basin indices from NOAA
Switch(05) = 0;         % Link model output indices to altimeter points

% Options
do_debug = 1;

% Global vars
load('global.mat');

%==========================================================================
if(Switch(01))
    
	% Loading reference data
	eqcrosses    = importdata([glbl.dir_info 'eqxl_j3.txt']);            % Jason-3
	reftrack_asc = importdata([glbl.dir_info 'reftrk_ascending.txt']);   % Ascending tracks
	reftrack_dsc = importdata([glbl.dir_info 'reftrk_descending.txt']);  % Descending tracks

	% Total number of tracks
	n_track = numel(eqcrosses);    % 254
  
	% Parse data
	lat_asc = reftrack_asc(:,1);
	lon_asc = reftrack_asc(:,2);
	lat_dsc = reftrack_dsc(:,1);
	lon_dsc = reftrack_dsc(:,2);

	% Total number of sample points per track
	n_point = numel(lat_asc);
  
	% Compute track longitudes and latitudes
    lon = zeros(n_track,n_point);
    lat = zeros(n_track,n_point);
    for it=1:n_track
    if(mod(it,2)~=0)
        lon(it,:) = lon_asc + eqcrosses(it);
        lat(it,:) = lat_asc;
    else
        lon(it,:) = lon_dsc + eqcrosses(it);
        lat(it,:)   = lat_dsc;
    end
    end
    clear it;
  
    % Limit points to 0-360
    lon(lon < 0  ) = lon(lon < 0  ) + 360;
    lon(lon > 360) = lon(lon > 360) - 360;
    
    % lon-lat to distance units
    %[x, y ] = grn2eqa(lat,lon,[-90, 0],referenceEllipsoid('wgs84'));
    
    % Along-track distance
    dist = zeros(n_track,n_point);
    for it=1:n_track
        xx = lon(it,:); xx = unwrap(xx.*(pi/180)).*(180/pi);    
        if(max(xx)>360);	xx = xx-360;    end
        if(min(xx)<0);      xx = xx+360;    end
        yy = lat(it,:);
        [x,y] = grn2eqa(yy,xx,[yy(1) xx(1)],referenceEllipsoid('wgs84'));
        dist(it,2:end) = sqrt( (x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2 );
    end
    clear it xx yy x y;
    dist = cumsum(dist,2);
  
    % Debug plot: map of distance-along-track's
    if(do_debug)
        figure('units','normalized','outerposition',[0 0 1 1]);
        	hold on;
                for it=1:n_track
                    scatter(lon(it,:),lat(it,:),2,dist(it,:),'o','filled'); 
                end
            hold off; axis tight; title('Altimeter sampling info');
            colormap(jet); cb=colorbar; set(get(cb,'ylabel'),'string','along-track distance (km)');
        print(gcf,'-dpng','Plots/Debug/Debug_Step01-01_AltimSamplingInfo.png');
        close all;  clear ia it cb;
    end
    
%     % Comparison to Jay's track points
%     tmp = load('F:\Datasets\PODAAC\altim_lonlat_jayshriver.txt');
%     xj  = tmp(:,1);
%     yj  = tmp(:,2);
%     nj  = numel(xj);
%     ij  = find( sign(diff(yj(2:end))) == -sign(diff(yj(1:end-1))) );
%     dj  = zeros(nj,1);
%     for i=2:nj
%         dj(i) = rE*circledist(xj(i-1),yj(i-1),xj(i),yj(i));
%     end
    
    % Save?
    save([glbl.dir_data '\track_info.mat'],'lon','lat','dist');
  
    % Clean-up
    clear eqcrosses* reftrack_* lon* lat* dist

end
%==========================================================================
if(Switch(02))
    
    % Read in HYCOM lon and lat
    tmp = load(glbl.file_mdlgrid,'lon','lat');
    lon_mdl = tmp.lon'; lon_mdl(lon_mdl>360) = lon_mdl(lon_mdl>360)-360;
    lat_mdl = tmp.lat';
    clear tmp

    % Read in HYCOM-bathymetry-derived-from-GEBCO file
    nt  = numel(lon_mdl);
    fid = fopen(glbl.file_bath,'r','b');
    tmp = fread(fid,[1,nt],'real*4');
    bath_gebco = reshape(tmp,size(lon_mdl)); 
    fclose(fid); clear fid tmp nt;
  
	% % Check
	% zz = bath_gebco(5:10:end,5:10:end);   zz(zz>1e10) = NaN;
	% surf(lon_mdl(5:10:end,5:10:end),lat_mdl(5:10:end,5:10:end),zz,'edgecolor','none'); view(2); colormap(jet); colorbar;

    % Generate the interpolant
    [x_mdl,y_mdl] = grn2eqa(lat_mdl,lon_mdl,[0 0],referenceEllipsoid('wgs84'));
    ntrplnt = scatteredInterpolant(double(x_mdl(:)),double(y_mdl(:)),bath_gebco(:),'linear');
    
    % Load info
    info  = load([glbl.dir_data '\track_info.mat']); 
    [x,y] = grn2eqa(info.lat,info.lon,[0 0],referenceEllipsoid('wgs84'));

    % Interpolate bathymetry to sampling points
    bath = ntrplnt(x,y);  
    
    % Set land values to NaN
    maxb = max(bath_gebco(bath_gebco < 1e10));
    bath(bath>maxb) = NaN;
    
    % Generate land-sea mask
    mask = ~isnan(bath);
    
    % Debug plot
    if(do_debug)
        figure('units','normalized','outerposition',[0 0 1 1]);
        	xx = info.lon(:);
            yy = info.lat(:);
            zz = bath(:);
            scatter(xx(:),yy(:),2,zz(:),'o','filled'); axis tight; colormap(jet); caxis([0 7000]); 
            title('Track Point Bathymetry');
            cb=colorbar; set(get(cb,'ylabel'),'string','bathymetry (m)');
        print(gcf,'-dpng','Plots/Debug/Debug_Step01-03_track-Bathymetry.png');
        close all; clear xx yy zz cb;
    end
    
    % Save
    save([glbl.dir_data '\track_info.mat'],'bath','mask','-append');
    
    % Clean-up
    clear info bath maxb mask ntrplnt lon_mdl lat_mdl x_mdl y_mdl x y bath_gebco;
    
end
%==========================================================================
if(Switch(03))
    
    %Load GEBCO distance-from-shore data
	lon_gebco    = ncread(glbl.file_dshore,'longitude');
	lat_gebco    = ncread(glbl.file_dshore,'latitude');
    dshore_gebco = ncread(glbl.file_dshore,'dist_to_coastline');
  
    % Convert dshore_gebco to kilometers
    dshore_gebco = dshore_gebco./1000;
  
    % Load sampling info
    info = load([glbl.dir_data '\track_info.mat'],'lon','lat');
    
    % lon-lat to x-y
    [YG,XG] = meshgrid(lat_gebco,lon_gebco);
    [x_gebco, y_gebco] = grn2eqa(YG,XG,[0 0],referenceEllipsoid('wgs84'));
    [x,y] = grn2eqa(info.lat,info.lon,[0 0],referenceEllipsoid('wgs84'));
    
    % Interpolate bathymetry to sampling points
    ii = find(YG<=66.15 & YG>=-66.15);
    ntrplnt = scatteredInterpolant(double(x_gebco(ii)),double(y_gebco(ii)),double(dshore_gebco(ii)));
    dshore  = ntrplnt(x,y);
    %dshore = interp2(y_gebco,x_gebco,dshore_gebco,y,x,'linear');

    % Debug plot
    if(do_debug)
        figure('units','normalized','outerposition',[0 0 1 1]);
            xx = info.lon(:); 
            yy = info.lat(:);
            zz = dshore(:);
            scatter(x(:),y(:),2,zz(:),'o','filled'); axis tight; colormap(jet); caxis([0 2000]);
            title('Track Point Distance from Shore');
            cb=colorbar; set(get(cb,'ylabel'),'string','bathymetry (m)');
        print(gcf,'-dpng','Plots/Debug/Debug_Step01-04_track-DistFromShore.png');
        close all; clear xx yy zz cb;
    end
    
    % Save
    save([glbl.dir_data '\track_info.mat'],'dshore','-append');
   
    % Clean-up
    clear info n_track n_point dshore *_gebco x XG y YG ii ntrplnt
    
end
%==========================================================================
if(Switch(04))
    
	% Get basin data
    basin_info  = importdata(glbl.file_basins);
    basin_lat   = basin_info.data(:,1);
    basin_lon   = basin_info.data(:,2);   basin_lon(basin_lon<0) = basin_lon(basin_lon<0)+360;
    basin_index = basin_info.data(:,3);
    clear basin_info basin_file;
  
    % Create interpolant
    ntrplnt = scatteredInterpolant(basin_lon,basin_lat,basin_index,'nearest');

    % Load track info
    load([glbl.dir_data '\track_info.mat'],'lon','lat','mask');
    
    % Interpolate
    basin          = NaN(size(lon));
    basin(mask==1) = ntrplnt(lon(mask==1),lat(mask==1));
    
    % Debug plot
    if(do_debug)
        figure('units','normalized','outerposition',[0 0 1 1]);
        	xx = lon(:);
            yy = lat(:);
            zz = basin(:); 
            scatter(xx(:),yy(:),2,zz(:),'o','filled'); axis tight; colormap(hsv); %caxis([0 2000]);
            title('Track Point Basin Indices');
            colorbar;
        print(gcf,'-dpng','Plots/Debug/Debug_Step01-04_track-BasinIndex.png');
        close all; clear xx yy zz cb;
    end
    
    % Save
    save([glbl.dir_data '\track_info.mat'],'basin','-append');
    
    % Clean-up
    clear lon lat basin basin_* ntrplnt
    
end
%==========================================================================
if(Switch(05))
    
  % Load track info
  info = load([glbl.dir_data '\track_info.mat'],'lon','lat');
  n_track    = size(info.lon,1);   
    
  % Loop through altimeter missions
  index_link = cell(n_track,1);

  % Variable to determine: index_link
  for it=1:n_track
  clc; disp(['On track # ' num2str(it) '...']);
    
      % Model file
      fmdl = [glbl.dir_mdl '\' glbl.str_mdl '\_Processed\Track' sprintf('%0.3d',it) '.nc']; 
        
      % Model points
      xm = ncread(fmdl,'x');    xm(xm>360) = xm(xm>360)-360;
      ym = ncread(fmdl,'y');
      
      % Find track point closest to each model point
      xt = info.lon(it,:);
      yt = info.lat(it,:);
      
      % Loop x if need be
      xm = unwrap(xm.*(pi/180)).*(180/pi);
      xt = unwrap(xt.*(pi/180)).*(180/pi);
      
      % Ensure x's unwrapped around the same boundary
      if(min(xm)>max(xt));      xm = xm-360;    end
      if(max(xm)<min(xt));      xm = xm+360;    end
      
      % Loop though model positions and link to track positions
      np = numel(xm);
      index_link{it} = NaN(np,1);
      for ip=1:np
          
        % Distances
        dd = sqrt( (xm(ip)-xt).^2 + (ym(ip)-yt).^2 );
        
        % Min distance
        index_link{it}(ip) = find(dd==min(dd),1);
        
        % Clean-up
        clear dd;
          
      end
      clear ip;
      
      % Clean-up
      clear xm ym xt yt pm pt np;
        
  end
  clear it;
    
    % Debug plot
    if(do_debug)
      figure;
        x = info.lon + sqrt(-1).*info.lat;     
        y = x;
        for it=1:n_track
          y(it,~ismember(1:size(y,2),index_link{it})) = NaN;
        end
        scatter(real(x(:)),imag(x(:)),2,'.k');
        hold on; scatter(real(y(:)),imag(y(:)),4,'.b'); hold off
      axis tight; box on; title('index link');
      print(gcf,'-dpng','Plots/Debug/Debug_step01-5_index-link.png');
      clear x y it; close all;
    end

      
  % Save (appends)
  save([glbl.dir_data '\track_info.mat'],'index_link','-append');    
      
end
%==========================================================================