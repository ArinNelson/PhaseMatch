% Collect baroclinic tide models and interpolate to track points
% Arin Nelson
% on 04/09/2021
% 
% Last updated by Arin Nelson on 09/12/2021
%==========================================================================
clc; close all;

% Switches
Switch     = zeros(9,1);
Switch(01) = 0;     % Read in Richard Ray's baroclinic tide model and interpolate to sample points
Switch(02) = 0;     % Read in Zaron 2019's baroclinic tide model and interpolate to sample points
Switch(03) = 1;     % Test Zaron 2019's baroclinic tide model by subtracting it from altimetry
    % SO FAR, CAN'T SEEM TO GET ZARON TIDE IN CORRECT PHASE... WHY IS THAT?

% Options
do_debug = 1;
load('global.mat');

% Notes
% Ray's tidal fit only exists for depths >1000m
% Ray's tidal fit is in units of cm? for amplitude, degrees for phase
% Ray's tidal fit is the FULL M2 tide, not the baroclinic tide!  
%   To get baroclinic tidedo along-track filtering. See analyze_gofs_maarten.m
% Zaron's tidal fit is in units of m, and empty values are 0's

%==========================================================================
if(Switch(01))
    
	% Richard Ray's baroclinic tide model
	load(glbl.str_raytide,'lon','lat','aM22','pM22');

	% Quick look
	%figure; scatter(lon(:),lat(:),2,log10(aM22(:)),'o','filled'); colormap(jet); colorbar;
  
    % Generate the interpolant
    RayTide = aM22 .* exp( sqrt(-1) .* (pM22.*(pi/180)) ); % Test: abs(RayTide)==aM22; angle(RayTide)==pM22;
    ntrplnt = scatteredInterpolant(double(lon(:)),double(lat(:)),RayTide(:),'linear');

    % Info
    info = load([glbl.dir_data 'track_info.mat']);
    n_track = size(info.lon,1);
    n_point = size(info.lon,2);
    
    % New fits
    tide_ray_fullm2 = NaN(n_track,n_point);
    
    % Valid points
    ii = find(info.mask == 1 & info.bath > glbl.bath_min);
      
    % Interpolate
    tide_ray_fullm2(ii) = ntrplnt(info.lon(ii),info.lat(ii));
    
    % Debug plot
    if(do_debug)
        figure('units','normalized','outerposition',[0 0 1 1]);
            zz = log10(abs(tide_ray_fullm2));
            scatter(info.lon(:),info.lat(:),2,zz(:),'o','filled'); axis tight; box on;
            colormap(jet); colorbar; caxis([0 2]); title('Ray M2 Amp. ( Log_{10}(cm) ) (Jason-3 Points)');
        print(gcf,'-dpng','Plots\Debug\Debug_Step03-1_RayTideModel.png');
        close all; clear i xx zz;
    end
      
    % Save
    save('Data\track_tidemodels.mat','tide_ray_fullm2');
      
    % Clean-up
    clear info n_track n_point tide_ray_fullm2;
    
end
%==========================================================================
if(Switch(02))

    % Misc info
    str_ed = {'M2','S2','K1','O1'};    
	n_ed   = numel(str_ed);
  
    % Load data, and convert from m to cm
    lon = ncread(glbl.str_ztide,'longitude');	nx = numel(lon);
    lat = ncread(glbl.str_ztide,'latitude');     ny = numel(lat);
    tide_zaron_re = NaN(nx,ny,n_ed);
    tide_zaron_im = NaN(nx,ny,n_ed);
    for i=1:n_ed
        tide_zaron_re(:,:,i) = ncread(glbl.str_ztide,[str_ed{i} 're']).*100; % Units of cm
        tide_zaron_im(:,:,i) = ncread(glbl.str_ztide,[str_ed{i} 'im']).*100; % Units of cm
    end
    clear i;
  
    % Zero's are empty
    %tide_zaron_re(tide_zaron_re==0) = NaN;
    %tide_zaron_im(tide_zaron_im==0) = NaN;

    % Meshgrid for easier interpolation
    [yz,xz] = meshgrid(lat,lon);
    zr = tide_zaron_re;     clear tide_zaron_re;
    zi = tide_zaron_im;     clear tide_zaron_im;

%     % First look
%     figure('units','normalized','outerposition',[0 0 1 1]);
%     for i=1:4
%         subplot(2,2,i);   
%             %imagesc(lon,lat,log10(tide_zaron_re(:,:,i)'.^2 + tide_zaron_im(:,:,i)'.^2)); 
%             surf(xz,yz,log10(zr(:,:,i).^2 + zi(:,:,i).^2),'edgecolor','none'); view(2);
%             colormap(jet); colorbar; caxis([-3 1]); 
%             title([str_ed{i} ' ( Log_{10}(cm) )']); set(gca,'ydir','normal','box','on'); axis tight; ylim([-60 60]);
%     end
%     clear i;
    
    % Sampling info
    info    = load('Data\track_info.mat');
    n_track = size(info.lon,1);
    n_point = size(info.lon,2);
    
    % New fits
    tide_zaron_re = NaN(n_track,n_point,n_ed);
    tide_zaron_im = NaN(n_track,n_point,n_ed);
      
	% Valid points
	ii = find(info.mask == 1);

    % Interpolate for each tidal constituent
	for ie=1:n_ed
          
        % RE  
        tmp = tide_zaron_re(:,:,ie) + sqrt(-1).*tide_zaron_im(:,:,ie);
        tmp(ii) = interp2(yz,xz,zr(:,:,ie),info.lat(ii),info.lon(ii),'linear');
        tide_zaron_re(:,:,ie) = real(tmp);
        tide_zaron_im(:,:,ie) = imag(tmp);

	end
    clear ie;
    
    % Debug plot
    if(do_debug)
      figure('units','normalized','outerposition',[0 0 1 1]);
        for ie=1:4
          subplot(2,2,ie);
            zz = log10( sqrt( tide_zaron_re(:,:,ie).^2 + tide_zaron_im(:,:,ie).^2 ) );
            scatter(info.lon(:),info.lat(:),2,zz(:),'o','filled'); axis tight; box on; set(gca,'color','k');
            colormap(jet); colorbar; caxis([-2 0]); title(['Zaron 2019 ' str_ed{ie} ' Amp. ( Log_{10}(cm) ) (Jason-3)']);
        end
      print(gcf,'-dpng','Plots\Debug\Debug_Step03-1_ZaronTideModel_Jason3.png');
      close all; clear i xx zz;
    end
    
    % Save
    save('Data\track_tidemodels.mat','tide_zaron_*','-append');
    
    % Clean-up
    clear info n_track n_point tide_zaron_*;

end
%==========================================================================
if(Switch(03))
    
    % If not loaded, load zaron tide data
    load('Data\track_tidemodels.mat','tide_zaron_re','tide_zaron_im');
    
    % Sampling info
    info    = load('Data\track_info.mat');
    n_track = size(info.lon,1);
    n_point = size(info.lon,2);
    
    % Altimeter data
    altm  = load('Data\track_altmdata.mat');
    anom  = cell(n_track,1);
    ztide = cell(n_track,1);
    vsub  = NaN(n_track,n_point);
    
    % Compute stationary tide time series, the anomaly, and the 
    for it=1:n_track
    
        % Time data
        nc = size(altm.time{it},2);
        
        % Stationary tide time series
        ztide{it} = zeros(n_point,nc,4);
        for ie=1:2
            
            % Relate tide index to Zaron tide
            iw = glbl.itide_zaron(ie);
            if(iw<3 | abs(info.lat)<35)
            
                % Trig stuff  
                trigarg = glbl.tide_omega(iw).*( (altm.time{it}-glbl.tide_tref) ) - (pi/180)*(glbl.tide_chi(iw)+glbl.tide_nu(iw));
                nt = size(trigarg,2);
                
                % Tide stuff
                zz = tide_zaron_re(it,:,ie) + sqrt(-1).*tide_zaron_im(it,:,ie);
                za = abs(zz(:));
                zp = angle(zz(:));
                
                % Resized
                za = repmat(za(:),[1 nt]);
                zp = repmat(zp(:),[1 nt]);
            
                % The tidal time series
                %ztide{it}(:,:,ie) = (repmat(za(:),[1 nt]).*glbl.tide_f(iw)).*exp( sqrt(-1).*(trigarg-repmat(zp(:),[1 nt])) );
                ztide{it}(:,:,ie) = (za.*glbl.tide_f(iw)) .* cos(trigarg - zp);
            
            end
            
        end
        
        % The anomaly
        tmp = altm.sla{it}; tmp(tmp>1e30) = NaN;
        anom{it} = tmp-nansum(ztide{it},3);
        
        % Variance reduction
        vsub(it,:) = nanvar(altm.sla{it},0,2) - nanvar(anom{it},0,2);
        plot(1:n_point,vsub(it,:)); title([sprintf('%0.3d',it) ': ' num2str(nanmean(vsub(it,:)))]);
        drawnow;
        
    end
    clear it;
        
end
%==========================================================================
