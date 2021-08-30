% Collect altimeter data and interpolate to standard and linearly-spaced tracks
% Arin Nelson
% Rewritten 12/12/2020
% V3 07/25/2021
%==========================================================================
clc; close all;

% Switches
Switch     = zeros(9,1);
Switch(01) = 0;     % Collect raw altimeter data
Switch(02) = 1;     % Interpolate to standard points
Switch(03) = 0;     % Check data

% Info
do_debug   = 1;
load('global.mat');

%==========================================================================
if(Switch(01))
    
	% Sampling info
	info = load([glbl.dir_data 'track_info.mat']);
    n_track = size(info.lon,1);

	% Initialize containers
	lon  = cell(n_track,1);     % Longitudes (^oE)
	lat  = cell(n_track,1);     % Latitudes (^oN)
	time = cell(n_track,1);     % Times (datenum)
	num  = cell(n_track,1);     % Total # of samples at this pt
	sla  = cell(n_track,1);     % sea level anomaly (cm)

	% Load data
	for it=1:n_track
	clc; disp(['Loading altimeter data track ' num2str(it) '...']);
    
        % Construct file name  
        path_data = [glbl.dir_altm '\' glbl.str_altm 'acolin\' glbl.str_altm 'acolin_p' sprintf('%0.4d',it) '.nc'];
  
        % Read in data
        lon{it}  = ncread(path_data,'lon' );  
        lat{it}  = ncread(path_data,'lat' );
        time{it} = ncread(path_data,'time');
        sla{it}  = ncread(path_data,'sla' );
        clear path_data;

        % Compute total number of cycles with valid observations
        num{it} = sum(~isnan(sla{it}),1);
      
        % Set time units to matlab datenum
        time{it} = (time{it}./86400) + glbl.altm_tref;
      
        % Convert sea level anomaly from units of m to cm
        sla{it} = sla{it}.*100;

    end
	clear it;
     
	% Debug check
	if(do_debug==1)
        figure('units','normalized','outerposition',[0.2 0 0.6 1]);
            xx = [];    for it=1:n_track;   xx = [xx(:); nanmean(lon{it},1)'];   end
            yy = [];    for it=1:n_track;   yy = [yy(:); nanmean(lat{it},1)'];   end
            nn = [num{:}];
            zm = [];    for it=1:n_track;   zm = [zm(:); nanmean(sla{it},1)'];   end
            zs = [];    for it=1:n_track;   zs = [zs(:); nanstd(sla{it},0,1)'];   end
            ii = 5:10:numel(xx);
            subplot(3,1,1); scatter(xx(ii),yy(ii),2,zm(ii),       'o','filled'); axis tight; box on; colormap(jet); colorbar; title('sla mean');
            subplot(3,1,2); scatter(xx(ii),yy(ii),2,log10(zs(ii)),'o','filled'); axis tight; box on; colormap(jet); colorbar; title('sla Log_{10}(std)'); 
            subplot(3,1,3); scatter(xx(ii),yy(ii),2,nn(ii),       'o','filled'); axis tight; box on; colormap(jet); colorbar; title('num');
        print(gcf,'-dpng',['Plots/Debug/Debug_Step02-01_' glbl.str_altm '_raw.png']);
        close all; clear xx yy nn zm zs ii
	end
    
    % Save
    save([glbl.dir_data '\raw_data_' glbl.str_altm '.mat'],'lon','lat','time','sla','num');
     
    %end
	%clear ia;
    
end
%==========================================================================
if(Switch(02))
    
%     % Load data
%     raw_data = load([glbl.dir_data '\raw_data_' glbl.str_altm '.mat']);
%     
%     % Sampling info
%     info = load([glbl.dir_data '\track_info.mat']);
%     n_track = size(info.lon,1);
%     n_point = size(info.lon,2);
%       
	% Init new variables
    n_cycle = 136;
    time    = cell( n_track,n_point,n_cycle);
    sla     = cell( n_track,n_point,n_cycle);
    num     = zeros(n_track,n_point);

	% Loop through tracks
	for it=1:n_track
	clc; disp(['On track #' num2str(it) '...']);
    
        % Model points on this track
        xm = info.lon(it,:)';
        ym = info.lat(it,:)';
        mm = info.mask(it,:)';   
        im = find(mm==1);
        
%         % lon-lat to kilometers
%         x_raw = zeros(size(lon_raw));
%         y_raw = zeros(size(lat_raw));
%         for ic=1:n_cycle
%             xx = nanfill(lon_raw(ic,:)');   xx = unwrap(xx.*(pi/180)).*(180/pi);
%             yy = nanfill(lat_raw(ic,:)');
%             if(max(xx)>360);	xx = xx-360;	end
%             if(min(xx)<0  );	xx = xx+360;    end
%             [x_raw(ic,:), y_raw(ic,:)] = grn2eqa(yy, xx,[0 info.lon(it,1)],referenceEllipsoid('wgs84'));
%         end

            
        % Loop through cycles & interpolate
        for ic=1:n_cycle
                
            % Gather data
            xx = raw_data.lon{it}(ic,:)';
            yy = raw_data.lat{it}(ic,:)';
            ss = raw_data.sla{it}(ic,:)';
            tt = raw_data.sla{it}(ic,:)';
            
            % Take only valid samples
            ii = find(~isnan(time_raw(ic,:)') & ~isnan(sla_raw(ic,:)'));
            if(~isempty(ii))

                    % Model points
                    om = unwrap(info.lon(it,im).*(pi/180)).*(180/pi);
                    if(max(om)>360);	om = om-360;    end
                    if(min(om)<0);      om = om+360;    end
                    pm = info.lat(it,im);
                    [xm,ym] = grn2eqa(pm,om,[0 info.lon(it,1)],referenceEllipsoid('wgs84'));
                    
                    % Data
                    ot = unwrap(lon_raw(ic,:).*(pi/180)).*(180/pi);  
                    if(max(ot)>360);	ot = ot-360;    end
                    if(min(ot)<0);      ot = ot+360;    end
                    pt = lat_raw(ic,:);
                    [xt,yt] = grn2eqa(pt,ot,[0 info.lon(it,1)],referenceEllipsoid('wgs84'));
                    tt = time_raw(ic,:)';
                    ss = sla_raw( ic,:)';
                    
                    % Fill in missing positions
                    xt = nanfill( xt(:) );    % linearly-interpolate missing lon's
                    yt = nanfill( yt(:) );    % linearly-interpolate missing lat's
                    
                    % Set NaNs to realmax
                    maxt = max(tt(ii));    mint = min(tt(ii));  tt(isnan(tt)) = realmax('single');
                    maxs = max(ss(ii));    mins = min(ss(ii));  ss(isnan(ss)) = realmax('single');
                    
                    % Some lon/lat at edges may still be missing
                    ii = find(~isnan(xt) & ~isnan(yt));
                    xt = xt(ii);
                    yt = yt(ii);
                    tt = tt(ii);
                    ss = ss(ii);
                    
                    % Interpolants
                    ntrplnt_t = scatteredInterpolant(xt(:),yt(:),tt(:),'linear');
                    ntrplnt_s = scatteredInterpolant(xt(:),yt(:),ss(:),'linear');
                    
                    % Interpolate
                    time{it}(im,ic) = ntrplnt_t(xm,ym);
                    sla{ it}(im,ic) = ntrplnt_s(xm,ym);
                    
                    % Remove erroneous values
                    time{it}(time{it}(:,ic)>maxt,ic) = NaN;
                    time{it}(time{it}(:,ic)<mint,ic) = NaN;
                    sla{it}( sla{it}( :,ic)>maxs,ic) = NaN;
                    sla{it}( sla{it}( :,ic)<mins,ic) = NaN;
                
                    % Lookit?
                    if(do_debug==1 & ic==n_cycle)
                        ss(ss>maxs) = NaN;
                        tt(tt>maxt) = NaN;
                        subplot(1,2,1); hold on
                            scatter(xt(:),yt(:),8,tt(:),'o','filled');
                            scatter(xm(:),ym(:),8,time{it}(im,ic),'s','filled');
                        hold off; axis tight; box on; colorbar;
                        subplot(1,2,2); hold on
                            scatter(xt(:),yt(:),8,ss(:),'o','filled');
                            scatter(xm(:),ym(:),8,sla{it}(im,ic),'s','filled');
                        hold off; axis tight; box on;   colorbar;
                        pause(1e-9);
                    end
                    
                    % Clean-up
                    clear ntrplnt_t ntrplnt_s;
                    
                end
                clear ii;
                
            end
            clear ic;
            
        % Num
        num(it,:) = sum(~isnan(sla{it}),2);
            
	end
	clear it;

    % Debug plot
    if(do_debug)
      figure('units','normalized','outerposition',[0 0 1 1]);
        xx = info.lon';     xx = xx(:);
        yy = info.lat';     yy = yy(:);
        ii = 1:numel(xx);
        subplot(2,3,1); % mask
          zz = info.mask';      
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('mask');
        subplot(2,3,2); % bath
          zz = info.bath';      
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('bath');
        subplot(2,3,3); % dshore
          zz = info.dshore';	
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('dshore');
        subplot(2,3,4); % num
          zz = num';  zz = zz(:);    
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('num');
        subplot(2,3,5); % sla mean
          zz = []; for it=1:n_track; zz = [zz(:)', nanmean(sla{it},2)']; end
          scatter(xx(ii),yy(ii),2,zz(ii),'o','filled'); axis tight; colormap(jet); colorbar; title('sla mean');
        subplot(2,3,6); % sla std
          zz = []; for it=1:n_track; zz = [zz(:); nanstd(sla{it},0,2)]; end
          scatter(xx(ii),yy(ii),2,log10(zz(ii)),'o','filled'); axis tight; colormap(jet); colorbar; title('sla Log_{10}(std)');
      print(gcf,'-dpng',['Plots/Debug/Debug_Step02-02_' str_mthd{id} '-' str_altm{ia} '.png']);
      clear xx yy zz ii; close all;
    %end
    end

    % Save
    save([dir_data '\track_altmdata.mat'],'time','sla','num');
    clear time sla num info n_track n_point
    
end
%==========================================================================
if(Switch(03))
    
    % Statistics
    info    = load([dir_data '\track_info.mat']);
    n_track = size(info.lon,1);
    n_point = size(info.lon,2);
    std     = zeros(size(info.lat));
    for i=1:n_track
        std(i,:) = nanstd(sla{i},0,2);
    end
    xx = info.lon(:,:,1);
    scatter(xx(:),info.lat(:),4,std(:),'o','filled'); colormap(jet); colorbar;
    
end