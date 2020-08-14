function scenario = drt_env(scenario)
    %drt_env: code to download relevant environmental forcings and set
    %these time series datasets needed to run erosion/accretion model
    %
    %Required Inputs: 'scenario' structure variable with the following
    %infromation,
    %       scenario.location.lat [value from -90 to 90]
    %       scenario.location.lon [value from -180 to 180]
    %       scenario.timing.start_date [model start date in datenum format]
    %       scenario.timing.end_date [model end date in datenum format]
    %       scenario.timing.dt [model time step in hrs]
    %       scenario.type [string variable indicating this is a 'hindcast' or forecast']
    %
    %Outputs:
    %       scenario.env.waves
    %       scenario.env.winds
    %       scenario.env.tides
   
    if strcmp(scenario.type, 'hindcast') == 1  %workflow for hindcast model setup

            %throw an error if the end date selected is past the period wheere WIS waves are present
            if scenario.timing.end_date > datenum(2017,12, 31)
                f = msgbox('Model simulation dates must be between 1980 and 2017. Select an end date thht falls in this range', 'Error','error');
                error('');
            end
            %develop scenario timing
            scenario.timing.times = scenario.timing.start_date:[scenario.timing.dt/24]:scenario.timing.end_date;         
            
            %acquire wave and wind data from the Wave Information Studies database
            scenario.wis = wis_determine_node(scenario); %find closest wis node to the relevant site    
            [scenario.env.waves, scenario.env.winds] = wis_download(scenario); %download data

            %acquire tide data for NOAA Tides and Currents
            scenario.noaa = noaa_determine_node(scenario);
                             scenario.env.tides = noaa_download_tides(scenario);

                             
            try %download the verified tide data, if avaialble
                 scenario.env.tides = noaa_download_tides(scenario);
            catch err %add in an error check an error can occur if no measured data exists during the time period wanted
                scenario.env.tides = noaa_download_tide_prediction(scenario);           
            end                
            
      elseif strcmp(scenario.type, 'forecast') == 1  %workflow for forecast model setup

            %first pull the wave data from the global WaveWatchIII model
            [scenario.env.waves, scenario.timing.times] = ww3_forecast_download(scenario);  

            %create timing information based on output from WaveWatchIII
            scenario.timing.start_date = min(scenario.timing.times);
            scenario.timing.end_date = max(scenario.timing.times);
            
            %acquire wind data from the Global Forecast System
            scenario.env.winds = gfs_forecast_download(scenario);
           
            %download tidal data from NOAA Tides and Currents
            scenario.noaa = noaa_determine_node(scenario);
            scenario.env.tides = noaa_download_tide_prediction(scenario);           

            %add in storm surge
            scenario.wis = wis_determine_node(scenario); %will not use wis data, but this easily tells us which coast the site is on
            zone = scenario.wis.closest_zone;
            
            try %the OpenDap ESTOFs site seems to have dependability issues, so a try catch is necessary here
                scenario.env.tides = download_ESTOFS(scenario, zone);
            catch err
                scenario.env.tides.surge = choosedialogSurge;
                scenario.env.tides.wl = scenario.env.tides.wl + scenario.env.tides.surge;                
            end   
            
    end
        
end

function wis = wis_determine_node(scenario)
    %wis_determine_node: finds the closest Wave Information Studies node to
    %the provided field site latitude and longitude
    %Required Inputs: 'scenario' structure variable with the following
    %infromation,
    %       scenario.location.lat [value from -90 to 90]
    %       scenario.location.lon [value from -180 to 180]   

    %load metadata file that provides all useable WIS nodes for this analysis
    wis.table = readtable('drt_env_station_list.xlsx', 'Sheet', 'WIS');

    %find closest node to given lat/lon
    distances = distAway(scenario.location.lat,scenario.location.lon,wis.table.Lat,wis.table.Lon);
    [minval, imin] = nanmin(distances);
    
    if minval > 5 %dont consider a node more than 5 degrees away
        f = msgbox('No environmental node close to the selected site', 'Error','error');
        error('');
    end
                
    %store relevant information
    wis.closest_node = wis.table.Station(imin);
    wis.closest_lat = wis.table.Lat(imin);
    wis.closest_lon = wis.table.Lon(imin);
    wis.closest_depth = wis.table.Depth_m_(imin);
    wis.closest_zone = wis.table.Region(imin);
    end

function distances = distAway(ptY, ptX, listY, listX)
    %ptX, ptY should be a single value
    distances = sqrt((listX(:)-ptX).^2 + (listY(:)-ptY).^2);
end
    
function [waves, winds] = wis_download(scenario)
    %wis_download: downloads WIS ONLNS file and loads data
    %Required Inputs: 'scenario' structure variable with the following
    %infromation,
    %       scenario.wis.closest_node 
    %       scenario.wis.closest_zone  


    %relevant information from scenario file
    station_num = scenario.wis.closest_node;   
    
    %check what part of the country the zone is in
    if strcmp(scenario.wis.closest_zone{1}, 'Atlantic') == 1
        zone = 'atl'; %note that there are different conventions depending on whether pulling onlns or thredds
        zone2 = 'Atlantic';
    elseif strcmp(scenario.wis.closest_zone{1}, 'Pacific') == 1
        zone = 'pac';
        zone2 = 'Pacific';
    elseif strcmp(scenario.wis.closest_zone{1}, 'GulfOfMexico') == 1
        zone = 'gom';
        zone2 = 'GulfOfMexico';
    end
    
    try %download wave data
        url = ['https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/', zone2, '/ST', num2str(station_num),'/ST', num2str(station_num),'.ncml#noprefetch'];

        %download and convert time
        time=ncread(url,'time'); % downloading time from server
        pause(1); %seems like need to give some time to thredds
        %tunit= ncreadatt(url,'time','units'); % reading attributes of variable time
        mtime=time/(3600.0*24)+datenum(1970,1,1);

        % finding index that corresponds to dates of interest
        itime=find(mtime >= scenario.timing.start_date & mtime<= scenario.timing.end_date); % indicies in netCDF record of data of interest

        % pulling data from server with itime index 
        time = mtime(itime); % record of time in matlab datetime 
        waveHs = ncread(url,'waveHs',min(itime),length(itime));
        %pause(1); %seems like need to give some time to thredds
        waveTp = ncread(url,'waveTp',min(itime),length(itime));
        %pause(1); %seems like need to give some time to thredds
        waveD = ncread(url,'waveMeanDirection',min(itime),length(itime));
        %pause(1); %seems like need to give some time to thredds
        windSpeed = ncread(url,'windSpeed',min(itime),length(itime));
        %pause(1); %seems like need to give some time to thredds
        windD = ncread(url,'windDirection',min(itime),length(itime));
        %pause(1); %seems like need to give some time to thredds 
                                    
    catch err %alternatively download text files instead using thredds server if down
        
        if strcmp(zone2, 'Pacific') == 1
            %last date ==
            if scenario.timing.end_date > datenum(2011, 12, 31)
                f = msgbox('WIS Thredds Download Failed: Pick a Date Before 2011 For the Pacific Basin to proceed with altnerative ONLNS File', 'Error','error');
                error('');                
            end
        elseif strcmp(zone2, 'Atlantic') == 1   
             if scenario.timing.end_date > datenum(2014, 12, 31)
                f = msgbox('WIS Thredds Download Failed: Pick a Date Before 2014 For the Atlantic Basin to proceed with altnerative ONLNS File', 'Error','error');
                error('');           
             end
        elseif strcmp(zone2, 'GulfOfMexico') == 1   
             if scenario.timing.end_date > datenum(2014, 12, 31)
                f = msgbox('WIS Thredds Download Failed: Pick a Date Before 2014 For the Atlantic Basin to proceed with altnerative ONLNS File', 'Error','error');
                error('');           
             end            
        end
        
        try %check to see fi file is already downloaded, speeds up process substantially
            file = dir(['ST', num2str(station_num), '*.onlns']);
            if numel(file) == 1
                data = dlmread(file.name);
            else
                data = dlmread(file(1).name);
            end            
        catch err %if file is not avaialble then download
            websave('wis.temp.onlns.zip', ['http://wis.usace.army.mil/data/', zone,'/onlns/raw/ST', num2str(station_num),'_ONLNS.zip']);
            unzip('wis.temp.onlns.zip');
            file = dir(['ST', num2str(station_num), '*.onlns']);
            if numel(file) == 1
                data = dlmread(file.name);
            else
                data = dlmread(file(1).name);
            end
        end
        
        
        %pull out relevant variables from the WIS ONLNS files
        time =datenum(num2str(data(:,1)), 'yyyymmddhhMMSS');
        [~, iunique] = unique(time); %get rid of repeat time issues here
        time = double(time(iunique));
        windSpeed = double(data(iunique,5));
        windD = double(data(iunique,6));
        waveHs = double(data(iunique,10));
        waveTp = double(data(iunique,12));
        waveD = double(data(iunique,16));
     end        
                
                
    %interpolate waves onto desired time interval
    waves.Hs_deepwater = interp1(time, waveHs, scenario.timing.times);
    waves.Tp = interp1(time, waveTp, scenario.timing.times);
    waves.D_deepwater = interp1(time, waveD, scenario.timing.times);
    waves.depth = scenario.wis.closest_depth;
    winds.windSpeed = interp1(time, windSpeed, scenario.timing.times);    
    winds.windDirection = interp1(time, windD, scenario.timing.times);   
    
    %shoal wave heights to 25 m for Stockdon and fix any other output issues
    waves = transform_waves(waves, scenario);
           
end

function noaa = noaa_determine_node(scenario)
    %noaa_determine_node: finds the closest NOAA tide gauge to
    %the provided field site latitude and longitude
    %Required Inputs: 'scenario' structure variable with the following
    %infromation,
    %       scenario.location.lat [value from -90 to 90]
    %       scenario.location.lon [value from -180 to 180]   


    %load metdata
    noaa.table = readtable('drt_env_station_list.xlsx', 'Sheet', 'NOAA_Tides');


    %find closest node to given lat/lon
    distances = distAway(scenario.location.lat,scenario.location.lon,noaa.table.Lat,noaa.table.Lon);
    [minval, imin] = nanmin(distances);
    
    if minval > 5 %dont consider a node more than 5 degrees away
        f = msgbox('No environmental node close to the selected site', 'Error','error');
        error('');
    end
      
    %store relevant information
    noaa.closest_node = noaa.table.Station(imin);
    noaa.closest_lat = noaa.table.Lat(imin);
    noaa.closest_lon = noaa.table.Lon(imin);
end

function tides = noaa_download_tides(scenario)
    %noaa_download_tides: function to download NOAA 
    %Required Inputs: 'scenario' structure variable with the following
    %infromation,
    %       scenario.timing.start_date
    %       scenario.timing.end_date     
    %       scenario.noaa.closest_node
    
        
    %initialize variables
    [startYear, ~,~, ~,~,~] = datevec(scenario.timing.start_date);
    [endYear, ~,~, ~,~,~] = datevec(scenario.timing.end_date);
    datum = 'NAVD'; %stations in the list should all be NAVD compatible
    gauge = num2str(scenario.noaa.closest_node);

    %loop through each year of vertified tide data for noaa erddap server
    wl = [];
    time = [];
    for yr = startYear:endYear
        
        %erdap server link
        website = ['https://opendap.co-ops.nos.noaa.gov/erddap/tabledap/IOOS_Hourly_Height_Verified_Water_Level.mat?STATION_ID%2CDATUM%2CBEGIN_DATE%2CEND_DATE%2Ctime%2CWL_VALUE&STATION_ID=%22', gauge,'%22&DATUM%3E=%22', datum,'%22&BEGIN_DATE%3E=%22', num2str(yr),'0101%2000%3A00%22&END_DATE%3E=%22', num2str(yr),'1231%2023%3A59%22'];
        
        %download data
        out_name = 'tides.mat'; %need to store temporary data to current folder
        options = weboptions;
        options.Timeout = 120;
        websave(out_name, website, options);
        
        %load and then clean up variables
        load(out_name);
        delete(out_name);

        %store data
        wltemp = IOOS_Hourly_Height_Verified_Wat.WL_VALUE;
        timetemp = double(IOOS_Hourly_Height_Verified_Wat.time/86400 + datenum(1970,1,1));
        wl = [wl; wltemp];
        time = [time; timetemp];

    end

%     %loop through each year of predicted tide data to fill in any data gaps in ther verified data
    yrs=startYear:endYear; % define all the year withing the start and end year limit
    dataout=[];    
    for i=1:length(yrs)

         st_date = datestr([yrs(i),1,1,0,0,0],'yyyymmdd');
         end_date   = datestr([yrs(i),12,31,23,59,59],'yyyymmdd');

         url=['http://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=',st_date,'&end_date=',end_date,'&datum=NAVD&station=',gauge,'&time_zone=GMT&units=metric&interval=h&format=CSV'];

         datatemp=webread(url);    

         dataout=[dataout;datatemp];
    end   
    datacell = table2cell(dataout); 
    for i=1:length(datacell)
        time_pred(i)=datenum(datacell{i,1});
        wl_pred(i)=datacell{i,2};
    end  

%  %loop through each year of vertified tide data for noaa erddap server
%     wl_pred = [];
%     time_pred = [];
%     for yr = startYear:endYear
%         
%         %erdap server link
%         website = ['https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsWLTP60.mat?STATION_ID%2CDATUM%2CBEGIN_DATE%2CEND_DATE%2Ctime%2CpredictedWL&STATION_ID=%22', gauge,'%22&DATUM%3E=%22MLLW%22&BEGIN_DATE%3E=%22', num2str(yr),'0101%2000%3A00%22&END_DATE%3E=%22', num2str(yr),'1231%2023%3A59%22'];
%         
%         %download data
%         out_name = 'tides.mat'; %need to store temporary data to current folder
%         options = weboptions;
%         options.Timeout = 120;
%         websave(out_name, website, options);
%         
%         %load and then clean up variables
%         load(out_name);
%         delete(out_name);
% 
%         %store data
%         wltemp = double(nosCoopsWLTP60.predictedWL);
%         timetemp = double(nosCoopsWLTP60.time/86400 + datenum(1970,1,1));
%         wl_pred = [wl_pred; wltemp];
%         time_pred = [time_pred; timetemp];
% 
%     end
    
    %interpolate data to output
    tides.wl = interp1gap(time, double(wl), scenario.timing.times, 6/24);  
    tides.pwl = interp1(time_pred, double(wl_pred), scenario.timing.times);  
    ibad = find(isnan(tides.wl) == 1);
    tides.wl(ibad) = tides.pwl(ibad);   
    
end


function [waves2, times] = ww3_forecast_download(scenario)
    %noaa_download_tides: function to download NOAA 
    %Required Inputs: 'scenario' structure variable with the following
    %infromation,
    %       scenario.timing.start_date
    %       scenario.timing.end_date     
    %       scenario.location.lat [value from -90 to 90]
    %       scenario.location.lon [value from -180 to 180]       

    %generate timing in format needed for erddap
    start_time = datestr(floor(now), 'yyyy-mm-ddT00:00:00Z');
    end_time = datestr(floor(now)+4, 'yyyy-mm-ddT00:00:00Z');

    %first just try the coordinates given
    waves = download_ww3(start_time, end_time, scenario.location.lat, scenario.location.lon);
    
    %next try the location of the closest wis node
    if numel(find(isnan(waves.Hs_deepwater) == 1)) > 0
        wis = wis_determine_node(scenario);
        lat= wis.closest_lat;
        lon = wis.closest_lon;
        waves = download_ww3(start_time, end_time, lat,lon);
    end
   
    %now try different locations if not
    %go west if in the west coast, east on the east coast, and south in the
    %gulf coast
    if numel(find(isnan(waves.Hs_deepwater) == 1)) > 0
        if strcmp(wis.closest_zone{1}, 'Atlantic') == 1
            lon = lon+1;
        elseif strcmp(wis.closest_zone{1}, 'Pacific') == 1
            lon = lon-1;
        elseif strcmp(wis.closest_zone{1}, 'GulfOfMexico') == 1
            lat = lat -1;
        end    
        waves = download_ww3(start_time, end_time, lat,lon);
    end    

    %try one last time
    if numel(find(isnan(waves.Hs_deepwater) == 1)) > 0
        if strcmp(wis.closest_zone{1}, 'Atlantic') == 1
            lon = lon+1;
        elseif strcmp(wis.closest_zone{1}, 'Pacific') == 1
            lon = lon-1;
        elseif strcmp(wis.closest_zone{1}, 'GulfOfMexico') == 1
            lat = lat -1;
        end    
        waves = download_ww3(start_time, end_time, lat,lon);
    end    
   
    %download local water depths at wave note since this is not part of the WW3 outout
    url = ['https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus.mat?z[(', num2str(waves.latitude),'):1:(', num2str(waves.latitude),')][(', num2str(waves.longitude-360),'):1:(', num2str(waves.longitude-360),')]'];
    websave('srtm.mat',url);
    
    %load water depth data
    load('srtm.mat');
    waves.depth = srtm15plus.z;
    delete('srtm.mat');
    
    %transform waves to shallow water
    waves = transform_waves(waves, scenario);
    
    %generate new timings
    scenario.timing.times = round(waves.time):[scenario.timing.dt/24]:floor(max(waves.time));
   
    %interpolate
    waves2.Hs_25m = interp1(waves.time, waves.Hs_25m, scenario.timing.times);
    waves2.Hs_deepwater = interp1(waves.time, waves.Hs_deepwater, scenario.timing.times);
    waves2.L_25m = interp1(waves.time, waves.L_25m, scenario.timing.times);
    waves2.D_25m = interp1(waves.time, waves.D_25m, scenario.timing.times);
    waves2.D_deepwater = interp1(waves.time, waves.D_deepwater, scenario.timing.times);
    waves2.Tp = interp1(waves.time, waves.Tp, scenario.timing.times);
    waves2.times = scenario.timing.times;  
    times = scenario.timing.times;
end

function waves = transform_waves(waves, scenario)

    %fix wave problems
        localD = abs(wrapTo180(waves.D_deepwater - scenario.grids.morphometrics.azimuth));
        
        %also set wave heights as zero if wave direction is headed offshore
        ifind = find(abs(localD)>= 90);     
        Hstemp = waves.Hs_deepwater;
        Hstemp(ifind) = 0;
        localD(ifind) = 0;
    
        %the refraction gets confused for really oblique waves, so set an
        %upper limit
        maxD = 60;
        ifind = find(abs(localD)>= maxD);     
        Hstemp = waves.Hs_deepwater;
        localD(ifind) = maxD;    
       
%    if nanmean(localD) < 45
        [waves.Hs_25m,waves.L_25m,waves.D_25m] = shoal_waves(double(Hstemp),abs(double(waves.depth)),double(localD),double(waves.Tp),25);                        
        waves.Hs_25m = real(waves.Hs_25m);
%     else %there may be something wrong with the azimuth so just use the waves as is
%             waves.Hs_25m = waves.Hs_deepwater;
%             waves.L_25m = 9.81*waves.Tp(:).*waves.Tp(:)./(2*pi);
%             waves.D_25m = zeros(size(waves.Hs_25m));
%     end
    %end
    
    if isfield(waves, 'Hs_25m') ~= 1
            waves.Hs_25m = waves.Hs_deepwater;
            waves.L_25m = 9.81*waves.Tp(:).*waves.Tp(:)./(2*pi);
            waves.D_25m = zeros(size(waves.Hs_25m));
    end
       
%     %also set wave heights as zero if wave direction is headed offshore
%     if nanmean(localD) < 45
%          waves.Hs_25m(abs(localD)>=65) = 0;
%     end
    
    %there are also other hiccups in the shoal code that leads to large
    %wave heights. a limiter is added to prevent this
    ibad = find(waves.Hs_25m > 20);
    waves.Hs_25m(ibad) = waves.Hs_deepwater(ibad);
    
    %lastly if there are nans, just use offshore condition
    ibad = find(isnan(waves.Hs_25m) ==  1);
    waves.Hs_25m(ibad) = waves.Hs_deepwater(ibad);  

end


function waves = download_ww3(start_time, end_time, lat, lon)

    %download waves from WaveWatchIII
    url = ['https://coastwatch.pfeg.noaa.gov/erddap/griddap/NWW3_Global_Best.mat?Tdir[(', start_time,'):1:(', end_time,')][(0.0):1:(0.0)][(', num2str(lat),'):1:(', num2str(lat),')][(', num2str(wrapTo360(lon)),'):1:(', num2str(wrapTo360(lon)),')],Tper[(', start_time,'):1:(', end_time,')][(0.0):1:(0.0)][(', num2str(lat),'):1:(', num2str(lat),')][(', num2str(wrapTo360(lon)),'):1:(', num2str(wrapTo360(lon)),')],Thgt[(', start_time,'):1:(', end_time,')][(0.0):1:(0.0)][(', num2str(lat),'):1:(', num2str(lat),')][(', num2str(wrapTo360(lon)),'):1:(', num2str(wrapTo360(lon)),')]'];
    websave('ww3.mat',url);
    
    %load data
    load('ww3.mat');    
    waves.time = double(ww3_global.time)/86400+datenum(1970,1,1);
    waves.Hs_deepwater = double(ww3_global.Thgt);
    waves.D_deepwater = double(ww3_global.Tdir);
    waves.Tp = double(ww3_global.Tper);
    waves.latitude = double(ww3_global.latitude);
    waves.longitude = double(ww3_global.longitude);          
end

function winds = gfs_forecast_download(scenario)
    start_time = datestr(floor(now), 'yyyy-mm-ddT00:00:00Z');
    end_time = datestr(floor(now)+4, 'yyyy-mm-ddT00:00:00Z');
    
    url = ['https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best.mat?ugrd10m[(', start_time,'):1:(', end_time,')][(', num2str(scenario.location.lat),'):1:(', num2str(scenario.location.lat),')][(', num2str(wrapTo360(scenario.location.lon)),'):1:(', num2str(wrapTo360(scenario.location.lon)),')],vgrd10m[(', start_time,'):1:(', end_time,')][(', num2str(scenario.location.lat),'):1:(', num2str(scenario.location.lat),')][(', num2str(wrapTo360(scenario.location.lon)),'):1:(', num2str(wrapTo360(scenario.location.lon)),')]'];
    websave('gfs.mat',url);
    load('gfs.mat');    
    
    time = ncep_global.time/(86400)+datenum(1970,1,1);
    u = ncep_global.ugrd10m;
    v = ncep_global.vgrd10m;
    [windSpeed,windD] = uv_to_wswd(u, v);
   
    %interpolate waves onto desired time interval
    winds.windSpeed = interp1(time, windSpeed, scenario.timing.times);    
    winds.windDirection = interp1(time, windD, scenario.timing.times);  
end


 function [ws,wd] = uv_to_wswd(u, v)
    % function that takes vectors (u,v) and converts to  
    % wind speed and wind direction (deg from North) vectors
    % This version copes with NaNs in vectors: 23 April 2007
    % function [ws,wd] = uv_to_wswd2(u, v)

    ws = NaN*ones(size(u));
    wd = NaN*ones(size(u));
    e = find(~isnan(u) & ~isnan(v));
    ws(e) = sqrt(u(e).*u(e) + v(e).*v(e));
    wd(e) = 270 - (180/pi)*atan2(v(e),u(e));
    for i = 1:length(wd)
        if (wd(i) > 360) 
            wd(i) = wd(i)-360; 
        end
    end

 end
  
function tides = noaa_download_tide_prediction(scenario)

    %initialize variables
    [startYear, ~,~, ~,~,~] = datevec(min(scenario.timing.times));
    [endYear, ~,~, ~,~,~] = datevec(max(scenario.timing.times));
    datum = 'NAVD'; %stations in the list should all be NAVD compatible
    gauge = num2str(scenario.noaa.closest_node);

    yrs=startYear:endYear; % define all the year withing the start and end year limit
    dataout=[];    
    for i=1:length(yrs)

         st_date = datestr([yrs(i),1,1,0,0,0],'yyyymmdd');
         end_date   = datestr([yrs(i),12,31,23,59,59],'yyyymmdd');

         url=['http://tidesandcurrents.noaa.gov/api/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date=',st_date,'&end_date=',end_date,'&datum=NAVD&station=',gauge,'&time_zone=GMT&units=metric&interval=h&format=CSV'];
         options = weboptions;
        options.Timeout = 120;
         datatemp=webread(url, options);    

         dataout=[dataout;datatemp];
    end 
    
    datacell = table2cell(dataout); % need to change the table to cell as unfortunately tabe datetime is not a string to apply datevec
    for i=1:length(datacell)
        time_pred(i)=datenum(datacell{i,1});
        wl_pred(i)=datacell{i,2};
    end  


    tides.wl = interp1(time_pred, double(wl_pred), scenario.timing.times);   
    
end

function [H1,L1,alpha1] = shoal_waves(H0,h0,alpha0,T,h1)
    %shoal.m shoals surface gravity waves to a different water depth
    %H0 [meters] given wave height
    %h0 [meters] given water depth
    %alpha0 [deg] given angle rel. to shore normal
    %T [sec] wave period
    %h1 [meters]  depth where characteristics are sought
    %
    [alpha1] = snells(alpha0,h0,T,h1);
    [k,n0,c0] = dispersion (2*pi./T,h0);
    [k,n1,c1] = dispersion (2*pi./T,h1);
    L1 = 2*pi./k;
    H1 = H0.*sqrt((c0.*n0)./(c1.*n1)).*sqrt(cos(alpha0*pi/180)./cos(alpha1*pi/180));

end

function [alpha] = snells(alpha0,h0,T,d)
    %snells.m usings Snells law to determine wave direction changes
    %Calculate the changing angle due to refraction for the 
    % case of straight parallel contours for monchromatic wave
    % alpha0[deg]
    % h0[m]
    % T[s]
    % d [m] matrix or array of depths from which to find angles
    % alpha [deg] is matrix or array of angles with size(alpha)=size(d) 
    % function [alpha] = snells(alpha0,h0,T,d)
    omega = 2*pi./T;
    [k0] = dispersion (omega,h0);c0 = omega./k0;
    [k] = dispersion (omega,d);c = omega./k;
    alpha = real(180/pi*asin(sin(alpha0*pi/180).*c./c0));

end
  
function [k,n,c] = dispersion (w,h,showflag)
    % dispersion solves the linear dispersion relation
    % given frequency, w = 2 pi/Period [1/sec], 
    % and water depth, h [meters] the dispersion routine returns the 
    % wavenumber, k = 2 pi/Wavelength [1/m],
    % n = ratio of group speed to phase speed, and 
    % c = phase speed [m/s]
    % [k n c] = dispersion (w,h)
    % function [k,n,c] = dispersion (w,h,showflag)

    if nargin<3;showflag = 1;end
    if find(h<=0)>0&showflag;
      %  disp('WARNING: water depth should be > 0');
      %  disp('A wavenumber of 10^10 will be returned for all depths <= 0  ');
    end
    ind = find(h<=0);h(ind)= 1;
    g = 9.81;
    k=w./sqrt(g*h);diff=max(max(w.^2-g*k.*tanh(k.*h)));
    while abs(diff) > 1*10^-8
      knew = k -(w.^2-g*k.*tanh(k.*h))./(-g*tanh(k.*h)-g*k.*h.*sech(k.*h).^2);
      k = knew;
      diff=max(max(w.^2-g*k.*tanh(k.*h)));
    end
    c=w./k;
    k(ind) = 10^10;
    n = .5*(1+(2*k.*h)./sinh(2*k.*h));
    n(ind) = 1;
    c(ind) = 0;
end


function tides = download_ESTOFS(scenario, zone)
    %Function to download ESTOFS surge data for input to the model

    input_lat = scenario.location.lat;
    input_lon = scenario.location.lon;
    
    %Time keeping
    time = now;
    timevec = datevec(time);
    if timevec(4)>=21
        datestring = datestr(ceil(now), 'yyyymmdd');
    else
        datestring = datestr(floor(now), 'yyyymmdd');
    end

    %Determine which geographic zone are in
    if strcmp(zone, 'Pacific') == 1
        estofs_zone = 'pac';
    elseif strcmp(zone, 'Atlantic') == 1
        estofs_zone = 'atl';        
    elseif strcmp(zone, 'GulfOfMexico') == 1
        estofs_zone = 'atl';        
    else
        error('No forecast data in this zone');
    end
        
    %Use only in forecast mode
    info = ncinfo(['https://nomads.ncep.noaa.gov:9090/dods/estofs_', estofs_zone,'/',datestring,'/estofs_', estofs_zone,'_conus_00z']);

    %find the closest node
    lon = double(info.Variables(1,3).Attributes(1,6).Value:info.Variables(1,3).Attributes(1,8).Value:info.Variables(1,3).Attributes(1,7).Value);
    lat = double(info.Variables(1,2).Attributes(1,6).Value:info.Variables(1,2).Attributes(1,8).Value:info.Variables(1,2).Attributes(1,7).Value);
    time = double(temptime1:1/24:temptime2);
    [LON, LAT] = meshgrid(lon, lat);
    [minval, ilatuse] = min(abs(lat - input_lat));
    [minval, ilonuse] = min(abs(lon - input_lon));
    
    %set up timing information
    temptime1 = info.Variables(1,1).Attributes(1,8).Value;
    temptime1 = datenum(temptime1(4:end), 'ddmmmyyyy')+str2num(temptime1(1:2))/24;
    temptime2 = info.Variables(1,1).Attributes(1,9).Value;
    temptime2 = datenum(temptime2(4:end), 'ddmmmyyyy')+str2num(temptime2(1:2))/24;
    dt = info.Variables(1,1).Attributes(1,10).Value;
    
    if dt> 0.04 & dt<0.045
        dt = 1/24;
    end
    
    %find longitude location closest to shore (but must be at or west of defined point) with surge values
    surge2 = NaN;
    count = 0;
    while sum(isnan(surge2))>0
        surge = ncread(['https://nomads.ncep.noaa.gov:9090/dods/estofs_', estofs_zone, '/',datestring,'/estofs_', estofs_zone,'_conus_00z'], 'etsrgsfc', [ilonuse+count ilatuse 1], [1 1 Inf]);
        clear surge2
        surge2(:,1) = surge(1,1,:);
        count = count-1;
        if count<-10
            surge2 = zeros(size(surge2));
        end
    end

    tides.surge = interp1(time, surge2, scenario.timing.times);  
    tides.wl = scenario.tides.wl+tides.surge;
    
end


function surge = choosedialogSurge

    d = dialog('Position',[300 300 250 150],'Name','Select One');
    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[20 80 210 40],...
           'String','Surge Forecast Download Failed. Manually add a sea level anomaly?');
       
    popup = uicontrol('Parent',d,...
           'Style','popup',...
           'Position',[75 70 100 25],...
           'String',{'No SLA    ';'0.1 m SLA';'0.25 m SLA'; '0.50 m SLA'; '1.0 m SLA '},...
           'Callback',@popup_callback);
       
    btn = uicontrol('Parent',d,...
           'Position',[89 20 70 25],...
           'String','Close',...
           'Callback','delete(gcf)');
       
    choice = 'No SLA';
       
    % Wait for d to close before running to completion
    uiwait(d);
    
       
   if strcmp(choice, 'No SLA')    
       surge = 0;
   elseif strcmp(choice, '0.10 m SLA')
       surge = 0.1;
   elseif strcmp(choice, '0.25 m SLA')
       surge = 0.25;       
   elseif strcmp(choice, '0.50 m SLA')
       surge = 0.5;   
   elseif strcmp(choice, '1.0 m SLA ')
       surge = 1;   
   end       
   
       function popup_callback(popup,event)
          idx = popup.Value;
          popup_items = popup.String;
          choice = char(popup_items(idx,:));
       end

end
