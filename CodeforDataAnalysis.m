
%% 1: make climate forcing
%%% make climate forcing for ORCHIDEE using ISIMIP
% create date: 2020/07/29
% modify date: 2020/08/02
% Yi

ISIMIP_old_dir = '/workdir/xiyi/DATA/ISIMIP/';
Outdir = '/workdir/xiyi/DATA/ISIMIP_for_ORCHIDEE/';

scens = {'obsclim'; 'historical'; 'ssp126'; 'ssp370'; 'ssp585'};
yrs = [1901 2016; 1850 2014; 2015 2100; 2015 2100; 2015 2100];
yrs_ORCHIDEE = [1960 2016; 1960 2014; 2015 2100; 2015 2100; 2015 2100];

gcms = {{'GSWP3-W5E5'}; {'GFDL-ESM4'; 'UKESM1-0-LL'; 'MPI-ESM1-2-HR'; 'IPSL-CM6A-LR'; 'MRI-ESM2-0'}; ...
    {'GFDL-ESM4'; 'UKESM1-0-LL'; 'MPI-ESM1-2-HR'; 'IPSL-CM6A-LR'; 'MRI-ESM2-0'}; ...
    {'GFDL-ESM4'; 'UKESM1-0-LL'; 'MPI-ESM1-2-HR'; 'IPSL-CM6A-LR'; 'MRI-ESM2-0'}; ...
    {'GFDL-ESM4'; 'UKESM1-0-LL'; 'MPI-ESM1-2-HR'; 'IPSL-CM6A-LR'; 'MRI-ESM2-0'};};

vars_ISIMIP = {'tasmax'; 'tasmin'; 'pr'; 'huss'; 'ps'; 'rsds'; 'rlds'; 'sfcwind'};
vars_ORCHIDEE = {'Tmax'; 'Tmin'; 'precip'; 'Qair'; 'PSurf'; 'SWdown'; 'LWdown'; 'Wind'};
vars_units = {'K'; 'K'; 'kg m-2 s-1'; 'kg kg-1'; 'Pa'; 'W m-2'; 'W m-2'; 'm s-1'};
vars_longname = {'Maximum Daily Air Temperature'; 'Minimum Daily Air Temperature'; ...
    'Precipitation'; 'Specific Humidity'; 'Surface Pressure'; 'Downward Shortwave Radiation'; ...
    'Downward Longwave Radiation'; 'Wind Speed'};

reso = 0.5;
[nb_lat, nb_lon] = deal(180/reso, 360/reso);

% example_file = 'F:\DATA\GSWP3\GSWP3_v1_halfdeg_mon_1964.nc';
example_file = '/workdir/data/others/xiyi/GSWP3/GSWP3_v1_halfdeg_3h_1964.nc';
nav_lat = ncread(example_file, 'nav_lat');
nav_lon = ncread(example_file, 'nav_lon');

occupation = ncread('/workdir/xiyi/DATA/ISIMIP/mask2D_67420.nc', 'occupation');
Areas = nan(size(occupation));
Areas(isnan(Areas)) = 1.000000020040877e+20;

mask = ncread('/workdir/xiyi/DATA/ISIMIP/mask2D_67420.nc', 'mask');
contfrac = nan(size(mask));
contfrac(mask>=0) = 1;
maskr = contfrac;
contfrac(isnan(contfrac)) = 1.000000020040877e+20;

for isc = 1:length(scens)

    gcms_isc = gcms{isc, 1};
    yrs_interval = [yrs(isc, 1) ceil(yrs(isc, 1)/10)*10+1:10:yrs(isc, 2); ceil(yrs(isc, 1)/10)*10:10:yrs(isc, 2)-1 yrs(isc, 2)]';

    for igcm = 1:length(gcms_isc)

        ISIMIP_old_subdir = [ISIMIP_old_dir, scens{isc}, '/', gcms_isc{igcm}, '/'];
        ISIMIP_new_subdir = [Outdir, scens{isc}, '/', gcms_isc{igcm}, '/'];
        if ~exist(ISIMIP_new_subdir, 'dir')
            mkdir(ISIMIP_new_subdir)
        end

        for iyr = yrs_ORCHIDEE(isc, 1):yrs_ORCHIDEE(isc, 2)

            disp(iyr)
            tic
            filename_new = [ISIMIP_new_subdir, gcms_isc{igcm}, '_halfdeg_daily_', num2str(iyr), '.nc'];

            days_iyr = datenum(iyr, 12, 31)-datenum(iyr, 1, 1)+1;
            days_begin = datenum(iyr, 1, 1)-datenum(yrs_interval(find(yrs_interval(:, 2)>=iyr, 1), 1), 1, 1)+1;
            tdaysecond = [0:86400:(days_iyr-1)*86400]';

            nb_days = length(tdaysecond);

            if exist(filename_new, 'file')
                delete(filename_new)
            end
            % 1: lat
            nccreate(filename_new, 'nav_lat', 'Dimensions', {'longitude', nb_lon, 'latitude', nb_lat}, ...
                'Format', 'netcdf4', 'Datatype', 'single');
            ncwriteatt(filename_new, 'nav_lat', 'long_name', 'Latitude');
            ncwriteatt(filename_new, 'nav_lat', 'units', 'degrees north');
            ncwrite(filename_new, 'nav_lat', nav_lat)

            % 2: lon
            nccreate(filename_new, 'nav_lon', 'Dimensions', {'longitude', nb_lon, 'latitude', nb_lat}, ...
                'Format', 'netcdf4', 'Datatype', 'single');
            ncwriteatt(filename_new, 'nav_lon', 'long_name', 'Longitude');
            ncwriteatt(filename_new, 'nav_lon', 'units', 'degrees east');
            ncwrite(filename_new, 'nav_lon', nav_lon);

            % 3: time
            nccreate(filename_new, 'time', 'Dimensions', {'time', nb_days, ...
                netcdf.getConstant('NC_UNLIMITED')}, 'Format', 'netcdf4', ...
                'Datatype', 'single', 'FillValue', 1.000000020040877e+20);
            ncwriteatt(filename_new, 'time', 'title', 'Time');
            ncwriteatt(filename_new, 'time', 'long_name', 'Time axis');
            ncwriteatt(filename_new, 'time', 'calendar', 'standard');
            ncwriteatt(filename_new, 'time', 'units', ['seconds since ', num2str(iyr), '-01-01 00:00:00']);
            ncwriteatt(filename_new, 'time', 'time_origin', [num2str(iyr), '-01-01 00:00:00']);
            ncwrite(filename_new, 'time', tdaysecond);

            % 4: contfrac
            nccreate(filename_new, 'contfrac', 'Dimensions', {'longitude', nb_lon, 'latitude', nb_lat}, ...
                'Format', 'netcdf4', 'Datatype', 'single', 'FillValue', 1.000000020040877e+20);
            ncwriteatt(filename_new, 'contfrac', 'long_name', 'Continental fraction');
            ncwriteatt(filename_new, 'contfrac', 'units', 1);
            ncwriteatt(filename_new, 'contfrac', 'coordinates', 'lon lat');
            ncwrite(filename_new, 'contfrac', contfrac);

            % 5: Areas
            nccreate(filename_new, 'Areas', 'Dimensions', {'longitude', nb_lon, 'latitude', nb_lat}, ...
                'Format', 'netcdf4', 'Datatype', 'single', 'FillValue', 1.000000020040877e+20);
            ncwriteatt(filename_new, 'Areas', 'long_name', 'Mesh areas');
            ncwriteatt(filename_new, 'Areas', 'units', 'm2');
            ncwriteatt(filename_new, 'Areas', 'coordinates', 'lon lat');
            ncwrite(filename_new, 'Areas', Areas);


            for ivar = 1:length(vars_ISIMIP)

                filename_old = dir([ISIMIP_old_subdir, '*_', vars_ISIMIP{ivar}, '_*',...
                    num2str(yrs_interval(find(yrs_interval(:, 2)>=iyr, 1), 1)), '*.nc']);

                % 6: other climate vars
                var_temp = ncread([filename_old.folder, '/', filename_old.name], ...
                    vars_ISIMIP{ivar}, [1, 1, days_begin], [inf inf days_iyr], [1 1 1]);
                var_temp = var_temp.*maskr;
                var_temp(isnan(var_temp)) = 1.000000020040877e+20;

                nccreate(filename_new, vars_ORCHIDEE{ivar}, 'Dimensions', ...
                    {'longitude', nb_lon, 'latitude', nb_lat, 'time', nb_days}, ...
                    'Format', 'netcdf4', 'Datatype', 'single', 'FillValue', 1.000000020040877e+20);
                ncwriteatt(filename_new, vars_ORCHIDEE{ivar}, 'long_name', vars_longname{ivar});
                ncwriteatt(filename_new, vars_ORCHIDEE{ivar}, 'units', vars_units{ivar});
                % ncwriteatt(filename_new, vars_ORCHIDEE{ivar}, 'cell_methods', cell_methods{ivar});
                ncwrite(filename_new, vars_ORCHIDEE{ivar}, var_temp);
                % clear var_temp

            end
            toc

            % 7: Global Attributes
            ncwriteatt(filename_new, '/', 'contact', 'yixi@pku.edu.cn');
            ncwriteatt(filename_new, '/', 'create date', '2020-07-29');

            % ncdisp(filename_new)

        end

    end

end



%% 2: extract vars from output file of ORCHIDEE
reso = 0.5;
[nb_lat, nb_lon] = deal(40/reso, 70/reso);

vars = {'rain'; 'evap'; 'evapor'; 'transpir'; 'runoff'; 'mrso'; 'drainage'; 'surface_runoff'; ...
    'mrsos'; 'TWBR'; 'humtot_low'; 'humtot_up'; 'npp'; 'gpp'; 'lai'; 'rivo'; 'evapot_corr'};

% mask
PROVmask = imread('G:\DATA\Tiff\China\province_reso05.tif');
PROVmask = double(PROVmask);

mask_cn = PROVmask;
mask_cn(mask_cn==-9999) = nan;
mask_cn(~isnan(mask_cn)) = 1;
clear PROVmask

s = area_weighted_any([15 55], [70 140], reso); % km2
s = s.*mask_cn;

depth_st = [0.0005 0.0020 0.0059 0.0137 0.0293 0.0606 ...
    0.1232 0.2483 0.4985 0.9990 1.7498 2.5005 3.5015 ...
    4.5525 5.6561 6.8148 8.0315 9.3091 10.6505 12.0589 ...
    13.5378 15.0907 16.7212 18.4332 20.2308 22.1183 ...
    24.1001 26.1811 28.3661 30.6604 33.0693 35.5988]; % node

load('F:\11Forest_Wetland\pro\CN\202010\md_sces.mat', 'md_sces');

for m = 1:size(md_sces, 1)
    
    disp(md_sces(m, 'periods'))
    disp(md_sces(m, 'models'))
    
    tic
    mddir = md_sces(m, 'directory2');
    mddir = mddir.directory2;
    mddir = mddir{1, 1};
    
    yb = md_sces(m, 'yrb');
    yb = yb.yrb;
    ye = md_sces(m, 'yre');
    ye = ye.yre;
    nb_yr = ye-yb+1;
    
    GWdir = ['G:\Simulations\Forest_Wetland\ORCHIDEE-GW\', mddir];
    SAVEdir = ['F:\11Forest_Wetland\data\CN\Sims\ORCHIDEE-GW\', mddir];
    
    for ivar = 1:length(vars)
        cprintf('*magenta', [vars{ivar}, '\n']);
        var_ts = nan(12, nb_yr);
        var_pattern = nan(nb_lat, nb_lon, 12, nb_yr);
        for iyr = yb:ye
            if mod(iyr,10)==0
                % disp(iyr)
            end
            filename = [GWdir, 'sechiba_history_igem_', num2str(iyr), '.nc'];
            if strcmp(vars{ivar}, 'rain')
                p = ncread(filename, 'pr');
                p = p*60*60*24; % mm/s--->mm/d
                p = flipud(rot90(p, 1));
                % mm/d--->mm/mon
                p = p.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(p.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = p; % mm/mon
                clear rain p
            elseif strcmp(vars{ivar}, 'evap')
                evap = ncread(filename, 'et');
                evap = evap*60*60*24; % mm/s--->mm/d
                evap = flipud(rot90(evap, 1));
                % mm/d--->mm/mon
                evap = evap.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(evap.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = evap; % mm/mon
                clear evap
            elseif strcmp(vars{ivar}, 'evapot_corr')
                evapot_corr = ncread(filename, 'evapot_corr');
                evapot_corr = evapot_corr*60*60*24; % mm/s--->mm/d
                evapot_corr = flipud(rot90(evapot_corr, 1));
                % mm/d--->mm/mon
                evapot_corr = evapot_corr.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(evapot_corr.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = evapot_corr; % mm/mon
                clear evap
            elseif strcmp(vars{ivar}, 'runoff')
                q = ncread(filename, 'mrro');
                q = q*60*60*24; % mm/s--->mm/d
                q = flipud(rot90(q, 1));
                % mm/d--->mm/mon
                q = q.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(q.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = q; % mm/mon
                clear q
            elseif strcmp(vars{ivar}, 'drainage')
                q = ncread(filename, 'mrrob');
                q = q*60*60*24; % mm/s--->mm/d
                q = flipud(rot90(q, 1));
                % mm/d--->mm/mon
                q = q.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(q.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = q; % mm/mon
                clear q
            elseif strcmp(vars{ivar}, 'surface_runoff')
                q = ncread(filename, 'mrros');
                q = q*60*60*24; % mm/s--->mm/d
                q = flipud(rot90(q, 1));
                % mm/d--->mm/mon
                q = q.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(q.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = q; % mm/mon
                clear q
            elseif strcmp(vars{ivar}, 'evapor')
                evap = ncread(filename, 'et');
                transpir = ncread(filename, 'tran');
                evapor = evap-transpir;
                evapor = flipud(rot90(evapor, 1));
                evapor = evapor*60*60*24; % mm/s--->mm/d
                % mm/d--->mm/mon
                evapor = evapor.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(evapor.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = evapor; % mm/mon
                clear evapnu
            elseif contains(vars{ivar}, 'transpir')
                transpir = ncread(filename, 'tran');
                transpir = flipud(rot90(transpir, 1));
                transpir = transpir*60*60*24; % mm/s--->mm/d
                % mm/d--->mm/mon
                transpir = transpir.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(transpir.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = transpir; % mm/mon
                clear transpir
            elseif strcmp(vars{ivar}, 'mrso')
                mrso = ncread(filename, 'humtot');
                mrso = flipud(rot90(mrso, 1)); % kg/m2
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(mrso.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = mrso; % mm/mon
                clear mrso
            elseif strcmp(vars{ivar}, 'TWBR')
                twbr = ncread(filename, 'TWBR');
                twbr = flipud(rot90(twbr, 1)); % kg/m2/s
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(twbr.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = twbr; % mm/mon
                clear mrsos
            elseif strcmp(vars{ivar}, 'mrsos')
                mrsos = ncread(filename, 'mrsos');
                mrsos = flipud(rot90(mrsos, 1)); % kg/m2
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(mrsos.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = mrsos; % mm/mon
                clear mrsos
            elseif strcmp(vars{ivar}, 'humtot_low')
                humtot_ns = ncread(filename, 'humtot_ns');
                humtot_ns = squeeze(humtot_ns(:, :, 4, :));
                humtot_ns = flipud(rot90(humtot_ns, 1)); % kg/m2
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(humtot_ns.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = humtot_ns; % mm/mon
                clear humtot_ns
            elseif strcmp(vars{ivar}, 'humtot_up')
                humtot_ns = ncread(filename, 'humtot_ns');
                soiltile = ncread(filename, 'soiltile');
                humtot_up = squeeze(nansum(humtot_ns(:, :, 1:3, :).*soiltile(:, :, 1:3), 3))./...
                    nansum(soiltile(:, :, 1:3), 3);
                humtot_up = flipud(rot90(humtot_up, 1)); % kg/m2
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(humtot_up.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = humtot_up; % mm/mon
                clear humtot_ns soiltile humtot_up
            elseif contains(vars{ivar}, 'gpp')
                gpp = ncread(filename, 'gpp');
                gpp = gpp*60*60*24*1e3; % kgC/m2/s-->gC/m2/d
                gpp = flipud(rot90(gpp, 1));
                % gC/m2/d--->gC/m2/mon
                gpp = gpp.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(gpp.*s, 1), 2))/1e15; % PgC/mon
                var_pattern(:, :, :, iyr-yb+1) = gpp; % gC/m2/mon
                clear gpp
            elseif strcmp(vars{ivar}, 'npp')
                npp = ncread(filename, 'npp');
                npp = npp*60*60*24*1e3; % kgC/m2/s-->gC/m2/d
                npp = flipud(rot90(npp, 1));
                % gC/m2/d--->gC/m2/mon
                npp = npp.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(npp.*s, 1), 2))/1e15; % PgC/mon
                var_pattern(:, :, :, iyr-yb+1) = npp; % gC/m2/mon
                clear npp
            elseif contains(vars{ivar}, 'lai')
                lai = ncread(filename, 'lai');
                lai = flipud(rot90(lai, 1));
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(lai.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = lai;
                clear lai maxvegetfrac
            elseif strcmp(vars{ivar}, 'rivo')
                rivo = ncread(filename, 'rivo');
                rivo = rivo*60*60*24/1e9; % m^3/s-->km^3/d
                rivo = flipud(rot90(rivo, 1));
                % km^3/d--->km^3/mon
                rivo = rivo.*permute(repmat(eomday(iyr, 1:12), [1, nb_lat, nb_lon]), [2 3 1]);
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(rivo.*s, 1), 2))/1e15; % km^3/mon
                var_pattern(:, :, :, iyr-yb+1) = rivo; % km^3/mon
                clear rivo
            elseif strcmp(vars{ivar}, 'tsl_70cm')
                tsl = ncread(filename, 'tsl'); % 140*80*32*12
                tsl = flipud(rot90(tsl, 1)); % K
                
                tsl_70cm = nan(nb_lat, nb_lon, 12);
                for mm = 1:12
                    for ilat = 1:nb_lat
                        for ilon = 1:nb_lon
                            if ~isempty(find(~isnan(squeeze(nanmean(tsl(ilat, ilon, :, mm), 3))), 1))
                                tsl_70cm(ilat, ilon, mm) = interp1(depth_st, squeeze(tsl(ilat, ilon, :, mm)), 0.7, 'linear');
                            end
                        end
                    end
                end
                var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(tsl_70cm.*s, 1), 2)./...
                    nansum(s(:)));
                var_pattern(:, :, :, iyr-yb+1) = tsl_70cm;
                clear tsl tsl_70cm
            end
            
        end
        % var_pattern = nanmean(var_pattern, 3); % mean annual
        if ~exist(SAVEdir, 'dir')
            mkdir(SAVEdir)
        end
        save([SAVEdir, vars{ivar}, '.mat'], 'var_ts', 'var_pattern');
    end
    toc    
end





%% 3: calculate WTD
reso = 0.5;
[nb_lat, nb_lon] = deal(40/reso, 70/reso);

vars = {'WTD'};

SMsat = ncread(['G:\Simulations\Forest_Wetland\ORCHIDEE-GW\', ...
    'sechiba_history_2000.nc'], 'tmcs');
SMsat = flipud(rot90(SMsat, 1)); % kg/m2

h = 2; % total soil depth

depth_st = [0.0005 0.0020 0.0059 0.0137 0.0293 0.0606 ...
    0.1232 0.2483 0.4985 0.9990 1.7498 2.5005 3.5015 ...
    4.5525 5.6561 6.8148 8.0315 9.3091 10.6505 12.0589 ...
    13.5378 15.0907 16.7212 18.4332 20.2308 22.1183 ...
    24.1001 26.1811 28.3661 30.6604 33.0693 35.5988]; % node

% mask
PROVmask = imread('G:\DATA\Tiff\China\province_reso05.tif');
PROVmask = double(PROVmask);

mask_cn = PROVmask;
mask_cn(mask_cn==-9999) = nan;
mask_cn(~isnan(mask_cn)) = 1;
clear PROVmask

s = area_weighted_CN(nb_lat, nb_lon); % m2
s = s.*mask_cn;

depth_st = [0.0005 0.0020 0.0059 0.0137 0.0293 0.0606 ...
    0.1232 0.2483 0.4985 0.9990 1.7498 2.5005 3.5015 ...
    4.5525 5.6561 6.8148 8.0315 9.3091 10.6505 12.0589 ...
    13.5378 15.0907 16.7212 18.4332 20.2308 22.1183 ...
    24.1001 26.1811 28.3661 30.6604 33.0693 35.5988]; % node

% mddir = 'obsclim_GSWP3-W5E5_S1_spinup2000_RFW_20210310\';
% mddir = 'obsclim_GSWP3-W5E5_S0_spinup2000_RFW_20210319\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_newLUC_RFW_20210323\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_newLUC_RFW_20210325\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_allwet_RFW_20210323\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_alldry_RFW_20210323\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_alldry_RFW_20210325\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_allwet_RFW_20210325\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_newLUC_RFW_20210329\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_newLUC2016_RFW_20210327\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_gt2_RFW_20210329\';
mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_1to2_RFW_20210329\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_lt1_RFW_20210329\';


yb = 2017;
ye = 2035;
nb_yr = ye-yb+1;

GWdir = ['G:\Simulations\Forest_Wetland\ORCHIDEE-GW\', mddir];
SAVEdir = ['F:\11Forest_Wetland\data\CN\Sims\ORCHIDEE-GW\', mddir];

for ivar = 1:length(vars)
    cprintf('*magenta', [vars{ivar}, '\n']);
    var_ts = nan(12, nb_yr);
    var_pattern = nan(nb_lat, nb_lon, 12, nb_yr);
    for iyr = yb:ye
        if mod(iyr,10)==0
            disp(iyr)
        end
        filename = [GWdir, 'sechiba_history_igem_', num2str(iyr), '.nc'];
        if strcmp(vars{ivar}, 'WTD')
            mrso = ncread(filename, 'humtot');
            mrso = flipud(rot90(mrso, 1)); % kg/m2
            
            tsl = ncread(filename, 'tsl'); % 140*80*32*12
            tsl = flipud(rot90(tsl, 1)); % K
            tsl_70cm = nan(nb_lat, nb_lon, 12);
            for mm = 1:12
                for ilat = 1:nb_lat
                    for ilon = 1:nb_lon
                        if ~isempty(find(~isnan(squeeze(nanmean(tsl(ilat, ilon, :, mm), 3))), 1))
                            tsl_70cm(ilat, ilon, mm) = ...
                                interp1(depth_st, squeeze(tsl(ilat, ilon, :, mm)), 0.7, 'linear');
                        end
                    end
                end
            end
            
            wtd = h*(mrso./SMsat-1);
            wtd_above = wtd.*(SMsat/1e3/h); % kg/m2-->m3/m3
            wtd(wtd>0) = wtd_above(wtd>0);
            wtd(tsl_70cm<=273.15) = nan;
            
            var_ts(:, iyr-yb+1) = squeeze(nansum(nansum(wtd.*s, 1), 2))/nansum(s(:)); % m
            var_pattern(:, :, :, iyr-yb+1) = wtd; % m
            clear mrso
        end
    end
    
    if ~exist(SAVEdir, 'dir')
        mkdir(SAVEdir)
    end
    save([SAVEdir, vars{ivar}, '.mat'], 'var_ts', 'var_pattern');
    
end






%% 4: parameter calibaration, using the obsclim S1
% the parameter should be the same for all scenarios
addpath('F:\11Forest_Wetland\pro\CN\202006\')

% 1: WTD-->Fwet
%%%: choose a set of benchmark data and                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                the method
benchmarkdata = 'RFW'; % {'RFW'; 'GIEMS2'}
method = 'max'; % {'max'; 'yrmax'}
% yrmax for GIEMS-2, period: 2000-2015

reso = 0.5;
[nb_lat, nb_lon] = deal(40/reso, 70/reso);

[ln, col] = extract_spatial([90 -90], [-180 180], [55 15], [70 140], reso);

%%%%%% 0.5 degree
BENCHdir = 'F:\1InundatedArea_180219\data\BenchmarkData\';
HYDEdir = 'G:\DATA\HYDE\';
if strcmp(benchmarkdata, 'RFW') == 1
    load([BENCHdir, 'FwetRFW_new.mat'], 'FwetRFW_reso05');
    %     load('F:\1InundatedArea_180219\data\BenchmarkData\MIRCA_rice_reso05.mat', ...
    %         'MIRCA_rice_2000');
    %     fwet_obs = FwetRFW_reso05-max(MIRCA_rice_2000, [], 3, 'omitnan');
    load([HYDEdir, 'rice_cn_2000_2017.mat'], 'rice_baseline_corr');
    rice_baseline_corr_max = max(rice_baseline_corr, [], 3, 'omitnan');
    fwet_obs = FwetRFW_reso05(ln(1):ln(2), col(1):col(2))-rice_baseline_corr_max;
    fwet_obs(fwet_obs<0) = 0;
    clear FwetRFW* MIRCA_rice_2000;
elseif strcmp(benchmarkdata, 'GIEMS2') == 1
    load([BENCHdir, 'FwetGIEMS2.mat'], 'FwetGIEMS2_reso05');
    % load([BENCHdir, 'MIRCA_rice_reso05.mat'], 'MIRCA_rice_2000');
    load([HYDEdir, 'rice_cn_2000_2017.mat'], 'rice_baseline_corr');
    FwetGIEMS2_reso05 = reshape(FwetGIEMS2_reso05, [360, 720 12 288/12]);
    FwetGIEMS2_reso05 = FwetGIEMS2_reso05(:, :, :, [2000:2015]-1992+1);
    rice_baseline_corr = rice_baseline_corr(:, :, [2000:2015]-2000+1);
    FwetGIEMS2_reso05 = squeeze(max(FwetGIEMS2_reso05, [], 3, 'omitnan'));
    FwetGIEMS2_reso05 = FwetGIEMS2_reso05(ln(1):ln(2), col(1):col(2), :);
    FwetGIEMS2_reso05 = FwetGIEMS2_reso05-rice_baseline_corr;
    FwetGIEMS2_reso05(FwetGIEMS2_reso05<0) = 0;
    fwet_obs = FwetGIEMS2_reso05;
    if strcmp(method, 'max') ==1
        fwet_obs = squeeze(max(FwetGIEMS2_reso05, [], 3, 'omitnan'));
    end
end

% fwet_obs = fwet_obs(ln(1):ln(2), col(1):col(2), :);
% figure; imagesc(fwet_obs)

% wtnames = {'obsclim_GSWP3-W5E5_S0_spinup2000_RFW_20210319'};
wtnames = {'obsclim_GSWP3-W5E5_S1_spinup2000_RFW_20210310'};

wtyb = 2017;
wtye = 2035;
WTDdir = ['F:\11Forest_Wetland\data\CN\Sims\ORCHIDEE-GW\', wtnames{1}, '\'];
PARAdir = 'F:\1InundatedArea_180219\data\Parameter\';
Outdir = 'F:\11Forest_Wetland\data\CN\Sims\TOPMODEL\Fwet\';


tic
%%% 1: input to optimize parameters
if strcmp(method, 'max') ==0
    if strcmp(benchmarkdata, 'SWAMPS') == 1
        obsyb0 = 2000; % !!!!!
        obsye0 = 2018; % !!!!! revise, 2019/09/06
    end
    if contains(benchmarkdata, 'GIEMS2') == 1
        obsyb0 = 2000; % !!!!!
        obsye0 = 2015; % !!!!! revise, 2020/02/21
    end
    if strcmp(benchmarkdata, 'GIEMS2part') == 1
        obsyb0 = 1992; % !!!!!
        obsye0 = 2005; % !!!!! revise, 2020/02/21
    end
    modyb0 = wtyb;
    modye0 = wtye;
    if modyb0 < obsyb0
        modyb0 = obsyb0;
    elseif modyb0 > obsye0
        error(message('No observation for Water Table begin year!'));
    end
    if modye0 > obsye0
        modye0 = obsye0;
    elseif modyb0 < obsyb0
        error(message('No observation for Water Table end year!'));
    end
    load([WTDdir, 'WTD.mat'], 'var_pattern');
    wtd = var_pattern;
    clear var_pattern
else
    modyb0 = wtyb;
    modye0 = wtye;
    load([WTDdir, 'WTD.mat'], 'var_pattern');
    wtd = var_pattern;
    clear var_pattern
end

wtd = wtd(:, :, :, (modyb0:modye0)-wtyb+1);

wtd_temp = reshape(wtd, [size(wtd, 1), size(wtd, 2), size(wtd, 3)*size(wtd, 4)]);
if strcmp(method, 'monmean') == 1
    wtd_new = reshape(wtd_temp, [size(wtd_temp, 1), size(wtd_temp, 2), 12, size(wtd_temp, 3)/12]);
    wtd_new = squeeze(nanmean(wtd_new, 4));   % mothly mean
elseif strcmp(method, 'yrmax') == 1
    wtd_new = reshape(wtd_temp, [size(wtd_temp, 1), size(wtd_temp, 2), 12, size(wtd_temp, 3)/12]);
    wtd_new = squeeze(max(wtd_new, [], 3, 'omitnan'));   % annual maximum
elseif strcmp(method, 'max') ==1
    [wtd_new, ind] = max(wtd_temp, [], 3, 'omitnan');
end
clear wt wtd wtd_temp;



%%% 2: observation: maximum, timeseries 180, yearly maximum, mean
%%% seasonal cycle
fwet_obs_temp = fwet_obs;




%%% 3: prepare all parameter sets
nb_param = 15;
para = nan(nb_lat, nb_lon, nb_param, 3);
index = nan(nb_param, 1);
n = 1;
for M = 1:15
    load([PARAdir, 'reso', num2str(reso), '_new_20191013\RFW_RFW-G17_SWAMPS\',...
        'gridPara_reso', num2str(reso), '_M', num2str(M), '.mat'], ...
        'vp', 'kp', 'qp');
    vp = rot90(vp, 1);
    qp = rot90(qp, 1);
    kp = rot90(kp, 1);
    para(:, :, n, 1) = vp(ln(1):ln(2), col(1):col(2));
    para(:, :, n, 2) = kp(ln(1):ln(2), col(1):col(2));
    para(:, :, n, 3) = qp(ln(1):ln(2), col(1):col(2));
    index(n, 1) = M;
    n = n+1;
end



%%% 4: compare mode and obs to calibrate parameters
vkq_fmax_opt = nan(nb_lat, nb_lon, 4);
RMSE_opt = nan(nb_lat, nb_lon);
M_opt = nan(nb_lat, nb_lon);
nb_min = nan(nb_lat, nb_lon);
% fwet_opt = nan([size(fwet_obs_temp), nb_param]);
for i=1:nb_lat
    for j=1:nb_lon
        wt_vec = squeeze(wtd_new(i, j, :)); %%% from 1
        fwet_obs_vec = squeeze(fwet_obs_temp(i, j, :)); %%% from 2
        vkq_p = squeeze(para(i, j, :, :)); %%% from 3
        ind = isnan(wt_vec);
        wt_vec(ind)=[];
        fwet_obs_vec(ind)=[];
        ind = isnan(fwet_obs_vec);
        % ind = fwet_obs_vec==0;
        wt_vec(ind)=[];
        fwet_obs_vec(ind)=[];
        if ~isempty(wt_vec)
            % ######
            % monthly mean, we need long-term max!!
            fwet_obs_max = max(squeeze(fwet_obs(i, j, :)));
            % ######
            
            [vkq_out, RMSE, ind2, num_min] = minParam_new(vkq_p, wt_vec, fwet_obs_vec, fwet_obs_max);
            
            vkq_fmax_opt(i, j, :) = vkq_out;
            RMSE_opt(i, j) = RMSE;
            M_opt(i, j) = index(ind2, 1);
            nb_min(i, j) = num_min;
        end
    end
    % disp(i)
end

if ~exist([Outdir, 'ParaOPT_new\'], 'dir')
    mkdir([Outdir, 'ParaOPT_new\']);
end
save([Outdir, 'ParaOPT_new\', wtnames{1}, '_', 'reso', num2str(reso), ...
    '_', benchmarkdata, '_', method, '_optPara_RMSE.mat'], ...
    'vkq_fmax_opt', 'RMSE_opt', 'M_opt', 'nb_min');
clear wtd* para fwet_obs_temp







%% 5: calculate wetland fraction
Outdir = 'F:\11Forest_Wetland\data\CN\Sims\TOPMODEL\Fwet\';

benchmarkdata = 'GIEMS2'; % {'RFW'; 'JRC_withpermannet'; 'RFW-G17'; 'GIEMS2'}
method = 'yrmax'; % {'max'; 'yrmax'}
reso = 0.5;
[nb_lat, nb_lon] = deal(40/reso, 70/reso);

[ln, col] = extract_spatial([90 -90], [-180 180], [55 15], [70 140], reso);

% soil freeze thaw
% SFT_Koppendir = 'F:\1InundatedArea_180219\data\Mask_SFT_Koppen\';
% sftyb = 1979;
% sftye = 2017;
% load([SFT_Koppendir, 'soil_freeze_thaw_reso', num2str(reso), ...
%     '_', num2str(sftyb), num2str(sftye), '_023.mat'], 'sft');
% sft_threshold = 5;
% sft = sft(ln(1):ln(2), col(1):col(2), :, :);

load([Outdir, 'ParaOPT_new\obsclim_GSWP3-W5E5_S1_spinup2000_RFW_20210310_', ...
    'reso', num2str(reso), '_', benchmarkdata, '_', method, '_optPara_RMSE.mat'], ...
    'vkq_fmax_opt', 'RMSE_opt', 'M_opt', 'nb_min');

% mddir = 'obsclim_GSWP3-W5E5_S1_spinup2000_RFW_20210310\';
% mddir = 'obsclim_GSWP3-W5E5_S0_spinup2000_RFW_20210319\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_newLUC_RFW_20210323\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_newLUC_RFW_20210325\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_allwet_RFW_20210323\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_alldry_RFW_20210323\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_alldry_RFW_20210325\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_allwet_RFW_20210325\';
mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_newLUC_RFW_20210329\';
% mddir = 'prc26_GSWP3-W5E5_S0_spinup2000_randomclimate_newLUC2016_RFW_20210327\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_gt2_RFW_20210329\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_1to2_RFW_20210329\';
% mddir = 'prc26_GSWP3-W5E5_S2_spinup2000_randomclimate_lt1_RFW_20210329\';

yb = 2017;
ye = 2035;
nb_yr = ye-yb+1;

WTDdir = ['F:\11Forest_Wetland\data\CN\Sims\ORCHIDEE-GW\', mddir];
md = extractBefore(mddir, '\');

load([WTDdir, 'WTD.mat'], 'var_pattern');
for yr = yb:ye
    for mm=1:12
        wt = var_pattern(:, :, mm, yr-yb+1);
        vp = vkq_fmax_opt(:, :, 1);
        kp = vkq_fmax_opt(:, :, 2);
        qp = vkq_fmax_opt(:, :, 3);
        fmax = vkq_fmax_opt(:, :, 4);
        x = (1 + vp .* exp(-kp .* (wt - qp)));
        x(x<0) = nan;
        x = x.^ (-1 ./ vp);
        fmax(isnan(x)) = NaN;
        fwet_opt = min(x, fmax);
        if ~exist([Outdir, md, '\', benchmarkdata, '\', method], 'dir')
            mkdir([Outdir, md, '\', benchmarkdata, '\', method])
        end
        save([Outdir, md, '\', benchmarkdata, '\', method, '\', md, '_', ...
            benchmarkdata, '_', method, '_fwet_', num2str(yr), sprintf('%.2d', mm), '.mat'], ...
            'fwet_opt');
    end
end



