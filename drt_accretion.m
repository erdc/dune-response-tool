function scenario = drt_accretion(scenario)
    %drt_accretion: code to run a simple aeolian sediment transport model
    %for calculating wind blown sediment fluxes into coastal dines
    %
    %Required Inputs: 'scenario' structure variable with the 'grids',
    %'erosion', and 'env' variables    
    %
    %Outputs:
    %       scenario.accretion
    
    %set up aeolian model
        twl = scenario.erosion.TWL;
        xprof = scenario.grids.XGrid;
        zprof = scenario.grids.ZGrid;
        dtoe = scenario.grids.morphometrics.dtoe;
        dhigh = scenario.grids.morphometrics.dhigh;
        [~, imax] = nanmax(scenario.grids.ZGrid);
        xtoe = linterp(zprof(1:imax),xprof(1:imax),dtoe);
        u_w = scenario.env.winds.windSpeed;
        windDir = scenario.env.winds.windDirection-scenario.grids.morphometrics.azimuth;
    
    %Utilizing Kawamura (1951) Approach for Wind-Driven Sediment Fluxes
        D50 = scenario.models.d50/1000; %grain size
        K = 0.4; %von karman constant
        z = 10; %assumed elevation of wind measurements
        zo = 2*D50/30;
        pa = 1.225; %air density
        ps = 2650; %sediment density
        g = 9.81; %gravity
        ustar = u_w.*K./log(z/zo); %shear velocity
        C = 1.87; %for typical sands
        M = 0; %assumed moisture content
        Ck = scenario.models.AeolianTransportCoefficient; %model coefficient now user input
        ustar_thresh = 0.1 * sqrt(g * D50 * (ps/pa)*(1 + C*M)); %threshold shear velocity
        Q = Ck*(pa/g)*(ustar - ustar_thresh).*(ustar + ustar_thresh).*(ustar + ustar_thresh);

    %Modify transport by the fetch effect per Delgado-Fernandez, Geomorphology, 2011
        %critical fetch length
        Fc = 4.38*u_w - 8.23; 

        %determine the beach width based on the total water level
        [maxval,imax] = nanmax(zprof);
        [minval,~] = nanmin(zprof);
        for itime = 1:numel(twl)
            if twl(itime)< maxval && twl(itime) > minval %if the twl is within the profile
                xwl = linterp(zprof(1:imax),xprof(1:imax), twl(itime));
                beachwidth(itime) = xtoe-xwl;
            elseif  twl(itime) <= minval %if the twl is lower than the profile even goes, set the beach width to be wide
                beachwidth(itime) = 100; %pick some high number if exceeded
            else
                beachwidth(itime) = 0;
            end
        end
        beachwidth(beachwidth<0) = 0;
        
        %determine the fetch for the specific conditions
        F = beachwidth.*tand(abs(windDir));
        F(F<0) = 0;
        F(find(isinf(F) ==1)) = 1000;
        Qtot = Q; %initialize variable
        for idx = 1:numel(Qtot)
            if F(idx) < Fc(idx)
                Qtot(idx) = Q(idx).*sin((pi/2)*F(idx)./Fc(idx));
            end
        end
        Qtot(Qtot <0) = 0;

        %lastly deal with flux to dunes based on angle to get flux to dune
            Qtot = Qtot.*cosd(abs(windDir));
            Qtot(Qtot <0) = 0;

    %convert to a volume flux (initially in kg/m/s)
        por = 0.4; %assumed porosity
        Qtot_m3m_dt = (Qtot/ps)*[scenario.timing.dt*60*60]/(1-por);

    %convert outputs
        scenario.accretion.dV = Qtot_m3m_dt;
    
end

function yy = linterp(x,y,xx)
    %linear interpolation routine
    nx = max(size(x));
    nxx = max(size(xx));
    yy = zeros(size(xx)); 
    j = 2;
    for i = 1:nxx
       while x(j) < xx(i)
             j = j+1;
       end
       alfa = (xx(i)-x(j-1))/(x(j)-x(j-1));
       yy(i) = y(j-1)+alfa*(y(j)-y(j-1));
    end
end