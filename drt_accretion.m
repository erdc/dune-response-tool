function scenario = drt_accretion(scenario)

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
    
    %Utilizing Kawamura Approach 
    D50 = scenario.models.d50/1000;
    K = 0.4;
    z = 10;
    zo = 2*D50/30;
    pa = 1.225; 
    ps = 2650;
    g = 9.81;
    ustar = u_w.*K./log(z/zo);
    C = 1.87; %for typical sands
    M = 0; %assumed moisture content
    Ck = 2.78;
    Ck = scenario.models.AeolianTransportCoefficient;
    ustar_thresh = 0.1 * sqrt(g * D50 * (ps/pa)*(1 + C*M));
    Q = Ck*(pa/g)*(ustar - ustar_thresh).*(ustar + ustar_thresh).*(ustar + ustar_thresh);

    %Modify by the fetch effect
    Fc = 4.38*u_w - 8.23; %per Delgado-Fernandez
    [maxval,imax] = nanmax(zprof);
    [minval,~] = nanmin(zprof);
    for itime = 1:numel(twl)
        if twl(itime)< maxval && twl(itime) > minval %if the twl is within the profile
            xwl = linterp(zprof(1:imax),xprof(1:imax), twl(itime));
            beachwidth(itime) = xtoe-xwl;
        elseif  twl(itime) <= minval %if the twl is lower than the profile even goes, set the beach width to be wide
            beachwidth(itime) = 500; %pick some high number
        else
            beachwidth(itime) = 0;
        end
    end
    beachwidth(beachwidth<0) = 0;
    F = beachwidth.*tand(windDir);
    F(F<0) = 0;
    F(find(isinf(F) ==1)) = 10000;
    Qtot = Q;
    for idx = 1:numel(Qtot)
        if F(idx) < Fc(idx)
            Qtot(idx) = Q(idx).*sin((pi/2)*F(idx)./Fc(idx));
        end
    end
    Qtot(Qtot <0) = 0;

    %lastly deal with flux to dunes based on angle to get flux to dune
    Qtot = Qtot.*cosd(windDir);
    Qtot(Qtot <0) = 0;

    %convert to a volume flux (initially in kg/m/s)
    por = 0.4;
    Qtot_m3m_dt = (Qtot/ps)*[scenario.timing.dt*60*60]/(1-por);

    %convert to outputs
    scenario.accretion.dV = Qtot_m3m_dt;
    

end