function scenario = drt_model(scenario)
    %drt_erosion: code to run the Palmsten and Holman (2012) dune retreat
    %model
      
    %initialize model in PH12 conventions
        xM = scenario.grids.XGrid+max(abs(scenario.grids.XGrid));
        xM = scenario.grids.XGrid;
        z = scenario.grids.ZGrid;
        time = [scenario.timing.times-scenario.timing.times(1)]*24*60*60; %converted to seconds
        WL = scenario.env.tides.wl;
        Ho = scenario.env.waves.Hs_25m;
        Lo = scenario.env.waves.L_25m;
        T = scenario.env.waves.Tp;
        Bo = scenario.grids.morphometrics.backshoreslope;
        dtoe = scenario.grids.morphometrics.dtoe;
        D50 = scenario.models.d50;
        WaveRunupFactor = scenario.models.WaveRunupFactor;
        DuneSlopeTrajectory = scenario.models.DuneSlopeTrajectory;
        DuneErodibility = scenario.models.DuneErodibility;   
       zShore = scenario.grids.morphometrics.zshore;
       shorechange = scenario.grids.morphometrics.shorechange;

    %run model
        output_times = linspace(time(1), time(end), 500);    
    
    %output_times 
        out = run_erosion_accretion_model(xM(:)',z(:)',time(:),WL(:),Ho(:),Lo(:),T(:),Bo,dtoe, output_times(:), D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility, zShore, shorechange, scenario);

    %store relevant model output and send back to main program
        scenario.output.dV_erosion = out.dV_erosion;
        scenario.output.dV_accretion = out.dV_accretion;
        scenario.output.dV_accretion2 = out.dV_accretion2;
        scenario.output.Z = out.zNew';
        scenario.output.TWL = out.TWL;
        scenario.output.xToe = out.xtoe;
        scenario.output.ztoe = out.ztoe;
        scenario.output.Nc = out.Nc;
        scenario.output.times = out.times;    
        scenario.output.zmat_times = [out.zmat_time./86400]+scenario.timing.times(1);
        scenario.output.zRunup = out.zRunup;
        scenario.output.beta = out.Beta;
        
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


function [zb] = avalanche(xprof,zb,tanalpha)
    dza=zeros(size(xprof));
    dx=abs(xprof(2)-xprof(1));

    for i=1:length(xprof)-1
        dz=zb(i+1)-zb(i);
        if abs(dz)>tanalpha*dx
            ddz=.5*(abs(dz)-tanalpha*dx)*sign(dz);
            dza(i)=dza(i)+ddz;
            dza(i+1)=dza(i+1)-ddz;
        end
    end
    zb=zb+dza;

end

function [scenario, zprof_update, Qtot_m3m_dt, dVacc2] = drt_accretion_instantaneous(scenario, zin, ztoe, tt)

    %drt_accretion: code to run a simple aeolian sediment transport model
    %for calculating wind blown sediment fluxes into coastal dines
    %
    %Required Inputs: 'scenario' structure variable with the 'grids',
    %'erosion', and 'env' variables    
    %
    %Outputs:
    %       scenario.accretion
    
    %set up aeolian model
        twl = scenario.erosion.TWL(tt);
        xprof = scenario.grids.XGrid;
        zprof = zin;
        dtoe = ztoe;
        dhigh = max(zin,[], 'omitnan');
        [~, imax] = max(zin,[], 'omitnan');
        if [imax+10]<numel(zin)
        imax = imax+10;
        end
        try
            xtoe = linterp(zprof,xprof,dtoe);          
        catch err
            try
              xtoe = interp1(zprof(1:imax),xprof(1:imax),dtoe);
            catch err
                try
                    xtoe = interp1(zprof(1:imax),xprof(1:imax),dtoe);
                catch err
                    ifind = find(zprof>=dtoe);
                    xtoe = xprof(ifind(1));
                end
            end
            error('Could Not Find Dune Toe')
        end       

    ifind = find(zprof<dtoe);
    xtoe = xprof(max(ifind)+1);
   
    %Utilizing Kawamura (1951) Approach for Wind-Driven Sediment Fluxes
        u_w = scenario.env.winds.windSpeed(tt);
        windDir = scenario.env.winds.windDirection(tt)-scenario.grids.morphometrics.azimuth;
        D50 = scenario.models.d50; %grain size
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
        [maxval,imax] = max(zprof,[], 'omitnan');
        [minval,~] = min(zprof,[],'omitnan');
           if twl< maxval && twl> minval %if the twl is within the profile
                xwl = linterp(zprof(1:imax),xprof(1:imax), twl);
                beachwidth = xtoe-xwl;
            elseif  twl <= minval %if the twl is lower than the profile even goes, set the beach width to be wide
                beachwidth = 100; %pick some high number if exceeded
            else
                beachwidth= 0;
            end
        beachwidth(beachwidth<0) = 0;

        %determine the fetch for the specific conditions
        F = beachwidth./cosd(abs(windDir));
        F(F<0) = 0;
        F(find(isinf(F) ==1)) = 1000;
        Qtot = Q; %initialize variable
            if F< Fc
                Qtot = Q.*sin((pi/2)*F./Fc);
            end
        Qtot(Qtot <0) = 0;

        %lastly deal with flux to dunes based on angle to get flux to dune
            Qtot = Qtot.*cosd(abs(windDir));
            Qtot(Qtot <0) = 0;

        %convert to a volume flux (initially in kg/m/s)
        por = 0.4; %assumed porosity
        Qtot_m3m_dt = (Qtot/ps)*[scenario.timing.dt*60*60]/(1-por);


        ifind = find([xprof >= xtoe] & [xprof < [xtoe + 15]]);
        xvals = xprof(ifind);
        try
            dx = abs(xvals(2)-xvals(1));
        catch err
            dx = abs(xprof(2)-xprof(1));
        end

    tot_x = abs(xvals(end)-xvals(1));
    %assuming a triangle that max deposition is 1/3 the distance of the
    %total deposition length
    if Qtot_m3m_dt > 0
        L = scenario.models.VegDepositionLengthScale;
        L = abs(L*cosd(abs(windDir)));

        if scenario.models.VegDepositionStyle == 1
    
            h = Qtot_m3m_dt/(L/2);
            theta1 = atan(h/(L/3));
            theta2 = atan(h/(L*2/3));
            xlocal = xprof-xtoe;
            xlocal = 0:0.01:L;
            ifind1 = find(xlocal>=0 & xlocal<L/3);
            xlocal1 = xlocal(ifind1);
            dz_local1 = xlocal1*tan(theta1);
            ifind2 = find(xlocal>=L/3 & xlocal<=L);
            xlocal2 = L-xlocal(ifind2);
            dz_local2 = xlocal2*tan(theta2);    
            %xlocalt = cosd(abs(windDir))*xlocal;
            xlocalt = xlocal + xtoe;
            dz_localcomb = [dz_local1 dz_local2];
            dz = zeros(size(xprof));
        elseif scenario.models.VegDepositionStyle == 2
            h = Qtot_m3m_dt/(L);
            theta1 = atan(h/(L));
            xlocal = 0:0.01:L;
            dz_local = xlocal*tan(theta1);
            xlocalt = xlocal + xtoe;
            dz_localcomb = fliplr(dz_local);
        else
              h = Qtot_m3m_dt/(L/2);
            theta1 = atan(h/(L/3));
            theta2 = atan(h/(L*2/3));
            xlocal = xprof-xtoe;
            xlocal = 0:0.01:L;
            ifind1 = find(xlocal>=0 & xlocal<L/3);
            xlocal1 = xlocal(ifind1);
            dz_local1 = xlocal1*tan(theta1);
            ifind2 = find(xlocal>=L/3 & xlocal<=L);
            xlocal2 = L-xlocal(ifind2);
            dz_local2 = xlocal2*tan(theta2);    
            %xlocalt = cosd(abs(windDir))*xlocal;
            xlocalt = xlocal + xtoe + scenario.models.VegDepositionCenter;
            dz_localcomb = [dz_local1 dz_local2];
            dz_localcomb = smoothdata(dz_localcomb, 'gaussian', 50);
            dz = zeros(size(xprof));          

        end

            %now need to scale to ensure that too much mass isnt added
            scaling = Qtot_m3m_dt / sum(dz_localcomb*0.01, 'omitnan');
            dz_localcomb = dz_localcomb*scaling;

            if numel(xlocalt) > 1
            zoff = interp1(xlocalt, dz_localcomb, xprof);
            elseif numel(xlocalt) == 1
    
                [~, imin] = min(abs(xprof - xlocalt), [], 'omitnan');
                zoff = zeros(size(xprof));
                zoff(imin) = dz_localcomb;
    
            else
                zoff = zeros(size(xprof));
            end      
            zoff(find(isnan(zoff) == 1)) = 0;
        
            dVacc2 = sum(zoff*dx, 'omitnan');
    
        else
            dVacc2 = 0;
            zoff = zeros(size(zprof));
            end

    zprof_update = zprof+zoff;


end



function out = run_erosion_accretion_model(xM,z,time,WL,Ho,Lo,T,Bo,dtoe,output_times, D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility, zShore, shorechange, scenario)

    %Set Model Coefficients
    d50 = D50/1000; %mm to m
    nsigma = 2; %in definition of R2, R16.. For R2, nsigma=2, R16 = nsigma=1; 
    g = 9.81; %gravity
    Kd = WaveRunupFactor; % coefficient to account for higher runup on dune
    Cs = DuneErodibility;
    Ac = 1.34*10^(-3);
    bc = 3.19*10^(-4);
    Btfac = DuneSlopeTrajectory;
    ifind = find(z>dtoe);
    val = z(ifind(1));
    st1 = ifind(1);
    deltax = abs(xM(2)-xM(1));
    
    %Model Initialization
    zbT = NaN(size(xM));
    zbT(st1:end)=dtoe;  %trajectory that dune toe receeds.
    dt = diff(time(1:2));%dt in seconds
    try
        xShore(1) = interp1(z, xM, zShore);
    catch err
        xShore(1) = linterp(z, xM, zShore);
    end

    %Main Program Loop
    output_num = 1;
    for tt=1:length(WL)
        current_output_time = output_times(output_num);
        if tt==1
            st = st1;
        else
            ifind = find(z>=dtoe);
            st = ifind(1);
        end

        zbT = NaN(size(xM));
        zbT(st:end)=dtoe;  %trajectory that dune toe receeds.    zbT = NaN(size(xM));
        zbT(1:st) = z(1:st);

        %dune toe position
        xToe(tt) = xM(st); 
 
        xShoreChange = time(tt) * -shorechange /(365*24*60*60);
        xShore(tt) = xShore(1) + xShoreChange;

        %dune volume
        V(tt) = sum(deltax.*(z(st:end)));    %measured in ref to z=0
        Vc = cumsum(deltax.*(z(st:end)-zbT(st:end)));  %cumulative volume above the dune trajectory

        Beta(tt) = tan(abs((zShore - dtoe)/(xShore(tt) - xToe(tt))));
        if abs(Beta(tt)) > 0.25
            Beta(tt) = 0.25;
        end
        Bt = Beta(tt)*Btfac;

        %stockdon for TWL
        etabar(tt) = 0.35.*Beta(tt).*sqrt(Ho(tt).*Lo(tt));
        sigma_s(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.004))./2.*nsigma./2;
        zR(tt) = 1.1.*(etabar(tt)+ sigma_s(tt));
        zTotal(tt) = zR(tt).*Kd + WL(tt);
        if zTotal(tt)>= [max(z)]
            zTotal(tt) = max(z);
        end 
        p(tt) = 0.5 * (1 - erf( (dtoe - (etabar(tt) + WL(tt))) / (sigma_s(tt) * sqrt(2)) ));
        Nc(tt) = p(tt).*(dt./T(tt));
        dVT(tt) = 4.*Cs.*(max(zTotal(tt)-dtoe,0)).^2.*Nc(tt);

        if dVT(tt)<=0
            dVT(tt) = 0; 
            ii=0;
        else
            Vc2 = Vc;
           dx = abs(xM(2)-xM(1));

            ifind = find(z>=dtoe);
            zdiff = z-dtoe;
            zdiff(zdiff<0) = 0;
            volcum = cumsum(zdiff)*dx;

            ifind = find(volcum <= dVT(tt) & z>= dtoe);
            z(ifind) = dtoe;

            if numel(ifind)>=1
                diffLast = dVT(tt)-volcum(ifind(end));
                dzNext = diffLast/dx;
                    if dzNext > 0 
                        try
                        z(ifind(end)+1) = z(ifind(end)+1)-dzNext;
                        catch err
                        end
                    end
            end
       
        end

        prof = z;

        %clean-up variables
        clear Vc
        scenario.erosion.TWL(tt) = zTotal(tt);

        %update with shoreline change
        ifind = find(z<=dtoe);
        xToe = xM(ifind(end));
        ifind = find(xM>=xShore(tt) & xM <= xToe);
        zbeach = linspace(zShore,z(ifind(end)), numel(ifind));
        prof(ifind) = zbeach;
        prof(1:ifind(1))= zShore;
        
        %Run Accretion Model
        [scenario, zprof_update, dVaccretion(tt), dVaccretion2(tt)] = drt_accretion_instantaneous(scenario, prof, dtoe, tt);
        z = zprof_update;
    
        if scenario.models.Avalanche == 1
            [z] = avalanche(xM,z,scenario.models.AvalancheAngle);
        end


        if time(tt)>= current_output_time
            zNew(output_num,:) = z; % assumes vertical cliff face; probably needs a slope adjustmnet
            %zNew(output_num,:) = [z(1:st1) zbT(st1+1:st) z(st+1:end)]; % assumes vertical cliff face; probably needs a slope adjustmnet
            actual_output_times(output_num) = time(tt);
            output_num = output_num + 1;
        end 

    end

    %xToe(tt+1) = xToe(tt)+dx(tt);
   % st = st+ii-1;
   zbOut = dtoe;
   % xToe = xToe(1:end-1);

    %save a structure of model results
    out.TWL = zTotal; %final profile assuming a vertical cliff face (can under-estimate volumes of erosion)
    out.zNew = zNew(1:end,:); %profiles at each time-step assuming a vertical cliff face (can under-estimate volumes of erosion)
    out.ztoe = zbOut; %elevation of dune toe through time at each time step
    out.xtoe = xToe; %cross-shore position of dune toe at each time step
    out.Nc = Nc; %number of dune collisions at each time step
    out.dV_erosion = dVT; %volume of eroded sand at each time step
    out.dV_accretion = dVaccretion;
    out.dV_accretion2 = dVaccretion2;
    out.zRunup = zR; 
    out.times = time;
    out.zmat_time = actual_output_times;
    out.Beta = Beta;

    end