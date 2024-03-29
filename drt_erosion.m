function scenario = drt_erosion(scenario)
    %drt_erosion: code to run the Palmsten and Holman (2012) dune retreat
    %model
      
    %initialize model in PH12 conventions
        xM = scenario.grids.XGrid+max(abs(scenario.grids.XGrid));
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

    %run model
        output_times = linspace(time(1), time(end), 10);    
    
    %output_times 
        PH12 = runPH12(xM(:)',z(:)',time(:),WL(:),Ho(:),Lo(:),T(:),Bo,dtoe, output_times(:), D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility);

    %store relevant model output and send back to main program
        scenario.erosion.dV = PH12.dV;
        scenario.erosion.Z = PH12.zNew';
        scenario.erosion.TWL = PH12.TWL;
        scenario.erosion.xToe = PH12.xtoe;
        scenario.erosion.ztoe = PH12.ztoe;
        scenario.erosion.Nc = PH12.Nc;
        scenario.erosion.dVR = PH12.dVResidual;
        scenario.erosion.dVT= PH12.dVT;
        scenario.erosion.times = PH12.times;    
        scenario.erosion.zmat_times = [PH12.zmat_time./86400]+scenario.timing.times(1);
        scenario.erosion.zRunup = PH12.zRunup;
end

function PH12 = runPH12(xM,z,time,WL,Ho,Lo,T,Bo,zb,output_times, D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility)
    %Function to run the Palmsten and Holman (2012) model for dune retreat
    %   inputs:
    %       xM - x vector (with origin offshore)
    %       z - elevation profile (same size as xM)
    %       time - time vector wave data in seconds
    %       WL - water level time-series (tide + surge)
    %       Ho - deep water wave height time-series 
    %       Lo - deep water wave length time-series
    %       T - peak wave period time-series
    %       Bo - initial beach slope
    %       zb - initial dune toe elevation 

    %Set Model Coefficients
    d50 = D50/1000; %mm to m
    nsigma = 2; %in definition of R2, R16.. For R2, nsigma=2, R16 = nsigma=1; 
    g = 9.81; %gravity
    Kd = WaveRunupFactor; % coefficient to account for higher runup on dune
    Cs = DuneErodibility;
    Ac = 1.34*10^(-3);
    bc = 3.19*10^(-4);
    Btfac = DuneSlopeTrajectory;

    %Model Initialization
    [val st1] = (min((z-zb(1)).^2)); %find grid point where initial dune toe is
    Bt = Bo*Btfac; %slope at which beta receeds. LEH04 = 1, PH11 = 0.54....
    zbT = [repmat(NaN,1,st1-1) Bt.*(xM(st1:end)-xM(st1)) + zb(1)];  %trajectory that dune toe receeds.
    ireplace = find(zbT>z);
    zbT(ireplace) = z(ireplace);
    dt = diff(time(1:2));%dt in seconds

    %Main Program Loop
    output_num = 1;

    for tt=1:length(WL)
        current_output_time = output_times(output_num);

        if tt==1
            st = st1;
        else
            st = st+ii-1;
        end

        %dune toe position
        xToe(tt) = xM(st);

        %dune volume
        V(tt) = sum(abs(diff(xM(1:2))).*(z(st:end)));    %measured in ref to z=0
        Vc = cumsum(abs(diff(xM(1:2))).*(z(st:end)-zbT(st:end)));  %cumulative volume above the dune trajectory
        Vc = Vc - Vc(1);

        Beta(tt) = Bo;   %initial dummy guess.

        %stockdon for TWL
        etabar(tt) = 0.35.*Beta(tt).*sqrt(Ho(tt).*Lo(tt)); %mean swash (setup)
        sigma_s(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.0004))./2.*nsigma./2;
        zR(tt) = 1.1.*(etabar(tt)+ sigma_s(tt));
        sigma_s2(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.0004))./2;
        zR(tt) = 1.1.*(etabar(tt)+ sigma_s(tt));

        zRLEH(tt) = 0.158.*sqrt(Ho(tt)./1.416.*Lo(tt));
        zTotal(tt) = zR(tt).*Kd + WL(tt);
        if zTotal(tt)>= [max(z)]
            zTotal(tt) = max(z);
        end 
        p(tt) =1-cdf('norm',zb(tt),etabar(tt)+WL(tt),sigma_s(tt));
        Nc(tt) = p(tt).*(dt./T(tt));
        if tt>1
            dV(tt) = 4.*Cs.*(max(zTotal(tt)-zb(tt),0)).^2.*Nc(tt);
            dVT(tt) = dV(tt) - dVResidual(tt-1);
        else
            dVT(tt) = 4.*Cs.*(max(zTotal(tt)-zb(tt),0)).^2.*Nc(tt);
        end
        if dVT(tt)<0
            ii=1;
        else
            [val ii] = (min((Vc-dVT(tt)).^2)); %find grid point where dune toe is
        end
        dx(tt) = xM(st+ii-1)-xToe(tt);
        dVResidual(tt) = Vc(ii)-dVT(tt);
        zb(tt+1) = Bt.*dx(tt) + zb(tt);  %trajectory that dune toe receeds.

        if time(tt)>= current_output_time
            zNew(output_num,:) = [z(1:st1) zbT(st1+1:st) z(st+1:end)]; % assumes vertical cliff face; probably needs a slope adjustmnet
            actual_output_times(output_num) = time(tt);
            output_num = output_num + 1;
        end

        %clean-up variables
        clear Vc
    end

    xToe(tt+1) = xToe(tt)+dx(tt);
    st = st+ii-1;
    zbOut = zb(1:end-1);
    xToe = xToe(1:end-1);

    %save a structure of model results
    PH12.TWL = zTotal; %final profile assuming a vertical cliff face (can under-estimate volumes of erosion)
    PH12.zNew = zNew(1:end,:); %profiles at each time-step assuming a vertical cliff face (can under-estimate volumes of erosion)
    PH12.ztoe = zbOut; %elevation of dune toe through time at each time step
    PH12.xtoe = xToe; %cross-shore position of dune toe at each time step
    PH12.Nc = Nc; %number of dune collisions at each time step
    PH12.dV = dV; %volume of eroded sand at each time step
    PH12.dVResidual = dVResidual;
    PH12.dVT = dVT; 
    PH12.zRunup = zR; 
    PH12.times = time;
    PH12.zmat_time = actual_output_times;

end