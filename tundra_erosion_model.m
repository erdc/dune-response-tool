function scenario = tundra_erosion_model(scenario)
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
        out = run_erosion_model(xM(:)',z(:)',time(:),WL(:),Ho(:),Lo(:),T(:),Bo,dtoe, output_times(:), D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility, zShore, shorechange, scenario);

    %store relevant model output and send back to main program
        scenario.output.dV_erosion = out.dV_erosion;
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


function out = run_erosion_model(xM,z,time,WL,Ho,Lo,T,Bo,dtoe,output_times, D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility, zShore, shorechange, scenario)

    %Set Model Coefficients
    d50 = D50/1000; %mm to m
    nsigma = 2; %in definition of R2, R16.. For R2, nsigma=2, R16 = nsigma=1; 
    g = 9.81; %gravity
    Kd = WaveRunupFactor; % coefficient to account for higher runup on dune
    Csmax = DuneErodibility;
    %erodibility = scenario.models.erodibility;
    thaw_depths = scenario.models.thaw_depth;
    Ac = 1.34*10^(-3);
    bc = 3.19*10^(-4);
    Btfac = DuneSlopeTrajectory;
    ifind = find(z>dtoe);
    val = z(ifind(1));
    st1 = ifind(1);
    deltax = abs(xM(2)-xM(1));
    maxrate = 10; %currently cant exceed 10 m3/m/hr of erosion
    dx = abs(xM(2)-xM(1));

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
    zbThawInit = z;

    for tt=1:length(WL)

        current_output_time = output_times(output_num);

        if tt==1
            st = st1;
        else
            itoe = find(z>=dtoe);
            st = itoe(1);
        end

        zbT = NaN(size(xM));
        zbT(st:end)=dtoe;  %trajectory that dune toe receeds.    zbT = NaN(size(xM));
        zbT(1:st) = z(1:st);

        if thaw_depths(tt) == 0 
            zbThawInit = zbT;
        end

        %dune toe position
        xToe(tt) = xM(st); 
 
        xShoreChange = time(tt) * -shorechange /(365*24*60*60);
        xShore(tt) = xShore(1) + xShoreChange;

        %dune volume
        V(tt) = sum(deltax.*(z(st:end)));    %measured in ref to z=0
        Vc = cumsum(deltax.*(z(st:end)-zbT(st:end)));  %cumulative volume above the dune trajectory
        %Vc(1)
        %Vc = Vc - Vc(1);

        Beta(tt) = tan(abs((zShore - dtoe)/(xShore(tt) - xToe(tt))));
        if abs(Beta(tt)) > 0.25
            Beta(tt) = 0.25;
        end
        Bt = Beta(tt)*Btfac;

        if Ho(tt) > 0
            %stockdon for TWL
            etabar(tt) = 0.35.*Beta(tt).*sqrt(Ho(tt).*Lo(tt)); %mean swash (setup)
            sigma_s(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.004))./2.*nsigma./2;
            zR(tt) = 1.1.*(etabar(tt)+ sigma_s(tt));
            sigma_s2(tt) = sqrt(Ho(tt).*Lo(tt).*(0.563.*(Beta(tt).^2)+0.004))./2;
            zR(tt) = 1.1.*(etabar(tt)+ sigma_s(tt));
    
            zRLEH(tt) = 0.158.*sqrt(Ho(tt)./1.416.*Lo(tt));
            zTotal(tt) = zR(tt).*Kd + WL(tt);
            if zTotal(tt)>= [max(z)]
                zTotal(tt) = max(z);
            end 
            p(tt) =1-cdf('norm',dtoe,etabar(tt)+WL(tt),sigma_s(tt));
            Nc(tt) = p(tt).*(dt./T(tt));

            %calculate erodibility as a function of thaw depth and local
            %topography

            %relative topography
            zrel = z - dtoe;
            zrel(zrel < 0) = 0;

            zbThawNow = zbThawInit - thaw_depths(tt);

            zbThawNowRel = zbThawInit - dtoe;
            ithawstart = find(zbThawNowRel > 0);

            itoe = find(zrel>0);
            
            thawDepthLocal = zbT-zbThawNow;
            thawDepthLocal(thawDepthLocal < 0) = 0;


            vert_erodible = [zrel-thawDepthLocal]./zrel; %looks at relative percentage above/below thaw depth
            vert_erodible = 1-vert_erodible;
            vert_erodible(vert_erodible>1) = 1;


            nxthaw = floor([thaw_depths(tt)*1.5]/dx);


            irel = find(zrel>=0);
            if nxthaw > 0
                try
                vert_erodible(ithawstart(1):[ithawstart(1)+nxthaw]) = 1;
                catch err
                end
            end

            Cs = Csmax;
            dVT(tt) = 4.*Cs.*(max(zTotal(tt)-dtoe,0)).^2.*Nc(tt);

            if dVT(tt)>maxrate
                dVT(tt) = maxrate;
            end





        if dVT(tt)<=0
            dVT(tt) = 0; 
            ii=0;
        elseif thaw_depths(tt) == 0
            dVT(tt) = 0; 
            ii=0;
        else
            Vc2 = Vc;

            %ifind = find(z>=dtoe);
            zdiff = z-dtoe;
            zdiff(zdiff<0) = 0;

            %
            relratios = vert_erodible; %vert erodible should be the relative height that is erodible
            %relratios(relratios == 0) = 0.005;
            %zapparent = zdiff./relratios;
            zapparent = zdiff.*relratios;


            volcum = cumsum(zapparent)*dx;
            volcumA = cumsum(zdiff)*dx;

            ifind2 = find(volcumA > dVT(tt));

            if numel(ifind2)>0

                    %ifind2(1)
                    %volcumA

                    volcumA(ifind2(1)) = 0;
                    if max(volcumA) < dVT(tt)
                        volcumA(ifind2(1)) = dVT(tt);
                        voldiff = volcumA(ifind2(1))-volcumA(ifind2(1)-1);
                        rat = [voldiff*dx]./zdiff(ifind2(1));
                        if rat > 1
                            rat = 1;
                        end
                        zapparent(ifind2(1)) = zapparent(ifind2(1))*rat;
                    end
        
                    %fix dz grid
                    ifind2 = find(volcumA == 0);
                    zapparent(ifind2) = 0;  
                    
                    %find cells where there is erodible sediment
                    ifind2 = find(zapparent > 0);
        
                    if numel(ifind2)> 0
        
                        proceed = 1;
                        count = 0;
                        dVTactual = 0;
                        while proceed == 1
                            if numel(ifind2)>count
                            count = count+1;
                                
                                if zapparent(ifind2(count)) > zdiff(ifind2(count))
                                    zapparent(ifind2(count)) = zdiff(ifind2(count));
                                end


                                if zapparent(ifind2(count)) == 0 %no volume change if there is nothing to erode at the seaward cell
                                    dVTactual = 0;
                                    proceed = 0;
                                elseif zapparent(ifind2(count)) < zdiff(ifind2(count)) %remove volume if not the full cell is removed
                                    z(ifind2(count)) = z(ifind2(count)) - zapparent(ifind2(count));
                                    dVTactual = dVTactual + [zapparent(ifind2(count))*dx];
                                    proceed = 0;
                                else %otherwise the full cell should be removed and should repeat loop
                                    z(ifind2(count)) = dtoe;
                                    proceed = 1;
                                    dVTactual = dVTactual + zdiff(ifind2(count))*dx;
                                end
        
                            else
                                proceed = 0;
                            end

                            if count > 20 %dont go more than 20 cells landward in a single time step
                                proceed = 0;
                            end
                        end
        
        
                       
                    dVT(tt) = dVTactual;

                    end

            else

                    dVT(tt) = 0;

            end
       
        end

                else
            dVT(tt) = 0;
        end

        prof = z;

        %clean-up variables
        clear Vc
        try
            scenario.erosion.TWL(tt) = zTotal(tt);
        catch err
            scenario.erosion.TWL(tt) = NaN;
            zR(tt) = NaN;
            zTotal(tt) = NaN;
            Nc(tt) = NaN;
        end
            
        %update with shoreline change
        ifind = find(z<=dtoe);
        xToe = xM(ifind(end));
        xToe2(tt) = xToe;
        ifind = find(xM>=xShore(tt) & xM <= xToe);

       % display(['ifind = ', num2str(ifind)])

        if numel(ifind)>2
            zbeach = linspace(zShore,z(ifind(end)), numel(ifind));
            prof(ifind) = zbeach;
            prof(1:ifind(1))= zShore;
            z = prof;
        else
            ifind = find(z>dtoe);
            ifind2 = find(xM>=xShore & xM <=xM(ifind(1)));
            zbeach = linspace(zShore,z(ifind(1)), numel(ifind2));
            prof(ifind2) = zbeach;
            prof(1:ifind2(1))= zShore;
            z = prof;        
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
    out.xtoe = xToe2; %cross-shore position of dune toe at each time step
    out.Nc = Nc; %number of dune collisions at each time step
    out.dV_erosion = dVT; %volume of eroded sand at each time step
    out.zRunup = zR; 
    out.times = time;
    out.zmat_time = actual_output_times;
    out.Beta = Beta;

    end