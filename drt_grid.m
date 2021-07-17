function scenario = drt_grid(scenario)
    %drt_erosion: code to set up grids for input to the dune response tool
    
    %Build Profile 
        %assumed geometric and gridding parameters for simplicity
        scenario.grids.morphometrics.dunewidth = 10; 
        scenario.grids.morphometrics.dx = 0.1;
        scenario.grids.morphometrics.smoothlev = 10;
        scenario.grids.morphometrics.shore = 0;

        %set up profile segments
        numdune = round(scenario.grids.morphometrics.dunewidth/scenario.grids.morphometrics.dx);
        Zbackdune = (scenario.grids.morphometrics.dhigh-2)*ones(ceil(20/scenario.grids.morphometrics.dx),1);
        Zbackduneslope = (scenario.grids.morphometrics.dhigh-0.5):(scenario.grids.morphometrics.duneslope)*scenario.grids.morphometrics.dx:scenario.grids.morphometrics.dhigh;
        Zdune1 = scenario.grids.morphometrics.dhigh*ones(numdune,1);
        Zdune2 = scenario.grids.morphometrics.dhigh:-(scenario.grids.morphometrics.duneslope)*scenario.grids.morphometrics.dx:scenario.grids.morphometrics.dtoe;
        Zbeach = scenario.grids.morphometrics.dtoe:-(scenario.grids.morphometrics.backshoreslope)*scenario.grids.morphometrics.dx:scenario.grids.morphometrics.shore;                      
        
        %combine all segments
        Zall = [ Zbackduneslope(:); Zdune1(:); Zdune2(:); Zbeach(:)];
        Zall = flipud(Zall);

    %Smooth Data
    Zall = smoothdata(Zall, 'gaussian', scenario.grids.morphometrics.smoothlev);
    
    %Grid Setup
    XGrid = -[0:1:[numel(Zall)-1]]*abs(scenario.grids.morphometrics.dx);
    XGrid = fliplr(XGrid)';
    
    %Output Data
    scenario.grids.XGrid = XGrid(:)';
    scenario.grids.ZGrid = Zall(:)';
    
end