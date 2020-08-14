function drt_plotting_ensemble(scenario, scenarioLow, scenarioHigh)
    %save temptemp.mat scenario
    %plot erosion data first
    xlims = [min(scenario.timing.times) max(scenario.timing.times)];
    
    figure('units','normalized','outerposition',[0 0 1 1])
    
    %Subplot for Wind Time Series
    ax1 = subplot(2,3,1);
    yyaxis left
    plot(scenario.timing.times, scenario.env.winds.windSpeed, 'k-', 'LineWidth', 3)
    xpts = linspace(min(xlims), max(xlims), 5);
    xpts = unique(floor(xpts));
    xlim(xlims)
    set(gca, 'XTick', xpts)
    datetick('x', 'keeplimits', 'keepticks')
    ylabel('Wind Speed (m/s)')
    grid on
    set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold')
    title('Winds')    
    yyaxis right
    plot(scenario.timing.times, scenario.env.winds.windDirection,'.','LineWidth', 0.5, 'Color', [0.5 0.5 0.5], 'MarkerSize', 12)
    ylabel('Wind Direction (^o)')
    ylim([0 360])
    
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = [0.5 0.5 0.5];
    
    %Subplot for Wave Time Series   
    ax2 = subplot(2,3,2);
    plot(scenario.timing.times, scenario.env.waves.Hs_deepwater,'k-','LineWidth', 3)
    xlim(xlims)
    set(gca, 'XTick', xpts)
    datetick('x', 'keeplimits', 'keepticks')
    ylabel('H_s (m)')
    grid on
    set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold')
    title('Waves')
    yyaxis right
    plot(scenario.timing.times, scenario.env.waves.D_deepwater,'.','LineWidth', 0.5, 'Color', [0.5 0.5 0.5], 'MarkerSize', 12)
    ylabel('Wave Direction (o)')
    ylim([0 360])
    
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = [0.5 0.5 0.5];
       
    %Subplot for Water Level Time Series
    subplot(2,3,3)
    h1 = plot(scenario.timing.times, scenario.env.tides.wl, 'Color', [0.5 0.5 0.8], 'LineWidth', 3)
    xlim(xlims)
    datetick('x', 'keeplimits')
    hold on
    h2 = plot(scenario.timing.times, scenario.erosion.TWL, 'k','LineWidth', 3)   
    plot(scenario.timing.times, scenarioLow.erosion.TWL, 'k-','LineWidth', 1)   
    plot(scenario.timing.times, scenarioHigh.erosion.TWL, 'k','LineWidth', 1)   
    
       
    max_twl = max([scenario.erosion.TWL(:) scenarioLow.erosion.TWL(:) scenarioHigh.erosion.TWL(:)]');
    min_twl = min([scenario.erosion.TWL(:) scenarioLow.erosion.TWL(:) scenarioHigh.erosion.TWL(:)]');
    twl_all = [max_twl(:); flipud(min_twl(:))];
    
    timing_all = [scenario.timing.times(:); flipud(scenario.timing.times(:))];
    ptwl = patch(timing_all, twl_all, 'k');
    ptwl.FaceAlpha = 0.1;
    
    
    xlim(xlims)    
    plot([min(scenario.timing.times) max(scenario.timing.times)], [scenario.grids.morphometrics.dtoe scenario.grids.morphometrics.dtoe], '--', 'Color', [0.9 0.3 0.3], 'LineWidth', 2)
    plot([min(scenario.timing.times) max(scenario.timing.times)], [scenario.grids.morphometrics.dhigh scenario.grids.morphometrics.dhigh], '--', 'Color', [0.3 0.1 0.1], 'LineWidth', 2)
    
   ylim([min([min(scenario.env.tides.wl) scenario.grids.morphometrics.dtoe])-1 max([max(scenario.erosion.TWL) scenario.grids.morphometrics.dhigh])+1])
   
   xo = xlims(1) + [xlims(2)-xlims(1)]*0.05;
   text(xo, [scenario.grids.morphometrics.dhigh+0.25], 'Dune Crest', 'FontWeight', 'bold', 'Color', [0.3 0.1 0.1])
   text(xo, [scenario.grids.morphometrics.dtoe+0.25], 'Dune Toe', 'FontWeight', 'bold', 'Color', [0.9 0.3 0.3])
    
    %ax1.YAxis(1).Color = 'k';
    %ax1.YAxis(2).Color = [0.5 0.5 0.5];
    
    
    %xpts = linspace(min(xlims), max(xlims), 4);
    set(gca, 'XTick', xpts)
    xlim(xlims)
    datetick('x', 'keeplimits', 'keepticks')
    ylabel('SWL & TWL(m)')    
    grid on
    set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold')
    title('Water Levels')
    
    legend([h1(1) h2(1) ptwl(1)], 'SWL', 'TWL', 'Ensemble TWL', 'Location', 'NorthEast')
    
    %Subplot for Morphology Change
    subplot(2,3,6)
    cmap = parula(numel(scenario.erosion.zmat_times)+1);
    cumDV_erosion = -cumsum(scenario.erosion.dV);
    cumDV_accretion = cumsum(scenario.accretion.dV);
    cumDV_net = cumDV_accretion+cumDV_erosion;
    
    cumDV_erosion_high = -cumsum(scenarioHigh.erosion.dV);
    cumDV_accretion_high = cumsum(scenarioHigh.accretion.dV);
    cumDV_net_high = cumDV_accretion_high+cumDV_erosion_high;    
    
    cumDV_erosion_low = -cumsum(scenarioLow.erosion.dV);
    cumDV_accretion_low = cumsum(scenarioLow.accretion.dV);
    cumDV_net_low = cumDV_accretion_low+cumDV_erosion_low;
    
    max_erosion = max([cumDV_erosion_low(:) cumDV_erosion(:) cumDV_erosion_high(:)]');
    min_erosion = min([cumDV_erosion_low(:) cumDV_erosion(:) cumDV_erosion_high(:)]');    
    erosion_all = [max_erosion(:); flipud(min_erosion(:))];

    max_accretion = max([cumDV_accretion_low(:) cumDV_accretion(:) cumDV_accretion_high(:)]');
    min_accretion = min([cumDV_accretion_low(:) cumDV_accretion(:) cumDV_accretion_high(:)]');    
    accretion_all = [max_accretion(:); flipud(min_accretion(:))];
    
    max_net = max([cumDV_net_low(:) cumDV_net(:) cumDV_net_high(:)]');
    min_net = min([cumDV_net_low(:) cumDV_net(:) cumDV_net_high(:)]');    
    net_all = [max_net(:); flipud(min_net(:))];
    
    timing_all = [scenario.timing.times(:); flipud(scenario.timing.times(:))];
    
    ylims = [min(net_all)-0.5 max(net_all)+0.5]; 
    
    hold on  
%     ha = plot(scenario.timing.times, cumDV_erosion, 'b-', 'LineWidth', 3);   
%     he = plot(scenario.timing.times, cumDV_accretion, 'r-', 'LineWidth', 3);   
    hn = plot(scenario.timing.times, cumDV_net, 'k-', 'LineWidth', 5);   
    hnlow = plot(scenario.timing.times, max_net, 'k-', 'LineWidth', 2);   
    hnhigh = plot(scenario.timing.times, min_net, 'k-', 'LineWidth', 2);   
    
    pn = patch(timing_all, net_all, 'k');
    pn.FaceAlpha = 0.1;

    for itime = 1:numel(scenario.erosion.zmat_times)        
         plot([scenario.erosion.zmat_times(itime) scenario.erosion.zmat_times(itime)], ylims, '--', 'LineWidth', 2, 'Color', cmap(itime, :))   
    end
    xlim(xlims)
    ylim(ylims)
    set(gca, 'XTick', xpts)
    datetick('x', 'keeplimits', 'keepticks')    
    ylabel('\Delta V_{dune} (m^3/m)')    
    grid on
    set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold')       
    title('Dune Volume Change')
    legend([pn(1)], 'Ensemble Envelope', 'Location', 'SouthWest')
    
    subplot(2,3,[4 5])
    hold on
    for itime = 1:numel(scenario.erosion.zmat_times)
        plot(scenario.grids.XGrid, scenario.erosion.Z(:,itime), '-', 'Color', cmap(itime,:), 'LineWidth', 4)
    end
    ylims = [nanmin(nanmin(scenario.erosion.Z)) nanmax(nanmax(scenario.erosion.Z))];
    ylims = [0 nanmax(nanmax(scenario.erosion.Z))+1];
    ylim(ylims)
    ylabel('Z (m, NAVD)')
    xlabel('Cross-Shore Distance (m)')
    grid on
    set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold')    
    title('Dune Profile Change (Base Case, Erosion Only)')    
    
    if cumDV_net(end)>-0.5
        col = [0.2 0.9 0.2];
        tex = 'No or Minimal Net Dune Erosion Predicted';
    elseif cumDV_net(end)>-5
        col = [1 0 1];
        tex = 'Minor Dune Erosion Predicted';
    else
        col = [1 0 0];
        tex = 'Substantial Dune erosion predicted';
    end
        
    hold on
    p = patch([min(scenario.grids.XGrid) max(scenario.grids.XGrid) max(scenario.grids.XGrid) min(scenario.grids.XGrid) min(scenario.grids.XGrid)], ...
        [ylims(1) ylims(1) ylims(2) ylims(2) ylims(1)], col);
    p.FaceAlpha = 0.075;
    if min(scenario.grids.ZGrid)>0
        xlim([min(scenario.grids.XGrid) max(scenario.grids.XGrid)])
        x0 = min(scenario.grids.XGrid);
    else
        x0 = linterp(scenario.grids.ZGrid, scenario.grids.XGrid, 0);
        xlim([x0 max(scenario.grids.XGrid)])
    end
    
    %Add additional info onto plot
    xlims = xlim;
    xo = xlims(1) + [xlims(2)-xlims(1)]*0.8;
    xo2 = xlims(1) + [xlims(2)-xlims(1)]*0.5;
    yo2 = ylims(1) + [ylims(2)-ylims(1)]*0.08;

    plot(xlims, [scenario.grids.morphometrics.dtoe scenario.grids.morphometrics.dtoe], '--', 'Color', [0.9 0.3 0.3], 'LineWidth', 2)
    plot(xlims, [scenario.grids.morphometrics.dhigh scenario.grids.morphometrics.dhigh], '--', 'Color', [0.3 0.1 0.1], 'LineWidth', 2)

    text(xo, [scenario.grids.morphometrics.dhigh+0.25], 'Dune Crest', 'FontWeight', 'bold', 'Color', [0.3 0.1 0.1])
    text(xo, [scenario.grids.morphometrics.dtoe+0.25], 'Dune Toe', 'FontWeight', 'bold', 'Color', [0.9 0.3 0.3])
    text(xo2, yo2, tex, 'FontWeight', 'bold', 'Color', col)
    legend(datestr(scenario.erosion.zmat_times(:)), 'Location', 'NorthWest')
   
end