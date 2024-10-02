
%load example intput file
load example_inputs.mat

%call model scropt
model_output = tundra_erosion_model(scenario);

%create example output plot
figure,
hold on
cmap = jet(500);
for itime = 1:500
    if itime == 1
        h1 = plot(model_output.grids.XGrid, model_output.output.Z(:,itime), 'Color', cmap(itime, :), 'LineWidth', 2);
    elseif itime == 500
        h2 = plot(model_output.grids.XGrid, model_output.output.Z(:,itime), 'Color', cmap(itime, :), 'LineWidth', 2);
    else
        plot(model_output.grids.XGrid, model_output.output.Z(:,itime), 'Color', cmap(itime, :), 'LineWidth', 2);
    end
end
xlabel('Cross-Shore Distance (m)')
ylabel('Elevation (m)')
legend([h1(1) h2(1)], 'Initial Morphology', 'Final Prediction')
xlim([150 260])
