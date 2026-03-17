%load structure
load example_input_and_output.mat

%modify shoreline change rate
scenario.grids.morphometrics.shorechange = -1; %update to -1 m/yr SCR

%run model
scenario = drt_model_instantaneous(scenario);

% for MATLAB R2025b compatibility, use this instead:
% scenario = drt_model_instantaneous_R2025b(scenario);

%plot model results
figure, 
hold on
cmap = jet(500);
for idx = 1:500
    plot(scenario.grids.XGrid, scenario.output.Z(:, idx), 'Color', cmap(idx,:))
end
xlabel('Cross Shore Distance (m)')
ylabel('Elevation (m)')

