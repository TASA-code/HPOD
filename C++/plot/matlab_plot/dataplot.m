clear; clc

data = importdata("../../build/bin/output.txt");

% Split each line into components
components = cellfun(@(x) strsplit(x), data, 'UniformOutput', false);

% Find the maximum number of components in any line
numComponents = cellfun(@length, components);
maxComponents = max(numComponents);

% Ensure all lines have the same number of components by padding with NaNs
paddedComponents = cellfun(@(x) [x, repmat({NaN}, 1, maxComponents - length(x))], components, 'UniformOutput', false);

% Convert the padded cell array of components into a matrix of doubles
componentsMatrix = cellfun(@(x) str2double(x), paddedComponents, 'UniformOutput', false);

% Concatenate the matrices into one large matrix
componentsMatrix = vertcat(componentsMatrix{:});

% Extract the required components using array slicing
ECI = componentsMatrix(:, [3 4 5]);
ECEF = componentsMatrix(:, [9 10 11]);
lon = componentsMatrix(:, 17);
lat = componentsMatrix(:, 18);


figure()
subplot(2,2,1)
plot3(ECI(:,1), ECI(:,2), ECI(:,3),'LineWidth',2,'color','k')
hold on
grid on
plot_globe()

subplot(2,2,2)
plot3(ECEF(:,1), ECEF(:,2), ECEF(:,3),'LineWidth',2,'color','k')
hold on
grid on
plot_globe()


subplot(2,2,3:4)
ground_track(lat,lon)
