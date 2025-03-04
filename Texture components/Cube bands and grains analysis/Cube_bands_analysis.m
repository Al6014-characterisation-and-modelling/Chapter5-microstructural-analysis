%% Plotting the texture components from 2D EBSD maps

%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98])};
% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Import the sample I want to plot and the direction of the cross section

alloy = 'NoBA'; % here I choose BA or NoBA

direction = 'RD'; %here I choose RD or TD cross section direction


if strcmp(alloy, 'BA') && strcmp(direction, 'TD')
    fname = 'BA_NDTD.ctf'
    rot1 = rotation('Euler', 90*degree, 90*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 1.1*degree); 

    region = [80 -2964 6800 986]; % large BA_TD

    figname = 'BA_TD_'
end

if strcmp(alloy, 'BA') && strcmp(direction, 'RD')
    fname = ['BA_NDRD.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree);
    rot2 = rotation('Euler', 0*degree, 0*degree, -0.2*degree); 

    region = [3500 -1176 6800 986] %[700 -1176 15000 986] large BA_RD

    figname = 'BA_RD_'

end


if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD')
    fname = ['NoBA_NDTD.ctf']

    rot1 = rotation('Euler', 90*degree, 90*degree, 0*degree); % Apply yhis to
    rot2 =rotation('Euler', 0*degree, 0*degree, -0.8*degree); 

    region = [0 -1200 6800 940]; % large NoBA_TD - use same area than BA map

    figname = 'NoBA_TD_'

end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD')
    fname = ['NoBA_NDRD.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 =rotation('Euler', 0*degree, 0*degree, 0.5*degree); 

    region = [1150 -980 6800 940]; %[250 -980 14000 940] larger area

    figname = 'NoBA_RD_'

end

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
    'convertSpatial2EulerReferenceFrame');
 % 'convertEuler2SpatialReferenceFrame');

 %%  Plotting convention and rotations

ebsd = rotate(ebsd,rot1,'keepXY'); % rotate the orientation data
ebsd = rotate(ebsd,rot2); % rotate the orientation data 

%% Mark and crop region defined before 
image;
plot(ebsd)
rectangle('position',region,'edgecolor','r','linewidth',2)
cropped_ebsd = ebsd(inpolygon(ebsd,region));
ebsd=cropped_ebsd

saveas(gcf, [figname 'marked_EBSD.png']);

%% Plot ebsd map of the selected region
IPF_map(ebsd, 'Aluminium', vector3d.Z)
saveas(gcf, [figname 'ebsd_ipfZ.png']);


%% Crop section and save ebds map

cropped_ebsd = ebsd(inpolygon(ebsd,region));
ebsd=cropped_ebsd
IPF_map(ebsd, 'Aluminium', vector3d.Z)
saveas(gcf, [figname 'ebsd_ipfZ.png']);

%% Visualise subsections

%region_surface = [80 -2214 6800 250]; % large BA_TD - surface
%region_surface = [0 -510 6800 250]; % large NoBA_TD - surface

%region_subsurface = [80 -2330 6800 120]; % large BA_TD - subsurface


%rectangle('position',region,'edgecolor','r','linewidth',2)
%cropped_ebsd = ebsd(inpolygon(ebsd,region));
%ebsd2=cropped_ebsd
%IPF_map(ebsd2, 'Aluminium', RD)

%region_bulk = [80 -2714 6800 500]; % large BA_TD - bulk
%region_bulk = [0 -1010 6800 500]; % large NoBA_TD - bulk

%Sample C
%region_centered = [300 -400 350 380]; % bent sample C 

%Sample D
%region_centered = [400 -430 350 380]; % bent sample D


%% Visualise subsections
%IPF_map(ebsd, 'Aluminium', RD)
%rectangle('position',region_surface,'edgecolor','black','linewidth',2)
%rectangle('position',region_subsurface,'edgecolor','black','linewidth',2)
%rectangle('position', region_bulk,'edgecolor','black','linewidth',2)
%rectangle('position', region_centered,'edgecolor','black','linewidth',2)

%% Get texture plots only from the cropped area
%cropped_ebsd = ebsd(inpolygon(ebsd,region));
%cropped_ebsd = ebsd(inpolygon(ebsd,region_subsurface));
%cropped_ebsd = ebsd(inpolygon(ebsd,region_bulk));
%cropped_ebsd = ebsd(inpolygon(ebsd,region_centered));

%IPF_map(cropped_ebsd, 'Aluminium', RD)
%ebsd=cropped_ebsd
%IPF_map(ebsd, 'Aluminium', RD)

%% Calculate a list of grains

[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, ...
    'angle',5*degree); 

%% Plot cube grains

CubicCS = ebsd('Aluminium').CS;
tol_angle = 10;
% Define grid size (adjust as needed)
%grid_size = [20 6];  % 20 points in x, 6 points in y
grid_size = [1 20]; 
% Initialize array to store cube percentage in each region
cube_fraction_map = zeros(grid_size(2), grid_size(1));

% Loop through each grid cell (based on point coordinates)
for i = 1:grid_size(1)
  for j = 1:grid_size(2)
    % Define x and y range for current grid cell
    x_min = min(ebsd.x) + (i-1) * (range(ebsd.x) / grid_size(1));
    x_max = x_min + (range(ebsd.x) / grid_size(1));
    y_min = min(ebsd.y) + (j-1) * (range(ebsd.y) / grid_size(2));
    y_max = y_min + (range(ebsd.y) / grid_size(2));
    
    % Select points within the grid cell
    ebsd_grid_cell = ebsd(ebsd.x >= x_min & ebsd.x <= x_max & ...
                           ebsd.y >= y_min & ebsd.y <= y_max);
    
    % Calculate cube fraction in this cell
     %ori_cube = orientation.byMiller([0 0 1], [1 0 0], CubicCS);
    ori_cube =  orientation.byEuler(0*degree, 0*degree, 0*degree, CubicCS);
    ebsd_cell_cube = ebsd_grid_cell.findByOrientation(ori_cube, tol_angle*degree);
    cube_fraction_map(j,i) = length(ebsd_cell_cube) / length(ebsd_grid_cell) * 100; % Convert to percentage
  end
end

% Plot the EBSD map
figure;
ebsd_cube = ebsd('Aluminium').findByOrientation(ori_cube, tol_angle*degree);
IPF_map(ebsd_cube, ebsd.orientations, vector3d.Z);
saveas(gcf, [figname 'cube_grains_ebsd.png']);



%% Plot cube bands

CubicCS = ebsd('Aluminium').CS;

% Define tolerance angle for cube orientation
tol_angle = 10;

% Define grid size (adjust as needed)
grid_size = [10 10]; %COLUMNS, ROW % 10 points in x, 10 points in y

% Initialize arrays to store cube percentage in each region
cube_fraction_map = zeros(grid_size(2), grid_size(1));

% Loop through each grid cell (based on point coordinates)
for i = 1:grid_size(1)
    for j = 1:grid_size(2)
        % Define x and y range for current grid cell
        x_min = min(ebsd.x) + (i-1) * (range(ebsd.x) / grid_size(1));
        x_max = x_min + (range(ebsd.x) / grid_size(1));
        y_min = min(ebsd.y) + (j-1) * (range(ebsd.y) / grid_size(2));
        y_max = y_min + (range(ebsd.y) / grid_size(2));

        % Select points within the grid cell
        ebsd_grid_cell = ebsd(ebsd.x >= x_min & ebsd.x <= x_max & ...
            ebsd.y >= y_min & ebsd.y <= y_max);

        % Calculate cube fraction in this cell for cube component 
        ori_cube = orientation.byEuler(0*degree, 0*degree, 0*degree, CubicCS);
        ebsd_cell_cube = ebsd_grid_cell.findByOrientation(ori_cube, tol_angle*degree);
        cube_fraction_cube = length(ebsd_cell_cube) / length(ebsd_grid_cell) * 100; % Convert to percentage

        % Calculate cube fraction in this cell for cube22ND component 
        ori_cube22ND = orientation.byEuler(22*degree, 0*degree, 0*degree, CubicCS);
        ebsd_cell_cube22ND = ebsd_grid_cell.findByOrientation(ori_cube22ND, tol_angle*degree);
        cube_fraction_cube22ND = length(ebsd_cell_cube22ND) / length(ebsd_grid_cell) * 100; % Convert to percentage

        % Calculate cube fraction in this cell for cube22RD component 
        ori_cube22RD = orientation.byEuler(0*degree, 22*degree, 0*degree, CubicCS);
        ebsd_cell_cube22RD = ebsd_grid_cell.findByOrientation(ori_cube22RD, tol_angle*degree);
        cube_fraction_cube22RD = length(ebsd_cell_cube22RD) / length(ebsd_grid_cell) * 100; % Convert to percentage

         % Calculate cube fraction in this cell for cube22RD component 
        ori_cube22TD = orientation.byEuler(90*degree, 22*degree, 0*degree, CubicCS);
        ebsd_cell_cube22TD = ebsd_grid_cell.findByOrientation(ori_cube22TD, tol_angle*degree);
        cube_fraction_cube22TD = length(ebsd_cell_cube22TD) / length(ebsd_grid_cell) * 100; % Convert to percentage

         % Calculate cube fraction in this cell for cube22RD component 
        ori_cube45ND = orientation.byMiller([1 0 0], [0 1 1], CubicCS);
        ebsd_cell_cube45ND = ebsd_grid_cell.findByOrientation(ori_cube45ND, tol_angle*degree);
        cube_fraction_cube45ND = length(ebsd_cell_cube45ND) / length(ebsd_grid_cell) * 100; % Convert to percentage


        % Combine cube and cube22ND fractions
        cube_fraction_map(j,i) = cube_fraction_cube; % + cube_fraction_cube22ND + cube_fraction_cube22RD + cube_fraction_cube22TD + cube_fraction_cube45ND;
    end
end

% Plot the EBSD map with IPF coloring
figure;
IPF_map(ebsd, ebsd.orientations, vector3d.Z);

% Overlay the grid and cube fraction values with colormap and transparency
hold on;

% Define the maximum value for red color
max_red_value = 12; % Cube fraction values above this will be colored as red

% Calculate the minimum and maximum cube fraction values
min_cube_fraction = 0; % Minimum value for the colormap (0%)
max_cube_fraction = 12; % Maximum value for the colormap (20%)

% Define a colormap from yellow to red with fewer colors
nColors = 16; % Number of colors in the colormap
color_map = zeros(nColors, 3);
color_map(:, 1) = linspace(1, 1, nColors); % Red channel stays at 1 (full red)
color_map(:, 2) = linspace(1, 0, nColors); % Green channel decreases from 1 to 0 (yellow to red)
color_map(:, 3) = linspace(0, 0, nColors); % Blue channel stays at 0 (no blue)

% Loop through each grid cell to plot
for i = 1:grid_size(1)
    for j = 1:grid_size(2)
        % Calculate the boundaries of each grid cell
        x_min = min(ebsd.x) + (i-1) * (range(ebsd.x) / grid_size(1));
        x_max = x_min + (range(ebsd.x) / grid_size(1));
        y_min = min(ebsd.y) + (j-1) * (range(ebsd.y) / grid_size(2));
        y_max = y_min + (range(ebsd.y) / grid_size(2));

        % Determine the color based on cube fraction
        cube_percent = cube_fraction_map(j,i);
        
        % Map cube fraction to the colormap index
        colormap_idx = round(interp1([min_cube_fraction, max_cube_fraction], [1, nColors], cube_percent));

        % Clip values outside the colormap range
        colormap_idx = max(1, min(colormap_idx, nColors));

        % Get the color value from the colormap
        color_val = color_map(colormap_idx, :);

        % Draw filled rectangle with interpolated color and transparency
        fill([x_min x_max x_max x_min], [y_min y_min y_max y_max], color_val, 'EdgeColor', 'none', 'FaceAlpha', 0.8);

        % Draw the grid cell border
        rectangle('Position', [x_min, y_min, (x_max-x_min), (y_max-y_min)], 'EdgeColor', 'black');

        % Calculate the center position of each grid cell
        x_pos = (x_min + x_max) / 2;
        y_pos = (y_min + y_max) / 2;

        % Overlay the cube fraction value
        text(x_pos, y_pos, num2str(cube_fraction_map(j,i), '%.2f'), ...
            'Color', 'black', 'FontSize', 10, 'HorizontalAlignment', 'center');
        hScale = findobj(gca,'Tag','mtexScale');
        uistack(hScale,'top');

    end
end

% Set the colormap to the custom yellow-to-red colormap
colormap(color_map);

% Add colorbar with appropriate ticks and labels
caxis([min_cube_fraction, max_cube_fraction]); % Set colorbar limits based on min and max values

% Generate colorbar ticks and labels (without decimals)
colorbar_ticks = linspace(min_cube_fraction, max_cube_fraction, 5);
colorbar_tick_labels = arrayfun(@(x) sprintf('%.0f%%', x), colorbar_ticks, 'UniformOutput', false);

colorbar('Ticks', colorbar_ticks, 'TickLabels', colorbar_tick_labels);

hold off;

saveas(gcf, [figname 'cube_bands.png']);

%%
