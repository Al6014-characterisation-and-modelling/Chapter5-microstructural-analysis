%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.


%% Import the sample I want to plot and the direction of the cross section

alloy = 'NoBA'; % here I choose BA or NoBA

direction = 'TD'; %here I choose RD or TD cross section direction

size = 'small' %here I choose between small or large EBSD

export_csv = 'yes' % keep this with 'yes' to export equivalent diameter and aspect ratio in a csv file

cube_threshold = 0.5
tol_angle = 10 %select the tolerance degree for cube orientation

if strcmp(alloy, 'BA') && strcmp(direction, 'TD') && strcmp(size, 'large')
    fname = 'BA_NDTD_large.ctf'
    rot1 = rotation('Euler', 90*degree, 90*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 1.1*degree); 

    region = [80 -2964 6800 986]; % large BA_TD

    figname = 'BA_TD_large_'
    cleaning = 0
end


if strcmp(alloy, 'BA') && strcmp(direction, 'TD') && strcmp(size, 'small')
    fname = 'BA_NDTD_small.ctf'
    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 0*degree); 

    %region = [0 -970 1380 920]; % BA_TD
    region = [0 -970 1380 920] %region for statistics

    figname = 'BA_TD_small_'
    cleaning = 10
end


if strcmp(alloy, 'BA') && strcmp(direction, 'RD') && strcmp(size, 'large')
    fname = ['BA_NDRD_large.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree);
    rot2 = rotation('Euler', 0*degree, 0*degree, -0.2*degree); 

    region = [3500 -1176 6800 986] %[700 -1176 15000 986] large BA_RD

    figname = 'BA_RD_large_'
    cleaning = 0

end

if strcmp(alloy, 'BA') && strcmp(direction, 'RD') && strcmp(size, 'small')
    fname = ['BA_NDRD_small.ctf']
    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    
    %region = [0 -1040 1500 950]; %[xmin ymin xmax-xmin, ymax-ymin] BA_RD
    region = [0 -1040 1380 950] %region for statistics

    figname = 'BA_RD_small'
    cleaning = 10

end


if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD') && strcmp(size, 'large')
    fname = ['NoBA_NDTD_large.ctf']

    rot1 = rotation('Euler', 90*degree, 90*degree, 0*degree); % Apply yhis to
    rot2 =rotation('Euler', 0*degree, 0*degree, -0.8*degree); 

    region = [100 -1190 6800 930] %[100 -1200 6800 940]; % large NoBA_TD - use same area than BA map

    figname = 'NoBA_TD_large'
    cleaning = 0

end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD') && strcmp(size, 'small')
    fname = ['NoBA_NDTD_small.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    %region = [0 -1010 1450 985]; % NoBA_TD
    region = [0 -1010 1380 985] %region for statistics

    figname = 'NoBA_TD_small'
    cleaning = 10

end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD') && strcmp(size, 'large')
    fname = ['NoBA_NDRD_large.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 =rotation('Euler', 0*degree, 0*degree, 0.5*degree); 

    region = [1150 -980 6800 940]; %[250 -980 14000 940] larger area

    figname = 'NoBA_RD_large'
    cleaning = 0

end


if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD') && strcmp(size, 'small')
    fname = ['NoBA_NDRD_small.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    %region = [0 -1000 1400 950]; % NoBA_RD
    region = [0 -1000 1380 950] %region for statistics

    figname = 'NoBA_RD_small'
    cleaning = 10

end


%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98])};
% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
    'convertSpatial2EulerReferenceFrame');
 % 'convertEuler2SpatialReferenceFrame');


%%  Plotting convention and rotations

ebsd = rotate(ebsd,rot1,'keepXY'); % rotate the orientation data
ebsd = rotate(ebsd,rot2); % rotate the orientation data 
%% Plot the complete EBSD map:

plot(ebsd)

ori = ebsd('Aluminium').orientations;

IPF_map(ebsd, 'Aluminium', vector3d.Z)
saveas(gcf, [figname 'ebsd_ipfZ.png']);

%% Mark region
image;
plot(ebsd)
rectangle('position',region,'edgecolor','r','linewidth',2)
cropped_ebsd = ebsd(inpolygon(ebsd,region));
ebsd=cropped_ebsd

saveas(gcf, [figname 'marked_EBSD.png']);

%% Plot new ebsd
image;

%IPF_map(cropped_ebsd, 'Aluminium', ND)
ebsd=cropped_ebsd
IPF_map(ebsd, 'Aluminium', vector3d.Z)

saveas(gcf, [figname 'ebsd_ipfZ.png']);

%% % Calculate a list of grains

[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, ...
    'angle',5*degree);

%% Optional Cleaning the data
% Removing small grains – denoise
clean_ebsd = clean_grains(ebsd, grains, cleaning);  % grains with "cleaning" pixels or less will be removed

[cleaned_grains, clean_ebsd.grainId, clean_ebsd.mis2mean] = calcGrains(clean_ebsd, 'angle',5*degree);   % 5° tolerance

image;
IPF_grains_map(cleaned_grains, 'Aluminium', vector3d.Z)
saveas(gcf, [figname 'ebsd_clean_grains.png']);

%% Threshold the cube orientation
CubicCS = ebsd('Aluminium').CS;
ori_cube = orientation.byEuler(0*degree, 0*degree, 0*degree, CubicCS);

% Extract grain IDs from the cleaned grains
grainIDs = cleaned_grains.id;
cubeFraction = zeros(length(grainIDs), 1);

% Loop over each grain to compute its cube fraction
for k = 1:length(grainIDs)
    % Get the indices of pixels belonging to grain k in the cleaned EBSD data
    idx = (clean_ebsd.grainId == grainIDs(k));
    ebsd_grain = clean_ebsd(idx);
    
    % Restrict to the Aluminium phase to avoid comparing different phases
    ebsd_grain = ebsd_grain('Aluminium');
    
    if isempty(ebsd_grain)
        cubeFraction(k) = 0;
    else
        % Find the pixels in this grain that are cube oriented
        ebsd_grain_cube = ebsd_grain.findByOrientation(ori_cube, tol_angle*degree);
        cubeFraction(k) = length(ebsd_grain_cube) / length(ebsd_grain);
    end
end

% Select only those grains that are predominantly cube oriented
isCubeGrain = cubeFraction >= cube_threshold;
cube_grains = cleaned_grains(isCubeGrain);
%%

image;
IPF_grains_map(cube_grains, 'Aluminium', vector3d.Z)
saveas(gcf, [figname 'cube_grains.png']);

%% Calculate Equivalent Diameter
equivDiameters = 2 * sqrt(cube_grains.area ./ pi);

% Display histogram of equivalent diameters
figure;
histogram(equivDiameters, 'BinWidth', 2); % Adjust bin width as needed
xlabel('Equivalent Diameter (\mum)');
ylabel('Frequency');
title('Histogram of Cube Grain Equivalent Diameters');

%% Plot grains colored by equivalent diameter
figure;
plot(cube_grains, equivDiameters);
cb = colorbar;
% Set colorbar limits (adjust max_value as needed)
max_value = 70;  % Example max value in microns (adjust as needed)
caxis([0 max_value]); % Limits the colormap scale

cb.Label.String = 'Equivalent Diameter (\mum)';
cb.Label.FontSize = 20; % Increase colorbar label font size
cb.FontSize = 25; % Increase colorbar tick labels font size

saveas(gcf, [figname 'Equivalent_diameter.png']);


%% Calculate aspect ratio
aspectRatio = cube_grains.aspectRatio % Major axis length


%% Plot grains with aspect ratio as color
figure;
plot(cube_grains, aspectRatio);

cb = colorbar;

% Set colorbar limits (adjust max_value as needed)
max_value = 2;  % Maximum value for color scaling
caxis([1 max_value]); % Limits the colormap scale

% Create a custom colormap: blue at min, red at max
custom_colormap = [0 0 1; 1 0 0]; % [Blue; Red] in RGB
colormap(interp1([1, 2], custom_colormap, linspace(1, 2, 256))); 

% Define tick values and labels
tick_values = linspace(1, max_value, 5); % Adjust for better spacing
tick_labels = string(tick_values); % Convert to string
tick_labels(end) = '>2'; % Replace last tick label

% Apply custom ticks and labels
cb.Ticks = tick_values;
cb.TickLabels = tick_labels;

% Improve colorbar labels
cb.Label.String = 'Aspect Ratio';
cb.Label.FontSize = 20; % Increase colorbar label font size
cb.FontSize = 25; % Increase colorbar tick labels font size

saveas(gcf, [figname 'aspect_ratio.png']);


%% Save equivalen diamters and aspect ratios into a csv file

if strcmp(export_csv, 'yes')
    T = table(cube_grains.id, equivDiameters, aspectRatio, ...
              'VariableNames', {'GrainID', 'EquivalentDiameter', 'AspectRatio'});
    writetable(T, [figname 'GrainsData.csv']);
end
