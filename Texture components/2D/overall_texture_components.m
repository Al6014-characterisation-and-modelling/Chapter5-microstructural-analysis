%% Plotting the texture components from 2D EBSD maps with top & bottom surfaces combined

%% Specify Crystal and Specimen Symmetries
% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98])};
% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Import the sample I want to plot and the direction of the cross section
alloy     = 'NoBA';    % here I choose BA or NoBA
direction = 'RD';    % here I choose RD or TD cross section direction

% initialize figname
figname = '';

if strcmp(alloy, 'BA') && strcmp(direction, 'TD')
    fname  = 'BA/BA_NDTD.ctf';
    rot1   = rotation('Euler', 90*degree, 90*degree,   0*degree);
    rot2   = rotation('Euler',  0*degree,   0*degree, 1.1*degree);
    figname = 'BA_TD_';
    delta= 25 % delta 25 for BA and delta 24 for No BA
    % Define regions for top & bottom surfaces
    region = [80 -2964 6800 986]; % large BA_TD
end

if strcmp(alloy, 'BA') && strcmp(direction, 'RD')
    fname   = 'BA/BA_NDRD.ctf';
    rot1    = rotation('Euler',   0*degree,   0*degree,   0*degree);
    rot2    = rotation('Euler',   0*degree,   0*degree,  -0.2*degree);
    delta= 25 % delta 25 for BA and delta 24 for No BA
    region = [3500 -1176 6800 986] %[700 -1176 15000 986] large BA_RD
    figname = 'BA_RD_';
end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD')
    fname   = 'NoBA/NoBA_NDTD.ctf';
    rot1    = rotation('Euler', 90*degree, 90*degree,   0*degree);
    rot2    = rotation('Euler',  0*degree,   0*degree, -0.8*degree);
    delta= 24 % delta 25 for BA and delta 24 for No BA
    region = [0 -1200 6800 940]; % large NoBA_TD - use same area than BA map
    figname = 'NoBA_TD_';
end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD')
    fname   = 'NoBA/NoBA_NDRD.ctf';
    rot1    = rotation('Euler',   0*degree,   0*degree,   0*degree);
    rot2    = rotation('Euler',   0*degree,   0*degree,   0.5*degree);
    delta= 24 % delta 25 for BA and delta 24 for No BA
    region = [1150 -980 6800 940]; %[250 -980 14000 940] larger area
    figname = 'NoBA_RD_';
end

%% Helper function to load, rotate and crop a region
function e = load_and_crop_region(fname, CS, rot1, rot2, region)
    e = EBSD.load(fname, CS, 'interface', 'ctf', 'convertSpatial2EulerReferenceFrame');
    e = rotate(e, rot1, 'keepXY');
    e = rotate(e, rot2);
    e = e(inpolygon(e, region));
end

%% Load, rotate and crop bulk ebsd
ebsd_bulk    = load_and_crop_region(fname, CS, rot1, rot2, region);

%% Plotting the complete EBSD map of combined surface
figure;
plot(ebsd_bulk);
IPF_map(ebsd_bulk, 'Aluminium', vector3d.Z);
saveas(gcf, [figname 'bulk_ebsd_ipfZ.png']);

%% Calculate grains on the combined surface
[grains_bulk, ebsd_bulk.grainId, ebsd_bulk.mis2mean] = calcGrains(ebsd_bulk, 'angle', 5*degree);
psi = calcKernel(grains_bulk('Aluminium').meanOrientation);

%% 5. Compute ODF
ori = ebsd_bulk('Aluminium').orientations;
ori.SS = specimenSymmetry('orthorhombic');
odf = calcDensity(ori, 'kernel', psi);
odf.SS = specimenSymmetry('222');

%% 6. Define Texture Components
prefs = {
    orientation('Euler',0,0,0,'degree'), ...      % Cube
    orientation('Euler',59,34,65,'degree'), ...   % S
    orientation('Euler',45,15,10,'degree'), ...   % Q
    orientation('Euler',65,45,0,'degree'), ...    % P
    orientation('Euler',90,30,45,'degree'), ...   % Copper
    orientation('Euler',53,36,60,'degree'), ...   % R
    orientation('Euler',22,0,0,'degree'), ...     % cube22ND
    orientation('Euler',0,22,0,'degree'), ...     % cube22RD
    orientation('Euler',90,22,0,'degree'), ...    % cube22TD
    orientation('Euler',45,0,0,'degree'), ...     % cube45ND
    orientation('Euler',0,55,45,'degree'), ...    % E1
    orientation('Euler',60,55,45,'degree'), ...   % E2
    orientation('Euler',30,55,45,'degree'), ...   % F1
    orientation('Euler',90,55,45,'degree'), ...   % F2
    orientation('Euler',90,27,45,'degree'), ...   % D
    orientation('Euler',35,45,0,'degree'), ...    % Brass
    orientation('Euler',0,45,0,'degree') ...     % Goss
};
names = {'Cube','S','Q','P','Copper','R','cube22ND','cube22RD','cube22TD','cube45ND','E1','E2','F1','F2','D','Brass','Goss'};

%% Calculation of texture components
disp('Total percentage of texture components: ');
%CalclComponents returns the prefered orientation and the percentage of 
% orientations that crawled to each of them.
[component, vol] = calcComponents(odf);
component , vol*100 
% this will give me the total percentage of texture components


%% 8. Calculate the percentage of each texture component


% Define the names of the variables
variableNames = {'Cube', 'S', 'Q', 'P', 'copper', 'R', 'cube22ND', 'cube22RD', 'cube22TD', 'cube45ND','E1','E2', 'F1', 'F2', 'D', 'Brass', 'Goss'};


% Calculate the percentages of each texture component
%pref_ori=[cube, S, Q, P, copper, R, cube22ND, cube22RD, cube22TD, cube45ND, E1, E2, F1, F2, D, Brass, Goss];
comp = volume(odf, [prefs{:}], delta/10 * degree)*100;

% Loop through 'comp' that displays values with their names
for i = 1:length(comp)
    variableName = variableNames{i};
    disp([variableName, ': ', num2str(comp(i), '%.2f')]);
end

% Display the total percentage
disp(['Total percentage: ', num2str(sum(comp), '%.1f')]);
%This total percentage should be similar to total percentage of texture
%components, if not, modify delta slightly


%% 9. Save texture component percentages to CSV file

% Create a table with texture component names and their corresponding percentages
textureComponentsTable = table(variableNames', comp', 'VariableNames', {'TextureComponent', 'Percentage'});

% Write the table to a CSV file
writetable(textureComponentsTable, [figname 'texture_components_overall.csv']);

