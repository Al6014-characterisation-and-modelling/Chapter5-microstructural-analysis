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
    delta= 22;
    % Define regions for top & bottom surfaces
    region_top    = [80  -2178 6800 200];   % BA_TD - top surface
    region_bottom = [80  -2964 6800 200];   % BA_TD - bottom surface
    figname = 'BA_TD_surface_';

    %delta= 21;
    %region_top    = [80  -2278 6800 100];   % BA_TD - top subsurface
    %region_bottom = [80  -2764 6800 100];   % BA_TD - bottom subsurface
    %figname = 'BA_TD_subsurface_';
end

if strcmp(alloy, 'BA') && strcmp(direction, 'RD')
    fname   = 'BA/BA_NDRD.ctf';
    rot1    = rotation('Euler',   0*degree,   0*degree,   0*degree);
    rot2    = rotation('Euler',   0*degree,   0*degree,  -0.2*degree);

    delta= 20; 
    region_top    = [3500 -390 6800 200];    % BA_RD - top (use same as large map)
    region_bottom = [3500 -1176 6800 200];   % BA_RD - bottom (use same as large map)
    figname = 'BA_RD_surface_';

    %delta= 18;
    %region_top    = [3500 -490 6800 100];    % BA_RD - top (use same as large map)
    %region_bottom = [3500 -976 6800 100];   % BA_RD - bottom (use same as large map)
    %figname = 'BA_RD_subsurface_';
end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD')
    fname   = 'NoBA/NoBA_NDTD.ctf';
    rot1    = rotation('Euler', 90*degree, 90*degree,   0*degree);
    rot2    = rotation('Euler',  0*degree,   0*degree, -0.8*degree);

    delta= 20; 
    region_top    = [0 -460 6800 200];    % NoBA_TD - top surface
    region_bottom = [0 -1200 6800 200];    % NoBA_TD - bottom surface
    figname = 'NoBA_TD_surface_';

    %delta= 21;
    %region_top    = [0 -560 6800 100];    % NoBA_TD - top surface
    %region_bottom = [0 -1000 6800 100];    % NoBA_TD - bottom surface
    %figname = 'NoBA_TD_subsurface_';
end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD')
    fname   = 'NoBA/NoBA_NDRD.ctf';
    rot1    = rotation('Euler',   0*degree,   0*degree,   0*degree);
    rot2    = rotation('Euler',   0*degree,   0*degree,   0.5*degree);

    delta= 22;
    region_top    = [1150 -240 6800 200];   % NoBA_RD - top surface
    region_bottom = [1150 -980 6800 200];   % NoBA_RD - bottom surface  
    figname = 'NoBA_RD_surface_';
        
    %delta= 21.5;
    %region_top    = [1150 -340 6800 100];   % NoBA_RD - top surface
    %region_bottom = [1150 -780 6800 100];   % NoBA_RD - bottom surface  
    %figname = 'NoBA_RD_subsurface_';
end

%% Helper function to load, rotate and crop a region
function e = load_and_crop_region(fname, CS, rot1, rot2, region)
    e = EBSD.load(fname, CS, 'interface', 'ctf', 'convertSpatial2EulerReferenceFrame');
    e = rotate(e, rot1, 'keepXY');
    e = rotate(e, rot2);
    e = e(inpolygon(e, region));
end

%% Load, rotate and crop top & bottom surfaces
ebsd_top    = load_and_crop_region(fname, CS, rot1, rot2, region_top);
ebsd_bottom = load_and_crop_region(fname, CS, rot1, rot2, region_bottom);

% Concatenate into one EBSD object for surfaces
ebsd_surface = [ebsd_top; ebsd_bottom];

%% Plotting the complete EBSD map of combined surface
figure;
plot(ebsd_surface);
IPF_map(ebsd_surface, 'Aluminium', vector3d.Z);
saveas(gcf, [figname 'surface_ebsd_ipfZ.png']);


%% 4. Compute Grains & Kernel Compute Grains & Kernel Compute Grains & Kernel
[grains, ebsd_surface.grainId, ebsd_surface.mis2mean] = calcGrains(ebsd_surface, 'angle', 5*degree);
psi = calcKernel(grains('Aluminium').meanOrientation);

%% 5. Extract Orientations from Combined Surface
ebsd_orient = ebsd_surface('Aluminium').orientations;

%% 6. Define Texture Components (Euler angles)
prefs = [
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
    orientation('Euler',0,45,0,'degree') ...      % Goss
];
names = {'Cube','S','Q','P','Copper','R','cube22ND','cube22RD','cube22TD', ...
         'cube45ND','E1','E2','F1','F2','D','Brass','Goss'};

%% 7. Calculate Texture Component Percentages Directly from EBSD Orientations
% Use MTex's calcComponents on the raw orientations
[componentIDs, vol] = calcComponents(ebsd_surface, prefs, delta*degree);
vol_perc = vol * 100;

%% 8. Display Results and Save to CSV
% Show each component percentage on the command window
disp('Surface Texture Component Percentages:');
for i = 1:numel(names)
    fprintf('%-10s: %6.2f%%', names{i}, vol_perc(i));
end
% Show total
fprintf('Total Surface %%: %6.2f%%', sum(vol_perc));

% Export results to CSV
T = table(names', vol_perc', 'VariableNames', {'TextureComponent','Percentage'});
writetable(T, [figname 'texture_components.csv']);
