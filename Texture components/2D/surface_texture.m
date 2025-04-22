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

    % Define regions for top & bottom surfaces
    region_top    = [80  -2178 6800 200];   % BA_TD - top surface
    region_bottom = [80  -2964 6800 200];   % BA_TD - bottom surface
end

if strcmp(alloy, 'BA') && strcmp(direction, 'RD')
    fname   = 'BA/BA_NDRD.ctf';
    rot1    = rotation('Euler',   0*degree,   0*degree,   0*degree);
    rot2    = rotation('Euler',   0*degree,   0*degree,  -0.2*degree);
    region_top    = [3500 -390 6800 200];    % BA_RD - top 
    region_bottom = [3500 -1176 6800 200];   % BA_RD - bottom
    figname = 'BA_RD_';
end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD')
    fname   = 'NoBA/NoBA_NDTD.ctf';
    rot1    = rotation('Euler', 90*degree, 90*degree,   0*degree);
    rot2    = rotation('Euler',  0*degree,   0*degree, -0.8*degree);
    region_top    = [0 -460 6800 200];    % NoBA_TD - top surface
    region_bottom = [0 -1200 6800 200];    % NoBA_TD - bottom surface
    figname = 'NoBA_TD_';
end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD')
    fname   = 'NoBA/NoBA_NDRD.ctf';
    rot1    = rotation('Euler',   0*degree,   0*degree,   0*degree);
    rot2    = rotation('Euler',   0*degree,   0*degree,   0.5*degree);
    region_top    = [1150 -240 6800 200];   % NoBA_RD - top surface
    region_bottom = [1150 -980 6800 200];   % NoBA_RD - bottom surface            
    figname = 'NoBA_RD_';
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

%% Calculate grains on the combined surface
[grains_surf, ebsd_surface.grainId, ebsd_surface.mis2mean] = calcGrains(ebsd_surface, 'angle', 5*degree);

%% Compute kernel & ODF for the surface orientations
psi_surf    = calcKernel(grains_surf('Aluminium').meanOrientation);
ori_surf    = ebsd_surface('Aluminium').orientations;
ori_surf.SS = specimenSymmetry('monoclinic');
odf_surf    = calcDensity(ori_surf, 'kernel', psi_surf);

%% Plot Surface ODF at phi2 = 0°, 45°, 65°
for phi2 = [0 45 65]*degree
    figure;
    plot(odf_surf, 'phi2', phi2, 'antipodal', ...
         'linewidth',1, 'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0 6]);
    mtexColorbar('location','eastoutside','title','mrd');
    colormap('sky');
    title(sprintf('Surface ODF, \phi_2 = %d°', round(phi2/degree)));
    saveas(gcf, [figname sprintf('surface_odf_%02d.png', round(phi2/degree))]);
end

%% Plot Surface Pole Figures
h1 = Miller(0,0,1, ori_surf.CS);
h2 = Miller(1,1,0, ori_surf.CS);
h3 = Miller(1,1,1, ori_surf.CS);

figure;
plotPDF(ori_surf, [h1 h2 h3], 'contourf','minmax','complete', ...
        'xAxisDirection','north','yAxisDirection','east');
mtexColorbar('location','eastoutside','title','mrd');
setColorRange([0 3]);
colormap('sky');
title('Surface Pole Figures');
saveas(gcf, [figname 'surface_pole_figures.png']);
