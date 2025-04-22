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
    fname = 'BA/BA_NDTD.ctf'
    rot1 = rotation('Euler', 90*degree, 90*degree, 0*degree); 
    rot2 = rotation('Euler', 0*degree, 0*degree, 1.1*degree); 

    %region = [80 -2964 6800 986]; % large BA_TD
    region = [80 -2764 6800 586]; % large BA_TD - bulk

    % find a way of adding top and bottom sufaces and plot the 111 pole figure

    figname = 'BA_TD_'
end

if strcmp(alloy, 'BA') && strcmp(direction, 'RD')
    fname = ['BA/BA_NDRD.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree);
    rot2 = rotation('Euler', 0*degree, 0*degree, -0.2*degree); 

    %region = [3500 -1176 6800 986] %[700 -1176 15000 986] large BA_RD
    region = [3500 -976 6800 586] % bulk BA_RD

    figname = 'BA_RD_'

end


if strcmp(alloy, 'NoBA') && strcmp(direction, 'TD')
    fname = ['NoBA/NoBA_NDTD.ctf']

    rot1 = rotation('Euler', 90*degree, 90*degree, 0*degree); % Apply this to
    rot2 =rotation('Euler', 0*degree, 0*degree, -0.8*degree); 

    %region = [0 -1200 6800 940]; % large NoBA_TD - use same area than BA map
    region = [0 -1000 6800 540]; % bulk NoBA TD 

    figname = 'NoBA_TD_'

end

if strcmp(alloy, 'NoBA') && strcmp(direction, 'RD')
    fname = ['NoBA/NoBA_NDRD.ctf']

    rot1 = rotation('Euler', 0*degree, 0*degree, 0*degree); 
    rot2 =rotation('Euler', 0*degree, 0*degree, 0.5*degree); 

    % region = [1150 -980 6800 940]; %[250 -980 14000 940] larger area
    region = [1150 -780 6800 540]; % bulk NoBA RD

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

%% Plotting the complete EBSD map:
image;
plot(ebsd)

IPF_map(ebsd, 'Aluminium',vector3d.Z)

%% Mark and crop region defined before 
image;
plot(ebsd)
rectangle('position',region,'edgecolor','r','linewidth',2)
cropped_ebsd = ebsd(inpolygon(ebsd,region));
ebsd=cropped_ebsd

saveas(gcf, [figname 'bulk_marked_EBSD.png']);
%% Plot ebsd map of the selected region
IPF_map(ebsd, 'Aluminium', vector3d.Z)
saveas(gcf, [figname 'bulk_ebsd_ipfZ.png']);


%% Mark and crop subsections

%Bands of texture
%region = [80 -2214 6800 250]; % large BA_TD - surface
%region = [80 -2714 6800 500]; % large BA_TD - bulk
%region = [0 -510 6800 250]; % large NoBA_TD - surface
%region = [0 -1010 6800 500]; % large NoBA_TD - bulk

%region_surface = [80 -2214 6800 250]; % large BA_TD - surface
%region_surface = [0 -510 6800 250]; % large NoBA_TD - surface

%region_subsurface = [80 -2330 6800 120]; % large BA_TD - subsurface
%region_bulk = [80 -2714 6800 500]; % large BA_TD - bulk
%region_bulk = [0 -1010 6800 500]; % large NoBA_TD - bulk

%rectangle('position',region,'edgecolor','r','linewidth',2)
%cropped_ebsd = ebsd(inpolygon(ebsd,region));
%ebsd2=cropped_ebsd
%IPF_map(ebsd2, 'Aluminium', ND)

%% Visualise subsections
%IPF_map(ebsd, 'Aluminium', ND)
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

%% Plotting the grains in the selected area 
% Calculate a list of grains

[grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, ...
    'angle',5*degree); 

%% Plotting ODF
psi=calcKernel(grains('Aluminium').meanOrientation);  %I'm using the raw data

ori = ebsd('Aluminium').orientations;
ori.SS = specimenSymmetry('orthorhombic'); 

odf = calcDensity(ori,'kernel',psi); 

figure;
plot(odf,'phi2',[0]*degree,'antipodal','linewidth',1,'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0,6]);
mtexColorbar ('location','eastoutside','title','mrd');
colormap("sky");
ax = gca;
%ax.YDir = 'reverse';
saveas(gca, [figname 'bulk_odf_0.png']);


figure;
plot(odf,'phi2',[45]*degree,'antipodal','linewidth',1,'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0,6]);
mtexColorbar ('location','eastoutside','title','mrd');
colormap("sky");
ax = gca;
%ax.YDir = 'reverse';
saveas(gca, [figname 'bulk_odf_45.png']);

figure;
plot(odf,'phi2',[65]*degree,'antipodal','linewidth',1,'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0,6]);
mtexColorbar ('location','eastoutside','title','mrd');
colormap("sky");
ax = gca;
%ax.YDir = 'reverse';
saveas(gca, [figname 'bulk_odf_65.png']);


%% Pole Figures

%cs=ebsd('Aluminium').CS;
%ss = specimenSymmetry('monoclinic');

ori.SS = specimenSymmetry('monoclinic'); %We use monoclinic to don't force the symmmetry
cs = ori.CS;

h1=Miller(0,0,1,cs); % (0,0,1) pole
h2=Miller(1,1,0,cs); % (1,1,0) pole
h3=Miller(1,1,1,cs); % (1,1,1) pole

figure
%plotPDF(ori,[h1, h2, h3],'antipodal', 'contourf', 'minmax','complete','xAxisDirection', 'north', 'yAxisDirection', 'east'); 
plotPDF(ori, [h1, h2, h3], 'contourf', 'minmax', 'complete', 'xAxisDirection', 'north', 'yAxisDirection', 'east'); 

%view(90, -90);
mtexColorbar ('location','eastoutside','title','mrd');
setColorRange([0 3])
colormap("sky");

saveas(gcf, [figname 'bulk_pole_figures.png']);
