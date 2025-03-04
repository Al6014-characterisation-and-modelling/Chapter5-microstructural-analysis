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

alloy = 'BA'; % here I choose BA or NoBA

direction = 'TD'; %here I choose RD or TD cross section direction


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



%% Clean data to plot KAM

figure;
F = halfQuadraticFilter;
F.alpha = 0.5;
ebsdS = smooth(ebsd,F,'fill',grains); % denoise the orientation map

%% Plot KAM
figure;
kam = ebsd.KAM('threshold', 5*degree) / degree;
%kam(kam < 1) = NaN; % Remove noisy data
set(gcf, 'Color', 'k');
% lets plot it
plot(ebsdS,ebsdS.KAM('threshold',5*degree) ./ degree,'micronbar','off')
plot(ebsd,kam,'micronbar','off')
setColorRange([0,5])
mtexColorbar
mtexColorMap LaboTeX
%hold on
%IPF_map(cropped_ebsd, 'Aluminium', vector3d.Z)
%plot(grains.boundary,'lineWidth',1.5)
%hold off


%% Plot the IPF map with transparency on top of the KAM image

figure;
kam = ebsd.KAM / degree;
set(gcf, 'Color', 'k');

% Plot the KAM
plot(ebsd, kam, 'micronbar', 'off')
setColorRange([0, 5])
mtexColorbar
mtexColorMap LaboTeX

hold on
% Create the IPF map
ipfKey = ipfHSVKey(ebsd('indexed'));
colorData = ipfKey.orientation2color(ebsd('indexed').orientations);

% Plot the IPF map as an image with transparency
imageHandle = plot(ebsd('indexed'), colorData);
set(imageHandle, 'FaceAlpha', 0.6); % Adjust transparency level here (0 to 1)

% Optionally, plot the grains boundaries
% plot(grains.boundary, 'lineWidth', 1.5)
hold off


%% Select a grain in the previous map and read the x and y coordinates


IPF_map(cropped_ebsd, 'Aluminium', vector3d.Z)
hold on


imageHandle = plot(ebsd, kam, 'micronbar', 'off');
setColorRange([0, 5])
mtexColorbar
mtexColorMap LaboTeX
set(imageHandle, 'FaceAlpha', 0.7); % Adjust transparency level here (0 to 1)

%x=336, y=-305;
x=693; y=-356;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
%x = 596; y = -286;
x=565; y=-398;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
%x = 340; y = -152;
x=536; y=-342;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
%x = 400; y = -319;
x=496; y=-399;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
%x=570; y=-226;
x=623; y=-260;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
%x=608; y=-233;
x=547; y=-299;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
%x=322; y=-167;
x=529; y=-305;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
x=520; y=-423;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
x=698; y=-373;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
x=569; y=-285;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')


x=551; y=-330;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')


x=618; y=-311;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')

x=595; y=-361;
%plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')

x=558; y=-324;
plot(grains(x,y).boundary,'linewidth',2,'linecolor','black')
hold off
%% Plot the orientation of the selected grain
selectedGrain = grains(x,y);
selectedGrainOrientation = selectedGrain.meanOrientation;
disp(['Selected Grain Orientation: ', char(selectedGrainOrientation)])
%% Plot texture components


% Define the crystal symmetry
CubicCS = ebsd('Aluminium').CS;

% Define all orientation components
ori_cube = orientation.byMiller([0 0 1], [1 0 0], CubicCS);
ori_brass = orientation.byMiller([1 1 0], [1 1 2], CubicCS);
ori_S = orientation.byEuler(59*degree, 34*degree, 65*degree, CubicCS);
ori_copper = orientation.byMiller([1 1 2], [1 1 1], CubicCS);
ori_goss = orientation.byMiller([0 1 1], [1 0 0], CubicCS);
ori_cube45 = orientation.byMiller([1 0 0], [0 1 1], CubicCS);
ori_cube22ND = orientation.byEuler(22*degree, 0*degree, 0*degree, CubicCS);
ori_cube22RD = orientation.byEuler(0*degree, 22*degree, 0*degree, CubicCS);
ori_cube22TD = orientation.byEuler(90*degree, 22*degree, 0*degree, CubicCS);
ori_rotGoss = orientation.byMiller([0 -1 1], [0 1 1], CubicCS);
ori_R = orientation.byMiller([1 2 4], [2 1 1], CubicCS);
ori_Q = orientation.byMiller([0 1 3], [2 -3 1], CubicCS);
ori_P = orientation.byMiller([0 1 1], [1 2 2], CubicCS);
ori_E = orientation.byMiller([1 1 1], [1 1 0], CubicCS);
ori_F = orientation.byMiller([1 1 1], [1 1 2], CubicCS);

% Define the tolerance angle
tol_angle = 15 * degree;

% Extract EBSD data for each orientation component
ebsd_cube = ebsd('Aluminium').findByOrientation(ori_cube, tol_angle);
ebsd_brass = ebsd('Aluminium').findByOrientation(ori_brass, tol_angle);
ebsd_S = ebsd('Aluminium').findByOrientation(ori_S, tol_angle);
ebsd_copper = ebsd('Aluminium').findByOrientation(ori_copper, tol_angle);
ebsd_goss = ebsd('Aluminium').findByOrientation(ori_goss, tol_angle);
ebsd_cube45 = ebsd('Aluminium').findByOrientation(ori_cube45, tol_angle);
ebsd_cube22ND = ebsd('Aluminium').findByOrientation(ori_cube22ND, tol_angle);
ebsd_cube22RD = ebsd('Aluminium').findByOrientation(ori_cube22RD, tol_angle);
ebsd_cube22TD = ebsd('Aluminium').findByOrientation(ori_cube22TD, tol_angle);
ebsd_rotGoss = ebsd('Aluminium').findByOrientation(ori_rotGoss, tol_angle);
ebsd_R = ebsd('Aluminium').findByOrientation(ori_R, tol_angle);
ebsd_Q = ebsd('Aluminium').findByOrientation(ori_Q, tol_angle);
ebsd_P = ebsd('Aluminium').findByOrientation(ori_P, tol_angle);
ebsd_E = ebsd('Aluminium').findByOrientation(ori_E, tol_angle);
ebsd_F = ebsd('Aluminium').findByOrientation(ori_F, tol_angle);

% Plot the global EBSD map
%figure;
%IPF_map(ebsd, 'Aluminium', RD);
%hold on

% Overlay each orientation component with legend
%IPF_map(ebsd_cube, 'Aluminium', RD); 
%IPF_map(ebsd_brass, 'Aluminium', RD); 
%IPF_map(ebsd_S, 'Aluminium', RD); 
%IPF_map(ebsd_copper, 'Aluminium', RD); 
%IPF_map(ebsd_goss, 'Aluminium', RD); 
%IPF_map(ebsd_cube45, 'Aluminium', RD); 
%IPF_map(ebsd_cube22ND, 'Aluminium', RD); 
%IPF_map(ebsd_cube22RD, 'Aluminium', RD); 
%IPF_map(ebsd_cube22TD, 'Aluminium', RD); 
%IPF_map(ebsd_rotGoss, 'Aluminium', RD); 
%IPF_map(ebsd_R, 'Aluminium', RD); 
%IPF_map(ebsd_Q, 'Aluminium', RD); 
%IPF_map(ebsd_P, 'Aluminium', RD); 
%IPF_map(ebsd_E, 'Aluminium', RD); 
%IPF_map(ebsd_F, 'Aluminium', RD); 

%hold off


%% Plot the grains where highest misorientations are found

ori1 = orientation.byEuler(68.5*degree, 3.6*degree, 123.7*degree, CubicCS);
ori2 = orientation.byEuler(210.3*degree, 35.3*degree, 73.8*degree, CubicCS);
ori3 = orientation.byEuler(253.5*degree, 48*degree, 84.4*degree, CubicCS);

ebsd_ori1 = ebsd('Aluminium').findByOrientation(ori1 ,20*degree)
ebsd_ori2 = ebsd('Aluminium').findByOrientation(ori2 ,20*degree)
ebsd_ori3 = ebsd('Aluminium').findByOrientation(ori3 ,20*degree)


IPF_map(ebsd_ori1, 'Aluminium', vector3d.Z)

IPF_map(ebsd_ori2, 'Aluminium', vector3d.Z)

IPF_map(ebsd_ori3, 'Aluminium', vector3d.Z)

%% From the known texture components, the highest misorientations lies
% on grains close to:

% Define the known texture component orientations
ori_cube = orientation.byMiller([0 0 1], [1 0 0], CubicCS);
ori_brass = orientation.byMiller([1 1 0], [1 1 2], CubicCS);
ori_S = orientation.byEuler(59*degree, 34*degree, 65*degree, CubicCS);
ori_copper = orientation.byMiller([1 1 2], [1 1 1], CubicCS);
ori_goss = orientation.byMiller([0 1 1], [1 0 0], CubicCS);
ori_cube45 = orientation.byMiller([1 0 0], [0 1 1], CubicCS);
ori_cube22ND = orientation.byEuler(22*degree, 0*degree, 0*degree, CubicCS);
ori_cube22RD = orientation.byEuler(0*degree, 22*degree, 0*degree, CubicCS);
ori_cube22TD = orientation.byEuler(90*degree, 22*degree, 0*degree, CubicCS);
ori_rotGoss = orientation.byMiller([0 -1 1], [0 1 1], CubicCS);
ori_R = orientation.byMiller([1 2 4], [2 1 1], CubicCS);
ori_Q = orientation.byMiller([0 1 3], [2 -3 1], CubicCS);
ori_P = orientation.byMiller([0 1 1], [1 2 2], CubicCS);
ori_E = orientation.byMiller([1 1 1], [1 1 0], CubicCS);
ori_F = orientation.byMiller([1 1 1], [1 1 2], CubicCS);

% Define the given orientations
%ori1 = orientation.byEuler(68.5*degree, 3.6*degree, 123.7*degree, CubicCS);
ori1 = orientation.byEuler(325.6*degree, 46.3*degree, 64.1*degree, CubicCS);
%ori2 = orientation.byEuler(201.3*degree, 35.3*degree, 73.8*degree, CubicCS);
ori2 = orientation.byEuler(287.8*degree, 18.8*degree, 11.5*degree, CubicCS);
ori3 = orientation.byEuler(253.5*degree, 48*degree, 84.4*degree, CubicCS);

% Calculate misorientations with each texture component
misori1 = angle(ori1, [ori_cube, ori_brass, ori_S, ori_copper, ori_goss, ...
                      ori_cube45, ori_cube22ND, ori_cube22RD, ori_cube22TD, ...
                      ori_rotGoss, ori_R, ori_Q, ori_P, ori_E, ori_F]);
misori2 = angle(ori2, [ori_cube, ori_brass, ori_S, ori_copper, ori_goss, ...
                      ori_cube45, ori_cube22ND, ori_cube22RD, ori_cube22TD, ...
                      ori_rotGoss, ori_R, ori_Q, ori_P, ori_E, ori_F]);
misori3 = angle(ori3, [ori_cube, ori_brass, ori_S, ori_copper, ori_goss, ...
                      ori_cube45, ori_cube22ND, ori_cube22RD, ori_cube22TD, ...
                      ori_rotGoss, ori_R, ori_Q, ori_P, ori_E, ori_F]);

% Find the minimum misorientation for each orientation
[min_misori1, idx1] = min(misori1);
[min_misori2, idx2] = min(misori2);
[min_misori3, idx3] = min(misori3);


% Define the list of orientation objects along with their names
ori_list = {'ori_cube', 'ori_brass', 'ori_S', 'ori_copper', 'ori_goss', ...
            'ori_cube45', 'ori_cube22ND', 'ori_cube22RD', 'ori_cube22TD', ...
            'ori_rotGoss', 'ori_R', 'ori_Q', 'ori_P', 'ori_E', 'ori_F'};


fprintf('Orientation 1 is closest to: %s\n', (ori_list{idx1}));
fprintf('Orientation 2 is closest to: %s\n', (ori_list{idx2}));
fprintf('Orientation 3 is closest to: %s\n', (ori_list{idx3}));


%% Plot these texture components


%CubicCS = ebsd('Aluminium').CS;
%ebsd_cube = ebsd('Aluminium').findByOrientation(ori_cube, 10*degree);
%ebsd_cube22ND = ebsd('Aluminium').findByOrientation(ori_cube22ND, 15*degree);
%ebsd_goss = ebsd('Aluminium').findByOrientation(ori_goss, 15*degree);
%ebsd_rotGoss = ebsd('Aluminium').findByOrientation(ori_rotGoss, 15*degree);
%ebsd_brass = ebsd('Aluminium').findByOrientation(ori_brass, 15*degree);

%IPF_map(ebsd_cube22ND, 'Aluminium', RD)
%IPF_map(ebsd_goss, 'Aluminium', RD)
%IPF_map(ebsd_rotGoss, 'Aluminium', RD)
%IPF_map(ebsd_brass, 'Aluminium', RD)
%IPF_map(ebsd_cube, 'Aluminium', RD)

%IPF_map(ebsd_cube, 'Aluminium', RD); 
%IPF_map(ebsd_brass, 'Aluminium', RD); 
%IPF_map(ebsd_copper, 'Aluminium', RD); 
%IPF_map(ebsd_goss, 'Aluminium', RD); 
%IPF_map(ebsd_cube22RD, 'Aluminium', RD); 
%IPF_map(ebsd_rotGoss, 'Aluminium', RD); 
%IPF_map(ebsd_R, 'Aluminium', RD); 

%IPF_map(ebsd_cube, 'Aluminium', RD)

