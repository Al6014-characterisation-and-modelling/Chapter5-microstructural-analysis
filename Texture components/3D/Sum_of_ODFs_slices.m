%% Specify Crystal and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98])};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

% path to files
pname = 'NoBA\SlicesXY';
matfiles = dir(fullfile(pname, '*.ctf'));
nfiles = length(matfiles);

%% Initialize ODF cell array
odf_cell = cell(1, nfiles);

% Apply the rotation
rot=rotation('Euler', 0*degree, 0*degree, 0*degree); 

RD = vector3d.X;
TD = vector3d.Y;
ND = vector3d.Z;

%% Loop over each file
for fileIdx = 1:nfiles
    fname = fullfile(pname, matfiles(fileIdx).name);
    
    %% Import the Data
    ebsd = EBSD.load(fname,CS,'interface','ctf',...
      'convertEuler2SpatialReferenceFrame');
    ebsd=rotate(ebsd,rot,'keepXY'); 
    ori = ebsd('Aluminium').orientations;

    [grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, ...
        'angle',5*degree); % 10*degree is a conversion radians-degrees

    psi=calcKernel(grains('Aluminium').meanOrientation);

    ori = ebsd('Aluminium').orientations;
    ori.SS = specimenSymmetry('orthorhombic');
    %% Your existing code for calculations
    
    % ... (rest of your script)
    % Calculate ODF for the current file
    odf_current = calcDensity(ori,'kernel',psi);
    
    % Store the ODF result in the cell array
    odf_cell{fileIdx} = odf_current;
end

%% Accumulate and normalize the ODFs
odf_new = zeros(size(odf_cell{1}));
for fileIdx = 1:nfiles
    odf_new = odf_new + odf_cell{fileIdx};
end
odf_new = odf_new / nfiles;


%% Plot odf_new

odf_new.SS = specimenSymmetry('222');

figure;
plot(odf_new,'phi2',[0]*degree,'antipodal','linewidth',1,'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0,6]);
mtexColorbar ('location','eastoutside','title','mrd');
colormap("sky");
ax = gca;
%ax.YDir = 'reverse';
saveas(gcf, fullfile(pname, 'odf_phi2_0.png'));


figure;
plot(odf_new,'phi2',[45]*degree,'antipodal','linewidth',1,'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0,6]);
mtexColorbar ('location','eastoutside','title','mrd');
colormap("sky");
ax = gca;
%ax.YDir = 'reverse';
saveas(gcf, fullfile(pname, 'odf_phi2_45.png'));

figure;
plot(odf_new,'phi2',[65]*degree,'antipodal','linewidth',1,'colorbar','cs','ss','contourf',0:1:6,'colorRange',[0,6]);
mtexColorbar ('location','eastoutside','title','mrd');
colormap("sky");
ax = gca;
%ax.YDir = 'reverse';
saveas(gcf, fullfile(pname, 'odf_phi2_65.png'));


%% Save the 'odf_new' variable to the .mat file
fname = fullfile(pname, 'odf.mat');
save(fname);


%% Plot pole figures
figure;
cs = crystalSymmetry('m-3m');
plotPDF(odf_new,Miller(1,1,1,cs),'antipodal','contourf',0:0.15:3,'minmax','colorRange',[0,3],'colorbar','on','complete','upper')
view(90, -90);
mtexColorbar;
mtexColorMap sky
title('1 1 1', 'Position', [0, 0], 'HorizontalAlignment', 'center');
% Add a title annotation at the top
annotation('textbox', [0.05, 0.85, 0.20, 0.14], 'String', '', ...
    'HorizontalAlignment', 'center', 'FontSize', 20, 'EdgeColor', 'none');

saveas(gcf, fullfile(pname, 'pole_111.png'));