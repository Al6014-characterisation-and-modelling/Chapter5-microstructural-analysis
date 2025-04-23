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
%pname = 'BA\SlicesXY';
%matfiles = dir(fullfile(pname, '*.ctf'));
%nfiles = length(matfiles);
%figname = 'BA_3D_';
%delta= 25 % delta 25 for BA and delta 24 for No BA

pname = 'NoBA\SlicesXY';
matfiles = dir(fullfile(pname, '*.ctf'));
nfiles = length(matfiles);
figname = 'NoBA_3D_';
delta= 24 % delta 25 for BA and delta 24 for No BA

%delta= 22 % delta 20 for NoBA and delta 25 for BA


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
    %ebsd=rotate(ebsd,rot,'keepXY'); 
    %ori = ebsd('Aluminium').orientations;

    [grains, ebsd.grainId, ebsd.mis2mean] = calcGrains(ebsd, ...
        'angle',5*degree); 

    psi=calcKernel(grains('Aluminium').meanOrientation);

    ori = ebsd('Aluminium').orientations;
    ori.SS = specimenSymmetry('orthorhombic');

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

%% Calculate the three highest values in the odf and the corresponding orientations
disp('Highest values and their corresponding orientations: ');

[value,orimax] = max(odf_new,'numLocal',3)


%% Calculation of texture components
disp('Total percentage of texture components: ');
%CalclComponents returns the prefered orientation and the percentage of 
% orientations that crawled to each of them.
[component, vol] = calcComponents(odf_new);
component , vol*100 
% this will give me the total percentage of texture components



%% Define texture components to find
cube=orientation('Euler',0*degree,0*degree,0*degree);
S=orientation('Euler',59*degree,34*degree,65*degree);
Q=orientation('Euler',45*degree,15*degree,10*degree);
P=orientation('Euler',65*degree,45*degree,0*degree);
copper=orientation('Euler',90*degree,30*degree,45*degree);
R=orientation('Euler',53*degree,36*degree,60*degree);
cube22ND=orientation('Euler',22*degree,0*degree,0*degree);
cube22RD=orientation('Euler',0*degree,22*degree,0*degree);
cube22TD=orientation('Euler',90*degree,22*degree,0*degree);
cube45ND=orientation('Euler',45*degree,0*degree,0*degree);
E1=orientation('Euler',0*degree,55*degree,45*degree);
E2=orientation('Euler',60*degree,55*degree,45*degree);
F1=orientation('Euler',30*degree,55*degree,45*degree);
F2=orientation('Euler',90*degree,55*degree,45*degree);
D=orientation('Euler',90*degree,27*degree,45*degree);
Brass=orientation('Euler',35*degree,45*degree,0*degree);
Goss=orientation('Euler',0*degree,45*degree,0*degree);

%% Calculate the percentage of each texture component


% Define the names of the variables
variableNames = {'Cube', 'S', 'Q', 'P', 'copper', 'R', 'cube22ND', 'cube22RD', 'cube22TD', 'cube45ND','E1','E2', 'F1', 'F2', 'D', 'Brass', 'Goss'};


% Calculate the percentages of each texture component
pref_ori=[cube, S, Q, P, copper, R, cube22ND, cube22RD, cube22TD, cube45ND, E1, E2, F1, F2, D, Brass, Goss];
comp = volume(odf_new, pref_ori, delta/10 * degree)*100;

% Loop through 'comp' that displays values with their names
for i = 1:length(comp)
    variableName = variableNames{i};
    disp([variableName, ': ', num2str(comp(i), '%.2f')]);
end

% Display the total percentage
disp(['Total percentage: ', num2str(sum(comp), '%.1f')]);
%This total percentage should be similar to total percentage of texture
%components, if not, modify delta slightly


%% Save texture component percentages to CSV file

% Create a table with texture component names and their corresponding percentages
textureComponentsTable = table(variableNames', comp', 'VariableNames', {'TextureComponent', 'Percentage'});

% Write the table to a CSV file
writetable(textureComponentsTable,  [figname 'texture_components.csv']);

%% save data in odf.mat

save([figname 'odf.mat'], 'odf_new');

