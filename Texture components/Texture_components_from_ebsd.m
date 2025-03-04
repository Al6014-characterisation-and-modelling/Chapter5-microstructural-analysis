%% CALCULATIONS OF VOLUME FRACTIONS FROM EBSD (better than doing it with ODF)

% Recommendation: it is better to calculated the volume fractions from

% imported orientations rather than doing it with the orientation

% distribution function (ODFs).

% Definition of common orientations

% definition of orientations through Euler angles as ori = orientation.byEuler(30*degree,50*degree,10*degree,cs)

% definition of orientations through Miller indices by the crystal directions point towards the specimen directions Z and X

% orim = orientation.byMiller([0 1 1],[1 0 0],cs)

%If you did not crop the microstructure you can leave ebsd as variable. Otherwise if you cropped it, you should use ebsdc.

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
pname = 'C:\Users\Laura\Dropbox (The University of Manchester)\4th year - Crystal Plasticity and Microstructure Modelling for Rolled Aluminium Sheet\Experimental\3D-EBSD\BA\SlicesXY';
matfiles = dir(fullfile(pname, '*.ctf'));
nfiles = length(matfiles);

%% Initialize ODF cell array
Vfractot_cell1 = cell(1, nfiles);
Vfractot_cell2 = cell(1, nfiles);

RD = vector3d.X;
TD = vector3d.Y;
ND = vector3d.Z;


%% Components

Bs=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([0 1 1],[2 1 1],CS)));

S=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byEuler(59*degree,34*degree,65*degree,CS)));

Co=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byEuler(90*degree,30*degree,45*degree,CS)));

Goss=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([0 1 1],[1 0 0],CS)));

Cube45=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') *orientation.byMiller([1 0 0],[0 1 1],CS)));

Cube=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([0 0 1],[1 0 0],CS)));

R=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([1 2 4],[2 1 1],CS)));

Cube22ND=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byEuler(22*degree,0*degree,0*degree,CS)));

Cube22RD=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byEuler(0*degree,22*degree,0*degree,CS)));

Cube22TD=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byEuler(90*degree,22*degree,0*degree,CS)));

rotGoss=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([0 -1 1],[0 1 1],CS)));

Q=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') *orientation.byMiller([0 1 3],[2 -3 1],CS)));

P=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([0 1 1],[1 2 2],CS)));

E=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([1 1 1],[1 1 0],CS)));

F=unique(specimenSymmetry('1')*(specimenSymmetry('mmm') * orientation.byMiller([1 1 1],[1 1 2],CS)));

TXTtolAngle=10*degree;

TXTtolAngle2=10*degree;

TXTtolAngle3=5*degree;

%% Loop over each file
for fileIdx = 1:nfiles
    fname = fullfile(pname, matfiles(fileIdx).name);
    
    %% Import the Data
    ebsd = EBSD.load(fname,CS,'interface','ctf',...
      'convertEuler2SpatialReferenceFrame');

    volBs = volume(ebsd.orientations, Bs(1), TXTtolAngle)*100+volume(ebsd.orientations, Bs(2), TXTtolAngle)*100;

    volS = volume(ebsd.orientations, S(1), TXTtolAngle)*100+volume(ebsd.orientations, S(2), TXTtolAngle)*100+volume(ebsd.orientations, S(3), TXTtolAngle)*100+volume(ebsd.orientations, S(4), TXTtolAngle)*100;

    volCo = volume(ebsd.orientations, Co(1), TXTtolAngle)*100+volume(ebsd.orientations, Co(2), TXTtolAngle)*100;

    volGoss = volume(ebsd.orientations, Goss(1), TXTtolAngle)*100;

    volCube = volume(ebsd.orientations, Cube(1), TXTtolAngle)*100;

    volCube45 = volume(ebsd.orientations, Cube45(1), TXTtolAngle)*100;

    volinvGoss = volume(ebsd.orientations, rotGoss(1), TXTtolAngle)*100;

    volQ = volume(ebsd.orientations, Q(1), TXTtolAngle)*100+volume(ebsd.orientations, Q(2), TXTtolAngle)*100+volume(ebsd.orientations, Q(3), TXTtolAngle)*100+volume(ebsd.orientations, Q(4), TXTtolAngle)*100;

    volP = volume(ebsd.orientations, P(1), TXTtolAngle)*100+volume(ebsd.orientations, P(2), TXTtolAngle)*100;

    volR = volume(ebsd.orientations, R(1), TXTtolAngle3)*100+volume(ebsd.orientations, R(2), TXTtolAngle3)*100+volume(ebsd.orientations, R(3), TXTtolAngle3)*100+volume(ebsd.orientations, R(4), TXTtolAngle3)*100;

    VolCube22ND = volume(ebsd.orientations, Cube22ND(1), TXTtolAngle)*100+volume(ebsd.orientations, Cube22ND(2), TXTtolAngle)*100;

    VolCube22RD = volume(ebsd.orientations, Cube22RD(1), TXTtolAngle)*100+volume(ebsd.orientations, Cube22RD(2), TXTtolAngle)*100;

    VolCube22TD = volume(ebsd.orientations, Cube22TD(1), TXTtolAngle)*100+volume(ebsd.orientations, Cube22TD(2), TXTtolAngle)*100;

    volE = volume(ebsd.orientations, E(1), TXTtolAngle)*100+volume(ebsd.orientations, E(2), TXTtolAngle)*100;

    volF = volume(ebsd.orientations, F(1), TXTtolAngle)*100+volume(ebsd.orientations, F(2), TXTtolAngle)*100;


    % Calculate components for the current file

    TextureComps1={'Cube','Cube 45','Cube22ND', 'Cube22RD', 'Cube22TD', 'tot Cube'};
    TextureComps2={'Brass','Copper','S', 'invGoss', 'Q', 'R', 'volE', 'volF'};

    Vfractot_current1=[volCube volCube45 VolCube22ND VolCube22RD VolCube22TD volCube+volCube45+VolCube22ND+VolCube22RD+VolCube22TD];
    Vfractot_current2=[volBs volCo volS volinvGoss volQ volP volE volF];

    % Store the components result in the cell array
    Vfractot_cell1{fileIdx} = Vfractot_current1;
    Vfractot_cell2{fileIdx} = Vfractot_current2;
end


%% Accumulate and normalize the components percentages
components_new1 = zeros(size(Vfractot_cell1{1}));
for fileIdx = 1:nfiles
    components_new1 = components_new1 + Vfractot_cell1{fileIdx};
end
components_new1 = components_new1 / nfiles

components_new2 = zeros(size(Vfractot_cell2{1}));
for fileIdx = 1:nfiles
    components_new2 = components_new2 + Vfractot_cell2{fileIdx};
end
components_new2 = components_new2 / nfiles

%%
figure;

bar(components_new1)

text(1:length(components_new1),components_new1,num2str(components_new1','%.2f'),'vert','bottom','horiz','center');

box off

title('No BA alloy')

ylabel('Volume fraction (%)')

axis([0 8 0 18])

set(gca, 'XTickLabel', TextureComps1)

saveas(figure(7),[pname '15deg.jpg']);

close;

%%
figure;

bar(components_new2)

text(1:length(components_new2),components_new2,num2str(components_new2','%.2f'),'vert','bottom','horiz','center');

box off

title('No BA alloy')

ylabel('Volume fraction (%)')

axis([0 8 0 10])

set(gca, 'XTickLabel', TextureComps2)

saveas(figure(7),[pname '15deg.jpg']);

close;

