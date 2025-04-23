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
% 3D-EBSD BA
pname = 'BA\slices';
matfiles = dir(fullfile(pname, '*.ctf'));
nfiles = length(matfiles);
figname = 'BA_3D_';

% 3D-EBSD No BA
%pname = 'NoBA\slices';
%matfiles = dir(fullfile(pname, '*.ctf'));
%nfiles = length(matfiles);
%figname = 'NoBA_3D_';

%% Initialize variables for texture components
totalBs = 0;
totalS = 0;
totalCo = 0;
totalGoss = 0;
totalCube = 0;
totalCube45 = 0;
totalinvGoss = 0;
totalQ = 0;
totalP = 0;
totalR = 0;
totalCube22ND = 0;
totalCube22RD = 0;
totalCube22TD = 0;
totalE = 0;
totalF = 0;

%% Define texture components

% % % % definition of symmetric equivalent

Bs = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([0 1 1], [2 1 1], CS)));

S = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byEuler(59*degree, 34*degree, 65*degree, CS)));

Co = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byEuler(90*degree, 30*degree, 45*degree, CS)));

Goss = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([0 1 1], [1 0 0], CS)));

Cube45 = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([1 0 0], [0 1 1], CS)));

Cube = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([0 0 1], [1 0 0], CS)));

R = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([1 2 4], [2 1 1], CS)));

Cube22ND = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byEuler(22*degree, 0*degree, 0*degree, CS)));

Cube22RD = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byEuler(0*degree, 22*degree, 0*degree, CS)));

Cube22TD = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byEuler(90*degree, 22*degree, 0*degree, CS)));

rotGoss = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([0 -1 1], [0 1 1], CS)));

Q = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([0 1 3], [2 -3 1], CS)));

P = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([0 1 1], [1 2 2], CS)));

E = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([1 1 1], [1 1 0], CS)));

F = unique(specimenSymmetry('1') * (specimenSymmetry('mmm') * orientation.byMiller([1 1 1], [1 1 2], CS)));

%% Calculation of texture components

TXTtolAngle = 10*degree;
TXTtolAngle3 = 5*degree;

%% Loop over each file

for fileIdx = 1:nfiles
    fname = fullfile(pname, matfiles(fileIdx).name);
    
    %% Import the Data
    ebsd = EBSD.load(fname, CS, 'interface', 'ctf', 'convertEuler2SpatialReferenceFrame');

    %% Calculate volumes for each texture component and accumulate
    totalBs = totalBs + volume(ebsd.orientations, Bs(1), TXTtolAngle) * 100 + volume(ebsd.orientations, Bs(2), TXTtolAngle) * 100;

    totalS = totalS + volume(ebsd.orientations, S(1), TXTtolAngle) * 100 + volume(ebsd.orientations, S(2), TXTtolAngle) * 100 + volume(ebsd.orientations, S(3), TXTtolAngle) * 100 + volume(ebsd.orientations, S(4), TXTtolAngle) * 100;

    totalCo = totalCo + volume(ebsd.orientations, Co(1), TXTtolAngle) * 100 + volume(ebsd.orientations, Co(2), TXTtolAngle) * 100;

    totalGoss = totalGoss + volume(ebsd.orientations, Goss(1), TXTtolAngle) * 100;

    totalCube = totalCube + volume(ebsd.orientations, Cube(1), TXTtolAngle) * 100;

    totalCube45 = totalCube45 + volume(ebsd.orientations, Cube45(1), TXTtolAngle) * 100;

    totalinvGoss = totalinvGoss + volume(ebsd.orientations, rotGoss(1), TXTtolAngle) * 100;

    totalQ = totalQ + volume(ebsd.orientations, Q(1), TXTtolAngle) * 100 + volume(ebsd.orientations, Q(2), TXTtolAngle) * 100 + volume(ebsd.orientations, Q(3), TXTtolAngle) * 100 + volume(ebsd.orientations, Q(4), TXTtolAngle) * 100;

    totalP = totalP + volume(ebsd.orientations, P(1), TXTtolAngle) * 100 + volume(ebsd.orientations, P(2), TXTtolAngle) * 100;

    totalR = totalR + volume(ebsd.orientations, R(1), TXTtolAngle3) * 100 + volume(ebsd.orientations, R(2), TXTtolAngle3) * 100 + volume(ebsd.orientations, R(3), TXTtolAngle3) * 100 + volume(ebsd.orientations, R(4), TXTtolAngle3) * 100;

    totalCube22ND = totalCube22ND + volume(ebsd.orientations, Cube22ND(1), TXTtolAngle) * 100 + volume(ebsd.orientations, Cube22ND(2), TXTtolAngle) * 100;

    totalCube22RD = totalCube22RD + volume(ebsd.orientations, Cube22RD(1), TXTtolAngle) * 100 + volume(ebsd.orientations, Cube22RD(2), TXTtolAngle) * 100;

    totalCube22TD = totalCube22TD + volume(ebsd.orientations, Cube22TD(1), TXTtolAngle) * 100 + volume(ebsd.orientations, Cube22TD(2), TXTtolAngle) * 100;

    totalE = totalE + volume(ebsd.orientations, E(1), TXTtolAngle) * 100 + volume(ebsd.orientations, E(2), TXTtolAngle) * 100;

    totalF = totalF + volume(ebsd.orientations, F(1), TXTtolAngle) * 100 + volume(ebsd.orientations, F(2), TXTtolAngle) * 100;
end

%% Calculate the average percentages
avgBs = totalBs / nfiles;
avgS = totalS / nfiles;
avgCo = totalCo / nfiles;
avgGoss = totalGoss / nfiles;
avgCube = totalCube / nfiles;
avgCube45 = totalCube45 / nfiles;
avginvGoss = totalinvGoss / nfiles;
avgQ = totalQ / nfiles;
avgP = totalP / nfiles;
avgR = totalR / nfiles;
avgCube22ND = totalCube22ND / nfiles;
avgCube22RD = totalCube22RD / nfiles;
avgCube22TD = totalCube22TD / nfiles;
avgE = totalE / nfiles;
avgF = totalF / nfiles;

%% Display the average percentages
disp(['Average Bs: ', num2str(avgBs), '%']);
disp(['Average S: ', num2str(avgS), '%']);
disp(['Average Co: ', num2str(avgCo), '%']);
disp(['Average Goss: ', num2str(avgGoss), '%']);
disp(['Average Cube: ', num2str(avgCube), '%']);
disp(['Average Cube45: ', num2str(avgCube45), '%']);
disp(['Average invGoss: ', num2str(avginvGoss), '%']);
disp(['Average Q: ', num2str(avgQ), '%']);
disp(['Average P: ', num2str(avgP), '%']);
disp(['Average R: ', num2str(avgR), '%']);
disp(['Average Cube22ND: ', num2str(avgCube22ND), '%']);
disp(['Average Cube22RD: ', num2str(avgCube22RD), '%']);
disp(['Average Cube22TD: ', num2str(avgCube22TD), '%']);
disp(['Average E: ', num2str(avgE), '%']);
disp(['Average F: ', num2str(avgF), '%']);

%% Plot the percentages

figure;

TextureComps={'Cube','Cube 45','Goss','Brass','Copper','S'}

Vfractot=[avgCube avgCube45 avgGoss avgBs avgCo avgS]

bar(Vfractot)

text(1:length(Vfractot),Vfractot,num2str(Vfractot','%.2f'),'vert','bottom','horiz','center');

box off

title('Rolling Components')

ylabel('Volume fraction (%)')

axis([0 8 0 18])

set(gca, 'XTickLabel', TextureComps)


figure;

TextureComps={'invGoss','Q','P', 'R','Cube22ND', 'Cube22RD', 'Cube22TD', 'volE', 'volF'}

Vfractot2=[avginvGoss avgQ avgP avgR avgCube22ND avgCube22RD avgCube22TD avgE avgF]

bar(Vfractot2)

text(1:length(Vfractot2),Vfractot2,num2str(Vfractot2','%.2f'),'vert','bottom','horiz','center');

box off

title('Recrystallization Components')

ylabel('Volume fraction (%)')

axis([0 12 0 18])

set(gca, 'XTickLabel', TextureComps)
