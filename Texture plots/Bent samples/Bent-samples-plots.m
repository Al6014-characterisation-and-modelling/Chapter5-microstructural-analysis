%% Specify Crystal and Specimen Symmetries
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98])};

% Plotting convention
setMTEXpref('xAxisDirection', 'east');
setMTEXpref('zAxisDirection', 'outOfPlane');

%% File Names and Paths
pname1 = 'C:\Users\GlassDesktop\Desktop\Bent samples\IA alloy\T-bending\with 15% pre-stretching\A_RD(xaxis)_ND(yaxis)';
fname1 = [pname1 '\A Data.ctf'];

pname2 = 'C:\Users\GlassDesktop\Desktop\Bent samples\No IA alloy\T-bending\With 15% pre-stretching\B_RD(xaxis)_ND(Yaxis)';
fname2 = [pname2 '\B_texture_highres.ctf'];

pname3 = 'C:\Users\GlassDesktop\Desktop\Bent samples\IA alloy\L-bending\With 15% pre-stretching\C_TD(xaxis)_ND(yaxis)';
fname3 = [pname3 '\C_texture Data.ctf'];

pname4 = 'C:\Users\GlassDesktop\Desktop\Bent samples\No IA alloy\L-bending\With 15% pre-stretching\D_TD(xaxis)_ND(yaxis)';
fname4 = [pname4 '\D_texture Data.ctf'];


pname5 = 'C:\Users\GlassDesktop\Desktop\Bent samples\IA alloy\L-bending\Without pre-stretching\G_sample_TD(xaxis)_ND(yaxis)\G_sample_TD_ND\ctfs';
fname5 = [pname5 '\G_sample Data.ctf'];

%sample H
pname6 = 'C:\Users\GlassDesktop\Desktop\Bent samples\No IA alloy\L-bending\Without pre-stretching\D_TD(xaxis)_ND(yaxis)';
fname6 = [pname6 '\EBSD Map Data 1 - EBSD Data Map_002.ctf'];

%% Import the Data

% create an EBSD variable containing the data
ebsd = EBSD.load(fname6,CS,'interface','ctf',...
    'convertSpatial2EulerReferenceFrame');

%rot=rotation('Euler', 0*degree, 90*degree, 0*degree); %This rotation for NDTD
rot=rotation('Euler', 0*degree, 0*degree, 0*degree); %This rotation for NDRD                                              % 
                                                      
ebsd=rotate(ebsd,rot,'keepXY'); 

RD = vector3d.X;
ND = vector3d.Y;
TD = vector3d.Z;


ori = ebsd('Aluminium').orientations;

IPF_map(ebsd, 'Aluminium', TD)


%% OPTIONAL STEP - Marking a small are to analyse 
plot(ebsd)

%region = [200 -500 1100 500]; %[xmin ymin xmax-xmin, ymax-ymin] B sample
%region = [0 -350 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] B sample small

%region = [50 -500 1100 500]; %[xmin ymin xmax-xmin, ymax-ymin] D sample
%region = [80 -350 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] D sample small
%region = [390 -500 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] D sample small2


%region = [0 -270 690 270]; %[xmin ymin xmax-xmin, ymax-ymin] C sample small
%region = [0 -350 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] C sample small2


%region = [350 -400 1180 280]; %[xmin ymin xmax-xmin, ymax-ymin] G sample
%region = [450 -470 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] G sample small

%region = [750 -470 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] G sample small2

%region = [150 -300 1180 280]; %[xmin ymin xmax-xmin, ymax-ymin] H sample
%region = [250 -350 900 350]; %[xmin ymin xmax-xmin, ymax-ymin] H sample small
region = [480 -380 500 350]; %[xmin ymin xmax-xmin, ymax-ymin] H sample small2


%region = [250 -760 1350 660]; %[xmin ymin xmax-xmin, ymax-ymin] G sample - 125microns

%rectangle('position',region,'edgecolor','r','linewidth',2)

rectangle('position',region,'edgecolor','k','linewidth',1);
cropped_ebsd = ebsd(inpolygon(ebsd,region));
%IPF_map(cropped_ebsd, 'Aluminium', ND)
IPF_map(cropped_ebsd, 'Aluminium', TD)

%ebsd=cropped_ebsd

%% Plotting KAM (Kernel Average Misorientation)

ebsd=cropped_ebsd

image;
kam = ebsd.KAM / degree;
plot(ebsd,kam,'micronbar','off')
setColorRange([0,5])
mtexColorbar
mtexColorMap LaboTeX
%hold on
%plot(grains.boundary,'lineWidth',1.5)
%hold off

%text('Parent',Parent1,'VerticalAlignment','baseline',...
%    'HorizontalAlignment','center',...
%    'FontSize',18,...
%    'Interpreter','latex',...
%    'String','\rm{\textbf{150 $\mu$m}}',...
%    'Position',[122.352773826458 -831.207681365576 0],...
%    'Color',[1 1 1]);


%%
%image;
%IPF_map(ebsd, 'Aluminium', TD);

%IPF_map(ebsd, 'Aluminium', RD);
%mtexColorbar('title', 'IPF Map - Aluminium');
%mtexColorMap LaboTeX;

%G sample
%hold on;
%createhgtransform(gca, [9 -852; 227 -852; 227 -808; 9 -808], [1 2 3 4], [36, -845; 195, -845; 195, -838; 36, -838]);
%hold off;


%H sample
%hold on;
%createhgtransform(gca, [9 -652; 227 -652; 227 -608; 9 -608], [1 2 3 4], [36, -645; 195, -645; 195, -638; 36, -638]);
%hold off;

IPF_map(ebsd, 'Aluminium', TD)
