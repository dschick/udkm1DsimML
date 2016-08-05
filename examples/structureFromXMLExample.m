%% Structure From XML Example
% In this example we show how to build a sample structure from a given XML
% file. The use of XML files allows to easily exchange sample structures
% and have a clean code for the upcoming simulations.
%
% Be sure to include all folders of the udkm1dsim toolbox to your MATLAB 
% search path. Especially the _/parameters_ and _/helpers_ folders with all 
% their subfolders.

%%
% here we clear the workspace and close all figures
clear all;
close all;
% these are some units and constants we can use later on
u            = units;
const        = constants;

%% Content of Structure XML File
% First of all we display the content of the XML file. You should read the
% comments in the file carefully. Generally, building structures in a
% MATLAB script or from an XML file is not much different.
type('sample.xml');
%% Load the Sample Structure
% The name of the structure is overwritten from the XML file. So the only
% necessary parameter is the path to the XML file we want to load.
S = structure('', './sample.xml');
%%
% display the properties of the structure and visualize
S.disp();
%%
S.visualize();

%% More Structure Methods
% In order to access the different unit cells in the structure we can first
% display the internal nummeration of the unit cells in the structure:
disp(S.getUniqueUnitCells);

%%
% Now we can also retrieve all positions per unique unit cell in the
% structure by simply asking for the cell array:
Pos = S.getAllPositionsPerUniqueUnitCell();

% Each cell of this array is linked to the above retrieved unique unit
% cells. Accordingly the first entry gives all positions of the PZT unit
% cell in the structure:
disp(Pos{1});

%%
% You can easily access the properties of a unit cell by its index by 
% simply typing:
disp(S.getUnitCellHandle(10).name);
disp(['mass: ' num2str(S.getUnitCellHandle(10).mass/u.kg) ' kg']);