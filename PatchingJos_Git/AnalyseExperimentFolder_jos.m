%this script analyses the files in a specified directory to populate the
%Cells.csv, MemTest.csv, ActionPotentials.csv files

%change the folder name, strain and animal register string 
%rename the first memtest files for each cell such that they have an '_01'
%at the end of their name and all subsequent memtest files if you want to 
%save the data from them. 

close all;
clear all;
%depends on file size to automatically detect AP, MemTest and RMP files
sFolder = 'P:\Patching\20211116';
sStrain = 'SHR';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sMemTestFile = 'P:\Patching\MemTestData.csv';
oCellData = GetCSVData('P:\Patching\Cells.csv');
oAPData = GetCSVData('P:\Patching\ActionPotentials.csv');
oAChSignal = GetCSVData('P:\Patching\AChResponse.csv');

% get new cells
%    GetNewCells(sFolder,sStrain,sAnimalRegister,oCellData,'P:\Patching');
% get mem test data
%   SaveMemTestData(sFolder, sMemTestFile, oCellData)
% % get new APs
%    GetNewAPs(sFolder, oCellData, oAPData);
%     GetNewAPsFixBridgeBalance(sFolder, oCellData, oAPData);
% %get Family data
%    GetICFamily(sFolder, oCellData);
%   GetICStimFamily(sFolder, oCellData);
% %need to look at and fix this rheobase code 
%  GetRheobase(sFolder, oCellData);
  GetAChSignal(sFolder, oCellData);
% GetAChSignal_Bt(sFolder, oCellData);
%adjust windows to search for peak if aqusition settings foe IVC have
%changed
% GetAChIVCSignal(sFolder);
% GetGlutSignal(sFolder, oCellData);



