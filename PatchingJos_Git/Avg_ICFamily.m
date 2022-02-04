clear all

sFolder = 'P:\Patching';
sStrain = 'SHR';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sICFamily = 'P:\Patching\ICFamily.csv';
cd(sFolder)

dStrainID = GetStrainID(sStrain);

data = xlsread(sICFamily);
cell_nums = unique(data(:,1));
VoltageData = zeros(13, length(cell_nums));
cell_data = zeros(12,1);

for i = 1:length(cell_nums)
    %get cell number
    this_data = cell_nums(i);
    %get data for this cell
    
    %change this line to avg diff parts f the signal 
    cell_data = data(data(:,1)==this_data,4);
    %if has less sweeps still insert 
   
 
    VoltageData(1, i) = cell_nums(i);
VoltageData(2:length(cell_data)+1, i) = cell_data(:,1);

end



% figure
% scatter(Mean_SHRData(:,1), AP_height)
% figure
% scatter(Mean_SHRData(:,1), Mean_SHRData(:,12))
dlmwrite([sFolder,'\JS_ICFamilyData.csv'],VoltageData,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');

