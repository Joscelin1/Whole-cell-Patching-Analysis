clear all

sFolder = 'P:\Patching';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sICFamily = 'P:\Patching\ICStimFamily.csv';
cd(sFolder)

% dStrainID = GetStrainID(sStrain);

data = xlsread(sICFamily);
cell_nums = unique(data(:,1));
data_to_save = zeros(9, length(cell_nums));
cell_data = zeros(12,1);

for i = 1:length(cell_nums)
    %get cell number
    this_data = cell_nums(i);
    %get data for this cell
    
    
    cell_data = data(data(:,1)==this_data,4);
    %if has less sweeps still insert 
   
 
    data_to_save(1, i) = cell_nums(i);
data_to_save(2:length(cell_data)+1, i) = cell_data(:,1);

end

% figure
% scatter(Mean_SHRData(:,1), AP_height)
% figure
% scatter(Mean_SHRData(:,1), Mean_SHRData(:,12))


dlmwrite([sFolder,'\JS_ICStimFamilyData.csv'],data_to_save,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');

