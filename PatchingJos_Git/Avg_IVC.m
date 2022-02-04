clear all

sFolder = 'P:\Patching';

sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sIVC = 'P:\Patching\AChIVC.csv';
cd(sFolder)

data = xlsread(sIVC);
cell_nums = unique(data(:,1));
VoltageData = zeros(8, length(cell_nums));
cell_data = zeros(7,1);

for i = 1:length(cell_nums)
    %get cell number
    this_data = cell_nums(i);
    %get data for this cell
    
    
    cell_data = data(data(:,1)==this_data,5);
    %if has less sweeps still insert 
   
 
    VoltageData(1, i) = cell_nums(i);
VoltageData(2:length(cell_data)+1, i) = cell_data(:,1);

end



% figure
% scatter(Mean_SHRData(:,1), AP_height)
% figure
% scatter(Mean_SHRData(:,1), Mean_SHRData(:,12))
dlmwrite([sFolder,'\JS_IVC.csv'],VoltageData,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');

