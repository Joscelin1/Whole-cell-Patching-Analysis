%this script averages first 3 responses/ sweeps and saves as a baseline
%value in excel - it then uses this to make rest of peaks a % of baseline
%and saves in seperate excel
%JS. 09.21


clear all

sFolder = 'P:\Patching';
sStrain = 'SHR';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sACh = 'P:\Patching\AChResponse.csv';
cd(sFolder)

data = xlsread(sACh);
cell_nums = unique(data(:,1));
baseline_matrix = zeros(length(cell_nums), 12);
AChHexData = [];
% cell_data = zeros(12,1);

for i = 1:length(cell_nums)
    %get cell number
    this_data = cell_nums(i);
    %get data for this cell
    cell_data = data(data(:,1)==this_data,:);
    %average first 3 sweeps of each cell as baseline 
    baseline_values = mean(cell_data(1:3, :));
    
    baseline_matrix(i, :) = baseline_values(1,:);
    
    %extract the amplitude data 
        baseline_amplitude = baseline_values(1, 8);
    %make respones size a % of baseline
    amp_data = (cell_data(:,8)/baseline_amplitude)*100;
    AChHexData(1, i) = cell_nums(i);
    AChHexData(2:length(amp_data)+1, i) = amp_data(:,1);

end

% figure
% scatter(Mean_SHRData(:,1), AP_height)
% figure
% scatter(Mean_SHRData(:,1), Mean_SHRData(:,12))


baseline_data = num2cell(baseline_matrix);
baseline_data(baseline_matrix == 0) = {[]};

xlswrite([sFolder,'\JS_AChHexData.xlsx'], AChHexData, 'New_Run');
xlswrite([sFolder,'\JS_AChData.xlsx'], baseline_data, 'New_Run');
