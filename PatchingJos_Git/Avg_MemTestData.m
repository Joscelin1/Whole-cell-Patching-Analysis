clear all

sFolder = 'P:\Patching';
sStrain = 'SHR';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sMemTestData = 'P:\Patching\MemTestData.csv';
cd(sFolder)

dStrainID = GetStrainID(sStrain);

Data = xlsread(sMemTestData);
%making a matrix for averaged sweeps (5x smaller rows)
mean_data = zeros(size(Data, 1), size(Data, 2));
cell_nums = unique(Data(:,1));
EarlyMemData = zeros(length(cell_nums), size(Data,2));
to_insert = zeros(1,  size(Data,2));
row_counter = 1;

for i = 1:length(cell_nums)
    %get cell number
    this_data = cell_nums(i);
    %get data for this cell
    cell_data = Data(Data(:,1)==this_data,:);
    %select first two mem tests and average - if onlt one has been take
    %still insert 
    if size(cell_data, 1) >1
        to_insert = mean(cell_data(1:2, :));
    else
        to_insert = cell_data(1, :);
    end
EarlyMemData(row_counter,:) = to_insert(1,:);
%    M = SHRdata(i:i+4, :);
%    M = mean(M);
%    M(1, 1:4) = SHRdata(i, 1:4);
%    mean_data(row_counter,:) = M(1,:);
   %increment row so that next time data is inserted into the next
   %row down
   row_counter = row_counter + 1;

end


% figure
% scatter(mean_data(:,1), AP_height)
% figure
% scatter(mean_data(:,1), mean_data(:,12))
dlmwrite([sFolder,'\JS_MeanMemTestData.csv'],EarlyMemData,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');

