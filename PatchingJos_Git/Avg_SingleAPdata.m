

sFolder = 'P:\Patching';
sStrain = 'SHR';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sActionPotentials = 'P:\Patching\ActionPotentials.csv';
cd(sFolder)

dStrainID = GetStrainID(sStrain);

SHRdata = xlsread(sActionPotentials);
%making a matrix for averaged sweeps (5x smaller rows)
Mean_SHRData = zeros(size(SHRdata, 1)/5, size(SHRdata, 2));
row_counter = 1;
for i = 1:5:length(SHRdata) 
    
   M = SHRdata(i:i+4, :);
   M = mean(M);
   M(1, 1:4) = SHRdata(i, 1:4);
   Mean_SHRData(row_counter,:) = M(1,:);
   %increment row so that next time data is inserted into the next
   %row down
   row_counter = row_counter + 1;

end


% figure
% scatter(Mean_SHRData(:,1), AP_height)
% figure
% scatter(Mean_SHRData(:,1), Mean_SHRData(:,12))
dlmwrite([sFolder,'\JS_MeanAPdata.csv'],Mean_SHRData,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');

