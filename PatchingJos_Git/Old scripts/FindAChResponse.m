close all;
clear all;
%depends on file size to automatically detect AP, MemTest and RMP files
sFolder = 'P:\Patching\Test\20210629\AChPuff';
sStrain = 'SHR';
sAnimalRegister = 'P:\Patching\002090_2021_AnimalRegister.xlsx';
sFileName = '2021_06_29_C01_0010.abf';
cd(sFolder)

%get the current animal ID
[aAnimalsnum,aAnimalstxt,aAnimalsraw] = xlsread(sAnimalRegister,sStrain);
sExperiment = [strtrim(sprintf('%2.0f',str2double(sFolder(end-1:end)))),'/',sFolder(end-3:end-2),'/',sFolder(end-7:end-4)];
switch (sStrain)
    case {'Wistar','NZWhite','SHR','WKY'}
        sColumn = 'doe';
    otherwise
        sColumn = 'doi';
end
index = strcmp(aAnimalsraw(:,strcmpi({aAnimalsraw{1,:}},sColumn)),sExperiment);
%delete first row as this corresponds to header
index(1)=[];
dAnimalID = aAnimalsnum(index,strcmpi({aAnimalsraw{1,:}},'animal_id'));

% % %get file list in directory
% % oFileList = GetFileList(sFolder);
% % %select the files that are the right size to be memtest data
% % aFileSize = [oFileList(:).bytes];
% % aAChPuffFiles = find(round(aFileSize,-3)>=36000000 );
close all
[aData,iSamplingInterval,oHeader] = abfload(sFileName);
Data = aData(:,1);
aTime = (1:1:size(Data))*iSamplingInterval/1000; %sampling interval in us


wo = 1/(iSamplingInterval/(10^6))/2; %sampling interval in us converted to sec 1/ for Hz
[z,p,k] = butter(3, 2000/wo, 'low'); % 3rd order filter with a cut off of 2000Hz for data sampled at "wo" Hz
[sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
oFilter = dfilt.df2sos(sos,g);

[z,p,k] = butter(3, 0.5/wo, 'low'); % 3rd order filter with a cut off of 0.5Hz for data sampled at "wo" Hz
[sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
oFilterTwo = dfilt.df2sos(sos,g);
%jj is file number might wan to use later 
jj = 1;
ii = 1;
iCellID = 1;
iStrainID = 1;
%loop through the sweeps in the file (kk)
aDataToWrite{ii}{jj} = cell(size(aData,3),12);
for kk = 1:size(aData,3)
    
            aDataToWrite{ii}{jj}{kk,1} = iCellID;
            aDataToWrite{ii}{jj}{kk,2} = iStrainID;
            aDataToWrite{ii}{jj}{kk,3} = jj; %recording number
            aDataToWrite{ii}{jj}{kk,4} = kk; %sweep number
            aDataToWrite{ii}{jj}{kk,5} = mean(aData(1:11000,1,kk)); %baseline (0sec to 11sec)
%find peak

            [aDataToWrite{ii}{jj}{kk,6}, iMin] = min(aData(:,1,kk)); %Maxpeak(so min of data)


% % %             [aDataToWrite{ii}{jj}{kk,7}, iMin] = min(aData(:,1,kk)); %APMin - dont need? 
% % %             %             aDataToWrite{ii}{jj}{kk,9} =
% % %             %             round(mean(aData(maxind+100:minind-100,2,kk))
% % %             %             - mean(aData(1:maxind-100,2,kk)),-2);
% % %             %             %stimcurrent %this was in column 6 before - dont
% % %             %             need for now 
% % %             %%find stim current
% % %             aDiffStim = diff(aData(:,2,kk));
% % %             [maxval maxind] = max(aDiffStim);
% % %             [minval minind] = min(aDiffStim);
% % %             if maxval < 1
% % %                 maxind = 25 * 1000 / iSamplingInterval;
% % %                 minind = 35 * 1000 / iSamplingInterval;
% % %             end
            
 Fig =    figure;
 axes1 = axes();
plot(aTime, aData(:,1,kk), 'b','parent',axes1)
hold on;
    FilteredData = filter(oFilter,padarray(aData(:,1,kk),200,'replicate','both'));
    FilteredData = flipud(FilteredData);
    FilteredData = filter(oFilter,FilteredData);
    FilteredData = flipud(FilteredData);
    FilteredData = FilteredData(201:end-200);
    plot(aTime, FilteredData, 'c')
    
    
    FilteredData = filter(oFilterTwo,padarray(FilteredData(:,1,kk),200,'replicate','both'));
    FilteredData = flipud(FilteredData);
    FilteredData = filter(oFilterTwo,FilteredData);
    FilteredData = flipud(FilteredData);
    FilteredData = FilteredData(201:end-200);


    plot(aTime, FilteredData, 'r', 'LineWidth',2,'parent',axes1)
    plot(aTime, aHalf,'parent',axes1)
    
    hold off;
    
     %compute threshold voltage
            [maxval maxAPind] = max(diff(diff(FilteredData)));
            aDataToWrite{ii}{jj}{kk,10} = aData(maxAPind,1,kk); %threshold            
             %compute half width
            %             dHalfHeight = (aDataToWrite{ii}{jj}{kk,7} - aDataToWrite{ii}{jj}{kk,8})/2;
            aHalf = aData(:,1,kk) > aDataToWrite{ii}{jj}{kk,10}; %(aDataToWrite{ii}{jj}{kk,7}-dHalfHeight);
            aHalf = find(aHalf);
            aDataToWrite{ii}{jj}{kk,7} = aTime(aHalf(end)) - aTime(aHalf(1)); %half-width
                        
            %compute AHP decay
            dDecayHalfHeight = (aDataToWrite{ii}{jj}{kk,5} - aDataToWrite{ii}{jj}{kk,8})/2;
            aDecayHalf = aData(iMin:end,1,kk) > (aDataToWrite{ii}{jj}{kk,5}-dDecayHalfHeight);
            aDecayHalf = find(aDecayHalf);
            aDataToWrite{ii}{jj}{kk,8} = aTime(aDecayHalf(1)+iMin-1) - aTime(iMin); %decay
            aDataToWrite{ii}{jj}{kk,9} = mean(aData(1:maxind-100,2,kk)); %holdingcurrent
            
            
end
