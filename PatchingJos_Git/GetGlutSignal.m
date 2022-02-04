function oGlutSignal = GetGlutSignal(sFolder,oCellData)
%this script analyses glut data frime glut files within experiment folders 
%JS 10.11.21
AChFolder = [sFolder, '\Glut'];
%select all the files in the ACh folder
oFileList = GetFileList(AChFolder);

%get their names
aFiles = {oFileList.name};
aFileTimes = {oFileList.date};
aCellNums = GetUniqueCellNumbers(aFiles);
puff_loc = 110000;

%%
close all; 
%loop through each cell
for ii = 1:size(aCellNums,1)
    %get the files for this cell
    aFile = strfind(aFiles,['C',num2str(aCellNums(ii),'%02d')]);
    aFileNames = aFiles(find(~cellfun(@isempty, aFile)));
    aFileTime = aFileTimes(find(~cellfun(@isempty, aFile)));
    
    %loop through all files of this cell
    for jj =  1:size(aFileNames,2)
        %read file
        [aData,iSamplingInterval,oHeader] = abfload([AChFolder,'\',aFileNames{jj}]);
        
        %get experiment and cellID
        dExperiment = str2num(sFolder(end-7:end));
        aCells = oCellData.Data(oCellData.Data(:,strcmp(oCellData.Header,'Experiment'))==dExperiment,:);
        iCellID = aCells(aCells(:,strcmp(oCellData.Header,'CellNumber'))==aCellNums(ii),...
            strcmp(oCellData.Header,'CellID'));

        dSamplingFreq = 1/(iSamplingInterval/1000/1000);
        aTime = (1:1:size(aData))*iSamplingInterval/1000; %sampling interval in ms?
        
        %% create filters for later
        wo = 1/(iSamplingInterval/(10^6))/2; %sampling interval in us converted to sec 1/ for Hz
        [z,p,k] = butter(3, 2000/wo, 'low'); % 3rd order filter with a cut off of 2000Hz for data sampled at "wo" Hz
        [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
        oFilter = dfilt.df2sos(sos,g);
        
        [z,p,k] = butter(3, 0.3/wo, 'low'); % 3rd order filter with a cut off of 0.5Hz for data sampled at "wo" Hz
        [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
        oFilterTwo = dfilt.df2sos(sos,g);
        
        %%
        
        %loop through the sweeps in the file (kk)
        aDataToWrite{ii}{jj} = cell(size(aData,3),8);
        for kk = 1:size(aData,3)
            aDataToWrite{ii}{jj}{kk,1} = iCellID;
            aDataToWrite{ii}{jj}{kk,2} = jj; %recording number
            Time = str2double(regexp(aFileTime{jj}(13:end),'\d*','Match'));%find all numbers in aFileTime and make an array
            aDataToWrite{ii}{jj}{kk,3} = Time; %time of recording
            aDataToWrite{ii}{jj}{kk,4} = kk; %sweep number
            baseline = mean(aData([1:11000 end-11000:end],1,kk)); %baseline (0sec to 11sec)+last 11 sec of each sweep
            aDataToWrite{ii}{jj}{kk,5} = baseline;
            stdev_baseline = std(aData(6000:11000,1,kk));
            Fig = figure;
            axes1 = axes();
%             plot(aData(:,1,kk), 'b','parent',axes1)
            hold on;
            FilteredData = filter(oFilter,padarray(aData(:,1,kk),200,'replicate','both'));
            FilteredData = flipud(FilteredData);
            FilteredData = filter(oFilter,FilteredData);
            FilteredData = flipud(FilteredData);
            FilteredData = FilteredData(201:end-200);
             plot(FilteredData, 'Color', [0.6, 0.6, 0.6])
            
            
            FilteredData = filter(oFilterTwo,padarray(FilteredData(:,1),200,'replicate','both'));
            FilteredData = flipud(FilteredData);
            FilteredData = filter(oFilterTwo,FilteredData);
            FilteredData = flipud(FilteredData);
            FilteredData = FilteredData(201:end-200);
            
            plot(FilteredData, 'r', 'LineWidth',2,'parent',axes1)

            %find peak
            [minval, minind] = min(FilteredData(100000:200000)); %min value or bottom of peak between expected time locations
            peak_depth = diff([baseline, minval]);
            %average of 3s where peak is expected 
            peak_window_avg = mean(FilteredData(130000:160000));
            stim_level = round(squeeze(mean(aData(:,2,kk),1)),-1);
                
                aDataToWrite{ii}{jj}{kk,6} = abs(peak_depth);
                 aDataToWrite{ii}{jj}{kk,7} = diff([baseline, peak_window_avg]);
                 aDataToWrite{ii}{jj}{kk,8} = stim_level; 
                %plot min peak
                plot(100000+minind, minval, 'ko', 'Markersize', 10, 'linewidth',2);  
                hold off;
                pause(0.2);

        end

%         %% rename the file
%         FromFile = [AChFolder,'\',aFileNames{jj}];
%         ToFile = [AChFolder,'\',aFileNames{jj}(1:19) sprintf('_VC_ACh_Hex_%02d.abf', jj)];
%         if ~strcmp(FromFile,ToFile);
%             movefile(FromFile,ToFile);
%         end
        
    end
    
    %% write out results
        aDataToWrite{ii} = vertcat(aDataToWrite{ii}{:});

end

 Data = vertcat(aDataToWrite{:});
%Data = cell2mat(Data);
dlmwrite('P:\Patching\GlutResponse.csv',Data,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');   
  
close all 
end