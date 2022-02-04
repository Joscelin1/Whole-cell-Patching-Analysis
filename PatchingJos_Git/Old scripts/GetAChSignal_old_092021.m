function oAChSignal = GetAChSignal(sFolder,oCellData)
AChFolder = [sFolder, '\ACh_Bt'];
%select all the files in the ACh folder
oFileList = GetFileList(AChFolder);

%get their names
aFiles = {oFileList.name};
aFileTimes = {oFileList.date};
aCellNums = GetUniqueCellNumbers(aFiles);

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
        aTime = (1:1:size(aData))*iSamplingInterval/1000; %sampling interval in us
        
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
        aDataToWrite{ii}{jj} = cell(size(aData,3),9);
        for kk = 1:size(aData,3)
            aDataToWrite{ii}{jj}{kk,1} = iCellID;
            aDataToWrite{ii}{jj}{kk,2} = jj; %recording number
            Time = str2double(regexp(aFileTime{jj}(13:end),'\d*','Match'));%find all numbers in aFileTime and make an array
            aDataToWrite{ii}{jj}{kk,3} = Time; %time of recording
            aDataToWrite{ii}{jj}{kk,4} = kk; %sweep number
            baseline = mean(aData([1:11000 end-11000:end],1,kk)); %baseline (0sec to 11sec)+last 11 sec of each sweep
            aDataToWrite{ii}{jj}{kk,5} = baseline;
            stdev_baseline = std(aData([1:11000 end-11000:end],1,kk));
            Fig = figure;
            axes1 = axes();
            plot(aData(:,1,kk), 'b','parent',axes1)
            hold on;
            FilteredData = filter(oFilter,padarray(aData(:,1,kk),200,'replicate','both'));
            FilteredData = flipud(FilteredData);
            FilteredData = filter(oFilter,FilteredData);
            FilteredData = flipud(FilteredData);
            FilteredData = FilteredData(201:end-200);
            plot(FilteredData, 'c')
            
            
            FilteredData = filter(oFilterTwo,padarray(FilteredData(:,1),200,'replicate','both'));
            FilteredData = flipud(FilteredData);
            FilteredData = filter(oFilterTwo,FilteredData);
            FilteredData = flipud(FilteredData);
            FilteredData = FilteredData(201:end-200);
            
            
            plot(FilteredData, 'r', 'LineWidth',2,'parent',axes1)
            
            
            
            %compute threshold voltage
            %find peak
            %             [maxval, maxind] = max(aData(:,1,kk));
            [minval, minind] = min(FilteredData(100000:200000)); %min value or bottom of peak between expected time locations
            
            
            peak_depth = diff([baseline, minval]);
            %if peak base is below baseline (weird baseline) skipp or peak
            %depth is not bigger than -5mV it is not real
            if  baseline < minval
                aDataToWrite{ii}{jj}{kk,6} = 0;
                aDataToWrite{ii}{jj}{kk,7} = 0;
                aDataToWrite{ii}{jj}{kk,8} = 0;
                aDataToWrite{ii}{jj}{kk,9} = 0;
            else
                
                
                
                
                %change this
                aDataToWrite{ii}{jj}{kk,6} = peak_depth;
                %plot min peak
                plot(100000+minind, minval, 'ko');
                %compute half height absolute value
                half_height = baseline + (peak_depth/2);
                %if half height is within 2SD of baseline variance dont
                %look for peak width or rise/ decay times
                if abs(diff([baseline, half_height])) < 2*stdev_baseline
                    aDataToWrite{ii}{jj}{kk,7} = 0;
                    aDataToWrite{ii}{jj}{kk,8} = 0;
                    aDataToWrite{ii}{jj}{kk,9} = 0;
                    plot(100000, 0, 'rx', 'MarkerSize', 40, 'linewidth',8);
                    pause(0.2);
                else
                
                
                
                %work out width at half height
                threshold = FilteredData(100000:300000) < half_height;
                peak_width_left = find(threshold, 1, 'first');
                peak_width_right = find(threshold, 1, 'last');
                peak_width = diff([peak_width_left, peak_width_right])/dSamplingFreq;%convert to seconds
                
                plot(100000+peak_width_left, half_height, 'k*');
                plot(100000+peak_width_right, half_height, 'k*');
                
                
                aDataToWrite{ii}{jj}{kk,7} = peak_width; %width at half height
                
                %compute rise time(drop in signal) between 20% - 90% of max
                %peak amplitude
                
                ten_percent_of_peak = FilteredData(100000:500000,1)<=(baseline + peak_depth*0.2); %excluding firsst and last bits as filter is wrong here
                ninty_percent_of_peak = FilteredData(100000:500000,1)<=(baseline + peak_depth*0.9);
                peak_rise = diff([find(ten_percent_of_peak, 1, 'first'), find(ninty_percent_of_peak, 1, 'first')])/dSamplingFreq; %time between 10% and 90% of peak depth + convert to seconds
                peak_decay = diff([find(ninty_percent_of_peak, 1, 'last'), find(ten_percent_of_peak, 1, 'last')])/dSamplingFreq;
                
                aDataToWrite{ii}{jj}{kk,8} = peak_rise; %time for signal to drop
                aDataToWrite{ii}{jj}{kk,9} = peak_decay; %time for signal to return to baseline
                
                plot_line = repmat((baseline + peak_depth*0.2), length(FilteredData(:,1)), 1);
                plot(plot_line, 'k'); %plot 20% decay threshold
                plot_line = repmat((baseline + peak_depth*0.9), length(FilteredData(:,1)), 1);
                plot(plot_line, 'k'); %plot 90% decay threshold
                
                hold off;
                pause(0.2);
                
                end
            end
        end
        
        
        
        
%         %% rename the file
%         FromFile = [AChFolder,'\',aFileNames{jj}];
%         ToFile = [AChFolder,'\',aFileNames{jj}(1:19) sprintf('_VC_ACh_Hex_%02d.abf', jj)];
%         if ~strcmp(FromFile,ToFile);
%             movefile(FromFile,ToFile);
%         end
        
    end
    
    %% write out results
%   if ~isempty(aFile)
%         %concatenate data 
        aDataToWrite{ii} = vertcat(aDataToWrite{ii}{:});
%     end   

end

 Data = vertcat(aDataToWrite{:});
%Data = cell2mat(Data);
dlmwrite('P:\Patching\AChResponse.csv',Data,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');   
  
close all 
end