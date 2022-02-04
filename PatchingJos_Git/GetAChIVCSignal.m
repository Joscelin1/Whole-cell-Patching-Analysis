function oAChSignal = GetAChIVCSignal(sFolder)
%for this script the AChIVCs must be labeled with the cell number at the
%beginining of the file. It takes all the files from an IVC_ACh folder
%which should be the file in the AnalyseExperimentFolder script 
AChFolder = [sFolder];
%select all the files in the ACh folder
oFileList = GetFileList(AChFolder);

%change these to adjust for altered protocol 
baseline_start = 70000;
baseline_end = 150000;




%get their names
aFiles = {oFileList.name};
aFileTimes = {oFileList.date};
aCellNums = GetUniqueCellNumbers(aFiles);

%%
close all;

%loop through all files of this cell
for jj =  1:size(aFiles,2)
    %read file
    [aData,iSamplingInterval,oHeader] = abfload([AChFolder,'\',aFiles{jj}]);
    
    %get experiment and cellID
    iCellID = str2num(aFiles{jj}(1:2));
    
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
    aDataToWrite{jj} = cell(size(aData,3),8);
    for kk = 1:size(aData,3)
        aDataToWrite{jj}{kk,1} = iCellID;
        aDataToWrite{jj}{kk,2} = jj; %recording number
        
        %stim level for this sweep
        aCurvature = diff(diff(aData(:,2,kk),1,1),1,1);
        [stimpks,stimlocs] = findpeaks(aCurvature(:,1),'MinPeakHeight',2.5);
        if isempty(stimpks)
            stim_level = round(squeeze(mean(aData(:,2,kk),1)),-1)';
        else
            stim_level = round(squeeze(mean(aData(stimlocs(2)-(20*1000)/iSamplingInterval:...
                stimlocs(2)-(10*1000)/iSamplingInterval,2,kk),1)),-1)';
        end
        aDataToWrite{jj}{kk,3} = stim_level; %stim_level
        
        
        
        Fig = figure;
        axes1 = axes();
        plot(aData(20000:end-50000,1,kk), 'b','parent',axes1)
        hold on;
        FilteredData = filter(oFilter,padarray(aData(20000:end-50000,1,kk),200,'replicate','both'));
        FilteredData = flipud(FilteredData);
        FilteredData = filter(oFilter,FilteredData);
        FilteredData = flipud(FilteredData);
        FilteredData = FilteredData(201:end-200);
        plot(FilteredData, 'c')
        
        %data is no cropped
        FilteredData = filter(oFilterTwo,padarray(FilteredData(:,1),200,'replicate','both'));
        FilteredData = flipud(FilteredData);
        FilteredData = filter(oFilterTwo,FilteredData);
        FilteredData = flipud(FilteredData);
        FilteredData = FilteredData(201:end-200);
        
        
        plot(FilteredData, 'r', 'LineWidth',2,'parent',axes1)
        baseline = mean(aData(baseline_start:baseline_end,1,kk)); %baseline (10sec to 20sec)
        aDataToWrite{jj}{kk,4} = baseline;
        stdev_baseline = std(aData(baseline_start:baseline_end,1,kk));
        
        %find peak
        %             [maxval, maxind] = max(aData(:,1,kk));
        [minval, minind] = min(FilteredData(baseline_end:(baseline_end+150000))); %min value or bottom of peak between expected time locations
        
        
        peak_depth = diff([baseline, minval]);
        %if peak base is below baseline (weird baseline) skipp or if    baseline < minval |
        %amplitude is within 2SD of baseline variance record 0
        if  abs(diff([baseline, minval])) < 2*stdev_baseline
            aDataToWrite{jj}{kk,5} = 0;
            aDataToWrite{jj}{kk,6} = 0;
            aDataToWrite{jj}{kk,7} = 0;
            aDataToWrite{jj}{kk,8} = 0;
            plot(100000, 0, 'rx', 'MarkerSize', 40, 'linewidth',8);
            pause(0.2);
        else
            
            aDataToWrite{jj}{kk,5} = peak_depth;
            %plot min peak
            plot(baseline_end+minind, minval, 'ko');
            %compute half height absolute value
            half_height = baseline + (peak_depth/2);
            %if half height is within 2SD of baseline variance dont
            %look for peak width or rise/ decay times
            
            
            
            
            %work out width at half height
            threshold = FilteredData((baseline_end-50000):end) < half_height;
            peak_width_left = find(threshold, 1, 'first');
            peak_width_right = find(threshold, 1, 'last');
            peak_width = diff([peak_width_left, peak_width_right])/dSamplingFreq;%convert to seconds
            
            plot((baseline_end-50000)+peak_width_left, half_height, 'k*');
            plot((baseline_end-50000)+peak_width_right, half_height, 'k*');
            
            
            aDataToWrite{jj}{kk,6} = peak_width; %width at half height
            
            %compute rise time(drop in signal) between 20% - 90% of max
            %peak amplitude
            
            ten_percent_of_peak = FilteredData(baseline_start:500000,1)<=(baseline + peak_depth*0.2); %excluding firsst and last bits as filter is wrong here
            ninty_percent_of_peak = FilteredData(baseline_start:500000,1)<=(baseline + peak_depth*0.9);
            peak_rise = diff([find(ten_percent_of_peak, 1, 'first'), find(ninty_percent_of_peak, 1, 'first')])/dSamplingFreq; %time between 10% and 90% of peak depth + convert to seconds
            peak_decay = diff([find(ninty_percent_of_peak, 1, 'last'), find(ten_percent_of_peak, 1, 'last')])/dSamplingFreq;
            
            aDataToWrite{jj}{kk,7} = peak_rise; %time for signal to drop
            aDataToWrite{jj}{kk,8} = peak_decay; %time for signal to return to baseline
            
            plot_line = repmat((baseline + peak_depth*0.3), length(FilteredData(:,1)), 1);
            plot(plot_line, 'k'); %plot 30% decay threshold
            plot_line = repmat((baseline + peak_depth*0.9), length(FilteredData(:,1)), 1);
            plot(plot_line, 'k'); %plot 90% decay threshold
            
            hold off;
            pause(0.2);
        end
    end
    
    
    
    %% write out results
    %   if ~isempty(aFile)
    %         %concatenate data
%     aDataToWrite{ii} = vertcat(aDataToWrite{ii}{:});
    %     end
    
end

Data = vertcat(aDataToWrite{:});
%Data = cell2mat(Data);
dlmwrite('P:\Patching\AChIVC.csv',Data,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');

close all
