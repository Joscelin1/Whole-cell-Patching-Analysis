function oICStimFamilyData = GetICStimFamily(sFolder,oCellData)
%this script takes stim family data and records no of AP and AHP properties
%(also fits tau)
%JS. 15.09.2021

%get the list of files in this folder to check for IC family files
%get file list in directory
oFileList = GetFileList(sFolder);
%select the files that are the right size to be stimfamily data
aFileSize = [oFileList(:).bytes];
allICStimFiles = find(round(aFileSize,-3)==8007000);

%get their names
aICStimFiles = {oFileList(allICStimFiles).name}';
%only keep ones that include stim files in name
iICStimFiles = strfind(aICStimFiles,'StimFamily');
%get file names
aCellFiles = aICStimFiles(find(~cellfun(@isempty, iICStimFiles)));

aTime = 1:250000;

%only proceed for folders where there is at least one IC family file
if ~isempty(aCellFiles(:,1))
    aCellNums = GetUniqueCellNumbers(aCellFiles);
    %initialise cell arrays to hold results
    aNumAPs = cell(size(aCellNums));
    aStimLevels = cell(size(aCellNums));
    aRecordingNum = cell(size(aCellNums));
    aAmpFirst = cell(size(aCellNums));
    aAmpLast = cell(size(aCellNums));
    aTrainAHP = cell(size(aCellNums));
    aTrainAHPduration = cell(size(aCellNums));
    aTrainAHPtau = cell(size(aCellNums));
    aIntervalChange = cell(size(aCellNums));
    apk_interval = [];
    bFirstFileLoaded = false;
    for ii = 1:size(aCellNums,1)
        %get the files for this cell
        aICStimFile = strfind(aCellFiles,['C',num2str(aCellNums(ii),'%02d')]);
        aFileNames = aCellFiles(find(~cellfun(@isempty, aICStimFile)));
        %only do this for first file of cell (this avoids taking other
        %files of the same size that are other recordings)
        for jj = 1
            %read file
            [aAbfData,iSamplingInterval,oHeader] = abfload([sFolder,'\',aFileNames{jj}]);
            dSamplingFreq = 1/(iSamplingInterval/1000/1000);
            %% get stim levels
            %get the stim channel index and make sure that the channels are
            %labelled properly..
            iStimChannel = strcmp(oHeader.recChUnits,'pA');
            if isempty(find(iStimChannel)) || length(find(iStimChannel)) > 1
                continue;
            end
            %find curvature peaks aka double diff of most rapidly changing bit (assume that all signals will have same
            %peak indices)
            aCurvature = diff(diff(aAbfData(:,iStimChannel,:),1,1),1,1);
            [stimpks,stimlocs] = findpeaks(aCurvature(:,1,2),'MinPeakHeight',2.5);
            %             plot(aCurvature(:,1));
            if length(stimlocs) > 2
                [B I] = sort(stimpks,'descend');
                [stimpks,stimlocs] = findpeaks(aCurvature(:,1,2),'MinPeakHeight',B(2)-1);
            end
            %get min and max stim current values
            aStimBaseline = squeeze(mean(aAbfData(stimlocs(1)-(10*1000)/iSamplingInterval:...
                stimlocs(1)-(5*1000)/iSamplingInterval,iStimChannel,:),1));
            aTheseStimLevel = round(squeeze(mean(aAbfData(stimlocs(2)-(20*1000)/iSamplingInterval:...
                stimlocs(2)-(10*1000)/iSamplingInterval,iStimChannel,:),1)),-1)';
            %% get steady state voltage response and number of APs
            %get the voltage channel
            iSignalChannel = strcmp(oHeader.recChUnits,'mV');
            aSignals = squeeze(aAbfData(:,iSignalChannel,:));
            %create filter to use to smooth data
            if bFirstFileLoaded
                if iOldSamplingInterval ~= iSamplingInterval
                    wo = 1/(iSamplingInterval/(10^6))/2;
                    [z p k] = butter(3, 1000/wo, 'low'); % 3rd order filter
                    [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
                    oFilter = dfilt.df2sos(sos,g);
                end
            else
                wo = 1/(iSamplingInterval/(10^6))/2;
                [z p k] = butter(3, 1000/wo, 'low'); % 3rd order filter
                [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
                oFilter = dfilt.df2sos(sos,g);
            end
            %find APs
            aTheseAPs = zeros(1,size(aSignals,2));
            aTheseAmpFirst = zeros(1,size(aSignals,2));
            aTheseAmpLast = zeros(1,size(aSignals,2));
            aTheseTrainAHP = zeros(1,size(aSignals,2));
            aTheseTrainAHPduration = zeros(1,size(aSignals,2));
            aTheseTrainAHPtau = zeros(1,size(aSignals,2));
            aTheseIntervalChange = zeros(1,size(aSignals,2));
            aThesepk_interval = [];
            %loop through sweeps
            for mm = 1:size(aSignals,2)
                %filter signal
                bFirstFileLoaded = true;
                FilteredData = filter(oFilter,padarray(aSignals(:,mm),200,'replicate','both'));
                FilteredData = flipud(FilteredData);
                FilteredData = filter(oFilter,FilteredData);
                FilteredData = flipud(FilteredData);
                FilteredData = FilteredData(201:end-200);
                %find the peaks to get number of APs
                [sigpks,siglocs] = findpeaks(FilteredData,'MinPeakHeight',20);
                if ~isempty(siglocs)
                    %save number of APs if there are any
                    aTheseAPs(mm) = numel(siglocs);
                    %for 200pA current only work out inter peak distances
                    %and save seperately 
                    if mm == 3
                    %measure distances between peaks
                    for j = 1:size(siglocs, 1)-1
                        %convert into ms
                        aThesepk_interval(1, j) = 1000*(siglocs(j+1, 1) - siglocs(j, 1))/dSamplingFreq;
                    end
                    end
                end
                
                %plot trace and peak locs
                figure
                hold on
                plot(FilteredData);
                plot(siglocs(:,1), sigpks(:,1), '*', 'MarkerFaceColor', 'g');

                %continue with analysis for files with more than 4 APs
                if length(siglocs) > 4 
                    baseline = mean(FilteredData(1:8000));
                    aTheseAmpFirst(mm) = diff([baseline, sigpks(2,1)]);
                    aTheseAmpLast(mm) = diff([baseline, sigpks(end,1)]);
                    aTheseIntervalChange(mm) = ((siglocs(2,1)-siglocs(1,1))/(siglocs(end,1)-siglocs(end-1,1))); %first interval/last interval 
                    %compute AHP decay
                    [min_AHP, min_loc] = min(FilteredData(siglocs(end):end));
                    %only looking for min after train has finished so last
                    %AP
                    min_loc = min_loc +siglocs(end);

                    %choose only min points that are close to end of train 
                    if min_loc < 75000
                        end_baseline = mean(FilteredData(end-8000:end));
                        aTheseTrainAHP(mm) = end_baseline - min_AHP;
                        AHPHalfHeight = min_AHP + (aTheseTrainAHP(mm)/2);
                        %search for points that cross the halfAHP depth from
                        %last AP onwards
                        half = FilteredData(1:end, 1) < AHPHalfHeight;
                        [outputlabel, regioncount] = bwlabel(half);
                        %find where signal crosses 1/2 AHP depth
                        a = outputlabel(min_loc, 1);
                        half_region = find(outputlabel==a);
                        AHP_width = aTime(half_region(end)) - aTime(half_region(1)); %half-width in miniutes
                        aTheseTrainAHPduration(mm) = AHP_width/dSamplingFreq;
                        
                        %continue plot if AHP exists
                        plot(min_loc, min_AHP,  '*', 'MarkerFaceColor', 'b');
                        plot(half_region(end), AHPHalfHeight, '.','MarkerFaceColor','r', 'MarkerSize', 10);
                        plot(half_region(1), AHPHalfHeight, '.','MarkerFaceColor','r', 'MarkerSize', 10);
                        
                        %% fit curve to 80% of AHP amplitude to work out Tau
                        AHP80Height = min_AHP + (aTheseTrainAHP(mm)*0.8);
                        %search for points that cross the halfAHP depth from
                        %last AP onwards
                        eighty = FilteredData(1:end, 1) < AHP80Height;
                        [output, region] = bwlabel(eighty);
                        %find where signal crosses 80% AHP depth
                        a = output(min_loc, 1);
                        eighty_region = find(output==a);
                        %take time from min peak onwards to fit exponential
                        x_time = aTime(min_loc+500:eighty_region(end)+1000);
                        %define the time points to evalulate the fit this
                        %should be more than x_time to be accurate
                        xval = linspace(min(x_time), max(x_time), 2000);
                        %get the signal to fit
                        y_signal = FilteredData(min_loc+500:eighty_region(end)+1000, 1);
                        %run exponential fit use 1/y_signal to make data an
                        %exponential decay rather than rise
                        OutData = FitExponential_AHP(x_time,1./y_signal, xval);
                        fit_x = OutData.xval;
                        %convert back to a rise exponential
                        fit_y = 1./OutData.yval;
                        
                        
                        plot(fit_x, fit_y, 'r')
                        hold off
%                         figure
%                         hold on
%                         plot(fit_x, fit_y, 'r')
%                         plot(x_time, y_signal, 'b')
%                         hold off
                        %calculate tau from fit parameters this is C in the exponential
                        %curve fit function
                        Tau = abs(1/OutData.est(3));
                        aTheseTrainAHPtau(mm) = Tau/dSamplingFreq;

                        %                     plot(min_loc, min_AHP,  '.', 'MarkerSize', 20);
                        
                    else
                        hold off
                        
                    end
                end
                
%             %plot example figure 
%             
%               
%                 figure
%                 hold on
%                 plot((aTime(1:150000)./dSamplingFreq), FilteredData(1:150000),'k', 'LineWidth', 2);
%                 plot((aTime(siglocs(:,1))./dSamplingFreq), sigpks(:,1), 'b*',  'MarkerSize', 8);
%                 
%                 plot((aTime(min_loc)./dSamplingFreq), min_AHP,  'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
%                 plot((aTime(1:150000)./dSamplingFreq), repmat(AHPHalfHeight, size(FilteredData(1:150000), 1), 1), 'r');
%                 plot((aTime(half_region(end))./dSamplingFreq), AHPHalfHeight, 'ro', 'MarkerSize', 8);
%                 plot((aTime(half_region(1))./dSamplingFreq), AHPHalfHeight, 'ro', 'MarkerSize', 8);
% %                 plot(fit_x, fit_y, 'b', 'LineWidth', 1)
%                 hold off
% %                 
%                  figure
%                         hold on
%                         plot(fit_x, fit_y, 'r','LineWidth', 1.5)
%                         plot(x_time, y_signal, 'k','LineWidth', 1.5)
%                         hold off
                
            end
            aRecordingNum{ii} = horzcat(aRecordingNum{ii},ones(1,numel(aTheseAPs))*jj);
            aStimLevels{ii} = horzcat(aStimLevels{ii},aTheseStimLevel);
            aNumAPs{ii} = horzcat(aNumAPs{ii},aTheseAPs);
            
            aAmpFirst{ii} = horzcat(aAmpFirst{ii}, aTheseAmpFirst);
            aAmpLast{ii} = horzcat(aAmpLast{ii}, aTheseAmpLast);
            aTrainAHP{ii} = horzcat(aTrainAHP{ii}, aTheseTrainAHP); %AHP amplitude 
            aTrainAHPduration{ii} = horzcat(aTrainAHPduration{ii}, aTheseTrainAHPduration); %AHP duration 
            aTrainAHPtau{ii} = horzcat(aTrainAHPtau{ii}, aTheseTrainAHPtau); %AHP tau 
           aIntervalChange{ii} = horzcat( aIntervalChange{ii}, aTheseIntervalChange); %first gap/last gap
             apk_interval = catpad(1, apk_interval, aThesepk_interval); 
            %save the sampling interval
            iOldSamplingInterval = iSamplingInterval;
            %%rename the StimFamily file
            %     FromFile = [sFolder,'\',aFileNames{jj}];
            %     ToFile = [sFolder,'\',aFileNames{jj}(1:19) '_IC_StimFamily_01.abf'];
            %     if ~strcmp(FromFile,ToFile);
            %         movefile(FromFile,ToFile);
            %         end
            
        end
        
        %% write out results
        %get experiment and cellID
        dExperiment = str2num(sFolder(end-7:end));
        aCells = oCellData.Data(oCellData.Data(:,strcmp(oCellData.Header,'Experiment'))==dExperiment,:);
        iCellID = aCells(aCells(:,strcmp(oCellData.Header,'CellNumber'))==aCellNums(ii),...
            strcmp(oCellData.Header,'CellID'));
        aDataToWrite = horzcat(ones(numel(aStimLevels{ii}),1)*iCellID,...
            aRecordingNum{ii}',...
            aStimLevels{ii}',...
            aNumAPs{ii}',...
            aAmpFirst{ii}',...
            aAmpLast{ii}',...
            aTrainAHP{ii}',...
            aTrainAHPduration{ii}', ...
            aTrainAHPtau{ii}',...
            aIntervalChange{ii}');
        dlmwrite('P:\Patching\ICStimFamily.csv',aDataToWrite,'-append',...
            'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');
        % write AP intervals to a seperate file 
         aDataToWrite_intervals = cat(2, iCellID, aThesepk_interval);
        dlmwrite('P:\Patching\ICStimFamily_intervals.csv', aDataToWrite_intervals,'-append',...
            'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');
    end
end

end