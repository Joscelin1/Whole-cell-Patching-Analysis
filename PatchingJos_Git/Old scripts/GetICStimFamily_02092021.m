function oICStimFamilyData = GetICStimFamily_02092021(sFolder,oCellData)
%get the list of files in this folder to check for IC family files
%get file list in directory
oFileList = GetFileList(sFolder);
%select the files that are the right size to be stimfamily data
aFileSize = [oFileList(:).bytes];
allICStimFiles = find(round(aFileSize,-3)==8007000);
%get their names
aICStimFiles = {oFileList(allICStimFiles).name}';
aTime = 1:250000;
%maybe not necisary?
aFileIndexes(:,1) = allICStimFiles;

%only proceed for folders where there is at least one IC family file
if ~isempty(aFileIndexes)
    aCellNums = GetUniqueCellNumbers(aICStimFiles);
    %initialise cell arrays to hold results
    aNumAPs = cell(size(aCellNums));
    aStimLevels = cell(size(aCellNums));
    aRecordingNum = cell(size(aCellNums));
    aIntervalFirst = cell(size(aCellNums));
    aIntervalLast = cell(size(aCellNums));
    aAmpFirst = cell(size(aCellNums));
    aAmpLast = cell(size(aCellNums));
    aTrainAHP = cell(size(aCellNums));
    aTrainAHPduration = cell(size(aCellNums));
    bFirstFileLoaded = false;
    for ii = 1:size(aCellNums,1)
        %get the files for this cell
        aICStimFile = strfind(aICStimFiles,['C',num2str(aCellNums(ii),'%02d')]);
        aFileNames = aICStimFiles(find(~cellfun(@isempty, aICStimFile)));
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
            aTheseIntervalFirst = zeros(1,size(aSignals,2));
            aTheseIntervalLast = zeros(1,size(aSignals,2));
            aTheseAmpFirst = zeros(1,size(aSignals,2));
            aTheseAmpLast = zeros(1,size(aSignals,2));
            aTheseTrainAHP = zeros(1,size(aSignals,2));
            aTheseTrainAHPduration = zeros(1,size(aSignals,2));
            datatowrite = zeros(7,size(aSignals,2));
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
                    
                end
                %plot trace and peak locs
                 figure
                    hold on
                    plot(FilteredData);
                    plot(siglocs(:,1), sigpks(:,1), '*', 'MarkerFaceColor', 'g');
                
                
                %continue with analysis for files with more than 4 APs
                if length(siglocs) > 4
                    aTheseIntervalFirst(mm) = (siglocs(2,1)-siglocs(1,1))/iSamplingInterval;
                    aTheseIntervalLast(mm) = (siglocs((end),1)-siglocs((end-1),1))/iSamplingInterval;
                    aTheseAmpFirst(mm) = sigpks(2,1);
                    aTheseAmpLast(mm) = sigpks(end,1);
                    %compute AHP decay
%                     baseline = mean(FilteredData(1:8000));
                    [min_AHP, min_loc] = min(FilteredData(55000:end));
                    %only looking for min after train has finished
                    min_loc = min_loc +55000;
                    
                    
                    
                    
                    if min_loc < 75000
                        end_baseline = mean(FilteredData(end-8000:end));    
                    aTheseTrainAHP(mm) = end_baseline - min_AHP;
                    %                     duration_AHP
                    AHPHalfHeight = min_AHP + (aTheseTrainAHP(mm)/2);
                    %search for points that cross the halfAHP depth from
                    %last AP onwards
                    half = FilteredData(1:end, 1) < AHPHalfHeight;
                    
                    [outputlabel, regioncount] = bwlabel(half);
                    %find where signal crosses 1/2 AHP depth
                    %                     a = find(outputlabel(60000:end, 1)~=0, 1, 'first')+60000;
                    a = outputlabel(min_loc, 1);
                    half_region = find(outputlabel==a);
                    AHP_width = aTime(half_region(end)) - aTime(half_region(1)); %half-width in miniutes
                    aTheseTrainAHPduration(mm) = AHP_width;
                    
                    %continue plot if AHP exists 
                    plot(min_loc, min_AHP,  '*', 'MarkerFaceColor', 'b');
                    plot(half_region(end), AHPHalfHeight, '.','MarkerFaceColor','r', 'MarkerSize', 10);
                    plot(half_region(1), AHPHalfHeight, '.','MarkerFaceColor','r', 'MarkerSize', 10);
                    
                    %find apropriate end for curve fitting (80% of AHP
                    %depth)
                    AHP90Height = min_AHP + (aTheseTrainAHP(mm)*0.8);
                    %search for points that cross the halfAHP depth from
                    %last AP onwards
                    ninety = FilteredData(1:end, 1) < AHP90Height;
                    
                    [output, region] = bwlabel(ninety);
                    %find where signal crosses 1/2 AHP depth
                    %                     a = find(outputlabel(60000:end, 1)~=0, 1, 'first')+60000;
                    a = output(min_loc, 1);
                    ninety_region = find(output==a);
                    
                    %take time from min peak onwards to fit exponential 
                    x_time = aTime(min_loc+500:ninety_region(end)+1000);
                    %define the time points to evalulate the fit this
                    %should be more than x_time to be accurate 
                    xval = linspace(min(x_time), max(x_time), 2000);
                    
                    %get the signal to fit
                    y_signal = FilteredData(min_loc+500:ninety_region(end)+1000, 1);
                    %run exponential fit use 1/y_signal to make data an
                    %exponential decay rather than rise
                    OutData = FitExponential_AHP(x_time,1./y_signal, xval);
                                     fit_x = OutData.xval;
                    %convert back to a rise exponential 
                    fit_y = 1./OutData.yval;
                    
                    
                    plot(fit_x, fit_y, 'r')
                    hold off
                    figure
                    hold on
                    plot(fit_x, fit_y, 'r')
                    plot(x_time, y_signal, 'b')
                    hold off
                    %calculate tau from fit parameters this is C in the exponential
                    %curve fit function
                    Tau = abs(1/OutData.est(3));
                    %calculate Q1 from area under the curve
                    %        Q1 = trapz(x_time,y_signal-oMemTest(jj).I1(ii));
%                     else
%                      aTheseTrainAHP(mm) = 0;
%                         aTheseTrainAHPduration(mm) = 0;
%                         
%                         
%                         %                     plot(min_loc, min_AHP,  '.', 'MarkerSize', 20);

                    else
                        hold off
                    
                    end    
                end
                
                %                  else
                %                     %get steady state value
                %                     aTheseSteadyState(mm) = round(mean(aSignals(stimlocs(2)-1000:stimlocs(2)-500,mm)),1);
                
                
            end
            aRecordingNum{ii} = horzcat(aRecordingNum{ii},ones(1,numel(aTheseAPs))*jj);
            aStimLevels{ii} = horzcat(aStimLevels{ii},aTheseStimLevel);
            aNumAPs{ii} = horzcat(aNumAPs{ii},aTheseAPs);
            aIntervalFirst{ii} = horzcat(aIntervalFirst{ii}, aTheseIntervalFirst);
            aIntervalLast{ii} = horzcat(aIntervalLast{ii}, aTheseIntervalLast);
            aAmpFirst{ii} = horzcat(aAmpFirst{ii}, aTheseAmpFirst);
            aAmpLast{ii} = horzcat(aAmpLast{ii}, aTheseAmpLast);
            aTrainAHP{ii} = horzcat(aTrainAHP{ii}, aTheseTrainAHP);
            aTrainAHPduration{ii} = horzcat(aTrainAHPduration{ii}, aTheseTrainAHPduration);
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
            aIntervalFirst{ii}',...
            aIntervalLast{ii}',...
            aAmpFirst{ii}',...
            aAmpLast{ii}',...
            aTrainAHP{ii}',...
            aTrainAHPduration{ii}');
        dlmwrite('P:\Patching\ICStimFamily.csv',aDataToWrite,'-append',...
            'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');
    end
end

end