function oNewAPData = GetNewAPsFixBridgeBalance_wip(sFolder,oCellData,oAPData)
%this function gets new AP data from the specified folder and fixed bridge
%ballance before measurements
%JS. 14/09/2021 - adapted from JA GetNewAP

% %get the cells for this experiment
dExperiment = str2double(sFolder(end-7:end));
aCellIndexes = oCellData.Data(find(oCellData.Data(:,strcmp(oCellData.Header,'Experiment')) == dExperiment),...
    strcmp(oCellData.Header,'CellID'));
aCellNums = oCellData.Data(find(oCellData.Data(:,strcmp(oCellData.Header,'Experiment')) == dExperiment),...
    strcmp(oCellData.Header,'CellNumber'));
% %get the AP file list
% %get file list in directory
oFileList = GetFileList(sFolder);
%select the files that are the right size to be memtest data
aFileSize = [oFileList(:).bytes];
aAPFiles = find(round(aFileSize,-3)==1007000 | round(aFileSize,-3)==1008000);
aFileNames = {oFileList(aAPFiles).name}';
aAPFileIndexes = strfind(aFileNames,'SingleAP');
aFirstAPFiles = find(~cellfun(@isempty,aAPFileIndexes));
aAPFiles = aAPFiles(aFirstAPFiles);
aFileNames = {oFileList(aAPFiles).name}';
aFirstFileIndices = strfind(aFileNames,'_01.abf');
aCellFileIndexes = aAPFiles(find(~cellfun(@isempty, aFirstFileIndices)));
%get cellIDs from file names
%find the cell numbers from the file name portion between 'C' and
%the next '_'
aCellIDs = cell(numel(aCellFileIndexes),1);
for mm = 1:numel(aCellIDs)
    startIndex = strfind(oFileList(aCellFileIndexes(mm)).name,'C');
    StringToCheck = oFileList(aCellFileIndexes(mm)).name(startIndex+1:end);
    endIndex = regexp(StringToCheck,'_');
    StringToCheck = StringToCheck(1:endIndex(1)-1);
    %single file so just take the cellid from the name
    aCellIDs{mm} = str2double(StringToCheck);
end
%concatenate numbers and look for unique cellids
aCellIDs = vertcat(aCellIDs{:});
aCellIDs = unique(aCellIDs);
%initialise cell array to hold data
aDataToWrite = cell(numel(aCellIDs),1);
bFirstFileLoaded = false;
%find cell indexes for given CellIDs
[C ia ib] = intersect(aCellNums,aCellIDs);
aCellIndexes = aCellIndexes(ia);
close all
% %loop through the cells
for ii = 1:numel(aCellIDs)
    %get cell details
    iExperiment = oCellData.Data(aCellIndexes(ii),strcmp(oCellData.Header,'Experiment'));
    iCellNumber = oCellData.Data(aCellIndexes(ii),strcmp(oCellData.Header,'CellNumber'));
    iCellID = oCellData.Data(aCellIndexes(ii),strcmp(oCellData.Header,'CellID'));
    iStrainID = oCellData.Data(aCellIndexes(ii),strcmp(oCellData.Header,'StrainID'));
    %find the AP files for this Cell
    if ii == numel(aCellFileIndexes)
        %this is the last cell for this experiment so just take the files
        %greater than the first memtest file for this cell
        APFileIndexes = find(aAPFiles >= aCellFileIndexes(ii));
    else
        APFileIndexes = find(aAPFiles < aCellFileIndexes(ii+1) & ...
            aAPFiles >= aCellFileIndexes(ii));
    end
    if isempty(APFileIndexes)
        aDataToWrite{ii} = [];
    else
        aDataToWrite{ii} = cell(size(APFileIndexes,1),1);
    end
    %loop through the AP files for this cell
    for jj = 1:numel(APFileIndexes)
        iFileIndex = aAPFiles(APFileIndexes(jj));
        %open the file using abfload
        sAbfFile = [sFolder,'\',oFileList(iFileIndex).name];
        [aAbfData,iSamplingInterval,oHeader] = abfload(sAbfFile);
        aTime = (1:1:size(aAbfData,1))*iSamplingInterval/1000;
        
        %create filter to use to smooth data
        if bFirstFileLoaded
            if iOldSamplingInterval ~= iSamplingInterval
                wo = 1/(iSamplingInterval/(10^6))/2;
                [z, p k] = butter(3, 1000/wo, 'low'); % 3rd order filter
                [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
                oFilter = dfilt.df2sos(sos,g);
            end
        else
            wo = 1/(iSamplingInterval/(10^6))/2;
            [z p k] = butter(3, 1000/wo, 'low'); % 3rd order filter
            [sos,g] = zp2sos(z,p,k); % Convert to 2nd order sections form
            oFilter = dfilt.df2sos(sos,g);
        end
        %% this bit inseted by JS 22.07.21 to compensated for bridge ballance
        % this is how much before and after the vertical step that is being
        % removed
        gapbefore = 2;
        gapafter = 15;
        figure
        %loop through the sweeps in the file
        aDataToWrite{ii}{jj} = cell(size(aAbfData,3),19);
        for kk = 1:size(aAbfData,3)
            aDataToWrite{ii}{jj}{kk,1} = iCellID;
            aDataToWrite{ii}{jj}{kk,2} = iStrainID;
            aDataToWrite{ii}{jj}{kk,3} = jj; %recording number
            aDataToWrite{ii}{jj}{kk,4} = kk; %sweep number
            %%find stim current
            %the derivative of channel2 aka aDiffStim to find the edges of the square wave
            aDiffStim = diff(aAbfData(:,2,kk));
            [maxval, maxind] = max(aDiffStim);
            [minval, minind] = min(aDiffStim);
            if maxval < 1
                maxind = 25 * 1000 / iSamplingInterval;
                minind = 35 * 1000 / iSamplingInterval;
            end
            
            channel1 = aAbfData(:,1,kk);
            %to find difference in height of signal increase
            offset = channel1(maxind+gapafter)- channel1(maxind-gapbefore);
            %drop the middel chunk
            channel1(maxind:minind) = channel1(maxind:minind) - offset;
            %remove sections arround shift and replace with linear gradient
            %section
            channel1(maxind-gapbefore:maxind+gapafter-1) = linspace(channel1(maxind-gapbefore), channel1(maxind+gapafter), gapbefore+gapafter);
            channel1(minind-gapbefore:minind+gapafter-1) = linspace(channel1(minind-gapbefore), channel1(minind+gapafter), gapbefore+gapafter);
            %convert back to aBf file so works with Jesses code
            aAbfData(:,1,kk) = channel1;
            
            %%
            aDataToWrite{ii}{jj}{kk,5} = mean(aAbfData(1:maxind-100,1,kk)); %baseline
            aDataToWrite{ii}{jj}{kk,6} = round(mean(aAbfData(maxind+100:minind-100,2,kk)) - mean(aAbfData(1:maxind-100,2,kk)),-2); %stimcurrent
            [aDataToWrite{ii}{jj}{kk,7}, iMax] = max(aAbfData(:,1,kk)); %Max
            [aDataToWrite{ii}{jj}{kk,8}, iMin] = min(aAbfData(:,1,kk)); %Min
            %filter in each direction because filter shifts signal (so is
            %shifted back by refiltering)
            bFirstFileLoaded = true;
            FilteredData = filter(oFilter,padarray(aAbfData(:,1,kk),200,'replicate','both'));
            FilteredData = flipud(FilteredData);
            FilteredData = filter(oFilter,FilteredData);
            FilteredData = flipud(FilteredData);
            FilteredData = FilteredData(201:end-200);
            
            %compute threshold voltage
            [maxval maxAPind] = max(diff(diff(FilteredData)));
            aDataToWrite{ii}{jj}{kk,10} = aAbfData(maxAPind,1,kk); %threshold this is same as maxval?
            
            %compute half width - width at half height from threshold
            dHalfHeight = (diff([aDataToWrite{ii}{jj}{kk,10}, aDataToWrite{ii}{jj}{kk,7}]))/2;
            aHalf = aAbfData(:,1,kk) > (aDataToWrite{ii}{jj}{kk,5} + dHalfHeight);
            aHalf = find(aHalf);
            aDataToWrite{ii}{jj}{kk,9} = aTime(aHalf(end)) - aTime(aHalf(1)); %half-width
            
            aDataToWrite{ii}{jj}{kk,11} = diff([aDataToWrite{ii}{jj}{kk,5}, aDataToWrite{ii}{jj}{kk,10}]); %AP Threshold from baseline
            aDataToWrite{ii}{jj}{kk,12} = mean(aAbfData(1:maxind-100,2,kk)); %holdingcurrent
            aDataToWrite{ii}{jj}{kk,13} = diff([aDataToWrite{ii}{jj}{kk,10}, aDataToWrite{ii}{jj}{kk,7}]); %AP amplitude
            aDataToWrite{ii}{jj}{kk,14} = diff([aDataToWrite{ii}{jj}{kk,8}, aDataToWrite{ii}{jj}{kk,5}]); %AHP amplitude (using baseline)
            aDataToWrite{ii}{jj}{kk,17} = aTime(iMin - iMax); %time from AP top to AHP bottom
            aDataToWrite{ii}{jj}{kk,18} = aTime(iMax - maxind); %time from baseline to threshold
            aDataToWrite{ii}{jj}{kk,19} = aTime(iMax - maxAPind); %time from threshold to AP
            %% only compute AHP decay if AHP is>2pA and min is <100ms after
            %stim
            if aDataToWrite{ii}{jj}{kk,14}> 2
                if iMin < 10000
                    dDecayHalfHeight = (aDataToWrite{ii}{jj}{kk,5} - aDataToWrite{ii}{jj}{kk,8})/2;
                    aDecayHalf = aAbfData(:,1,kk) < (aDataToWrite{ii}{jj}{kk,5}-dDecayHalfHeight);
                    %if want to compute 50%AHP from end baseline not
                    %begining 
%                     AHPHalfHeight = aDataToWrite{ii}{jj}{kk,8} + (baseline - aDataToWrite{ii}{jj}{kk,8})*0.8;
%                     aDecayHalf = aAbfData(:,1,kk) < AHPHalfHeight;
                    
                    [output, region] = bwlabel(aDecayHalf);
                    %find where signal crosses 80% AHP depth
                    a = output(iMin, 1);
                    aDecayHalfRegion = find(output==a);
                    
                    aDataToWrite{ii}{jj}{kk,15} = aTime(aDecayHalfRegion(end) - aDecayHalfRegion(1)); %AHP decay
                    
                    baseline = mean(aAbfData(3000:end,1,kk));
                    AHP80Height = aDataToWrite{ii}{jj}{kk,8} + (baseline - aDataToWrite{ii}{jj}{kk,8})*0.8;
                    %search for points that cross the 80%AHP depth from
                    %last AP onwards
                    aDecay80 = aAbfData(:,1,kk) < AHP80Height;
                    [output, region] = bwlabel(aDecay80);
                    %find where signal crosses 80% AHP depth
                    a = output(iMin, 1);
                    aDecay80Region = find(output==a);
                    
                    aDataToWrite{ii}{jj}{kk,16} = aTime(aDecay80Region(end) - aDecay80Region(1)); %AHP 80% decay
                    % to plot and check
                   %  plot min, decay time, tau fit to check
                    hold on
                    plot(aTime(1:8000), smooth(aAbfData(1:8000,1,kk)),'k', 'LineWidth', 3);
                    plot(aTime(1:8000), repmat(aDataToWrite{ii}{jj}{kk,5}-dDecayHalfHeight, size(aAbfData(1:8000,1,kk), 1), 1), 'r', 'LineWidth', 1); %1/2 height
                    plot(aTime(1, aDecayHalfRegion(1)), aAbfData(aDecayHalfRegion(1),1,kk), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); %AHP decay
                    plot(aTime(1, aDecayHalfRegion(end)), aAbfData(aDecayHalfRegion(end),1,kk), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); %AHP decay
                    plot(aTime(1:8000), repmat(AHP80Height, size(aAbfData(1:8000,1,kk), 1), 1), 'r', 'LineWidth', 1); %80% height
                    plot(aTime(1, aDecay80Region(1)), aAbfData(aDecay80Region(1),1,kk), 'ro',  'MarkerSize', 8, 'LineWidth', 1.5); %AHP decay 80
                    plot(aTime(1, aDecay80Region(end)), aAbfData(aDecay80Region(end),1,kk), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); %AHP decay 80
                    plot(aTime(1, iMin), aAbfData(iMin,1,kk), 'ro', 'MarkerFaceColor', 'r','MarkerSize', 10); %min
                    plot(aTime(1, maxAPind), aAbfData(maxAPind,1,kk), 'b*', 'MarkerSize',10, 'LineWidth', 1.5); %plot threshold
                    plot(aTime(1, iMax), aAbfData(iMax,1,kk), 'bo', 'MarkerSize',10, 'LineWidth', 2); %max
                    hold off
                    
                else
                    aDataToWrite{ii}{jj}{kk,15} = NaN;
                    aDataToWrite{ii}{jj}{kk,16} = NaN;
                    hold on
                    plot(aAbfData(:,1,kk),'k');
                    plot(iMin, aAbfData(iMin,1,kk), 'r+'); %min
                    hold off
                end
            else
                aDataToWrite{ii}{jj}{kk,15} = NaN;
                aDataToWrite{ii}{jj}{kk,16} = NaN;
                hold on
                plot(aAbfData(:,1,kk),'k');
                plot(iMin, aAbfData(iMin,1,kk), 'r+'); %min
                hold off
            end
            
            
            
            pause(0.2);
            
        end
        close all

        %save the sampling interval
        iOldSamplingInterval = iSamplingInterval;
        %% rename the file if necessary
        %         FromFile = sAbfFile;
        %         ToFile = [sFolder,'\',num2str(iExperiment),'Cell',...
        %             num2str(iCellNumber),'_IC_SingleAP_',sprintf('%02.0f',jj),'.abf'];
        %         if ~strcmp(FromFile,ToFile)
        %             fprintf('Renaming file %s to %s\n',FromFile,ToFile);
        % %             movefile(FromFile,ToFile);
        %         end
    end
    if ~isempty(APFileIndexes)
        %concatenate data
        aDataToWrite{ii} = vertcat(aDataToWrite{ii}{:});
    end
end
oNewAPData.Data = vertcat(aDataToWrite{:});
oNewAPData.Header = oAPData.Header;
dlmwrite('P:\Patching\ActionPotentials.csv',oNewAPData.Data,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');


