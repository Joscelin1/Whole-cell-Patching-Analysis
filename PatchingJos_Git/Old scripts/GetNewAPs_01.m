function oNewAPData = GetNewAPs(sFolder,oCellData,oAPData)
%this function gets new AP data from the specified folder

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
% % % % % %Work out here how to only take the first AP file 
% % % % % aAPFileIndexes = zeros(size(aCellNums));
% % % % % for i = 1: length(aCellNums)
% % % % % aAPFileIndexes = strfind(aFileNames, ['_C0', num2str(aCellNums(i))]);
% % % % % end

aFirstAPFiles = find(~cellfun(@isempty,aAPFileIndexes));
aAPFiles = aAPFiles(aFirstAPFiles);
aFileNames = {oFileList(aAPFiles).name}';
aFirstFileIndices = strfind(aFileNames,'_01.abf');
aCellFileIndexes = aAPFiles(find(~cellfun(@isempty, aFirstFileIndices)));
%get cellIDs from file names
%find the cell numbers from the file name portion between 'Cell' and
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
% oFig = figure();
% set(oFig,'position',[179   567   560   420]);
% oAxes = axes('parent',oFig);
% oOverlay = axes('parent',oFig,'position',get(oAxes,'position'));
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
        %loop through the sweeps in the file
        aDataToWrite{ii}{jj} = cell(size(aAbfData,3),12);
        for kk = 1:size(aAbfData,3)
            aDataToWrite{ii}{jj}{kk,1} = iCellID;
            aDataToWrite{ii}{jj}{kk,2} = iStrainID;
            aDataToWrite{ii}{jj}{kk,3} = jj; %recording number
            aDataToWrite{ii}{jj}{kk,4} = kk; %sweep number
            %%find stim current
            aDiffStim = diff(aAbfData(:,2,kk));
            [maxval maxind] = max(aDiffStim);
            [minval minind] = min(aDiffStim);
            if maxval < 1
                maxind = 25 * 1000 / iSamplingInterval;
                minind = 35 * 1000 / iSamplingInterval;
            end
            aDataToWrite{ii}{jj}{kk,5} = mean(aAbfData(1:maxind-100,1,kk)); %baseline
            aDataToWrite{ii}{jj}{kk,6} = round(mean(aAbfData(maxind+100:minind-100,2,kk)) - mean(aAbfData(1:maxind-100,2,kk)),-2); %stimcurrent
            [aDataToWrite{ii}{jj}{kk,7}, iMax] = max(aAbfData(:,1,kk)); %APMax
            [aDataToWrite{ii}{jj}{kk,8}, iMin] = min(aAbfData(:,1,kk)); %APMin
            
            bFirstFileLoaded = true;
            FilteredData = filter(oFilter,padarray(aAbfData(:,1,kk),200,'replicate','both'));
            FilteredData = flipud(FilteredData);
            FilteredData = filter(oFilter,FilteredData);
            FilteredData = flipud(FilteredData);
            FilteredData = FilteredData(201:end-200);
            %             iWindowSize = 15;
            %             oWindow = bartlett(iWindowSize);
            %             aFilteredDiff = diff(diff(filter(oWindow,sum(oWindow),aAbfData(:,1,kk))));
            
            %compute threshold voltage
            [maxval maxAPind] = max(diff(diff(FilteredData)));
            aDataToWrite{ii}{jj}{kk,10} = aAbfData(maxAPind,1,kk); %threshold
            
            %compute half width
            %             dHalfHeight = (aDataToWrite{ii}{jj}{kk,7} - aDataToWrite{ii}{jj}{kk,8})/2;
            aHalf = aAbfData(:,1,kk) > aDataToWrite{ii}{jj}{kk,10}; %(aDataToWrite{ii}{jj}{kk,7}-dHalfHeight);
            aHalf = find(aHalf);
            aDataToWrite{ii}{jj}{kk,9} = aTime(aHalf(end)) - aTime(aHalf(1)); %half-width
                        
            %compute AHP decay
            dDecayHalfHeight = (aDataToWrite{ii}{jj}{kk,5} - aDataToWrite{ii}{jj}{kk,8})/2;
            aDecayHalf = aAbfData(iMin:end,1,kk) > (aDataToWrite{ii}{jj}{kk,5}-dDecayHalfHeight);
            aDecayHalf = find(aDecayHalf);
            aDataToWrite{ii}{jj}{kk,11} = aTime(aDecayHalf(1)+iMin-1) - aTime(iMin); %decay
            aDataToWrite{ii}{jj}{kk,12} = mean(aAbfData(1:maxind-100,2,kk)); %holdingcurrent

            
%             plot(oAxes,aTime,aAbfData(:,1,kk),'k');
%             plot(oOverlay,aTime,vertcat(NaN,NaN,diff(diff(FilteredData))),'r');
%             set(oOverlay,'color','none');
%             axis(oOverlay,'off');
%             hold(oAxes,'on');
%             set(oAxes,'xlim',[20 50]);
%             set(oOverlay,'xlim',[20 50]);
%             scatter(oAxes,aTime(aHalf(end)),aAbfData(aHalf(end),1,kk),'k+');
%             scatter(oAxes,aTime(aHalf(1)),aAbfData(aHalf(1),1,kk),'k+');
%             scatter(oAxes,aTime(maxAPind),aAbfData(maxAPind,1,kk),'r+');
%             text(25,20,num2str(iCellID),'parent',oAxes);
%             hold(oAxes,'off');
        end
        %save the sampling interval
        iOldSamplingInterval = iSamplingInterval;
        %         %rename the file if necessary
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

