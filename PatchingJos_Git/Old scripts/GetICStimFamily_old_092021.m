function oICStimFamilyData = GetICStimFamily(sFolder,oCellData)
%get the list of files in this folder to check for IC family files
%get file list in directory
oFileList = GetFileList(sFolder);
%select the files that are the right size to be stimfamily data
aFileSize = [oFileList(:).bytes];
allICStimFiles = find(round(aFileSize,-3)==8007000);
%get their names
aICStimFiles = {oFileList(allICStimFiles).name}';

%maybe not necisary?
aFileIndexes(:,1) = allICStimFiles;

%only proceed for folders where there is at least one IC family file
if ~isempty(aFileIndexes)
    aCellNums = GetUniqueCellNumbers(aICStimFiles);
    %initialise cell arrays to hold results
    aNumAPs = cell(size(aCellNums));
    aSteadyState = cell(size(aCellNums));
    aStimLevels = cell(size(aCellNums));
    aRecordingNum = cell(size(aCellNums));
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
            %find curvature peaks (assume that all signals will have same
            %peak indices)
            aCurvature = diff(diff(aAbfData(:,iStimChannel,:),1,1),1,1);
            [stimpks,stimlocs] = findpeaks(aCurvature(:,1,2),'MinPeakHeight',2.5);
            plot(aCurvature(:,1));
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
            %find APs or steady state
            aTheseAPs = zeros(1,size(aSignals,2));
            aTheseSteadyState = zeros(1,size(aSignals,2));
            for mm = 1:size(aSignals,2)
                %filter signal
                bFirstFileLoaded = true;
                FilteredData = filter(oFilter,padarray(aSignals(:,mm),200,'replicate','both'));
                FilteredData = flipud(FilteredData);
                FilteredData = filter(oFilter,FilteredData);
                FilteredData = flipud(FilteredData);
                FilteredData = FilteredData(201:end-200);
                %find the peaks to get number of APs
                [sigpks,siglocs] = findpeaks(FilteredData,'MinPeakHeight',0);
                if ~isempty(siglocs)
                    %save number of APs
                    aTheseAPs(mm) = numel(siglocs);
                else
                    %get steady state value
                    aTheseSteadyState(mm) = round(mean(aSignals(stimlocs(2)-1000:stimlocs(2)-500,mm)),1);
                end
            end
            aRecordingNum{ii} = horzcat(aRecordingNum{ii},ones(1,numel(aTheseAPs))*jj);
            aStimLevels{ii} = horzcat(aStimLevels{ii},aTheseStimLevel);
            aNumAPs{ii} = horzcat(aNumAPs{ii},aTheseAPs);
            aSteadyState{ii} = horzcat(aSteadyState{ii},aTheseSteadyState);
            %save the sampling interval
            iOldSamplingInterval = iSamplingInterval;
             %%rename the StimFamily file
    FromFile = [sFolder,'\',aFileNames{jj}];
    ToFile = [sFolder,'\',aFileNames{jj}(1:19) '_IC_StimFamily_01.abf'];
    if ~strcmp(FromFile,ToFile);
        movefile(FromFile,ToFile);
        end
       
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
%             aSteadyState{ii}',...
            aNumAPs{ii}');
        dlmwrite('P:\Patching\ICStimFamily.csv',aDataToWrite,'-append',...
            'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');
    end
end

end