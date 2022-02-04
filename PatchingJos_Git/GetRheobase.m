function oGetRheobase = GetRheobase(sFolder,oCellData)
%this script takes rheobase files, renames them and saves the files that
%have 0 or 1 AP in acending order of stim input 
%saves to Rheobase excel 
%JS. 15/09/2021
%get the list of files in this folder to check for rheobase files
aDataToWrite = [];
%get file list in directory
oFileList = GetFileList(sFolder);
%select the files that are the right size to be rheobase
aFileSize = [oFileList(:).bytes];
aRheoFiles = find(aFileSize==87552);
%get their names
aFiles = {oFileList(aRheoFiles).name}';

%maybe not necisary?
aFileIndexes(:,1) = aRheoFiles;

%only proceed for folders where there is at least one rheobase file
if ~isempty(aFileIndexes)
    aCellNums = GetUniqueCellNumbers(aFiles);
    %initialise cell arrays to hold results
    aNumAPs = cell(size(aCellNums));
    aSteadyState = cell(size(aCellNums));
    aStimLevels = cell(size(aCellNums));
    aRecordingNum = cell(size(aCellNums));
    bFirstFileLoaded = false;
    
   
        
    for ii = 1:size(aCellNums,1)
        %get the files for this cell
        aRheoFile = strfind(aFiles,['C',num2str(aCellNums(ii),'%02d')]);
        aFileNames = aFiles(find(~cellfun(@isempty, aRheoFile)));
        %only do this for first file of cell (this avoids taking other
        %files of the same size that are other recordings)
        
          %get experiment and cellID
        dExperiment = str2num(sFolder(end-7:end));
        aCells = oCellData.Data(oCellData.Data(:,strcmp(oCellData.Header,'Experiment'))==dExperiment,:);
        iCellID = aCells(aCells(:,strcmp(oCellData.Header,'CellNumber'))==aCellNums(ii),...
            strcmp(oCellData.Header,'CellID'));
        
        aTheseAPs = []; %rheobase data for this cell
        for jj = 1:size(aFileNames,1)
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
            
            %the derivative of channel2 aka aDiffStim to find the edges of the square wave
            aDiffStim = diff(aAbfData(:,iStimChannel));
            [maxval, maxind] = max(aDiffStim);
            [minval, minind] = min(aDiffStim);
            % find stim current
            stim = round(mean(aAbfData(maxind+100:minind-100,2)) - mean(aAbfData(1:maxind-100,2)), 0); %stimcurrent
            
            %find the peaks to get number of APs 30mV threshold
            [sigpks,siglocs] = findpeaks(aAbfData(:,1),'MinPeakHeight',30);
%             hold on
%             plot(aAbfData(:,1))
            
            if numel(siglocs) <= 1
                aTheseAPs(jj, 1) = iCellID;
                aTheseAPs(jj, 2) = jj; %recoding num
                aTheseAPs(jj, 3) = stim; %stim level
                aTheseAPs(jj, 4) = numel(siglocs); %no of APs
                
            end
        end
%         hold off
        aTheseAPs= aTheseAPs(any(aTheseAPs,2),:);
        aTheseAPs = sortrows(aTheseAPs, 3);
        %save the sampling interval
        iOldSamplingInterval = iSamplingInterval;
        %             %%rename the StimFamily file
        %             FromFile = [sFolder,'\',aFileNames{jj}];
        %             ToFile = [sFolder,'\',aFileNames{jj}(1:19) sprintf('_IC_Rheobase_%02d.abf', jj)];
        %             if ~strcmp(FromFile,ToFile);
        %                 movefile(FromFile,ToFile);
        %             end
        
        
        
        %% write out results
       
%         aData = horzcat(iCellID, aTheseAPs(:,:));
        aDataToWrite = cat(1, aDataToWrite, aTheseAPs);
        
    end
end
dlmwrite('P:\Patching\Rheobase.csv',aDataToWrite,'-append',...
    'roffset',0,'coffset',0,'delimiter',',','precision','%10.3f');
end