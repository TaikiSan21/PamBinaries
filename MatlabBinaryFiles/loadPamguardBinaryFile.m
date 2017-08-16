function [dataSet, fileInfo] = loadPamguardBinaryFile(fileName)
% Program to load a Pamguard binary file

% open binary file and read data
try
    fid = fopen(fileName, 'r', 'ieee-be.l64');

    % initialize variables
    prevPos = -1;
    dataSet=[];
    clear fileInfo;
    fileInfo.readModuleHeader=@readStdModuleHeader;
    fileInfo.readModuleFooter=@readStdModuleFooter;
    
    % main loop
    while (1)
        
        % if for some reason we're stuck at one byte, warn the user and
        % abort
        pos = ftell(fid);
        if (pos == prevPos)
            disp(fprintf('File stuck at byte %d', pos));
            break;
        end
        prevPos = pos;
        
        % read in the length of the object and the type.  Move the file
        % pointer back to the beginning of the length variable, and switch
        % on the type.  If we couldn't read nextLen or nextType, assume
        % that means we've hit the end of the file and break out of loop
        [nextLen, nL] = fread(fid, 1, 'int32');
        [nextType, nT] = fread(fid, 1, 'int32');
        if (nL == 0 || nT == 0)
            break;
        end
        fseek(fid, -8, 'cof');
        switch nextType
            
            % Case 1: File Header.  Read in the file header, and then set 
            % the function handles depending on what type of data the
            % binary file holds.  The module type is typically specified in
            % the package class that extends PamControlledUnit, and is a
            % string unique to that module.
            case -1
                fileInfo.fileHeader = readFileHeader(fid);
                switch fileInfo.fileHeader.moduleType
                    
                    % AIS Processing Module
                    case 'AIS Processing'
                        fileInfo.objectType=0;
                        fileInfo.readModuleData=@readAISData;
                        
                    % Click Detector Module
                    case 'Click Detector'
                        fileInfo.objectType=1000;
                        fileInfo.readModuleData=@readClickData;
                        fileInfo.readModuleFooter=@readClickFooter;
                     
                    % Clip Generator Module
                    case 'Clip Generator'
                        fileInfo.objectType=[1 2];
                        fileInfo.readModuleData=@readClipData;
                       
                    % DbHt Module
                    % NOT TESTED YET
                    case 'DbHt'
                        fileInfo.objectType=1;
                        fileInfo.readModuleData=@readDbHtData;
                   
                    % Difar Module
                    % NOT TESTED YET
                    case 'DIFAR Processing'
                        fileInfo.objectType=0;
                        fileInfo.readModuleData=@readDifarData;
                        
                    % LTSA Module
                    case 'LTSA'
                        fileInfo.objectType=1;
                        fileInfo.readModuleHeader=@readLTSAHeader;
                        fileInfo.readModuleData=@readLTSAData;
                        
                    % Noise Monitor Module and Noise Band Monitor Module
                    % Note: The two modules have different types, but both
                    % use noiseMonitor.NoiseBinaryDataSource class to save
                    % data
                    case {'Noise Monitor', 'Noise Band'}
                        fileInfo.objectType=1;
                        fileInfo.readModuleHeader=@readNoiseMonHeader;
                        fileInfo.readModuleData=@readNoiseMonData;
                        
                    % Filtered Noise Measurement Module (Noise One Band)
                    case 'NoiseBand'
                        fileInfo.objectType=1;
                        fileInfo.readModuleData=@readNoiseBandData;
                        
                    % Right Whale Edge Detector Module
                    case 'RW Edge Detector'
                        fileInfo.objectType=0;
                        fileInfo.readModuleData=@readRWEDetectorData;
                        disp('Right Whale Edge Detector binary file detected');
                        disp('Note that the low, high and peak frequencies are actually');
                        disp('saved as FFT slices.  In order to convert values to Hz, they');
                        disp('must be multiplied by (sampleRate/fftLength)');
                        
                    % Whistle And Moans Module
                    case 'WhistlesMoans'
                        fileInfo.objectType=2000;
                        fileInfo.readModuleHeader=@readWMDHeader;
                        fileInfo.readModuleData=@readWMDData;
                        
                    % Note: PamRawDataBlock has it's own Binary Store (RawDataBinarySource),
                    % but it is created by multiple different processes so doesn't have one
                    % particular type, and may have a type shared with a different binary
                    % store (e.g. the DbHt module creates both a DbHtDataSource and a
                    % RawDataBinarySource, and they will both have type 'DbHt').
                    % The comments in the class indicate that the binary store should never
                    % be used for raw data and that the sole purpose of the class is to
                    % enable sending raw data over the network.  If this is ever changed
                    % and raw data is allowed to be stored in the binary files, we will
                    % have to come up with a way of determining where the raw data came
                    % from besides querying the type.

                    otherwise
                        disp(['Don''t recognize type ' fileInfo.fileHeader.moduleType '.  Aborting load']);
                        break;
                end
                
            % Case 2: File Footer.  The file version should have been set
            % when we read the file header.  If the file header is empty,
            % something has gone wrong so warn the user and exit
            case -2
                if (isempty(fileInfo.fileHeader))
                    disp('Error: found file footer before file header.  Aborting load');
                    break;
                end
                fileInfo.fileFooter = readFileFooterInfo(fid, fileInfo.fileHeader.fileFormat);
                
            % Case 3: Module Header.  The correct function handle should
            % have been set when we read the file header.  If the file
            % header is empty, something has gone wrong so warn the user
            % and exit
            case -3
                if (isempty(fileInfo.fileHeader))
                    disp('Error: found module header before file header.  Aborting load');
                    break;
                end
                fileInfo.moduleHeader = fileInfo.readModuleHeader(fid);
                
            % Case 4: Module Footer.  The correct function handle should
            % have been set when we read the file header.  If the file
            % header is empty, something has gone wrong so warn the user
            % and exit
            case -4
                if (isempty(fileInfo.fileHeader))
                    disp('Error: found module footer before file header.  Aborting load');
                    break;
                end
                fileInfo.moduleFooter = fileInfo.readModuleFooter(fid);
                
            % Case 5: Data.  The correct function handle should have been
            % set when we read in the file header.  If the file header is
            % empty, something has gone wrong so warn the user and exit
            otherwise
                if (isempty(fileInfo.fileHeader))
                    disp('Error: found data before file header.  Aborting load');
                    break;
                end
                dataPoint=readPamData(fid, fileInfo);
                
                % TO DO - Preallocation
                % at this stage we know the size of a single structure.
                % Therefore, we might be able to preallocate the dataSet
                % structure array to increase speed.  Check out the
                % accepted answer here:
                % http://stackoverflow.com/questions/13664090/how-to-initialize-an-array-of-structs-in-matlab
                % maybe when we get close to filling the array, we can just
                % preallocate another block of structures and append them
                % onto dataSet the way it's done below...
                dataSet = [dataSet dataPoint];
        end
    end

catch mError
    disp('Error reading file');
    disp(getReport(mError));
end

% close the file and return to the calling function
fclose(fid);
return;
