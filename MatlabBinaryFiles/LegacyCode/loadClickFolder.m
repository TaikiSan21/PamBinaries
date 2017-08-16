function [ clicksAll ] = loadClickFolder(folderName, species, channel)
%creates an array of from all binary files within folder and sub folder. In
%order to preserve memory the waveform data is deleted from each click.
% folderName='E:\Orkney2012\Analysis\Vert Survey';%count the number of clicks per interval. (seconds)
%click code
nFiles=1;

%get all the possible .pgdf in the folder. 
binaryFiles=findBinaryFiles(folderName);

%pre allocate array
clicksAll=[];

aclickFiltered.millis=0;
aclickFiltered.date=0;
aclickFiltered.startSample=0;
aclickFiltered.channelMap=0;
aclickFiltered.triggerMap=0;
aclickFiltered.type=0;
aclickFiltered.flags=0;
aclickFiltered.angleErrors=0;
aclickFiltered.duration=0;
aclickFiltered.nChan=0;
aclickFiltered.wave=[];
tot=1;
for j=1:length(binaryFiles);
    
     filename=binaryFiles(1,j);
      
     clicks= loadClickFile(filename{:});
     clear clickFiltered ;

     clickFiltered=repmat(aclickFiltered, 1, length(clicks));
      
     %now get the time out for just the selected clicks
     N=1;
     for k=1:length(clicks)
        if (clicks(k).type==species || species==-1)
            if (ismember(channel,getChannels(clicks(k).channelMap))==1 || channel==-1)
%                disp(['Porp clicks   ' num2str(N)]);
                clickFiltered(N)= clicks(k);
                clickFiltered(N).wave=[];
                N=N+1;
            end
       end
     end
     
     %%%%need to filter empty click structures%%%%
     cut=-1;
     for i=1:length(clickFiltered)
         if (clickFiltered(i).date<1)
             cut=i;
             break
         end
         if (cut~=-1)
             break;
         end
     end
     
     output=['Click loaded: %: ' num2str((j*100)/length(binaryFiles)), '%'];
     disp(output)
     nFiles=nFiles+1;

     if (cut==1 || cut==-1)
         continue;
     else
         clickFiltered=clickFiltered(1:cut); 
         clicksAll=[clicksAll clickFiltered];  
     end
     
     tot=tot+N;
     output=['Click loaded: thisFile: ' num2str(N), ' total:  '  num2str(tot)];
     disp(output)
        
end

end

