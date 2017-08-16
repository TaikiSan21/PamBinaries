%search a folder for files types

function [filenames]= findBinaryFiles(filePath)
filenames={};
subFiles=dir(filePath);
for i=1:length(subFiles);
    
   if (strcmp(subFiles(i).name,'.')==1 || strcmp(subFiles(i).name,'..')==1)
       continue; 
   end
    
   if (subFiles(i).isdir==1)
      subFolderName=[filePath,'\', subFiles(i).name];
      filenames=cat(2, filenames,findBinaryFiles(subFolderName));
   
   else
    binaryFileName=[filePath,'\',subFiles(i).name];
    binaryChar=char(binaryFileName);
    fileEnd=binaryChar(length(binaryChar)-3: length(binaryChar));
     %add to binary file list
     if (strcmp(fileEnd,'pgdf')==1)
       filenames=cat(2,filenames ,binaryFileName(1,:));
     end
   end

end

end
   