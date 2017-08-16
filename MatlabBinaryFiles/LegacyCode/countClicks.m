clear

%Count the number of clicks in specified intervals in a binary file.
%Iterates through all the subfolders in a typical 'binary' folder. 

folderName='C:\Users\spn1\Desktop\Vert Survey\Binary';
%count the number of clicks per interval. (seconds)
interval= 60*5;
%click code
species=1;

subFolders = dir(folderName);

for i=1:length(subFolders);
    
   if (strcmp( subFolders(i).name,'.')==1 || strcmp( subFolders(i).name,'..')==1)
       continue; 
   end
  
   subFolderName=[folderName,'\', subFolders(i).name];
   
   binaryFiles=dir(subFolderName);
   
   for j=1:length(binaryFiles);
       
       if (strcmp(binaryFiles(j).name,'.')==1 || strcmp(binaryFiles(j).name,'..')==1)
      	continue; 
       end
       
       binaryFileName=[subFolderName,'\',binaryFiles(j).name];
       disp(binaryFileName);
       binaryChar=char(binaryFileName);
       fileEnd=binaryChar(length(binaryChar)-3: length(binaryChar));
       %read the binary click file
       if (strcmp(fileEnd,'pgdf')==1)
       click= loadClickFile(binaryFileName);
       end
       
   end
    
end



