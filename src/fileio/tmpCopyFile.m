function [tmpCopyFilePath] = tmpCopyFile(filePath,tmpDir)
% Will make a temporary copy of fileName in tmpDir and unzip it if it ends in .gz or .zip
% If file is not compressed them tmpCopyFilePath is the same as filePath
% function [tmpCopyFilePath] = tmpCopyFile(filePath,tmpDir)
% --------------------------------------------------------------------------------------------------
% INPUT ARGUMENTS
% ---------------
% filePath: full path of file whose temp copy is needed
% tmpDir (OPTIONAL): target directory. If not supplied, it is obtained using tempdir()
% ---------------
% OUTPUT ARGUMENTS
% ---------------
% tmpCopyFilePath: full path of temp copy of (decompressed) fileName

[~,fileName,fileExt] = fileparts(filePath);
isCompressed = find( strcmpi(fileExt , {'.gz','.zip'} ) ); % []: not compress, 1:gz, 2:zip
if ~exist('tmpDir','var')
    tmpDir = tempdir();
end
tmpPrefix = tempname(tmpDir); % temp file prefix

if isCompressed
    tmpCopyFilePathGz = [tmpPrefix,fileName,fileExt];
    tmpCopyFilePath = [tmpPrefix,fileName];
    % copy the file
    [cpStatus,cpMsg,cpMsgId] = copyfile(filePath,tmpCopyFilePathGz,'f'); 
    assert(cpStatus,cpMsgId,'ERROR: Unable to copy file %s to temp directory %s: %s',filePath,tmpDir,cpMsg);
    % decompress the file
    switch isCompressed
        case 1
            [status,~] = system( sprintf( 'gunzip -c %s 1> %s 2> /dev/null' , tmpCopyFilePathGz, tmpCopyFilePath ) );
        case 2
            [status,~] = system( sprintf( 'unzip -p %s 1> %s 2> /dev/null' , tmpCopyFilePathGz, tmpCopyFilePath ) );
    end
    assert( ~status , 'ERROR: Unable to decompress file %s' , filePath );
    delete(tmpCopyFilePathGz); % remove temporary compressed file
else
    tmpCopyFilePath = filePath;
end

end