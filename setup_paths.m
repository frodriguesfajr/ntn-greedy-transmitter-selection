function setup_paths()
%SETUP_PATHS Configure paths for the reproducibility package.

rootDir = fileparts(mfilename('fullpath'));

srcDir = fullfile(rootDir,'src');
valladoDir = fullfile(rootDir,'third_party','vallado_sgp4');

if ~isfolder(srcDir)
    error('Source directory not found: %s',srcDir);
end

if ~isfolder(valladoDir)
    error('Vallado directory not found: %s',valladoDir);
end

addpath(srcDir,'-begin');
addpath(valladoDir,'-begin');

end
