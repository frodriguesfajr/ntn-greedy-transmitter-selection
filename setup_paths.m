function setup_paths()
% SETUP_PATHS Configure repository paths required by the MATLAB scripts.
%
% The function adds the project source directory and the included Vallado
% SGP4 routines to the beginning of the MATLAB search path.
%
% Required directories:
%   src/
%   third_party/vallado_sgp4/

rootDir = fileparts(mfilename('fullpath'));

srcDir = fullfile(rootDir,'src');
valladoDir = fullfile(rootDir,'third_party','vallado_sgp4');

if ~isfolder(srcDir)
    error('Source directory not found: %s',srcDir);
end

if ~isfolder(valladoDir)
    error('Vallado SGP4 directory not found: %s',valladoDir);
end

addpath(srcDir,'-begin');
addpath(valladoDir,'-begin');

end
