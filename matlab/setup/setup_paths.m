function repoRoot = setup_paths()
%SETUP_PATHS Configure paths for the NTN transmitter-selection repository.
%
% Canonical location:
%   <repo>/matlab/setup/setup_paths.m
%
% Usage:
%   repoRoot = setup_paths();
%
% The function adds the MATLAB source tree recursively. Data/output
% directories are accessed through explicit paths and are not added.

    thisFile = mfilename('fullpath');

    if isempty(thisFile)
        error('setup_paths.m must be saved before it is executed.');
    end

    setupDir = fileparts(thisFile);      % <repo>/matlab/setup
    matlabDir = fileparts(setupDir);     % <repo>/matlab
    repoRoot = fileparts(matlabDir);     % <repo>

    if ~isfolder(fullfile(repoRoot,'.git'))
        error([ ...
            'Could not confirm repository root from setup_paths.m.\n' ...
            'Expected .git under:\n%s'], ...
            repoRoot);
    end

    addpath(matlabDir,'-begin');
    addpath(genpath(matlabDir),'-begin');

    valladoDir = fullfile( ...
        repoRoot,'third_party','vallado_sgp4');

    if ~isfolder(valladoDir)
        error('Vallado SGP4 directory not found:\n%s',valladoDir);
    end

    addpath(valladoDir,'-begin');

    fprintf('Repository paths configured from:\n  %s\n',repoRoot);
end
