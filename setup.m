% SETUP - Add all project folders to the MATLAB path
%
% Run this script ONCE from the project root directory before using
% any other scripts. It adds the project root and all code subdirectories
% to the MATLAB path so that config(), export functions, analysis
% functions, and utilities can all find each other.
%
% Usage:
%   >> cd '/path/to/CDC Analysis'
%   >> setup
%
% After running setup, you can call any script from any working directory.

    root = fileparts(mfilename('fullpath'));

    addpath(root);
    addpath(fullfile(root, 'code', 'preprocessing'));
    addpath(fullfile(root, 'code', 'analysis'));
    addpath(fullfile(root, 'code', 'visualization'));
    addpath(fullfile(root, 'code', 'utils'));

    fprintf('CDC Analysis project paths added.\n');
    fprintf('  Root: %s\n', root);
