function F = grp_housekeeping(F)
% This function loads folder specifications and places everything into a
% single structure 'F'

fs          = filesep;
F.scripts   = [F.base fs '00_Scripts'];
F.data      = [F.base fs '01_Data'];
F.analysis  = [F.base fs '02_Analysis']; 

% These need editing <- probably include in the repo
%--------------------------------------------------------------------------
F.pwdfile   = ['~/Documents/IEEG specs' fs 'ros_ieeglogin.bin'];   %% < This needs changing to your own file
if isfield(F, 'ieeg'), addpath(genpath(F.ieeg)); end
addpath(genpath(F.scripts)); 
