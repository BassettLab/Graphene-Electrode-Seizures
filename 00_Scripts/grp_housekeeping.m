function F = grp_housekeeping
% This function loads folder specifications and places everything into a
% single structure 'F'

fs          = filesep;
F.base      = '/Volumes/GoogleDrive/My Drive/Research/1903 Mouse Calcium Data';
F.scripts   = [F.base fs '00_Scripts'];
F.data      = [F.base fs '01_Data'];
F.analysis  = [F.base fs '02_Analysis']; 

% These need editing <- probably include in the repo
%--------------------------------------------------------------------------
uname       = char(java.lang.System.getProperty('user.name'));
F.ieeg      = ['/Users/' uname '/Dropbox/Research/0002 Tools/IEEG'];
F.addpaths  = {['/Users/' uname '/Dropbox/Research/0002 Tools/Tools/cbrewer'] 
               ['/Users/' uname '/Dropbox/Research/0002 Tools/Tools/cfc']
               ['/Users/' uname '/Dropbox/Research/0002 Tools/nmfv1_4'] 
               ['/Users/' uname '/Dropbox/Research/0002 Tools/spm']};

F.pwdfile   = ['~/Documents/IEEG specs' fs 'ros_ieeglogin.bin'];   %% < This needs changing to your own file
addpath(genpath(F.ieeg))
addpath(genpath(F.scripts)); 
for a = 1:length(F.addpaths),   addpath(F.addpaths{a});     end 

spm('defaults', 'eeg'); 