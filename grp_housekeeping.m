function F = grp_housekeeping
% This function loads folder specifications and places everything into a
% single structure 'F'

fs          = filesep;
F.base      = '/Volumes/GoogleDrive/My Drive/Research/1903 Mouse Calcium Data';
F.scripts   = [F.base fs '00 - Scripts'];
F.data      = [F.base fs '01 - Data'];
F.analysis  = [F.base fs '02 - Analysis']; 
F.ieeg      = '/Users/roschkoenig/Dropbox/Research/0002 Tools/IEEG';
F.addpaths  = {'/Users/roschkoenig/Dropbox/Research/0002 Tools/Tools/cbrewer', 
               '/Users/roschkoenig/Dropbox/Research/0002 Tools/Tools/cfc'
                '/Users/roschkoenig/Dropbox/Research/0002 Tools/nmfv1_4'};

F.pwdfile   = ['~/Documents/IEEG specs' fs 'ros_ieeglogin.bin'];   %% < This needs changing to your own file
addpath(genpath(F.ieeg))
addpath(genpath(F.scripts)); 
for a = 1:length(F.addpaths),   addpath(F.addpaths{a});     end 

spm('defaults', 'eeg'); 