function F = grp_housekeeping

fs          = filesep;
F.base      = '/Volumes/GoogleDrive/My Drive/Research/1903 Mouse Calcium Data';
F.scripts   = [F.base fs '00 - Scripts'];
F.data      = [F.base fs '01 - Data'];
F.analysis  = [F.base fs '02 - Analysis']; 
F.ieeg      = '/Users/roschkoenig/Dropbox/Research/0002 Tools/IEEG';
F.addpaths  = {'/Users/roschkoenig/Dropbox/Research/0002 Tools/IEEG specs', 
               '/Users/roschkoenig/Dropbox/Research/0002 Tools/Tools/cbrewer'};

F.pwdfile   = [F.ieeg fs 'ros_ieeglogin.bin'];
addpath(genpath(F.ieeg))
addpath(genpath(F.scripts)); 
for a = 1:length(F.addpaths),   addpath(F.addpaths{a});     end 

spm('defaults', 'eeg'); 