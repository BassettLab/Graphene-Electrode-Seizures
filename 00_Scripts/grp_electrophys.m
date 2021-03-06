function C = grp_electrophys(F)
% This function takes the folder-location structure 'F' and returns a
% (large) single variable 'C' that contains all data used for subsequent
% analysis; 
% It will also save this dataset in the '01_Data' folder

% Housekeeping and accessing IEEG portal data
%--------------------------------------------------------------------------
fs      = filesep; 

% This needs some generalising
sess    = IEEGSession('10_18_18_bath1', 'roschkoenig', F.pwdfile);
sesb    = IEEGSession('10_18_18_baseline', 'roschkoenig', F.pwdfile); 

% Specify experiment details
%--------------------------------------------------------------------------
Fs     = sess.data.sampleRate; 

chlabs = sess.data.channelLabels;
chans  = [ 2, 1, 16, 15;        % channel id mapped onto grid layout
           5, 8, 11, 14; 
           4, 7, 10, 13;
           3, 6, 9, 12 ]; 
       
chanx  = [115, 182, 252, 320;   % x position in pixels
          120, 188, 256, 324; 
          124, 192, 260, 328;
          129, 197, 265, 333];
      
chany  = [148, 144, 140, 135;   % y position in pixels
          217, 212, 207, 203;
          284, 280, 275, 271; 
          352, 348, 343, 339];
       
nchans = length(chans(:));      % Number of channels
nsamps = 60 * 14 * Fs;          % Number of samples in 14min window of interest

% Specify channel details
%--------------------------------------------------------------------------
clear C             % Structure in which data will  be collated
k = 0;              % Counting variable 

for r = 1:size(chans,1)     % Looping through rows of channels
for c = 1:size(chans,2)     % Looping through columns of channels 
    
    k = k + 1;
    % Find channel label that matches channel id given in the 'chans' layout
    id       = find(~cellfun(@isempty, regexp(chlabs(:,1), num2str(chans(r,c), '%02.f'))));
    C(k).row = r;   % keeps track of row position from 'chans' layout matrix
    C(k).col = c;   % keeps track of col position from 'chans' layout matrix 
    C(k).oi  = id(1);       % original id (i.e. index in chanlist)
    C(k).id  = chans(r,c);   
    C(k).x   = chanx(r,c);  % x position in pixels from layout 
    C(k).y   = chany(r,c);  % y position in pixels from layout
    
    % Load two minutes baseline and derive z-score
    %----------------------------------------------------------------------
    C(k).bas = sesb.data.getvalues(2*60*Fs : 4*60*Fs-1, C(k).oi); 
end
end

% Find onset trigger - this will generate a time vector in minutes
%--------------------------------------------------------------------------
trig = find(strcmp(chlabs(:,1), 'Trigger'));    % Trigger signal corresponds to imaging onset time 
tdat = sess.data.getvalues(1:nsamps, trig);     
[val loc] = max(tdat); 
tvec = ([1:length(tdat)] - loc) / (Fs * 60);  

toload = find(tvec >= 0);                       % Only load samples that are coming after the trigger 
tvec   = tvec(toload);                          

% Load data into the channel spec variable
%--------------------------------------------------------------------------
for c = 1:length(C)
    C(c).dat = sess.data.getvalues(toload, C(c).oi); 
    C(c).Fs  = Fs; 
end

fbp{1} = [1 40];        % 'LFP' range
fbp{2} = [80 250];      % High gamma range

dat     = horzcat(C.dat)';
bdat    = horzcat(C.bas)'; 

% Notch (or narrow bandstop) filter to remove 60Hz line noise
%--------------------------------------------------------------------------
clear fdat

disp('Applying notch filtering to the data')             
ndat = dat;     nbdat = bdat;   % ndat - notch filtered data;
for m = 1:5     % filtering for multiple harmonics of line noise
    disp(['Harmonic ' num2str(m) ': ' num2str(60*m) 'Hz'])
    ndat    = bandstop(ndat', [60*m-1, 60*m+1], Fs, 'Steepness', 0.75)';
    nbdat   = bandstop(nbdat', [60*m-1, 60*m+1], Fs, 'Steepness', 0.75)';
end

% Bandpass filter into set frequency bands
%--------------------------------------------------------------------------
for f = 1:length(fbp)
    disp(['Applying bandpass filter: ' num2str(fbp{f}(1)) ' to ' num2str(fbp{f}(2)) 'Hz'])  
    fdat{f}     = bandpass(ndat', fbp{f}, Fs, 'Steepness', 0.75)'; 
    fbdat{f}    = bandpass(nbdat', fbp{f}, Fs, 'Steepness', 0.75)'; 
end

% Reattach filtered data to the Channel object
%--------------------------------------------------------------------------

for f = 1:length(fdat)
for c = 1:length(C)
    C(c).filt{f} = fdat{f}(c,:); 
    C(c).filb{f} = fbdat{f}(c,:);
end
end

% Save everything in one massive file
%--------------------------------------------------------------------------
save([F.data fs 'Electrophysiology' fs 'Channel_Data.mat'], 'C', '-v7.3'); 