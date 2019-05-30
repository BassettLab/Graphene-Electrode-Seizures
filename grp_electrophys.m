function C = grp_electrophys(F)
% Housekeeping and accessing IEEG portal data
%--------------------------------------------------------------------------
fs      = filesep; 
sess    = IEEGSession('10_18_18_bath1', 'roschkoenig', F.pwdfile);
sesb    = IEEGSession('10_18_18_baseline', 'roschkoenig', F.pwdfile); 

% Specify experiment details
%--------------------------------------------------------------------------
Fs     = 5000; 

chlabs = sess.data.channelLabels;
chans  = [ 2, 1, 16, 15; 
           5, 8, 11, 14; 
           4, 7, 10, 13;
           3, 6, 9, 12 ]; 
       
chanx  = [115, 182, 252, 320;
          120, 188, 256, 324; 
          124, 192, 260, 328;
          129, 197, 265, 333];
      
chany  = [148, 144, 140, 135;
          217, 212, 207, 203;
          284, 280, 275, 271; 
          352, 348, 343, 339];
       
nchans = length(chans(:)); 
nsamps = 60 * 14 * Fs;       

% Specify channel details
%--------------------------------------------------------------------------
clear C
k = 0; 
for r = 1:size(chans,1)
for c = 1:size(chans,2)
    
    k = k + 1;
    id       = find(~cellfun(@isempty, regexp(chlabs(:,1), num2str(chans(r,c), '%02.f'))));
    C(k).row = r;
    C(k).col = c;
    C(k).oi  = id(1);       % original id (i.e. index in chanlist)
    C(k).id  = chans(r,c); 
    C(k).x   = chanx(r,c); 
    C(k).y   = chany(r,c); 
    
    % Load two minutes baseline and derive z-score
    %----------------------------------------------------------------------
    C(k).bas = sesb.data.getvalues(10000:20000, C(k).oi); 
end
end

% Find onset trigger
%--------------------------------------------------------------------------
trig = find(strcmp(chlabs(:,1), 'Trigger'));
tdat = sess.data.getvalues(1:nsamps, trig); 
[val loc] = max(tdat); 
tvec = ([1:length(tdat)] - loc) / (Fs * 60); 

toload = find(tvec >= 0); 
tvec   = tvec(toload); 

%% Load data into the channel spec variable
%--------------------------------------------------------------------------
for c = 1:length(C)
    C(c).dat = sess.data.getvalues(1:nsamps, C(c).oi); 
end

fbp{1} = [1 40];
fbp{2} = [80 250];
fbp{3} = [500 2450]; 

dat     = horzcat(C.dat)';
bdat    = horzcat(C.bas)'; 

% Notch (or narrow bandstop) filter to remove 60Hz line noise
%--------------------------------------------------------------------------
clear fdat
ndat = dat;     nbdat = bdat;
for m = 1:5
    ndat    = ft_preproc_bandstopfilter(ndat, Fs, [60*m-1 60*m+1], 16000, 'firws');
    nbdat   = ft_preproc_bandstopfilter(nbdat, Fs, [60*m-1 60*m+1], 16000, 'firws'); 
end

% Bandpass filter into set frequency bands
%--------------------------------------------------------------------------
for f = 1:length(fbp)
    fdat{f}     = ft_preproc_bandpassfilter(ndat, Fs, fbp{f}, 16000, 'firws'); 
    fbdat{f}    = ft_preproc_bandpassfilter(nbdat, Fs, fbp{f}, 16000, 'firws'); 
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