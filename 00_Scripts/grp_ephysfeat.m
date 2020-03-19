function [F, Fcor] = grp_ephysfeat(C)

% Electrophys feature calculator
%==========================================================================
% 1) Phase locked high gamma
% 2) Ictal spike amplitude and frequency
% 3) Bandpower in different frequency bands
% 4) Correlation structure in different frequency bands 

% Everything will be stored in a Feature structure 'F' with approx 1Hz
% sampling rate / interpolation

% Set up timing parameters
%--------------------------------------------------------------------------
Fs      = C(1).Fs;  
rnge    = 1:length(C(1).dat);

% Identify starting points of each imaging frame (sampled at 10Hz)
srts    = fix([rnge(1)/Fs : 0.1 : rnge(end)/Fs] * Fs); 
dlen    = length(srts)-1; 
doplt   = 0;

%% Phase locked high gamma
%--------------------------------------------------------------------------
% Calculating phase locked high gamma - this code is based on work by Dr
% Shennan Aibel-Weiss (Ictal  Map v0.1, developed 2011-14)

disp('Calculating phase locked high gamma'); 
clear PL 
for c = 1:length(C)

    % Define filters and sampling rates 
    %----------------------------------------------------------------------
    bp{1}   = [4 30]; 
    bp{2}   = [80 250]; 
    td      = C(c).dat; 

    % Calculate hilbert transform instantaneous amplitude
    %----------------------------------------------------------------------
    for b = 1:2
    ft{b}   = fir1(1000, [bp{b}(1)/(Fs/2), bp{b}(2)/(Fs/2)]); 
    df{b}   = filtfilt(ft{b}, 1, td); 
    hlb{b}  = hilbert(df{b});
    amplo{b}  = abs(hlb{b});  
    end
    
    % Calculate phase of low frequencies and envelope phase of high freq
    %----------------------------------------------------------------------
    lphase = atan2(imag(hlb{1}), real(hlb{1})); % LFP signal phase
    envhlb = hilbert(amplo{2});                 % Hilbert transform on high gamma envelope (instantaneous amplitude)
    phsamp = atan2(imag(envhlb), real(envhlb)); % High Gamma envelope phase

    % Sliding window to calculate phase locked high gamma (normalised by
    % baseline)
    %----------------------------------------------------------------------
    % Richard to change this into a higher frequency sampling (akin to
    % imaging)
    
    win     = 3 * Fs;                       % 3 second window 
    nwin = fix(length(df{1}) / win); 
    clear plv plhg
    for n = 1:nwin
        sid     = [1:win] + (n-1)*win;      % sample ids that belong to the window 
        i       = sqrt(-1);
        plv(n)  = abs(mean(exp(i*(lphase(sid) - phsamp(sid)))));
        plhg(n) = abs(mean(amplo{2}(sid).*exp(i*(lphase(sid) - phsamp(sid)))));
    end

    % Interpolate to full dataset and save
    %----------------------------------------------------------------------
    PL(c).val = plv;
    PL(c).hg  = plhg; 
    F(c).plhg = interp1(linspace(1,dlen,length(PL(c).hg)), PL(c).hg, 1:dlen);

end

% % Ictal spikes
% %--------------------------------------------------------------------------
% % Ictal spike frequency and amplitude
% %--------------------------------------------------------------------------
% disp('Calculating ictal spika amplitude and frequency'); 
% clear spks ampl smspk smamp
% for c = 1:16
% [z , mu, sig]   = zscore(C(c).filt{1}); 
% [pks loc]       = findpeaks(mu-C(c).filt{1}(rnge));
% loc             = loc(pks > 5 * sig); 
% pks             = pks(pks > 5 * sig); 
% 
% for t = 1:length(srts)-1
%     td      = C(c).filt{1}(srts(t):srts(t+1)-1);
%     spks(c,t) = length(intersect(find(loc >= srts(t)), find(loc < srts(t+1))));
%     ampl(c,t) = max(td) - min(td); 
% end
% 
% smamp(c,:) = smooth(ampl(c,:), 50); 
% smspk(c,:) = smooth(spks(c,:), 50); 
% 
% F(c).spamp = smamp(c,:);
% F(c).spcnt = smspk(c,:); 
% 
% end

% Bandpower
%--------------------------------------------------------------------------
disp('Calculating Bandpower'); 
for c = 1:length(C)     % Looping through channels 
    amplo = abs(hilbert(C(c).filt{1}));     % LFP range amplitude derived from Hilberr transform
    amphi = abs(hilbert(C(c).filt{2}));     % High gamma range ampl derived from Hilbert transform
    
    F(c).lopow = 0;
    F(c).hipow = 0; 
    
    win     = 3 * Fs;                       % Richard to make more steps as above
    nwin = fix(length(df{1}) / win);                   
    for n = 1:nwin                          % Loop through all windows 
        sid  = [1:win] + (n-1)*win;         
        ld = amplo(sid); 
        hd = amphi(sid); 
        F(c).lopow(n) = mean(ld); 
        F(c).hipow(n) = mean(hd);
    end
    
	F(c).lopow = interp1(linspace(1,dlen,length(F(c).lopow)), F(c).lopow, 1:dlen);
	F(c).hipow = interp1(linspace(1,dlen,length(F(c).hipow)), F(c).hipow, 1:dlen);
    
end

% Functional network
%--------------------------------------------------------------------------
disp('Calculating functional network edges'); 

% Dummy to find indices for unique channel pairings / edges 
%--------------------------------------------------------------------------
dum     = ones(length(C)); dum = triu(dum,1); dum = find(dum); 
win     = 3 * Fs;   % 3 second sliding window - Richard to adapt
nwin = fix(length(df{1}) / win);
clear Fcor

for f = 1:2     % Loop through two frequency bands 
    dat = []; fcor    = []; 
    % loop over channels to pull out data
    for c = 1:length(C),    dat = [dat; C(c).filt{f}];      end;    dat = dat';

    % Sliding window 
    for n = 1:nwin
        sid     = [1:win] + (n-1)*win;  % Sample IDs
        cd      = corr(dat(sid, :));    % Correlation matrix between channels 
        fcor    = [fcor, cd(dum)];      % Save unique edges (upper triangle) as single column in fcor matrix
    end
    
    clear Fc
    for r = 1:size(fcor,1)
        Fc(r,:) = interp1(linspace(1,dlen,size(fcor,2)), fcor(r,:), 1:dlen);
    end

    Fcor(f).dat = Fc; 
end

% Plot everything together
%--------------------------------------------------------------------------
if doplt
    colormap gray
    pn = 7;     % plot number

    subplot(pn,1,1);    imagesc(vertcat(F.plhg));   colorbar;   title('Phase locked high gamma'); 
    subplot(pn,1,2);    imagesc(vertcat(F.spamp));  colorbar;   title('Ictal spike amplitude'); 
    subplot(pn,1,3);    imagesc(vertcat(F.spcnt), [0 0.15]);    colorbar;   title('Ictal spike frequency');
    subplot(pn,1,4);    imagesc(vertcat(F.lopow), [0 1000]);    colorbar;   title('LFP broadband power'); 
    subplot(pn,1,5);    imagesc(vertcat(F.hipow), [0 20]);      colorbar;   title('High gamma power'); 
    subplot(pn,1,6);    imagesc(Fcor(1).dat, [.8 1]);           colorbar,   title('LFP correlation');
    subplot(pn,1,7);    imagesc(Fcor(2).dat, [0 .6]);           colorbar,   title('High gamma correlation'); 
end



