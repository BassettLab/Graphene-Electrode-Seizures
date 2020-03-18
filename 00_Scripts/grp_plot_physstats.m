function grp_plot_physstats(C,H,Fs,fastrange,wsort,toplot,win)

% Define which plots should be generated in the 'toplot' variable
% This should be a cell array containing any of the following
% 'spikes', 'isi', 'spike correlation', 'lags', 'all'
% The 'win' variable defines the sliding window size in seconds
%
% Examples: 
% grp_plot_physstats(C,H,Fs,fastrange,wsort,{'isi', 'spike correlation'},3); 
% grp_plot_physstats(C,H,Fs,fastrange,wsort,'all',3); 

% Transform toplot to cell string of individual plots
%--------------------------------------------------------------------------
toplot = cellstr(toplot); 
if ~isempty(find(strcmp(toplot, 'all')))
    toplot = {'spikes', 'isi', 'spike correlation', 'lags'}; 
end

ts      = [500, 950, 1250, 2200, 3500, 4500]; 
plotall = 0;
zcutoff = 3.9;    % z-score cut off 
pkwidth = 130;  % minimum peak width to be considered 
mincnt  = 3;    % minimum peaks per window to be counted  
pthrld  = 0.05; 

%==========================================================================
% Plot example spikes extracted from around peak factor times
%==========================================================================
% Find 10s segments around the peak ts
%--------------------------------------------------------------------------
if ~isempty(find(strcmp(toplot, 'spikes')))
figure(1)
segrange = [];
for t = 1:length(ts)
    segrange = [segrange;fastrange(ts(t))-5*Fs : fastrange(ts(t))+5*Fs-1];
    subplot(6,1,t)
    plot(C(2).filt{1}(segrange(t,:)))
end

% Calculate baseline z-scoring
%--------------------------------------------------------------------------
c       = 1;
mbase   = mean(C(c).filt{1}(segrange(1,:)));
sdbase  = std(C(c).filt{1}(segrange(1,:))); 

% Plot example IEDs
%--------------------------------------------------------------------------
for t = 1:length(ts)
    
    % Pull relevant data segment, normalise and detect peaks
    %----------------------------------------------------------------------
    dat = (C(c).filt{1}(segrange(t,:)) - mbase) / sdbase;
    [val loc]   = findpeaks(-dat, 'minpeakwidth', pkwidth);
    id          = intersect(intersect(find(val > zcutoff), find(loc>1000)), find(loc<49000));
    
    % Plot actual data
    %--------------------------------------------------------------------------
    clear ds
    cl = cbrewer('qual', 'Dark2', 12);
    for d = 1:length(id)
        subplot(1,6,t)
        ds(:,d) = loc(id(d))-floor(0.05*Fs) : ...
                  loc(id(d))+floor(0.15*Fs)-1;
        plot(dat(ds(:,d)), 'color', [cl(wsort(t),:) 0.2]), hold on
        ylim([-35,10])
    end
    plot(mean(dat(ds)'), 'color', cl(wsort(t),:), 'linewidth', 1)
end
end

%==========================================================================
% Quantify complementary neurophysiology features using sliding window
%==========================================================================
if ~isempty(find(strcmp(toplot,'isi'))) || ...
   ~isempty(find(strcmp(toplot, 'spike correlation')))

[val loc] = min(abs(fastrange - (fastrange(1)+Fs*(win/2))));   first = loc + 1; 
[val loc] = min(abs(fastrange - (fastrange(end)-Fs*(win/2)))); last  = loc - 1; 

for c = 1:16
% Prepare from baseline recording for z-scoring
%--------------------------------------------------------------------------
mbase   = mean(C(c).filt{1}(fastrange(1):fastrange(200)));
sdbase  = std(C(c).filt{1}(fastrange(1):fastrange(200))); 

i       = 0; 
clear mtrough isi fact spcorr

for k = first:win*10:last
    i   = i + 1;
    seg = fastrange(k-(win*10/2)): fastrange(k+(win*10/2))-1; 
    dat = C(c).filt{1}(seg);
    dat = (dat - mbase) / sdbase;
    [val loc]   = findpeaks(-dat, 'minpeakwidth', pkwidth);
    id          = find(val > zcutoff);
    
    % Calculate amplitude and inter-stim interval for the time window
    %----------------------------------------------------------------------
    mtrough(i) = mean(dat(loc(id)));
    if length(id) > mincnt, isi(i) = mean(diff(loc(id)));
    else                    isi(i) = NaN;   end 
    
    % Calculate between spike correlation for the time window
    %----------------------------------------------------------------------
    spks  = [];
    for j = 1:length(id)
        try spks = [spks; dat(loc(id(j))-0.05*Fs : loc(id(j))+0.15*Fs)]; end
    end
    if size(spks,1) > mincnt,   spcorr(i) = mean(mean(corr(spks')));
    else,                       spcorr(i) = NaN;    end

    % Find most closely associated factor
    %----------------------------------------------------------------------
    [val loc] = min(abs(fastrange-seg(1)));     slostart = loc; 
    [val loc] = min(abs(fastrange-seg(end)));   slostop  = loc; 
    
    h = H(:,slostart:slostop);
    [val, loc] = max(mean(h'));     fact(i) = loc; 
    
end

C(c).stats.spcorr = spcorr;
C(c).stats.ampl   = mtrough; 
C(c).stats.fact   = fact;
C(c).stats.isi    = isi; 
end

% Calculate group means
%--------------------------------------------------------------------------
clear amp isi spcorr
cl = cbrewer('qual', 'Dark2', 12);
cl = cl((C(1).stats.fact),:); 
clear spcorr isi
for c = 1:length(C)
    spcorr(c,:) = C(c).stats.spcorr;
    isi(c,:) = C(c).stats.isi;
    amp(c,:) = C(c).stats.ampl;
end
spcorr = nanmean(spcorr);  isi = nanmean(isi);  amp = nanmean(amp);
end

%==========================================================================
% Calculate one group mean lag matrix
%==========================================================================
if ~isempty(find(strcmp(toplot,'lags')))
Dat = [];
for c = 1:length(C), Dat = [Dat; C(c).filt{1}];  end
mDat = mean(Dat); 

% Prepare from baseline recording for z-scoring
%--------------------------------------------------------------------------
mbase   = mean(mDat(fastrange(1):fastrange(200)));
sdbase  = std(mDat(fastrange(1):fastrange(200))); 

% Loop through time windows
%--------------------------------------------------------------------------
clear lagsum
i = 0;
for k = first:win*10:last
    
    % Extract spike times from mean trace
    %----------------------------------------------------------------------
    i   = i + 1;
    seg = fastrange(k-(win*10/2)): fastrange(k+(win*10/2))-1; 
    dat = mDat(seg);
    dat = (dat - mbase) / sdbase;
    [val loc]   = findpeaks(-dat, 'minpeakwidth', pkwidth);
    id          = find(val > zcutoff);
    
    % Calculate crosscorrelation between channels and find lags 
    %----------------------------------------------------------------------
    lagmat = [];
    for j = 1:length(id)
    if loc(id(j)) > 0.05*Fs
        spkseg          = seg(loc(id(j)))-0.05*Fs : seg(loc(id(j)))+0.1*Fs;
        spks            = Dat(:,spkseg);
        [cmat, lags]    = xcorr(spks');
        [val, lag]      = max(cmat); 
        lagmat(j,:,:)   = reshape(lags(lag), [length(C), length(C)]); 
    end    
    end
    
    % Average lags across spikes, and sum total lags across channels
    %----------------------------------------------------------------------
    if ~isempty(lagmat)
        lagmat = squeeze(mean(lagmat)); 
        lagsum(i) = sum(sum(abs(triu(lagmat))))/(length(C)-1)^2;
    else 
        lagsum(i) = 0;
    end
end
end

%==========================================================================
% Plot results
%==========================================================================
cl      = cbrewer('qual', 'Dark2', 12);
pltid   = find(strcmp(toplot,'isi')+strcmp(toplot,'spike correlation')+strcmp(toplot,'lags')); 
plts    = toplot(pltid); 

for f = 1:6
    figure(2)
    if f > 1 
        if ~isempty(find(strcmp(toplot, 'isi'))),                 oldisi = isi(id);   end
        if ~isempty(find(strcmp(toplot, 'spike correlation'))),   oldspcorr = spcorr(id); end
        if ~isempty(find(strcmp(toplot, 'lags'))),                oldlagsum = lagsum(id); end
    end
    
    id = find(C(1).stats.fact == wsort(f)); 
    ns = length(id); 
    
    if ~isempty(find(strcmp(toplot, 'isi')))
    subplot(length(pltid),1,find(strcmp(plts,'isi')))
    jit = randn(ns,1)/10;
    scatter(ones(ns,1)*f+jit, (1./isi(id)), [], cl(wsort(f),:), 'filled'), hold on
    plot([f-0.2 f+0.2], [nanmedian(1./isi(id)) nanmedian(1./isi(id))], 'k', 'linewidth', 1); 
    xlabel('Factors'), ylabel('Frequency of ictal discharges'); 
    xlim([0 7])
    end
    
    
    if ~isempty(find(strcmp(toplot, 'spike correlation')))
    subplot(length(pltid),1,find(strcmp(plts,'spike correlation')))
    jit = randn(ns,1)/10;
    scatter(ones(ns,1)*f+jit, spcorr(id), [], cl(wsort(f),:), 'filled'), hold on
    plot([f-0.2 f+0.2], [nanmedian(spcorr(id)) nanmedian(spcorr(id))], 'k', 'linewidth', 1); 
    xlabel('Factors'), ylabel('Stereotypy of ictal discharges'); 
    xlim([0 7])
    end
    
    if ~isempty(find(strcmp(toplot, 'lags')))
    subplot(length(pltid),1,find(strcmp(plts,'lags')))
    jit = randn(ns,1)/10;
    scatter(ones(ns,1)*f+jit, lagsum(id), [], cl(wsort(f),:), 'filled'), hold on
    plot([f-0.2 f+0.2], [nanmedian(lagsum(id)) nanmedian(lagsum(id))], 'k', 'linewidth', 1); 
    xlabel('Factors'), ylabel('Sum of between-channel lag times'); 
    xlim([0 7])
    end
    
    % Simple minded stats testing
    %----------------------------------------------------------------------
    if f > 1, try
        if ranksum(isi(id),oldisi) < pthrld && ~isempty(find(strcmp(toplot, 'isi')))
            subplot(length(pltid),1,find(strcmp(plts,'isi')))
            yval = prctile([1./isi(id),1./oldisi], 90);
            text(f-0.55, yval, '*',  'fontweight', 'bold', 'fontsize', 16, ...
                 'horizontalalignment', 'center')
        end
        if ranksum(spcorr(id),oldspcorr) < pthrld && ~isempty(find(strcmp(toplot, 'spike correlation')))
            subplot(length(pltid),1,find(strcmp(plts,'spike correlation')))
            yval = prctile([spcorr(id),oldspcorr], 90);
            text(f-0.55, yval, '*', 'fontweight', 'bold', 'fontsize', 16, ...
                 'horizontalalignment', 'center')
        end
        
        if ~isempty(find(strcmp(toplot, 'lags'))) && ranksum(lagsum(id),oldlagsum) < pthrld
            subplot(length(pltid),1,find(strcmp(plts,'lags')))
            yval = prctile([lagsum(id),oldlagsum], 90);
            text(f-0.55, yval, '*', 'fontweight', 'bold', 'fontsize', 16, ...
                 'horizontalalignment', 'center')
        end
        end, end
    
    % Plot it all together
    %----------------------------------------------------------------------
    if plotall, figure(2)
    scatter3(log(1./isi(id)), spcorr(id), lagsum(id), [], cl(wsort(f),:), 'filled'); hold on
    end
end
