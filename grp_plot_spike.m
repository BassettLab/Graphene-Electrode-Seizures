function grp_plot_spike(C, H, wsort, fastrange, Fs,c)

ts      = 1:50:4500; 
cols    = cbrewer('div', 'Spectral', length(ts));
ccols   = cbrewer('qual', 'Dark2', 12); 

clear std
spikes = [];
[val, Hclass] = max(H(wsort,:));
spikeH = [];
spike = []; 

for ti = 1:length(ts)
    disp(['Running segment ' num2str(ti) ' of ' num2str(length(ts))]); 
    t           = ts(ti); 
    tr          = fastrange(t): (fastrange(t)+10*Fs); 
    dt          = C(c).filt{1}(tr); 
    [pks locs]  = findpeaks(-dt);
    locs        = locs(find(pks > mean(dt) + 5*std(dt)));

    for l = 1:length(locs)
        srt = fix(locs(l) - 0.5*Fs);
        stp = fix(locs(l) + 0.5*Fs); 
        if srt > 0 && stp <= length(dt) 
            spike = dt(srt : stp);
            base  = dt((locs(l)-500):(locs(l)-100)); 
            spike = spike - mean(base);
            if ~isempty(spikes), spikes = vertcat(spikes, spike);
            else,                spikes = spike;
            end
            
            spikeH = [spikeH, Hclass(t)];
            
        end
        
    end
end

cols = flip(cbrewer('div', 'Spectral', size(spikes,1)));
cols = cbrewer('qual', 'Dark2', 12);
cols = cols(wsort,:); 

for s = 1:size(spikes,1)
    plot(spikes(s,:) + 0*s, 'color', cols(spikeH(s),:)),   hold on
    xlim([2000,4000]);
end
