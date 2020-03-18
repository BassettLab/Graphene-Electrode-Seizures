function grp_plot_tsrs(C, IMG, trange, fastrange, Fs, dolong)
figure
% Plot example time series
%--------------------------------------------------------------------------
srtmin = 1.75; 
stpmin = 3;
shortr = [srtmin*60*Fs:100:stpmin*60*Fs];
shorti = [srtmin*60*10:stpmin*60*10]; 

longr  = fastrange(1):100:fastrange(end); 
longi  = trange;  

if dolong  
    thisr = longr; 
    thisi = longi; 
else
    thisr = shortr;
    thisi = shorti;
end

seq = flip(cbrewer('div', 'PuOr', 6)); 
seq = seq([1 2 5 6],:); 
gry = [.5 .5 .5]; 
cols = [gry; gry; gry; gry; ...
        seq(1,:); seq(2,:); seq(3,:); seq(4,:); ...
        gry; gry; gry; gry; ...
        gry; gry; gry; gry];
    
for c = 1:length(C)
    subplot(2,1,1)
        if dolong 
            t = linspace(0,(thisr(end)-thisr(1))/5000/60,length(thisr)); 
            for s = [1 length(shortr)] 
                subplot(2,1,1)
%                 [val ind] = min(abs(longr-shortr(s)));
%                 plot([t(ind),t(ind)], [-16000, 0], 'k');  hold on    
            end
            
        else,       t = linspace(srtmin,stpmin,length(thisr));          end
        plot(t,C(c).filt{1}(thisr) - 1000*c, 'color', cols(c,:)), hold on
        xlim([-Inf Inf]); 
%         
    subplot(2,1,2)
        if dolong,  t = linspace(0,(thisi(end)-thisi(1))/10/60,length(thisi)); 
        else,       t = linspace(srtmin,stpmin,length(thisi));          end
        plot(t, IMG(c).ampl(thisi), 'color', cols(c,:)), hold on
        xlim([-Inf Inf]); 
        
end

if dolong
set(gcf, 'Position', [50, 300, 1900, 800])
% print(gcf, '-dpng', '~/Desktop/long_seizure', '-r500')
end
