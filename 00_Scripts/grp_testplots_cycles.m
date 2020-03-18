%% Set up
t_sec = ((472:5000)-472)/10;
sampl = (472*500:5000*500);

srt = sampl(1);
seg = fix(length(sampl)/200);
cols = flip(cbrewer('div', 'Spectral', 200)); 

a = C(13).filt{1};
b = C(8).filt{1}; 

%% Go through loadings and find segments 
mis = ones(1, size(H,2)); 
for h = 1:size(H,2)
    [mv mi] = max(H(:,h)); 
    mis(1,h) = mi; 
end
%%
    
for s = 1:fix(length(sampl)/seg)
%     ts = seg+(s-1)*20000;
%     thislc = intersect(plc(plc >= ts(1)), plc(plc <= ts(end)));
%     ta = zeros(length(thislc), 1000);
%     tb = zeros(length(thislc), 1000); 
%     for t = 1:length(thislc)
%         ta(t,:) = a(thislc(t)-299:thislc(t)+700);
%         tb(t,:) = b(thislc(t)-299:thislc(t)+700); 
%     end
%     plot3(ones(1,1000)*s, mean(ta), mean(tb), 'color', cols(s,:)); hold on
%     
    thisseg = (srt:srt+seg-1)+(s-1)*seg;
    plot3(b(thisseg), ones(1,seg)*s, a(thisseg), 'color', cols(s,:))
    hold on
    
end
% grid on

%% --------------------------------------------------------------------------

clear rc fc

fc = fcE(1).dat; 
for r = 1:length(IMG),  rc(:,r) = IMG(r).ampl;  end

plid = 472:5000; 

mrc = mean(rc(plid,:)');
mfc = mean(fc(:,plid));

subplot(3,1,1)
scatter(mfc, mrc, 10, 'filled', 'markerfacealpha', 0.2);

% within and between recruitment correlation
%--------------------------------------------------------------------------
rc = (rc > prctile(rc(:), 40));
fs = mean(fc(:,plid));
rs = sum(rc(plid,:)'); 
rs = mrc; 

clear mf
rrange = min(rs):max(rs);
for r = 1:length(rrange)
    thisi = find(rs == rrange(r));
    if isempty(thisi), disp(r); end
    mf(r) = mean((fs(thisi)));
end

colormap(flip(cbrewer('div', 'Spectral', 100))); 
subplot(3,1,2)
scatter(rs+0.2*randn(length(rs),1)', log(fs), 10, plid, 'filled', 'markerfacealpha', 0.2); hold on
plot(rrange, mf)

subplot(3,1,3)
scatter(rs+0.2*randn(length(rs),1)', (fs), 10, plid, 'filled'); hold on
plot(rrange, mf)
xlim([100, Inf])

figure 
scatter3(rs, (fs), plid)