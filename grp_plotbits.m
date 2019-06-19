%% Plot W matrix features
%--------------------------------------------------------------------------
colormap(flip(cbrewer('div', 'RdBu', '100'))); 
figure(1)
for f = 1:7
    subplot(7,1,f)
    imagesc(W(find(fid == f),wsort), [0 0.1])
end

figure(2)
colormap(flip(cbrewer('div', 'RdBu', '100'))); 
imagesc(H(wsort,:))


