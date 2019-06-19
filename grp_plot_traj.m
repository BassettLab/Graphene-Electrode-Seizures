function grp_plot_traj(F,trange,H,wsort)
figure

% Plot trajectories
%==========================================================================
thisset = {'Extent of seizure'};
% cfeat   = grp_featmatmaker(trange, thisset, EPH, IMG, fcE); 
% cols    = log(mean(cfeat)); 
colormap(flip(cbrewer('div', 'Spectral', 200))); 
if ~exist('binm'), load([F.analysis filesep 'Wave map' filesep 'Binary.mat']); end
cols = (sum(sum(binm,2),3)); 
cols = cols(trange,:); 

P       = {[1,2,3],[4,5,6]};
for k = 1:length(P)
    p = P{k}; 
    subplot(1,length(P),k)
    scatter3(H(wsort(p(1)),:), H(wsort(p(2)),:), H(wsort(p(3)),:), [], cols, 'filled', ...
             'markerfacealpha', 1); hold on

    view(gca, 140, 20);
    title(['Trajectory color coded by ' thisset]);
    xlabel(['Factor ' num2str(p(1))]);
    ylabel(['Factor ' num2str(p(2))]); 
    zlabel(['Factor ' num2str(p(3))]); 
    
    % To check position of specific point in plot
    %----------------------------------------------------------------------
    pid = [];
    if ~isempty(pid)
    scatter3(H(wsort(p(1)),pid), H(wsort(p(2)),pid), H(wsort(p(3)),pid), 200, 'r', 'filled')
    end
    
end

set(gcf, 'Position', [100,100,1600,1000])
