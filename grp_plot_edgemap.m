function grp_plot_edgemap(F, trange) 

fs = filesep;

% Load detected edges
%--------------------------------------------------------------------------
if ~exist('edgm'), load([F.analysis fs 'Wave map' fs 'Edge.mat']); end


% Plot edges on top of each other
%--------------------------------------------------------------------------
cols = flip(cbrewer('div', 'RdYlBu', length(trange))); 

for ei = 1:length(trange)           % Run through each time window
    e = trange(ei);                 % e encodes actual time frame 
    
    if mod(ei, 1000) == 0,   disp(['Iteration ' num2str(ei)]);    end
    if mod(ei,50) == 0              % Only plot ever n-th edge 
        
    edgi = find(squeeze(edgm(e,:,:)));  
    [r, c] = ind2sub([size(edgm,2), size(edgm,3)], edgi);
    
    if ~isempty(r)                  % Plot only when you have found an edge
        scatter(c,r, 40, cols(ei,:), 'filled', 'markerfacealpha', 0.5);   hold on
        set(gca, 'ydir', 'reverse'); 
        xlim([0 512]); 
        ylim([0 512]); 
        axis square
        axis off
    end
    end
end

% % Plot locations of electrodes on top of the plot
% %--------------------------------------------------------------------------
% for c = 1:length(C)
%     scatter(C(c).x, C(c).y, 200, 'k', 'filled')
%     text(C(c).x+15, C(c).y + 5, num2str(c))
% end