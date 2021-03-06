% Graphene electrode code - runs analysis of Graphene electrode data
%==========================================================================
% This code contains analysis code to run the analysis on graphene, steps
% of this can be run by toggling on or off the following variables



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Manual Definitions Required Below %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Base folder                           <- Change to your local paths
%--------------------------------------------------------------------------
F.base = '/Volumes/GoogleDrive/My Drive/Research/1903 Mouse Calcium Data/Graphene_Seizure_Mapping';
F.ieeg = ['/Users/roschkoenig/Dropbox/Research/0002 Tools/IEEG']; 
%   ^-- only required if IEEG portal download (pullieeg) required

% Getting EEG Datafeatures              <- Define electrophys processing  
%--------------------------------------------------------------------------
pullieeg = 0;       % Downloads EEG data from IEEG portal
physfeat = 1;       % Estimates neurophysiology derived data features

% Getting Calcium Imaging Datafeatures  <- Define imaging processing
%--------------------------------------------------------------------------
estwave = 0;        % Identify seizure core and travelling wave front  
getfeat = 1;        % Translate wave dynamics to features for further analysis

% Run NMF
%--------------------------------------------------------------------------
runnmf  = 1;        % run the combined non-negative matrix factorisation

% Plotting functions - add to the vector which ones you want
%--------------------------------------------------------------------------
% This can be toggled on or off further down
% 1) grp_plot_nmf - example weights and loadings for NMF calculated above
% 2) grp_plot_traj - plots state-space trajectories 
% 3) grp_plot_spike - plots binned ictal spike averages
% 4) grp_plot_tsrs - plots example time series
% 5) grp_plot_edgemap - plots map of changing ictal wavefront map over time
% 6) grp_plot_snapshot - plots different features at different timepoints
% 7) grp_plot_physstats - plots exploratory analysis of electrophys

plots = [7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
%==========================================================================
fs          = filesep;              
F           = grp_housekeeping(F);

% Running actual analysis
%==========================================================================
% Get EEG Datafeatures
%--------------------------------------------------------------------------
if pullieeg,    C   = grp_electrophys(F);   end         % Pull data from IEEG portal 
if ~exist('C'), load([F.data fs 'Electrophysiology' fs 'Channel_Data.mat']); end
if physfeat,    [EPH, fcE] = grp_ephysfeat(C);     end  % Calculate Ephys features from data stored in C

% Get Imaging Datafeatures
%--------------------------------------------------------------------------
if estwave,     F   = grp_imaging(F);       end
if getfeat,     IMG = grp_imgfeat(C, F);    end

% Matrix decomposition
%--------------------------------------------------------------------------
if runnmf 
    [W,H,FMT,fid,fset,trange,fastrange,wsort] = grp_nmf(EPH,IMG,fcE,C(1).Fs,0,6,'all');     
end

% Plotting functions (use the if statements as toggles)
%==========================================================================
for p = plots
switch p
    case 1, grp_plot_nmf(W,H,C,trange,fastrange,fset,fid,wsort);    
    case 2, grp_plot_traj(F,trange,H,wsort);                    
    case 3, grp_plot_spike(C, H, wsort, fastrange, C(1).Fs, 9);          
    case 4, grp_plot_tsrs(C, IMG, trange, fastrange, C(1).Fs, 0);   
    case 5, grp_plot_edgemap(F, trange);     
    case 6, grp_plot_snapshot(F,C,trange, fastrange, C(1).Fs, fset, W, H, fid, wsort);
    case 7, grp_plot_physstats(C,H,C(1).Fs,fastrange,wsort,{'isi', 'spike correlation'}, 3); 
end
end

