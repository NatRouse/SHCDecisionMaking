close all; clear all; clc

matchcount_red = struct2array(load('matchcount_red.mat','match_count'));
matchcount_green = struct2array(load('matchcount_green.mat','match_count'));
matchcount_blueblack = struct2array(load('matchcount_blueblack.mat','match_count'));

act = [.1 .2 .3 .4 .5 .6 .7 .8 .9];
duration = [1 10 20 30 40 50]*0.1;
biasval = [1e-5 1e-4 1e-3 1e-2 1e-1];


figure(1)
tiledlayout(3,6)
for i = 1:size(matchcount_red,3)
    nexttile
    heatmap(biasval,act,squeeze(matchcount_red(:,:,i))', ...
        'ColorLimits',[0 100],'CellLabelColor','none'),...
        % 'XDisplayLabels',nan*ones(1,6),'YDisplayLabels',nan*ones(1,9))
    colorbar off
end

for i = 1:size(matchcount_blueblack,3)
    nexttile
    heatmap(biasval,act,squeeze(matchcount_blueblack(:,:,i))', ...
        'ColorLimits',[0 100],'CellLabelColor','none'),...
        % 'XDisplayLabels',nan*ones(1,6),'YDisplayLabels',nan*ones(1,9))
    if mod(i,size(matchcount_blueblack,3))~=0
        colorbar off
    end
end

for i = 1:size(matchcount_green,3)
    nexttile
    heatmap(biasval,act,squeeze(matchcount_green(:,:,i))', ...
        'ColorLimits',[0 100],'CellLabelColor','none'),...
        % 'XDisplayLabels',nan*ones(1,6),'YDisplayLabels',nan*ones(1,9))
    colorbar off
end








