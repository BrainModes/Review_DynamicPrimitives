%% compute simulated Kong2021 FCDs
clear
close all
clc

load('Kong2021_fMRI.mat')

fc_len = round(60/0.72); % ~60 s
regions = 68;
ts_len = 1200;
Isubdiag = find(tril(ones(regions),-1));

edges = -1.0:0.02:1.0;

ts=TC';

fcs = zeros(length(Isubdiag), length(1:ts_len-fc_len));

for ii = 1:ts_len-fc_len
    tmp = corr(ts(ii:ii+fc_len-1,1:regions));
    fcs(:,ii) = tmp(Isubdiag);
end

fcd=corr(fcs);

IsubdiagFCD = find(tril(ones(length(fcd)),-1));
figure,imagesc(fcd)
figure;histogram(fcd(IsubdiagFCD),edges)
hist_fcd = histcounts(fcd(IsubdiagFCD),edges);
save('simfcd_Kong2021.mat','hist_fcd', 'edges', 'fcd')


%%
% compute fast simulated Kong2021 power spectra
clear
close all
clc

regions = 68;



res_pxx = zeros(257,1);
c = 0;

load('Kong2021_Syn.mat')

ts = Syn';

for ii = 1:regions
    [pxx, f]=pwelch(zscore(ts(:,ii)),256/2,[],256*2,100);
    res_pxx = res_pxx + pxx;
    c = c+1;
end


res_pxx = res_pxx ./ c;

save('simpsd_fast_Kong2021.mat','res_pxx','f')


%%
% compute simulated Kong2021 fmri PSD
clear
close all
clc


res_pxx = zeros(257,1);
c = 0;

load('Kong2021_fMRI.mat')

regions = 68;
ts = TC';

for ii = 1:regions
    [pxx, f]=pwelch(zscore(ts(:,ii)),256/2,[],256*2,1/0.72);
    res_pxx = res_pxx + pxx;
    c = c+1;
end


res_pxx = res_pxx ./ c;

save('simpsd_slow_Kong2021.mat','res_pxx','f')

%%
% compute simulated Hopf fmri power spectra
clear
close all
clc

regions = 360;
res_pxx = zeros(257,1);
c = 0;

load('Hopf_optworkpnt_sim_G0.22025_a-0.008.mat')

ts = xs;

for ii = 1:regions
    [pxx, f]=pwelch(zscore(ts(:,ii)),256/2,[],256*2,1);
    res_pxx = res_pxx + pxx;
    c = c+1;
end

res_pxx = res_pxx ./ c;

save('simpsd_slow_Deco2017.mat','res_pxx','f')

%% compute simulated Hopf FCD
clear
close all
clc

load('Hopf_optworkpnt_sim_G0.22025_a-0.008.mat')

fc_len = 60; % ~60 s
regions = 360;
ts_len = 1200;
Isubdiag = find(tril(ones(regions),-1));

edges = -1.0:0.02:1.0;

ts=squeeze(bold_data);
ts=ts(end-ts_len+1:end,:);

fcs = zeros(length(Isubdiag), length(1:ts_len-fc_len));

for ii = 1:ts_len-fc_len
    tmp = corr(ts(ii:ii+fc_len-1,1:regions));
    fcs(:,ii) = tmp(Isubdiag);
end

fcd=corr(fcs);

IsubdiagFCD = find(tril(ones(length(fcd)),-1));
hist_fcd = histcounts(fcd(IsubdiagFCD),edges);
save('simfcd_Deco2017.mat','hist_fcd', 'edges', 'fcd')
