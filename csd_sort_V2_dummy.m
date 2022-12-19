%% Load CSD and other data of the animal
tic % timer starts
clear variables
close all
clc
filepath =('/Users/vishalkannan/Documents/CortXplorer/Ebru project/DATA'); % main file path
cd(filepath); 
fileinfo = dir('a*'); filenames = {fileinfo.name}'; %Select all the files that starts with 'a'
n_files = length(filenames); 
save_flag = 0; % 1 = saves the data structure into mat files 
%%
for file_id = 1:1 % change it to n_files to run all the files (now it runs animal a003 only)
file_id  % Select the animal file to run
% load(filenames{file_id});
DATA = importdata(filenames{file_id}).(filenames{file_id}(1:end-4)); % importing the animal data structure

% DATA = a005; % XXXXXX IMPORTANT XXXX % change animal name according to file id 
phase = {DATA.phase}';
phase_name = unique(phase(cellfun('isclass',phase,'char')));
phase_length = length(phase_name);

cont = {DATA.contingency}';
cr = {DATA.cr}';
CSD = {DATA.CSD}';
rt = {DATA.rt}';
ord_stim = {DATA.order_stim}';
hit_count = 0;
miss_count =0;
fa_count = 0;
correj_count = 0;
fa_plus6s = 0;
  csd_go = cell(5,1);  csd_nogo = cell(5,1);
    csd_hit = cell(5,1);csd_miss = cell(5,1);csd_fa = cell(5,1);csd_correj = cell(5,1);
%% Segregation of CSD as per hit,miss,false alarm and correct rejections

for i = 1:phase_length % loop for each phase
    idx_ss = find(strcmp(phase,phase_name{i}));
    length_ss(i) = length(idx_ss); % no.of sessions per phase
    % initialize empty csd array per phase
    csd_go_complete = [];csd_nogo_complete=[];
    csd_hit_complete = []; csd_miss_complete = [] ;
    csd_fa_complete = [] ; csd_correj_complete = [];
    for j = 1:length_ss(i) % loop for each session in a phase
        CSD_ss = CSD{idx_ss(j)};
        cont_ss = cont{idx_ss(j)};
        cr_ss = cr{idx_ss(j)};
        rt_ss = rt{idx_ss(j)};
        ord_stim_ss = ord_stim{idx_ss(j)};
        % Initialize all the counts to zero for each session
        go_count = 0;nogo_count = 0;
        hit_count=0;miss_count=0;fa_count=0;correj_count=0; fa_plus6s = 0;
        clear go_csd nogo_csd hit_csd miss_csd fa_csd correj_csd % clear the csd arrays
        for k = 1:size(CSD_ss,3) % loop for each trial in one session of each phase

            if ord_stim_ss(k) == 3 || ord_stim_ss(k) == 4 % Condition for GO stimulus
                go_count = go_count+1;
                go_csd(:,:,go_count) = CSD_ss(:,:,k);
            end

            if ord_stim_ss(k) == 34 || ord_stim_ss(k) == 33 % Condition for NOGO stimulus
                nogo_count = nogo_count+1;
                nogo_csd(:,:,nogo_count) = CSD_ss(:,:,k);
            end

            if cr_ss(k) ==1 && cont_ss(k) ==1 % hit condition
                hit_count = hit_count+1;
                hit_csd(:,:,hit_count) = CSD_ss(:,:,k);
            end

            if cr_ss(k) ==0 && cont_ss(k) ==1 % miss condition
                miss_count = miss_count+1;
                miss_csd(:,:,miss_count) = CSD_ss(:,:,k);

            end
            %             cr_ss(k) ==1
            if  rt_ss(k) <= 6 && cont_ss(k) ==-1 % false alarm condition (rt < 6s) ( using reaction time and not cr)
                fa_count = fa_count+1;
                fa_csd(:,:,fa_count) = CSD_ss(:,:,k);
            end
            if cr_ss(k) ==0 && cont_ss(k) ==-1 % correct rejection condition
                correj_count = correj_count+1;
                correj_csd(:,:,correj_count) = CSD_ss(:,:,k);
            end
%             if cr_ss(k) ==2 && cont_ss(k) ==-1 % false alarm condition (rt > 6s)
%                 fa_plus6s = fa_plus6s+1;
%                 csd_fa_plus6s(:,:,fa_plus6s) = CSD_ss(:,:,k);
%             end
        end
        if hit_count==0, hit_csd = nan(32,12400);end % dummy array
        if miss_count==0, miss_csd = nan(32,12400);end % dummy array
        if fa_count==0, fa_csd = nan(32,12400);end % dummy array
        if correj_count==0, correj_csd = nan(32,12400);end % dummy array


%         % storing the go and nogo counts per session
        go_sess_count(i,j) = go_count;
        nogo_sess_count(i,j) = nogo_count;
        hit_sess_count(i,j) = hit_count;
        miss_sess_count(i,j) = miss_count;
        fa_sess_count(i,j) = fa_count;
        correj_sess_count(i,j) = correj_count;

% Appending the csds of each session without averaging
        csd_go_complete = cat(3,csd_go_complete,go_csd);
        csd_nogo_complete = cat(3,csd_nogo_complete,nogo_csd);
        csd_hit_complete = cat(3,csd_hit_complete,hit_csd);
        csd_miss_complete = cat(3,csd_miss_complete,miss_csd);
        csd_fa_complete = cat(3,csd_fa_complete,fa_csd);
        csd_correj_complete = cat(3,csd_correj_complete,correj_csd);
    end
    sess_count(i) = j;

    % Averaging per phase
    csd_go_phaseavg(:,:,i) = nanmean(csd_go_complete,3);
    csd_nogo_phaseavg(:,:,i) = nanmean(csd_nogo_complete,3);
    csd_hit_phaseavg(:,:,i) = nanmean(csd_hit_complete,3);
    csd_miss_phaseavg(:,:,i) = nanmean(csd_miss_complete,3);
    csd_fa_phaseavg(:,:,i) = nanmean(csd_fa_complete,3);
    csd_correj_phaseavg(:,:,i) = nanmean(csd_correj_complete,3);
%     csd_fa_plus6s_phases(:,:,i) = nanmean(csd_fa_plus6s_ss,3);
  
    csd_go{i,1} = csd_go_complete;
    csd_nogo{i,1} = csd_nogo_complete;
    csd_hit{i,1} = csd_hit_complete;
    csd_miss{i,1} = csd_miss_complete;
    csd_fa{i,1} = csd_fa_complete;
    csd_correj{i,1} = csd_correj_complete;
end

% storing the go and nogo counts per phase
go_phase_count = sum(go_sess_count,2);
nogo_phase_count = sum(nogo_sess_count,2);
hit_phase_count = sum(hit_sess_count,2);
miss_phase_count = sum(miss_sess_count,2);
fa_phase_count = sum(fa_sess_count,2);
correj_phase_count= sum(correj_sess_count,2);





%% plotting
stimulus = 1; % choose from 1 to 4; 5 = all stimulus 
all_stim = 1:12400; 
stim_arr = [200:1200;3200:4200;6200:7200;9200:10200];
if stimulus <= 4, stim = stim_arr(stimulus,:);end
if stimulus == 5, stim = all_stim;end

 % choose this to include all stimulus

for i = 1:phase_length
if stimulus == 5
fig_title = sprintf('Animal : %s, Phase : %s, Stimulus : all',filenames{file_id}(1:4),phase_name{i});
else
fig_title = sprintf('Animal : %s, Phase : %s, Stimulus : %.0f',filenames{file_id}(1:4),phase_name{i},stimulus);
end
figure
sgtitle(fig_title,'fontweight', 'bold')

subplot(3,2,1) %plot for go_csd
imagesc(csd_go_phaseavg(:,stim,i))
% x_ticks = [1:1000];
% str_x = stim; % time from -1 to +5 seconds
% set(gca,'XTick',x_ticks); set(gca,'xTickLabel',str_x,'fontsize',8);
title('Go-csd','fontweight', 'bold')
colormap('jet')
colorbar
caxis([-2.5e-7 2.5e-7])

subplot(3,2,2) %plot for nogo_csd
imagesc(csd_nogo_phaseavg(:,stim,i))
title('Nogo-csd','fontweight', 'bold')
colormap('jet')
colorbar
caxis([-2.5e-7 2.5e-7])

subplot(3,2,3) %plot for hit_csd
imagesc(csd_hit_phaseavg(:,stim,i))
title('Hit','fontweight', 'bold')
colormap('jet')
colorbar
caxis([-2.5e-7 2.5e-7])

subplot(3,2,4) %plot for miss_csd
imagesc(csd_miss_phaseavg(:,stim,i))
title('Miss','fontweight', 'bold')
colormap('jet')
colorbar
caxis([-2.5e-7 2.5e-7])

subplot(3,2,5) %plot for false allarm_csd
imagesc(csd_fa_phaseavg(:,stim,i))
title('False alarm','fontweight', 'bold')
colormap('jet')
colorbar
caxis([-2.5e-7 2.5e-7])

subplot(3,2,6) %plot for correct rejction_csd
caxis([-2.5e-7 2.5e-7])
imagesc(csd_correj_phaseavg(:,stim,i))
title('Correct Rejection','fontweight', 'bold')
colormap('jet')
colorbar
caxis([-2.5e-7 2.5e-7])

end

%% Storing the csd in a data structure

clear data
cond_names = ["Go";"Nogo";"Hit";"Miss";"False alarm";"Correct rejection"];
for i = 1:phase_length
    clear data
    csd_cond = {csd_go{i};csd_nogo{i};csd_hit{i};csd_miss{i};csd_fa{i};csd_correj{i}};
for j = 1:length(cond_names)
% data(i).phase = phase_name{i};
data(j).condition = cond_names(j);
data(j).csd = csd_cond{j};
data(j).GS1 = [4;5;6;7;8;9;10;11];
data(j).SGS1 = [2,3,4,5];
data(j).layer6 = [26;27;28;29;30;31];
data(j).layer5a = [16;17;18;19;20];
data(j).layer5b = [19;20;21;22;23;24;25];
end
if save_flag == 1 
 save_filename = sprintf('%s_%s_sorted.mat',filenames{file_id}(1:end-4),phase_name{i});%animal_phase_sorted
 save(save_filename,'data','-v7.3'); % save the struct as .mat file
 dis = sprintf('Data saved in %s',save_filename);
 disp(dis);
end
end

toc % timer ends
end
% t = struct2table(data);