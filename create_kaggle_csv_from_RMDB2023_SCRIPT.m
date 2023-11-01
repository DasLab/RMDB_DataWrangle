%% create_kaggle_csv_from_RMDB2023_SCRIPT.m
filedir = 'data_sets/published_rdat_10_29_2023';
d = dir([filedir,'/*.rdat']);
for i =  104; %810:length(d)
    fprintf('Reading in %d of %d...\n',i,length(d));
    dataset_name{i} = d(i).name;
    r{i} = read_rdat_file(sprintf('%s/%s',filedir,d(i).name));
end

%%
dataset_name = cell(length(d),1);
modifier = cell(length(d),1);
num_profiles = zeros(length(d),1);
num_res = zeros(length(d),1);
num_vals = zeros(length(d),1);
for i = 1:length(r)
    dataset_name{i} = d(i).name;
    modifier{i} = strip(get_tag(r{i}.annotations,'modifier'));
    if isempty(modifier{i}) | strcmp(modifier,' '); modifier{i} = ''; end;
    num_profiles(i) = size(r{i}.reactivity,2);
    num_res(i) = size(r{i}.reactivity,1);
    num_vals(i) = num_profiles(i) * num_res(i);
end
t = table(dataset_name,modifier,num_profiles,num_res,num_vals)


%%
mods = unique(modifier);
fprintf('\n\n%20s %5s %9s %9s\n','modifier','Ndata','num_prof','num_res')
for i = 1:length(mods)
    idx = find(strcmp(modifier,mods{i}));
    numvals_for_mod(i) = sum(num_vals(idx));
end
[~,mod_sort_idx] = sort(numvals_for_mod,'descend');

for n = 1:length(mods)
    i = mod_sort_idx(n);
    idx = find(strcmp(modifier,mods{i}));
    fprintf('%20s %5d %9d %9d\n',mods{i},length(idx),sum(num_profiles(idx)),sum(num_vals(idx)));
end

fprintf('%20s %5d %9d %9d\n','TOTAL',length(d),sum(num_profiles),sum(num_vals));

%% BAD RDATs
bad_rdats = {'HC16M2R_1M7_0001.rdat','HC16M2R_1M7_0002.rdat','HC16M2R_1M7_0003.rdat'}; %,'ETERNA_R82_0000.rdat','ETERNA_R82_0001.rdat','ETERNA_R83_0000.rdat','ETERNA_R83_0001.rdat'}; % still unresolved.
for i = 1:length(dataset_name)
    ok(i) = ~any(contains(bad_rdats,dataset_name{i}));
end

%% Just grab 1M7 for initial tests.
mod = '1M7';
LENGTH_CUTOFF = 457; % max size in Ribonanza test set
idx = find(ok' & strcmp(modifier,mod) & num_res <= LENGTH_CUTOFF); % there are some long ones, including HoxA9, skip those
outdir_csv = 'data_sets/published_rdat_10_29_2023_CSV/rmdb_v1/1M7/';
for i = idx'
    outfile = [outdir_csv,'/train_',dataset_name{i},'.csv'];
    condition = mod;
    experiment_type_out = mod;
    dataset_name_out = dataset_name{i};
    d = struct();
    d.conditions = {condition};
    [r_norm,d.sequences,d.BLANK_OUT5,d.BLANK_OUT3,r_norm_err] = get_r_norm_from_rdat( r{i} );
    nprof = size(r_norm,1);
    output_idx = [1:nprof];

    % special treatment of Lucks data.
    datatype = get_tag(r{i},'datatype');
    if ~isempty(datatype) & ~isempty(datatype{1}); 
        r_norm_err = nan + 0*r_norm;
        output_idx = find(contains(datatype,'REACTIVITY:')); 
        rel_err = 1./sqrt(r_norm(output_idx+1,:)); % actual counts for 1M7 are in line after RHO values
        r_norm_err(output_idx,:) = r_norm(output_idx,:) .* rel_err;
    end;

    [d.r_norm, d.r_norm_err] = normalize_reactivity(r_norm,r_norm_err,output_idx,d.BLANK_OUT5, d.BLANK_OUT3, d.conditions );
    if isempty(d.r_norm_err); d.r_norm_err = NaN + 0*d.r_norm;end
    d.reads = NaN*ones(nprof,1);
    d.signal_to_noise = NaN*ones(nprof,1);
    for k = 1:size(d.r_norm,1)
        d.signal_to_noise(k) = ubr_estimate_signal_to_noise_ratio(d.r_norm(k,:)',d.r_norm_err(k,:)');
    end
    output_kaggle_csv(outfile, d, output_idx, condition, experiment_type_out, dataset_name_out);
end


%% Prepare final table
infiles = dir( [outdir_csv,'/train*.csv'] );
t_train = concatenate_tables(infiles);
t_train = renamevars(t_train,'id','sequence_id')
check_dataset_stats(t_train);

%%
data_v1_version = 'v1.1.0'
outdir = sprintf('data_sets/published_rdat_10_29_2023_CSV/rmdb_v1/%s/',data_v1_version);
if ~exist(outdir,'dir'); mkdir(outdir); end;
outfile = sprintf('%s/rmdb_data.%s.csv',outdir,data_v1_version);
fprintf('Outputting %d rows to %s.\n',height(t_train),outfile);
writetable(t_train,outfile);
strip_nan(outfile);

%%
histogram(t_train.signal_to_noise,[0:0.2:100])
xlabel('signal_to_noise','Interpreter','none');
ylabel('counts');
title(outfile,'interpreter','none');
fprintf('Number of profiles with SN_filter=1: %d/%d (%.1f %%)\n',... ...
    sum(t_train.SN_filter),height(t_train),100*sum(t_train.SN_filter)/height(t_train));

%%


%%
names = t_train.Properties.VariableNames;
cols = find(contains(names,'reactivity_') & ~contains(names,'reactivity_error'));
imagesc(cell2mat(table2cell(t_train(:,cols))),[-1,1]);
colormap([0.7 0.7 0.7; redwhiteblue(-0.1,0.1)])
colorbar
set(gcf,'color','white');
xlabel('Position');ylabel('Profile')
