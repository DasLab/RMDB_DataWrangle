%% create_kaggle_csv_from_RMDB2023_SCRIPT.m
filedir = 'data_sets/published_rdat_11_09_2023';
rdats = dir([filedir,'/*.rdat']);
for i = 1:length(rdats); % %758; % 21 %%247; %find(contains(dataset_name,'TTYH2'))'
    fprintf('Reading in %d of %d...\n',i,length(rdats));
    dataset_name{i} = rdats(i).name;
    r{i} = read_rdat_file(sprintf('%s/%s',filedir,rdats(i).name));
end

%%
dataset_name = cell(length(rdats),1);
modifier = cell(length(rdats),1);
num_profiles = zeros(length(rdats),1);
num_res = zeros(length(rdats),1);
num_vals = zeros(length(rdats),1);
% TITR and MGTI were Joe's experiments -- DMS on miniTTR's
% ETFMN and ETE3D date to 2014 when we still used NMIA, not 1M7
filename_tags = {'DMS','MCA',  '1M7','NOMOD','NMD',  'MGPH',       'MG50'      ,'PH10',        '50C','MGTI','TITR','GLX',    'SYN41','ETFMN','ETE3D'};
mods_for_tags = {'DMS','MOHCA','1M7','nomod','nomod','deg_Mg_pH10','deg_Mg_50C','deg_pH10','deg_50C','DMS', 'DMS' ,'glyoxal','nomod','NMIA' ,'NMIA'};
for i = 1:length(r)
    dataset_name{i} = rdats(i).name;
    modifier{i} = strip(get_tag(r{i}.annotations,'modifier'));
    if isempty(modifier{i}) | strcmp(modifier,' '); 
        modifier{i} = ''; 
        for k = 1:length(filename_tags)
            if contains(dataset_name{i},filename_tags{k})
                modifier{i} = mods_for_tags{k};
                break;
            end
        end
    end;
    modifier{i} = strrep( modifier{i},'none','nomod'); 
    modifier{i} = strrep( modifier{i},'None','nomod'); 
    if strcmp(modifier{i},'SHAPE')
        if contains(dataset_name{i},'1M7')
            modifier{i} = '1M7';
        else
            modifier{i} = 'NMIA';
        end
    end
    if any(contains(r{i}.annotations,'cotranscription')) modifier{i} = [modifier{i},'_cotx']; end
    if any(contains(r{i}.annotations,'TECprobe-ML')) modifier{i} = [modifier{i},'_cotx']; end
    if any(contains(r{i}.annotations,'MOHCA')) modifier{i} = 'MOHCA'; end;
    if length(r{i}.sequences)>10 & all(contains(r{i}.sequences(2:end),'X')) & ~contains(r{i}.sequence,'X'); modifier{i} = 'DMS_M2_seq'; end;
    seq_lens = cellfun(@length,r{i}.sequences);
    if all( (seq_lens(2:end)-seq_lens(1:end-1))>0 ) & size(r{i}.reactivity,2) > 10 & ~strcmp(modifier{i},'MOHCA')
        fprintf('%s modifier:%s experiment:%s\n',dataset_name{i},modifier{i},get_tag(r{i}.annotations,'experiment'));
    end
    num_profiles(i) = size(r{i}.reactivity,2);
    num_res(i) = size(r{i}.reactivity,1);
    num_vals(i) = num_profiles(i) * num_res(i);
    seq_length(i) = length(r{i}.sequence);
end
%t = table(dataset_name,modifier,num_profiles,num_res,num_vals)

% how many RDATs are singly devoted each kind of modification? 
% Data sets with mixed modifiers will not be categorized here.
mods = unique(modifier); numvals_for_mod = [];
for i = 1:length(mods)
    idx = find(strcmp(modifier,mods{i}));
    numvals_for_mod(i) = sum(num_vals(idx));
end
[~,mod_sort_idx] = sort(numvals_for_mod,'descend');

fprintf('\n\n%20s %5s %9s %9s\n','modifier','Ndata','num_prof','num_res')
for n = 1:length(mods)
    i = mod_sort_idx(n);
    idx = find(strcmp(modifier,mods{i}));
    fprintf('%20s %5d %9d %9d\n',mods{i},length(idx),sum(num_profiles(idx)),sum(num_vals(idx)));
end

fprintf('%20s %5d %9d %9d\n','TOTAL',length(modifier),sum(num_profiles),sum(num_vals));

%% how many profiles have each kind of modification? Look profile by profile for modifier DATA_ANNOTATION tags.
mods_per_profile = {}; num_vals_per_profile = [];
for i = 1:length(r)
    if length(modifier{i})>1 & ~strcmp(modifier{i},' ');
        data_annotations_mods = repmat({modifier{i}},1,length(r{i}.data_annotations));
        check_mods = get_tag(r{i}.data_annotations,'modifier');
        mod_specified = find(cellfun(@length,check_mods));
        data_annotations_mods(mod_specified) = check_mods(mod_specified);
    else
        data_annotations_mods = get_tag(r{i}.data_annotations,'modifier');
        modlen = cellfun(@length,data_annotations_mods);
        if max(modlen) == 0
            fprintf('All blank modifiers (%3d): %s\n',i,dataset_name{i});
        elseif any(modlen==0)
            %fprintf('Some blank modifiers (%3d): %s\n',i,dataset_name{i});
            data_annotations_mods( find(modlen==0)) = {'nomod'};
        elseif any(modlen==1)
            %fprintf('Some single-length-name mystery modifiers (%3d): %s\n',i,dataset_name{i});
            %data_annotations_mods( find(modlen==1)) = {'nomod'};
        end
    end
    data_annotations_mods = strrep(data_annotations_mods,'Nomod','nomod');
    data_annotations_mods = strrep(data_annotations_mods,'none','nomod');
    data_annotations_mods = strrep(data_annotations_mods,'Glyoxal','glyoxal');
    data_annotations_mods = strrep(data_annotations_mods,'Terbium','terbium');
    data_annotations_mods = strrep(data_annotations_mods,'SHAPE','NMIA');
    if any(contains(r{i}.annotations,'cotranscription')) | any(contains(r{i}.annotations,'TECprobe-ML')) 
        original_modifier = strip(get_tag(r{i}.annotations,'modifier'));
        data_annotations_mods = strrep(data_annotations_mods,original_modifier,[original_modifier,'_cotx']); 
    end
    data_annotations_mods(find(contains(data_annotations_mods,'HRF'))) = {'hydroxyl_radical'};
    data_annotations_mods(find(contains(data_annotations_mods,'CMC'))) = {'CMCT'};
    data_annotations_mods(find(contains(data_annotations_mods,'DMS:'))) = {'DMS'};
    data_annotations_mods(find(contains(data_annotations_mods,'1M7:'))) = {'1M7'};
    mods_for_rdat{i} = data_annotations_mods;
    mods_per_profile = [mods_per_profile,data_annotations_mods];
    if any(strcmp(mods_per_profile,''))
        fprintf('Some blank modifiers (%3d): %s\n',i,dataset_name{i});
        break
    end
    num_vals_per_profile = [num_vals_per_profile, repmat(num_res(i),1,length(r{i}.data_annotations))];
end

mods = unique(mods_per_profile); 
numvals_for_mod = [];
for i = 1:length(mods)
    idx = find(strcmp(mods_per_profile,mods{i}));
    numvals_for_mod(i) = sum(num_vals_per_profile(idx));
end
[~,mod_sort_idx] = sort(numvals_for_mod,'descend');

fprintf('\n%30s %9s %9s\n','modifier','num_prof','num_res')
for n = 1:length(mod_sort_idx)
    i = mod_sort_idx(n);
    idx = find(strcmp(mods_per_profile,mods{i}));
    fprintf('%30s %9d %9d\n',mods{i},length(idx),sum(num_vals_per_profile(idx)));
end

fprintf('%30s %9d %9d\n\n\n','TOTAL',length(mods_per_profile),sum(num_vals_per_profile));

%% BAD RDATs
bad_rdats = {'HC16M2R_1M7_0002.rdat','HC16M2R_1M7_0003.rdat'}; 
for i = 1:length(dataset_name)
    ok(i) = ~any(contains(bad_rdats,dataset_name{i}));
end

%% Just grab 1M7 for initial tests.
output_mods = {'1M7','DMS','NMIA','BzCN','CMCT','DMS_M2_seq','BzCN_cotx','DMS_cotx','deg_Mg_50C','deg_Mg_pH10','deg_50C','deg_pH10'};
LENGTH_CUTOFF = 512; % max size in Ribonanza test set is 457; go up to 512 to make power of 2.
complement = struct('A','U','U','A','C','G','G','C');
for q = 1:length(output_mods)
    mod = output_mods{q};
    outdir_csv = sprintf('data_sets/rmdb_v1/%s/',mod);
    for i = 1:length(r); % 758
        if ~ok(i) continue; end;
        if (seq_length(i) > LENGTH_CUTOFF); continue; end;
        mod_idx = find(strcmp(mods_for_rdat{i},mod));
        if length(mod_idx) == 0; continue; end;
        if any(contains(r{i}.sequences,'X'))
            if contains(r{i}.sequence,'X') % may be M2-seq.
                warning(sprintf('Sequences with X in %s',dataset_name{i}));
                continue;
            end
        end
        outfile = [outdir_csv,'/train_',dataset_name{i},'.csv'];
        condition = mod;
        experiment_type_out = mod;
        dataset_name_out = dataset_name{i};
        d = struct();
        d.conditions = {condition};
        [r_norm,d.sequences,d.BLANK_OUT5,d.BLANK_OUT3,r_norm_err] = get_r_norm_from_rdat( r{i} );
        if all(isnan(r_norm(:))); continue; end; % NEIL1_DMS_0020 edge case.
        if isempty(r_norm_err) | all((r_norm_err(:)==0) | isnan(r_norm_err(:))); % rough guess for error
            r_norm_err = 0.1*mean(r_norm( ~isnan(r_norm) ))+0.2*r_norm;
        end

        output_idx = mod_idx;
        datatype = get_tag(r{i},'datatype');
        if  ~isempty(datatype) & any(strcmp(datatype,'READS')); output_idx = find( ~strcmp(datatype,'READS') ); end
        if length(output_idx) == 0; continue; end;

        % special treatment of old Lucks data -- moved into read_rdat_file.m
        [d.r_norm, d.r_norm_err] = normalize_reactivity(r_norm,r_norm_err,output_idx,d.BLANK_OUT5, d.BLANK_OUT3, d.conditions );
        if isempty(d.r_norm_err); d.r_norm_err = NaN + 0*d.r_norm;end
        nprof = size(r_norm,1);
        d.reads = NaN*ones(nprof,1);
        d.signal_to_noise = NaN*ones(nprof,1);
        for k = 1:size(d.r_norm,1)
            d.signal_to_noise(k) = ubr_estimate_signal_to_noise_ratio(d.r_norm(k,:)',d.r_norm_err(k,:)');
        end

        % Handle 'X' for M2-seq
        if any(contains(d.sequences,'X'));
            for k = 1:nprof
                seq = d.sequences{k};
                pos = strfind(seq,'X');
                assert(length(pos)<2);
                if isempty(pos); continue; end;
                seq(pos) = complement.(upper(r{i}.sequence(pos)));
                d.sequences{k} = seq;
            end
            assert(~any(contains(d.sequences,'X')));
        end

        output_kaggle_csv(outfile, d, output_idx, condition, experiment_type_out, dataset_name_out);
    end
end

%% Prepare final table
for q = 1:length(output_mods)
    mod = output_mods{q};
    outdir_csv = sprintf('data_sets/rmdb_v1/%s/',mod);
    infiles = dir( [outdir_csv,'/train*.csv'] );
    t_train_all{q} = concatenate_tables(infiles);
    %t_train_all{q} = concatenate_tables(infiles(1:min(20,length(infiles))));
end
t_train = concatenate_tables(t_train_all);
t_train = renamevars(t_train,'id','sequence_id')
t_train.sequence = strrep(upper(t_train.sequence),'T','U');
%check_dataset_stats(t_train);
clear t_train_all

%%
%data_v1_version = 'v1.0.0'; % original output with just three RDAT entries for 1M7
%data_v1_version = 'v1.1.0'; % all available 1M7 entries except HC16M2R
%data_v1_version = 'v1.2.0'; % all available 1M7, DMS, CMCT,... entries
data_v1_version = 'v1.3.0'; % enforce sequence length cutoff

outdir = sprintf('data_sets/rmdb_v1/%s/',data_v1_version);
if ~exist(outdir,'dir'); mkdir(outdir); end;
outfile = sprintf('%s/rmdb_data.%s.csv',outdir,data_v1_version);
fprintf('Outputting %d rows to %s.\n',height(t_train),outfile);
writetable(t_train,outfile);
strip_nan(outfile);

%% Check signal-to-noise
histogram(t_train.signal_to_noise,[0:0.2:100])
xlabel('signal_to_noise','Interpreter','none');
ylabel('counts');
title(outfile,'interpreter','none');
fprintf('Number of profiles with SN_filter=1: %d/%d (%.1f %%)\n',... ...
    sum(t_train.SN_filter),height(t_train),100*sum(t_train.SN_filter)/height(t_train));

%% Check sequences
for i = 1:length(t_train.sequence)
    seq = t_train.sequence{i}; %upper(t_train.sequence{i});
    seq_ok = all(seq=='A'|seq=='C'|seq=='G'|seq=='U');
    if ~seq_ok; fprintf('Problem sequence in %s: %s\n',t_train.dataset_name{i},seq); end;
end

%% Visualize on heatmap
names = t_train.Properties.VariableNames;
cols = find(contains(names,'reactivity_') & ~contains(names,'reactivity_error'));
good_idx =find(t_train.SN_filter);
%good_idx =find(t_train.SN_filter & ~strcmp(t_train.experiment_type,'M2_seq' ));
%imagesc(cell2mat(table2cell(t_train(good_idx(100:100:end),cols))),[-1,1]);
clf
chunk_size = 10000;
Nchunks = ceil(length(good_idx)/chunk_size);
for q = 1:Nchunks
    fprintf('Showing image chunk %d of %d\n',q,Nchunks);
    axes('Position',[1/(2*(Nchunks+1))+(q-1)/(Nchunks+1),0.1, 1/(Nchunks+1),0.8])
    chunk_idx = good_idx( 1+(q-1)*chunk_size: min(q*chunk_size,end));
    r_show = cell2mat(table2cell(t_train(chunk_idx,cols)));
    imagesc(r_show,[-1,1]);
    colormap([0.7 0.7 0.7; redwhiteblue(-0.1,0.1)])
    %colorbar
    set(gcf,'color','white');
    xlabel('Position'); %ylabel('Profile')
    ylim([0,chunk_size])
    %set(gca,'ytick',[]);
    set(gca,'ytick',[100:100:12000],'yticklabel',(q-1)*10000+[100:100:12000]);
end

%% Check an individual profile.
plot(cell2mat(table2cell(t_train(good_idx(110000),cols))))

%% 
clear t_train;
save workspaces/save_create_kaggle_csv_from_RMDB2023.mat % -v7.3

%% check S/N
dsname = strcat(t_train.dataset_name,{' '}, t_train.experiment_type);
[dsnames,IA,IC] = unique(dsname,'stable');
for n = 1:length(dsnames)
    idx = find(strcmp(dsname,dsnames{n}));
    nprof_per_dataset(n) = length(idx);
end

%% check S/N
dsname = strcat(t_train.dataset_name,{' '}, t_train.experiment_type);
[dsnames,IA,IC] = unique(dsname,'stable');
[~,sortidx] = sort(nprof_per_dataset);
for n = sortidx
    idx = find(strcmp(dsname,dsnames{n}));
    num_SN_filter(n)=sum(t_train.SN_filter(idx));
    num_profiles(n)=length(idx);
    fprintf('%20s %4d/%4d (%.1f%%)\n',dsnames{n},sum(t_train.SN_filter(idx)),length(idx),100*sum(t_train.SN_filter(idx))/length(idx));
end
dsnames(find(num_SN_filter==0))


