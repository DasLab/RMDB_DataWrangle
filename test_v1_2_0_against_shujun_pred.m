%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stress test by comparison to Shujun model predictions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_train = readtable('~/Desktop/rmdb_data.v1.2.0.csv');
%%
parquetwrite('data_sets/predictions/rmdb_data.v1.2.0.csv.parquet', t_train); % not synced to github
%%
names = t_train.Properties.VariableNames;
cols = find(contains(names,'reactivity_')& ~contains(names,'reactivity_error'));
tic; r = cell2mat(table2cell(t_train(:,cols))); toc
%%
tic
t_pred_DMS = readtable('data_sets/predictions/rmdb_pred_DMS_sequence_only.csv');
t_pred_2A3 = readtable('data_sets/predictions/rmdb_pred_2A3_sequence_only.csv');
toc
%%
tic
assert(height(t_train)==height(t_pred_DMS));
assert(height(t_train)==height(t_pred_2A3));
names = t_pred_DMS.Properties.VariableNames;
cols = find(contains(names,'reactivity_')& ~contains(names,'reactivity_error'));
r_pred_DMS = cell2mat(table2cell(t_pred_DMS(:,cols)));
r_pred_2A3 = cell2mat(table2cell(t_pred_2A3(:,cols)));
toc
%%
r_data = r(:,1:size(r_pred_DMS,2));
%% check MAE
dsname = strcat(t_train.dataset_name,{' '}, t_train.experiment_type);
[dsnames,IA,IC] = unique(dsname,'stable');
for i = 1:length(dsnames)
    if mod(i,10)==0; fprintf('Doing %d out of %d [%s]...\n',i,length(dsnames),dsnames{i});end;
    idx = find(strcmp(dsname,dsnames{i}));
    idx_in_pred = [];
    for n = 1:length(idx)
        ids = find(strcmp(t_pred_DMS.sequence_id,t_train.sequence_id{idx(n)}));
        assert(length(ids)>0);
        idx_in_pred(n) = ids(1);
    end
    data = r_data(idx,:);
    pred_DMS = r_pred_DMS(idx_in_pred,:);
    pred_2A3 = r_pred_2A3(idx_in_pred,:);

    vals = abs(pred_DMS-data); 
    pos = find(~isnan(vals));
    mae_DMS(i) = mean(vals(pos));
    if length(pos)>0
        corr_DMS(i)=corr(reshape(pred_DMS(pos),[],1),reshape(data(pos),[],1));
        corr_spearman_DMS(i)=corr(reshape(pred_DMS(pos),[],1),reshape(data(pos),[],1),'type','Spearman');
    end

    vals = abs(pred_2A3-data); 
    pos = find(~isnan(vals));
    mae_2A3(i) = mean(vals(pos));
    if length(pos)>0
        corr_2A3(i)=corr(reshape(pred_2A3(pos),[],1),reshape(data(pos),[],1));
        corr_spearman_2A3(i)=corr(reshape(pred_2A3(pos),[],1),reshape(data(pos),[],1),'type','Spearman');
    end
end

%%
nts = 'ACGU';
for i = 1:length(dsnames)
    if mod(i,10)==0; fprintf('Doing %d out of %d [%s]...\n',i,length(dsnames),dsnames{i});end;
    idx = find(strcmp(dsname,dsnames{i}));
    data = r_data(idx,:);
    for n = 1:length(nts)
        vals = [];
        for k = 1:length(idx)
            pos = strfind(t_train.sequence{idx(k)},nts(n));
            vals = [vals, r_data(idx(k),pos)];            
        end
        meanval(i,n) = mean(vals( ~isnan(vals) ));
    end
end


%%
set(gcf,'color','white');
exptypes = t_train.experiment_type(IA(1:end));
x = find(~strcmp(exptypes(1:end-1),exptypes(2:end)));
x = [0;x;length(exptypes)];

subplot(4,1,1)
plot([corr_2A3;corr_DMS]'); hold on
plot(x*[1 1],[-0.5 1],'k'); hold off
legend('2A3','DMS'); title('Correlation coefficient, Pearson, to Shujun prediction');

subplot(4,1,2)
plot([corr_spearman_2A3;corr_spearman_DMS]'); hold on
plot(x*[1 1],[-0.5 1],'k'); hold off
legend('2A3','DMS'); title('Correlation coefficient, Spearman, to Shujun prediction');

subplot(4,1,3)
plot([mae_2A3;mae_DMS]'); hold on
plot(x*[1 1],[0 4],'k'); hold off
legend('2A3','DMS'); title('Mean absolute error, to Shujun prediction');

subplot(4,1,4)
nt_colors = [1 0.5 0; 0 0.5 0; 1 0 0;0 0 0.8];
for n = 1:4
    plot(meanval(:,n),'color',nt_colors(n,:)); hold on
end
plot(x*[1 1],[0 4],'k'); hold off
legend(nts'); title('Mean value at A,C,G,U');
ylim([0 2])
for i = 1:5; text( (x(i)+x(i+1))/2, 0, exptypes{x(i+1)-1}, 'horizontalalign','center','verticalalign','top','interpreter','none'); end;