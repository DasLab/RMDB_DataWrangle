r_old = read_rdat_file('../HC16M2R_1M7_0001.rdat');

%% From Kalli's script that she dug out:
reorder = [1, 8,9,10,1,11,12,13,1,14,15,16,1,17,18,19,1,2,3,4,1,5,6,7,1,44,45,46,1,47,48,49,1,50,51,52,1,8,53,54,1,11,55,56,1,14,57,58,1,17,59,60,1,2,61,62,1,5,63,64,1,65,66,67,1,68,69,70,1,38,39,40,1,41,42,43,1,20,21,22,1,23,24,25,1,26,27,28,1,29,30,31,1,32,33,34,1,35,36,37,1,71,24,72,1,73,27,74,1,75,30,76,1,77,33,78,1,79,36,80,1,81,82,83];
assert(length(reorder)==size(r_old.reactivity,2));
t = readtable('round1_mutations.csv');
t = t(reorder,:);

%% Check quartets
for i = 1:height(t)/4
    idx = (i-1)*4;
    assert(strcmp(t.mut1{idx+1},'WT')); % WT
    assert(strcmp(t.mut1{idx+2},t.mut1{idx+4})); % double mut1 = single mut1
    assert(strcmp(t.mut1{idx+3},t.mut2{idx+4})); % double mut2 = single mut2
end

%% Create remediated file.
r = r_old;
seq_shift = 51-106;
%r.sequence = strrep(r_old.sequence,'X','');
r.comments = {
    'hc16 RNA ligase mutation rescue testing stems P4, alt-P4, P7, and alt-P7.',...
    'Mutate-map-rescue read out by capillary electrophoresis.'...
    'Note that WT and some single mutants are repeated so that the profiles show quartets of WT, MutA, MutB, and MutAB.',...
    'No attenuation correction of background subtraction applied.',...
    'Experiment 031619_hc16m2r_1M7_Elim_578931 by W. Kladwang, K. Kappel.'};
for i = 1:height(t)
    r.data_annotations{i} = {};
    if ~strcmp(t.mut1{i},'WT')
        [mutpos,mut_seq,start_seq] = get_mutation_info_from_tag( strrep(t.mut1{i},'T','U')  );
        mutpos = mutpos+seq_shift; % the r.offset is because get_mutation_info_from_tag subtracts r.offset
        r.data_annotations{i} = [r.data_annotations{i},{sprintf('mutation:%s%d%s',start_seq,mutpos,mut_seq)}];
    end
    if length(t.mut2{i})>0
        [mutpos,mut_seq,start_seq] = get_mutation_info_from_tag( strrep(t.mut2{i},'T','U'));
        mutpos = mutpos+seq_shift;
        r.data_annotations{i} = [r.data_annotations{i},{sprintf('mutation:%s%d%s',start_seq,mutpos,mut_seq)}];
    end
end

check_rdat(r);

%% 
if ~exist('NEW','dir'); mkdir('NEW'); end;
outfile = 'NEW/HC16M2R_1M7_0001.rdat';
output_rdat_to_file(outfile,r);
r_new = read_rdat_file(outfile);

