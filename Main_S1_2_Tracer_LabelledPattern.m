met = 'glc_tracer';
Size = 6;
label_atoms = ['1';'2'];
label_ratio = 1;

emu_labels = cell(Size,2);
for k = 1:Size
    AtomIndex_k = nchoosek(1:Size,k);
    emu_ids_k = char([]);
    for i = 1:k
        emu_ids_k(:,i) = mydec2hex(AtomIndex_k(:,i));
        mid = strcat('m',num2str((0:k)'));
    end
    emu_mid = [];
    label_mid = [];    
    for j = 1:height(emu_ids_k)
        label_mid_j = zeros(height(mid),1);
        emu_ids_k_j = emu_ids_k(j,:);
        emu_mid = [emu_mid;strcat(string(emu_ids_k_j),'_',mid)];
        label_index = sum(ismember(emu_ids_k_j',label_atoms));
        label_mid_j(sum(ismember(emu_ids_k_j',label_atoms))+1) = label_ratio;
        label_mid = [label_mid;label_mid_j];
    end
    emu_labels{k,1} = emu_mid;
    emu_labels{k,2} = label_mid;
end
emu_ids_total = [];
emu_label_total = [];
for i = 1:Size
    emu_ids_total = [emu_ids_total;emu_labels{i,1}];
    emu_label_total = [emu_label_total;emu_labels{i,2}];
end

emu_tracer_mid = strcat(met,'_',emu_ids_total);
emu_tracer_var = emu_label_total;

LabelledPattern = table(emu_tracer_mid,emu_tracer_var);
writetable(LabelledPattern,strcat('Tracer_',regexprep(met,'\_tracer',''),'_',strjoin(string(label_atoms),''),'/Data_raw/LabelledPattern','.xlsx'))
