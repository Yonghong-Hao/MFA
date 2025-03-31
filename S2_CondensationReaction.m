function midTerms = S2_CondensationReaction(reactant_r,product_r,flux_r,Size_r)
%% emu reactant
% 将emu展开为emu_mid
emuSize_r = zeros(numel(reactant_r),1);
for i = 1:numel(reactant_r)
    emu = reactant_r{i};
    emuSplit = regexp(emu,'_','split');
    midR = emuSplit(end);
    emuSize_r(i) = length(midR{:});
end
emu_mid_r = cell(numel(reactant_r),1);
for i = 1:numel(reactant_r)
    emu = reactant_r{i};
    emu_mid_r(i) = {strcat(emu, '_m', cellfun(@num2str, num2cell(0:emuSize_r(i)),'UniformOutput',false))};
end

% 车轮式进行卷积运算
emu_midR1 = emu_mid_r{1};
for rr = 2:numel(reactant_r)
    emu_midR2 = emu_mid_r{rr};
    emu_midR = cell(numel(emu_midR1)+numel(emu_midR2)-1,1);
    for i = 1:numel(emu_midR1)
        for j = 1:numel(emu_midR2)
            iterm_i = emu_midR1{i};
            iterm_j = emu_midR2{j};
            % 分配律
            if ~ismember('+',iterm_i)
                iterm_ij = strcat(iterm_i,'*',iterm_j);
            else
                iterm_i_split = regexp(iterm_i,' \+','split');
                iterm_i_split_multiple = strcat(iterm_i_split,'*',iterm_j);
                iterm_i_split_multiple_sum = strjoin(iterm_i_split_multiple,' +');
                iterm_ij = iterm_i_split_multiple_sum;
            end
            emu_midR{i+j-1} = strcat(emu_midR{i+j-1},' +',iterm_ij);
        end
    end
    for k = 1:numel(emu_midR)
        str = emu_midR{k};
        str(1:2) = '';
        emu_midR{k} = str;
    end
    emu_midR1 = emu_midR;
end

% 合并同类项
emu_midR1_combine = emu_midR1;
for i = 1:numel(emu_midR1_combine)
    iterm_i = emu_midR1_combine{i};
    if ~ismember('+',iterm_i)
        iterm_i_multiple = regexp(iterm_i,'\*','split');
        [iterm_i_multiple_unique,~,ic] = unique(iterm_i_multiple');
        a_counts = accumarray(ic,1);
        for t = 1:numel(a_counts)
            if a_counts(t)~=1
                iterm_i_multiple_unique{t} = strcat(iterm_i_multiple_unique{t},'^',num2str(a_counts(t)));
            end
        end
        if numel(a_counts)>1
            iterm_i = strjoin(iterm_i_multiple_unique,'*');
        else
            iterm_i = iterm_i_multiple_unique;
        end
        emu_midR1_combine{i} = iterm_i;
    else
        iterm_i_sum = regexp(iterm_i,' \+','split');
        for s = 1:numel(iterm_i_sum)
            iterm_i_sum_multiple = regexp(iterm_i_sum{s},'\*','split');
            [iterm_i_sum_multiple_unique,~,ic] = unique(iterm_i_sum_multiple');
            a_counts = accumarray(ic,1);
            for t = 1:numel(a_counts)
                if a_counts(t)~=1
                    iterm_i_sum_multiple_unique{t} = strcat(iterm_i_sum_multiple_unique{t},'^',num2str(a_counts(t)));
                end
            end
            if numel(a_counts)>1
                iterm_i_sum{s} = strjoin(iterm_i_sum_multiple_unique,'*');
            else
                iterm_i_sum{s} = iterm_i_sum_multiple_unique{1,1};
            end
        end
        [iterm_i_sum_unique,~,ic] = unique(iterm_i_sum');
        a_counts = accumarray(ic,1);
        for t = 1:numel(a_counts)
            if a_counts(t)~=1
                iterm_i_sum_unique{t} = strcat(num2str(a_counts(t)),'*',iterm_i_sum_unique{t});
            end
        end
        if numel(a_counts)>1
            iterm_i = strjoin(iterm_i_sum_unique,'+');
        else
            iterm_i = iterm_i_sum_unique;
        end
    end
    emu_midR1_combine{i} = iterm_i;
end

%% emu product
emu_midP = strcat(product_r,'_m',cellfun(@num2str, num2cell(0:Size_r),'UniformOutput',false));

%% emu flux
% 编辑flux格式：含有+的项左右添加（）
if ismember('+',flux_r{:})
    flux_r = strcat('(',flux_r,')');
end
emu_flux = repmat(flux_r,Size_r+1,1);

%% emu size
emu_size = repmat(Size_r,Size_r+1,1);

%% combine emu reactant and product
midTerms = table(emu_midR1_combine,emu_midP',emu_flux,emu_size,'VariableNames',{'ReactantMid','ProductMid','Flux','Size'});
end