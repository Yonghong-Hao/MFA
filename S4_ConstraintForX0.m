function T = S4_ConstraintForX0(VariablesRefernce)
emuAll = VariablesRefernce.emu_emu;
emuUnique = unique(emuAll);
% constraint matrix T
T = zeros(numel(emuUnique),numel(emuAll));
for i = 1:numel(emuUnique)
    emuIndex = strcmp(emuAll,emuUnique(i));
    emuCount = sum(emuIndex);
    T(i,emuIndex) = ones(1,emuCount);
end
end