%function dMdV = MatrixGradient(M,NumFluxes)
[sizeM_1, sizeM_2] = size(M);
dMdV = zeros(sizeM_1,sizeM_2,NumFluxes);
VVal = eye(NumFluxes);

for m = 1:sizeM_1
    for n = 1:sizeM_2
        M_mn = M{m,n};
        if ~isempty(M_mn)
            for t = 1:NumFluxes
                V = VVal(t,:);
                dMdV(m,n,t) = eval(M_mn);
            end
        end
    end
end
%end