function [c, ceq] = S4_mycon_mapping1(U,S)
ceq = S*exp(U');
c = [];
end