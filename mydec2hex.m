% 仅限16及16以内的数值
function d2h = mydec2hex(num)
[s1, s2] = size(num);
d2h = '';
for i = 1:s1
    for j = 1:s2
        num_ij = num(i,j);
        switch num_ij
            case 1
                d2h(i,j) = '1';
            case 2
                d2h(i,j) = '2';
            case 3
                d2h(i,j) = '3';
            case 4
                d2h(i,j) = '4';
            case 5
                d2h(i,j) = '5';
            case 6
                d2h(i,j) = '6';
            case 7
                d2h(i,j) = '7';
            case 8
                d2h(i,j) = '8';
            case 9
                d2h(i,j) = '9';
            case 10
                d2h(i,j) = 'A';
            case 11
                d2h(i,j) = 'B';
            case 12
                d2h(i,j) = 'C';
            case 13
                d2h(i,j) = 'D';
            case 14
                d2h(i,j) = 'E';
            case 15
                d2h(i,j) = 'F';
            case 16
                d2h(i,j) = 'G';
        end
    end
end
d2h = reshape(d2h,[],1);
end