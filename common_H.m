function [common_H] = common_H(H, f0)

[I, J] = size(H);
common_H = zeros(I, J);
for i = 1 : I % row
    for j = i : J % colume
        trans = 8 / f0; % basic freq -> 8

        if i == 1 || common_H(i - 1, j) == 0
            common_H(i, j) = H(ceil(trans * i), j);
        end
        
    end
end

end