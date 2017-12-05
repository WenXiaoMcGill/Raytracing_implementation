function[source] = rmdev(source)
% remove the deviation during the calculation
if isempty(find(abs(source)<1e-15, 1)) == 0
    source(find(abs(source)<1e-15,1)) = 0;
end