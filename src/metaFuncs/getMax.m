function y = getMax(x)
% y(:,1) = maximum value in x, y(:,2) = position of max
[m,i] = nanmax(x,[],2);
y = [m,i];
end