function y = getMin(x)
% y(:,1) = minimum value in x, y(:,2) = position of min
[m,i] = nanmin(x,[],2);
y = [m,i];
end
