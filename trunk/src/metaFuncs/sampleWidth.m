function y = sampleWidth(x,w)
% sample 'x' uniformly to get 'w' points
xLen = size(x,2);
if (w==0) || (w > xLen)
    y = x;
else
    inc = floor(xLen/w);
    y = x( : , (1:inc:w*inc) );
end
end