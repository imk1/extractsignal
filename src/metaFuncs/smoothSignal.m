function y = smoothSignal(x,w)
% smooth rows in x using window size w
[nSig,xLen] = size(x);
y = x;
if (w==1) || (w > xLen) || (xLen == 1) || all(all(isnan(x)))
    return;
end
for i = 1:nSig
    y(i,:) = kdeSmooth( x(i,:) , w/2 , 'triw' , 1 );
end
end