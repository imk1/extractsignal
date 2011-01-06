function y = samplePos(x,pos)
% sample the positions specified in 'pos' from x
xLen = size(x,2);
y = nan( size(x,1) , numel(pos) );
if all(pos < 1) 
    pos = round(pos*xLen); % convert fractional positions to actual positions
else
    pos = round(pos);
end
validPos = find( (pos > 0) & (pos <= xLen) ); % find valid positions
y(:,validPos) = x(:,pos(validPos)); % extract positional data
end