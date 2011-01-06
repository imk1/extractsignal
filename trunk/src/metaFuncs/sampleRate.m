function y = sampleRate(x,rate)
% Will sample a signal every 'rate' bp
y = x( : , (1:rate:end) );
end