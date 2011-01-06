function intervalData = slopIntervalData(intervalData , slopLeft , slopRight, slopStrand, slopUseSummit)
% adjusts intervals depending on slopParameters
% function intervalData = slopIntervalData(intervalData , slopLeft , slopRight, slopStrand, slopUseSummit)
% ================
% INPUT ARGUMENTS
% ================
% intervalData<dataset>: chr[nominal],start[double],stop[double],strand[nominal], OPTIONAL:summit[double]
% slopUseSummit('false','midpoint','peak'): whether to slop wrt ends or midpoint or peak/summit
% slopLeft<double>: left slop value (can be negative)
% slopRight<double>: right slop value (can be negative)
% slopStrand(true/false): if set to true left refers to 5' prime end

if (slopLeft==0 && slopRight==0)
    return;
end

switch slopUseSummit    
    case 'false' % do not use summit information i.e. slop from ends of the intervals
        
        if slopStrand % Use strand info
            plusStrand = (intervalData.strand == '+');
            intervalData.start = intervalData.start - (plusStrand * slopLeft + ~plusStrand * slopRight);
            intervalData.stop = intervalData.stop + (plusStrand * slopRight + ~plusStrand * slopLeft);
        else
            intervalData.start = intervalData.start - slopLeft;
            intervalData.stop = intervalData.stop + slopRight;
        end
        
    case 'midpoint'
        
        if slopStrand % Use strand info
            plusStrand = (intervalData.strand == '+');
            midPoint = floor( (intervalData.start + intervalData.stop) / 2 );
            intervalData.start = midPoint - (plusStrand * slopLeft + ~plusStrand * slopRight);
            intervalData.stop =  midPoint + (plusStrand * slopRight + ~plusStrand * slopLeft);            
        else
            intervalData.start = intervalData.start - slopLeft;
            intervalData.stop = intervalData.stop + slopRight;            
        end
        
    case 'peak'
        
        if slopStrand % Use strand info
            plusStrand = (intervalData.strand == '+');
            intervalData.start = intervalData.summit - (plusStrand * slopLeft + ~plusStrand * slopRight);
            intervalData.stop = intervalData.summit + (plusStrand * slopRight + ~plusStrand * slopLeft);            
        else
            intervalData.start = intervalData.summit - slopLeft;
            intervalData.stop = intervalData.summit + slopRight;            
        end        
        
    otherwise
        error('ERROR: Unsupported option for slopUseSummit');
        
end

% Check for invalid intervals
if (slopLeft < 0 || slopRight < 0) && nnz(intervalData.stop - intervalData.start < 0)
    fprintf( 2 ,'WARNING: You are using negative slop values (interval contraction) which is causing some intervals to have start > stop. These will be automatically flipped.' );
end

end
