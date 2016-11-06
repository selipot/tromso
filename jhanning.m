function [window] = jhanning(n)
window = frac(1,2)*(1-cos(frac(2*pi*(1:ceil(n/2)),n+1)))';

if iseven(n)
    window = [window;flipud(window)];
else
    window = [window(1:end-1);flipud(window)];
end
