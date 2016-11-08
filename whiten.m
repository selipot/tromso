function [rgbout] = whiten(rgbin,varargin)
% WHITEN(RGBIN,P) - generate a p% whitened RGB triplet compared to the input
% WHITEN(RGBIN) - generate a 50% whitened RGB triplet compared to the input
  if isempty(varargin)
    p = 50;
  else
    p = varargin{1};
  end
  
  rgbout = rgbin+(1-rgbin)*p/100;
  