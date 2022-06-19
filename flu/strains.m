% lookup table
global Virus;
% note that the strains and ordering are made to correspond to the data in cohen21 PLOS, see file exp/cohen21flu-21.dat
% also should correspond to the sequence alignment in msamat
Virus={'C09','V04', 'J57', 'HK99', 'AI68', 'SH13', 'JX13', 'HB09', 'P10', 'WA79'};
%gind1= @(x) find(ismember(upper(Virus(:)), upper(x))); % does not work correctly for vector x because find sorts indices
