% 8/7/21 : plot exprimental data from  Alex Cohen (Science 21), fig 3
%

fname='cohen21-14rbd.dat';

fid=fopen(fname,'r');

d=textscan(fid,'%s', 'delimiter',',\n');
d=reshape(d{1},4,[]);
d=d(1:end,2:end);

id=str2num(char(d(1,1:end)));
iggname=d(2,1:end); % strain & vaccine name
iggexp=str2num(char(d(end-1,1:end))); % mean igg signal
iggexpe=abs(str2num(char(d(end,1:end)))); % error bar

if ~exist('noxplot')
 noxplot=0;
end;
if noxplot ; return ; end

% plot
b=bar(id, iggexp, 'facecolor',0.75*[1 0 0]) ; hold on ;
errorbar(id, iggexp, 0, iggexpe,  'k.');

set(gca, 'xtick',id);
xlim([min(id)-1, max(id)+1]);

