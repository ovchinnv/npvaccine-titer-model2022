% plot data that corresponds to the sera after the prime ; Cohen 21, fig 3
%

if (~exist('nptype'))
 nptype='mosaic';
end

if (strcmp(nptype,'mosaic'))
 fname='cohen21flu-21.dat'; % mosaic
elseif (strcmp(nptype,'admix'))
 fname='cohen21flu-21ad.dat'; % admix
end

fid=fopen(fname,'r');

d=textscan(fid,'%s', 'delimiter',',\n');
d=reshape(d{1},4,[]);
d=d(1:end,2:end);

iggscale=0.1; % scale down for consistency with data in covid paper Cohen 21 (Science)
id=str2num(char(d(1,:)));
iggname=d(2,:); % strain+vaccine name
iggexp=str2num(char(d(end-1,:))) * iggscale ; % mean igg signal
iggexpe=abs(str2num(char(d(end,:)))) * iggscale ; % error bar

if ~exist('noxplot')
 noxplot=0;
end;
if noxplot ; return ; end

% plot
b=bar(id, iggexp, 'facecolor',0.75*[1 0 0]) ; hold on ;
errorbar(id, iggexp, 0, iggexpe,  'k.');

set(gca, 'xtick',id);

xlim([min(id)-1, max(id)+1]);
