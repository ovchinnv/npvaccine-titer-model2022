%
if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot') ;
end
%
strains;
qlabel=1;

if (~exist('model'))
 model='dist2ave';
end
%
if (~exist('tag'))
 tag='headstem';
end
%
if (~exist('enc'))
 enc='gr'; % default is grantham encoding
end
%
if (~exist('modver')) % model version should be defined in mod22
 modver=1;
end
%
if (~exist('noplot')) % whether to plot or just to check fit
 noplot=0;
end
%
load([model, '_', tag, '_', enc, '.mat'])
%noplot=0
close all;
lims=[-0.25 1.75];
lims=[-0.25 1.25];

cols=[ 255 92 103 ; 0 255 0 ; 255 165 0 ; 0 0 250 ]/255;
%
% create test sample :
% take all those not present in fitting set
iall_mat=ones(nallags,numel(vacs));
itrainsample_mat=sparse(itrainsample(:,1), itrainsample(:,2), 1);
itestsample_mat=iall_mat-itrainsample_mat ;
[ias,ivs]=find(itestsample_mat);
itestsample=[ias,ivs];
%
%
if (isempty(itestsample))
 itestsample=itrainsample ; % use all test strains
end
agtest=unique(itestsample(:,1)) ;
vactest=unique(itestsample(:,2)) ;
%
cols=repmat(cols, numel(agtest), 1);
%
if qnorm == 1
% scale :
 scale=0;
 iscale=0;
 for itest=1:size(itestsample,1) % test strains
  ia=itestsample(itest,1); % antigen
  iv=itestsample(itest,2); % vaccine
  scale = scale + iggemat( ia, iv ) ;
  iscale = iscale + 1 ;
 end
 enmat = iggemat / scale * iscale ;
else
 enmat = ones(size(iggmat));
end
oenmat=1./enmat;
%
clear iggexp1 iggexpe1 iggmod
wgt2=bestwgt.^2; % squared weights
% to weight coordinates :
ind=0;

for itest=1:size(itestsample,1) % this is a marix with two columns, ag index in first column, vaccine index in second
  ia=itestsample(itest,1); % antigen
  iv=itestsample(itest,2); % vaccine
%
  if (strcmp(model,'dist2ave')) % distance to average model
   dcoor=reshape(vcoor(iv,:)-coor(ia,:), ndim, []);
   ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor
%
   d2 = sum( wgt2 .* ndcoor2 );  % squared distance
%
   ind=ind+1;
%
   if (modver==1)
    iggmod(ind) = 1./(xint+d2^xp) ;  % model igg signal
   elseif (modver==2)
    iggmod(ind) = 1./(xint+d2)^xp;  % model igg signal
   end
%
  elseif (strcmp(model, 'avedist')) % average distance model
%
   dave=0 ;% average distance
   for k=vaccines{iv}

    dcoor=reshape(coor(k,:)-coor(ia,:), ndim, []);
    ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor

%    d = sqrt( sum( wgt2 .* ndcoor2(:,k,ia)' ) );  % squared distance to this strain
    d = sqrt( sum( wgt2 .* ndcoor2 ) );  % distance to this strain
    dave = dave + d ;
   end % for
   nstr=numel(vaccines{iv}); % normalization
   dave = dave / nstr ; % mean squared distance between the train strain and all vaccine strains
%
   ind=ind+1;
%
% model value :
   iggmod(ind) = 1./(xint+dave^xp);  % model igg signal
  end % if
%
  iggexp1(ind) = iggmat(ia,iv) ;
  iggexpe1(ind) = iggemat(ia,iv) ;
  oenorm(ind) = oenmat(ia,iv) ;
  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error
end

ibeg=1; % start at this row
iggmodt=iggmod ;
c=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
if exist('OCTAVE_VERSION')
 cs=spearman(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
else % matlab
 cs=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)', 'type', 'spearman')
end
err2=(iggmodt(:) - iggexp1(:)).^2;
e2=sum(err2(ibeg:end))
%
if (noplot)
 return
end
%
% note : plotting only "works" if the test samples are grouped in a certain way (see paper)
figure('position', [100, 100, 1000, 300]) ; hold off ;

%for i=1+isempty(tag):size(cols,1) % hack to throw away first exp datapoint depending on model
for i=1:size(cols,1)
inds=i ;
%%%%%%%%%%% plot
b=bar(id(inds), iggmodt(inds), 'facecolor',0.99*[1 1 1], 'barwidth', 0.3) ; hold on ;
set(b,'linewidth',2)
b=bar(id(inds), iggexp1(inds), 'facecolor',cols(i,:), 'barwidth',0.2) ; hold on ;
errorbar(id(inds), iggexp1(inds), 0, iggexpe1(inds), 'k.' );
set(gca, 'xtick',id);
set(gca, 'xticklabel',[]);

end
xlim([min(id)-1, max(id)+1]);

%return

%write out test strains
for i=1:numel(agtest)
 xx=size(cols,1)/numel(agtest) * (i-1) + 0.6;
 yy=1.9;
 text(xx,yy,upper(allags(i)), 'fontsize', 9);
% compute error and correlation _only_ for this test strain
 i1=1+(i-1)*numel(vactest);
 i2=i1+numel(vactest)-1;
 if (i1==1)
  i1=i1+isempty(tag)
 end
 cc=corr(iggmodt(i1:i2), iggexp1(i1:i2));
% ee=sum(err2(i1:i2))
 ee=sqrt(mean(err2(i1:i2)));
 text(xx,yy-0.1,['rmse=',num2str(ee)], 'fontsize', 9);
%
 plot([xx+4-0. xx+4-0.], [0 2], 'k--')
end

% overall error :
text(0.6,yy-0.25,['rmse-all=',num2str(sqrt(mean(err2)))], 'fontsize', 9);

if (qlabel)
 text(-2, 1.9, 'A', 'fontsize', 20) ;
end

set(gca, 'fontsize', 12)
set(gcf, 'paperpositionmode','auto')
xlim([min(id)-1, max(id)+1]);
set(gca, 'tickdir','out')
print(gcf, '-dtiff', [model, '-bar-', tag,'-',enc,'.tif'])
print(gcf, '-depsc2', [model, '-bar-', tag,'-',enc,'.eps'])

%return
figure
lw=2 ;
errorbar(iggmodt(ibeg:end), iggexp1(ibeg:end), iggexpe1(ibeg:end),  'k.'); hold on ; plot(lims, lims, 'k--', 'linewidth', lw)
plot(iggmodt(ibeg:end), iggexp1(ibeg:end),'ko', 'linewidth', lw); hold on ; plot(lims, lims, 'k--', 'linewidth', lw)

xlabel('\it Model IgG');
ylabel('\it Experimental IgG');

set(gca, 'fontsize', 14, 'linewidth', lw)

text(0,lims(2)-0.03*(diff(lims)),['R^2=', num2str(c.^2), ', c_P=', num2str(c), ', c_S=', num2str(cs)], 'fontsize', 15);

if (qlabel)
 text(-0.8, 1.5, 'B', 'fontsize', 22) ;
end

set(gcf, 'paperpositionmode','auto')
print(gcf, '-dtiff', [model,'-sc-',tag,'-',enc,'.tif'])
print(gcf, '-depsc2', [model,'-sc-',tag,'-',enc,'.eps'])

