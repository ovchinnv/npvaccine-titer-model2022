if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot') ;
 qoct=1;
else
 qoct=0;
end
%
fns={'dist2ave-gr-mosaic', 'avedist-gr-mosaic', 'dist2ave-at-mosaic'} ;
%
for imod=1:numel(fns)
fn=fns{imod}
%
% plot correlation coefficients for all models over vaccine splits, i.e. 4 choose 2 = 6 models :
if ~exist('qvacsplit') ; qvacsplit=0 ; end
if qvacsplit f
 n=[fn,'-vac']
else 
 fn=[fn,'-ags']
 if ~exist('qholdvac') ; qholdvac=1 ; end
 if (qholdvac)
  if ~exist('qtestholdvac') ; qtestholdvac=1 ; end
  if (qtestholdvac)
   fn=[fn,'-htest'];
  else
   fn=[fn,'-htrain'];
  end
 end
end
fn
load([fn,'.mat']);
%
qcsp=1; %whether to plot spearman
if ~exist('leg')
 leg={};
end
%%% compute average correlation coefficient & average error using nrep samples
d=reshape(allc,nrep,[]);
allca=mean(d,1);
allce=std(d,1);

d=reshape(allct,nrep,[]);
allcta=mean(d,1);
allcte=std(d,1);

d=reshape(allcs,nrep,[]);
allcsa=mean(d,1);
allcse=std(d,1);

d=reshape(allcst,nrep,[]);
allcsta=mean(d,1);
allcste=std(d,1);

d=reshape(alle2a,nrep,[]);
alle2aa=mean(d,1);
alle2ae=std(d,1);

d=reshape(alle2at,nrep,[]);
alle2ata=mean(d,1);
alle2ate=std(d,1);

if (nrep==1)
 allce=zeros(1,vind);
 allcte=allce;
 allcse=allce;
 allcste=allce;
 alle2ae=allce;
 alle2ate=allce;
end
%
% as a check, recompute any of the above from allerr :
nagtest=size(itestsample,1) ;
for i=1:vind
 alliggmod(i,1:nagtest) = allerrt(i,1:nagtest) ./ oenorm(1:nagtest) + iggexp1(1:nagtest);
 allct2(i)=corr(alliggmod(i,1:nagtest)', iggexp1(:));
end
assert(norm(allct2(:)-allct(:))<1e-7) % stop if they do not match
% compute mean prediction :
% 1) the average model has good correlation, but diminished titer magnitude due to smoothing/averaging;
% we can scale up the titers uniformly, which will preserve the correlation, but improve the error
%scale=mean(iggexp1)/mean(alliggmod(:))
scale=1
if (imod==1)
 iggmoda=mean(alliggmod,1)*scale ;
else
 iggmoda=[iggmoda ; mean(alliggmod,1)*scale] ;
end
% rmse
testrmse(imod) = sqrt(mean( (iggmoda(imod,:)-iggexp1).^2 ) );
end % models

% bar plot of models vs experiment :
f=figure(2);
set(f,'position', [100, 100, 1000, 300]) ; hold off ;
%
%%%%%%%%%%% plot
% to append exp data :
%iggmoda=[iggmoda ; iggexp1(:)'] ;
%b=bar(iggmoda', 'barwidth', 0.3) ; hold on ; set(b,'linewidth',1.5)
%b=bar(iggmoda') ; hold on ; set(b,'linewidth',1.5)
b=bar(iggmoda', 'edgecolor','none') ; hold on ; %set(b,'linewidth',1.5)
%b=bar(iggexp1', 'facecolor','k', 'barwidth',0.075) ; hold on ;
errorbar(1:numel(agtest), iggexp1, 1*iggexpe1, iggexpe1, 'ko--', 'linewidth',1. );
% show RMSE
leg={};
for i=1:numel(fns)
 leg=[leg {[fns{i},' (rmse=',num2str(testrmse(i)),')']}];
end
leg=[leg, {'experiment'}];
%
l=legend(leg);
legend boxoff; set(l, 'fontsize',10)
set(gca, 'xtick',id);
set(gca, 'xticklabel',allags);
set(gca, 'tickdir','out');
xlim([0,numel(agtest)]+0.5);
ylim([0,1.25]);
%
set(gcf, 'paperpositionmode','auto');
print(gcf,'-depsc2','test-flu.eps');

