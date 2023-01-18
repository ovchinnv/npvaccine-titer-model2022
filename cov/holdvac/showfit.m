if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot') ;
 qoct=1;
else
 qoct=0;
end
%
if (~exist('model'))
 model='dist2ave' ; % name of matlab script file that performs the fit
 model='avedist' ;
end
%
if (~exist('enc'))
 enc='gr' ; % name of matlab script file that performs the fit
end
%
fn=[model,'-',enc] ;
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
iggmoda=mean(alliggmod,1) ;

% bar plot of model vs experiment :
f=figure(2);
set(f,'position', [100, 100, 1000, 300]) ; hold off ;
%
for i=1:numel(agtest)
inds=i ;
%%%%%%%%%%% plot
b=bar(id(inds), iggmoda(inds), 'facecolor',0.99*[1 1 1], 'barwidth', 0.3) ; hold on ; set(b,'linewidth',2)
b=bar(id(inds), iggexp1(inds), 'facecolor',cols(itestsample(i,2),:), 'barwidth',0.2) ; hold on ;
errorbar(id(inds), iggexp1(inds), 0, iggexpe1(inds), 'k.' );
set(gca, 'xtick',id);
set(gca, 'xticklabel',[]);
set(gca, 'tickdir','out');
end
xlim([-1,numel(agtest)]);
%return
%
% correlation coefficient of average model
cta=corr(iggmoda(:),iggexp1(:))
if (qoct)
 csta=spearman(iggmoda(:),iggexp1(:))
else
 csta=corr(iggmoda(:),iggexp1(:), 'type', 'spearman')
end
%
okinds=1:vind ; % vaccination indices ; take all
%
if (1)
f=figure(1) ; hold on ;
%set(f,'position',[ 200, 100, 700, 525 ]) ;
set(f,'position',[ 200, 100, 750, 550 ]) ;
set(gca, 'fontsize', 14) ;
st='k.';
st='ro';
%st='gv';
%errorbar(allca(okinds), allcta(okinds), allcte(okinds),st)
if (qoct)
 scatter(allca(okinds), allcta(okinds),st(1),st(2))
else
 scatter(allca(okinds), allcta(okinds),st)
end
xlabel('c_P^{train}')
ylabel('c_P^{test}')
xlim([0 1])
ylim([0 1])
%
% average corr coef
cptrave=mean(allca(okinds))
cptrerr=std(allca(okinds))

cptave=mean(allcta(okinds))
cpterr=std(allcta(okinds))

leg={['\it C_P^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_P^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2), ';\it C_P^{test,ave}=',num2str(cta) ]};
%legend(leg, 'location', 'northwest'); legend boxoff;
legend(leg, 'location', 'southwest'); legend boxoff;

box on ;
set(gcf, 'paperpositionmode','auto')
%print(gcf, '-dpng', [model,'-cp.png']);
%print(gcf, '-depsc2', [fn,'-cp.eps']);

if (qcsp) % spearman
 st='gv';
% errorbar(allcsa(okinds), allcsta(okinds), allcste(okinds),st)
 if (qoct)
  scatter(allcsa(okinds), allcsta(okinds),st(1),st(2))
 else
  scatter(allcsa(okinds), allcsta(okinds),st)
 end
 xlabel('\it c_{P,S}^{train}')
 ylabel('\it c_{P,S}^{test}')
 xlim([0 1])
 ylim([0 1])
% average model value :
plot([0 1], [cta cta],['r--']) ;
plot([0 1], [csta csta],['g--'])
%
 cptrave=mean(allcsa(okinds))
 cptrerr=std(allcsa(okinds))
%
 cptave=mean(allcsta(okinds))
 cpterr=std(allcsta(okinds))
%
 leg=[leg {['\it C_S^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), ';\it C_S^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2), ';\it C_S^{test,ave}=',num2str(csta)]}];
% legend(leg, 'location', 'northwest'); legend boxoff;
 legend(leg, 'location', 'southwest'); legend boxoff;
%
 if exist('lbl')
  text(-0.12,0.99,lbl,'fontsize',18);
 end

 box on ;
 set(gca, 'fontsize',14)
 set(gcf, 'paperpositionmode','auto')
 print(gcf, '-dpng', [fn,'-cs.png']);
 print(gcf, '-depsc2', [fn,'-cs.eps']);
end

%errorbar(alle2aa(okinds), alle2ata(okinds), alle2ae(okinds),'k.')
%xlabel('E^{train}')
%ylabel('E^{test}')

end
