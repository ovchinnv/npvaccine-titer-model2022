% plot correlation coefficients for all models
if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot') ;
 qoct=1;
else
 qoct=0;
end
%
if (~exist('model'))
 model='dist2ave' ; % name of matlab script file that performs the fit
end
%
if (~exist('enc'))
 enc='gr' ;
end
%
if (~exist('nptype'))
 nptype='mosaic' ;
end
%
fn=[model,'-',enc,'-',nptype] ;
% plot correlation coefficients for all models over vaccine splits, i.e. 4 choose 2 = 6 models :
if ~exist('qvacsplit') ; qvacsplit=0 ; end
if qvacsplit
 fn=[fn,'-vac'];
else ; 
 fn=[fn,'-ags'];
 if ~exist('qholdvac') ; qholdvac=0 ; end
 if (qholdvac)
  if ~exist('qtestholdvac') ; qtestholdvac=0 ; end
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
%%% compute average correlation coefficient & average error using nrep samples
%
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
% as a check, recompute any of the above from alltest data :
assert(vind==size(alltestdata,1)) ;
nag= max(max(alltestdata(:,:,1)));
nvac=max(max(alltestdata(:,:,2)));
testmat=sparse(nag, nvac);
testmatexp=sparse(nag, nvac);
itestmat=sparse(nag, nvac); % count matrix
%
for i=1:vind
 allct2(i)=corr(alltestdata(i,:,3)', alltestdata(i,:,4)' ) ;
% compute average prediction :
 itestsample=squeeze(alltestdata(i,:,:)) ;
 testmat = testmat + sparse(itestsample(:,1),itestsample(:,2),itestsample(:,3),nag,nvac);
 testmatexp = testmatexp + sparse(itestsample(:,1),itestsample(:,2),itestsample(:,4),nag,nvac);
 itestmat = itestmat + sparse(itestsample(:,1),itestsample(:,2),1,nag,nvac);
end
testmat=testmat./itestmat ; % average
testmatexp=testmatexp./itestmat ; % average
iok=~isnan(testmat); % valid inds

assert(norm(allct2(:)-allct(:))<1e-7) % stop if they do not match
% compute mean prediction :
if (qoct)
 cta=corr(testmat(iok(:)),testmatexp(iok(:)))
 csta=spearman(testmat(iok(:)),testmatexp(iok(:)))
else
 [cta,pval]=corr(full(testmat(iok(:))),full(testmatexp(iok(:))))
 [csta,spval]=corr(full(testmat(iok(:))),full(testmatexp(iok(:))), 'type', 'spearman')
end
%
okinds=1:vind ; % vaccination indices ; take all
%
if (1)
f=figure(1) ; hold on ;
%set(f,'position',[ 200, 100, 700, 525 ]) ;
set(f,'position',[ 200, 100, 750, 550 ]) ;
set(gca, 'fontsize', 14)
st='k.';
st='ro';
%st='bx';
%st='gv';

%errorbar(allca(okinds), allcta(okinds), allcte(okinds),st)
if qoct
 scatter(allca(okinds), allcta(okinds),st(1),st(2))
else
 scatter(allca(okinds), allcta(okinds),st)
end
xlabel('\it c_P^{train}')
ylabel('\it c_P^{test}')
xlim([0 1])
ylim([0 1])

% average corr coef
cptrave=mean(allca(okinds));
cptrerr=std(allca(okinds));

cptave=mean(allcta(okinds));
cpterr=std(allcta(okinds));
% pvalues :
minctp=min(allctpval)
maxctp=max(allctpval)
avectp=mean(allctpval)

%l=sprintf('C_P^t %4.2f +/- %5.3f', cptrave, cptrerr')
leg={['\it C_P^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_P^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2), ';\it C_P^{test,ave}=',num2str(cta) ]};
%legend(leg, 'location', 'northwest'); legend boxoff;
legend(leg, 'location', 'southwest'); legend boxoff;

box on ;
set(gcf, 'paperpositionmode','auto')
%print(gcf, '-dpng', [fn,'-cp.png']);
%print(gcf, '-depsc2', [fn,'-cp.eps']);

if (1) % spearman
 st='gv';
% errorbar(allcsa(okinds), allcsta(okinds), allcste(okinds),st)
 if qoct
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
 cptrave=mean(allcsa(okinds));
 cptrerr=std(allcsa(okinds));

 cptave=mean(allcsta(okinds));
 cpterr=std(allcsta(okinds));
% pvalues :
 mincstp=min(allcstpval)
 maxcstp=max(allcstpval)
 avecstp=mean(allcstpval)
%
 leg=[leg {['\it C_S^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), ';\it C_S^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2), ';\it C_S^{test,ave}=',num2str(csta)]}];
% legend(leg, 'location', 'northwest'); legend boxoff;
 legend(leg, 'location', 'southwest'); legend boxoff;%

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
