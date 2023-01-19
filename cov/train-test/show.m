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
 enc='gr' ; % name of matlab script file that performs the fit
end
%
fn=[model,'-',enc] ;
% plot correlation coefficients for all models over vaccine splits, i.e. 4 choose 2 = 6 models :
if ~exist('qvacsplit') ; qvacsplit=0 ; end
if qvacsplit ; fn=[fn,'-vac']; else ; fn=[fn,'-ags']; end
%
load([fn,'.mat']);
%
qcsp=1; %whether to plot spearman
if ~exist('leg')
 leg={};
end

if ~exist('qtex') % plot latex table
 qtex=1;
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
% to include certain models
qany=1 ; % otherwise qall ; means that any of the strains must be present (vs all)
if ~exist('trains')
 trains={}
%trains={'sars-2'} ; % only models that have these strains in training set
%trains={'ratg13'} ;
%trains={'sars-2', 'ratg13'} ;
%trains={'sars'} ;

%trains=alltest(1) %
end
%%%%%
checks={} ;
%checks={'sars-2'} ; % only models that have these strains in testing set
%checks={'sars-2', 'ratg13'} ;

tinds=[];
for train = trains
 [~,tind]=ismember(train, allags);
 if (tind>0) ; tinds=[tinds tind] ; end
end
cinds=[];
for check = checks
 [~,cind]=ismember(check, allags);
 if (cind>0) ; cinds=[cinds cind] ; end
end
% go over all strains and pick those that contain tinds and do not contain cinds, whcih are part of a check set
okinds=[]
for i=1:vind
 if (qany)
  if ((isempty(tinds)|any(ismember(tinds,indags(i,:)))) & (isempty(cinds)|~all(ismember(cinds,indags(i,:)))) ) ; okinds=[okinds i]; end
 else
  if ((isempty(tinds)|all(ismember(tinds,indags(i,:)))) & (isempty(cinds)|~any(ismember(cinds,indags(i,:)))) ) ; okinds=[okinds i]; end
 end
end

f=figure(1) ; hold on ;
%set(f,'position',[ 200, 100, 700, 525 ]) ;
set(f,'position',[ 200, 100, 750, 550 ]) ;
set(gca, 'fontsize', 14) ;
st='k.';
st='ro';
if ~isempty(trains)
 st='bx';
end
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

% average corr coef
cptrave=mean(allca(okinds))
cptrerr=std(allca(okinds))
%
cptave=mean(allcta(okinds))
cpterr=std(allcta(okinds))
% pvalues :
minctp=min(allctpval)
maxctp=max(allctpval)
avectp=mean(allctpval)

if (qtex)
% output for table in latex :
 ftex=fopen([fn,'.tex'],'w')
%fprintf(ftex, 'model & $c_P^{train}$(m$\\pm$std) & $c_P^{test}$(m$\\pm$std) & $p^{test}_{val}$(min/max/avg)\n') ;
%fprintf(ftex, '%s & %3.2f$\\pm$%3.2f & %3.2f$\\pm$%3.2f & %2.1e/%2.1e/%2.1e', fn,cptrave,cptrerr,cptave,cpterr,minctp,maxctp,avectp);
 fprintf(ftex, 'model & $c_P^{train}$(m$\\pm$std) & $c_P^{test}$(m$\\pm$std) & $c_S^{train}$(m$\\pm$std) & $c_S^{test}$(m$\\pm$std) & $p^{test}_{val}$(min/max/avg)\\\\\n') ;
 fprintf(ftex, '%s & %3.2f$\\pm$%3.2f & %3.2f$\\pm$%3.2f & %2.1e/%2.1e/%2.1e',fn,cptrave,cptrerr,cptave,cpterr);
end

%l=sprintf('C_P^t %4.2f +/- %5.3f', cptrave, cptrerr')
if isempty(trains)
 leg={['\it C_P^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_P^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2) ]};
else
 leg=[leg {['\it C_P^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_P^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2) ]}];
end

legend(leg, 'location', 'northwest'); legend boxoff;

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
%
 cptrave=mean(allcsa(okinds))
 cptrerr=std(allcsa(okinds))
%
 cptave=mean(allcsta(okinds))
 cpterr=std(allcsta(okinds))
%
% pvalues :
 mincstp=min(allcstpval)
 maxcstp=max(allcstpval)
 avecstp=mean(allcstpval)
%

 leg=[leg {['\it C_S^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_S^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2) ]}];
 if (isempty(trains))
  legend(leg, 'location', 'northwest'); legend boxoff;
 end
%
 if exist('lbl')
  text(-0.12,0.99,lbl,'fontsize',18);
 end
%
 box on ;
 set(gcf, 'paperpositionmode','auto')
% print(gcf, '-dpng', [fn,'-cs.png']);
 print(gcf, '-depsc2', [fn,'-cs.eps']);
 if (qtex)
% output for table in latex :
% fprintf(ftex, '%3.2f$\\pm$%3.2f & %3.2f$\\pm$%3.2f & %3.2e/%3.2e/%3.2e\\\n', cptrave,cptrerr,cptave,cpterr,mincstp,maxcstp,avecstp);
% output pvalue for pearson correlation, only, to save line space
  fprintf(ftex, '& %3.2f$\\pm$%3.2f & %3.2f$\\pm$%3.2f & %2.1e/%2.1e/%2.1e\\\\\n', cptrave,cptrerr,cptave,cpterr,minctp,maxctp,avectp);
 end
end % qcsp
if (qtex) ; fclose(ftex) ; end

