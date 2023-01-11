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
 enc='gr' ; % name of matlab script file that performs the fit
end
%
fn=[model,'-',enc] ;
% plot correlation coefficients for all models over vaccine splits, i.e. 4 choose 2 = 6 models :
if ~exist('qvacsplit') ; qvacsplit=0 ; end
if qvacsplit ; fn=[fn,'-vac']; else ; fn=[fn,'-ags']; end

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
% to include certain models
%
qany=0 ; % otherwise qall ; means that any of the strains must be present (vs all)
trains={}
%trains={'sars-2'} ; % only models that have these strains in training set
%trains={'ratg13'} ;
%trains={'sars-2', 'ratg13'} ;
%trains={'sars'} ;

%trains=alltest(1) %
%%%%%
%checks={'sars-2'} ; % only models that have these strains in testing set
%checks={'sars-2', 'ratg13'} ;
checks={} ;

tinds=[];
for train = trains
 [~,tind]=ismember(train, allstrains);
 if (tind>0) ; tinds=[tinds tind] ; end
end
cinds=[];
for check = checks
 [~,cind]=ismember(check, allstrains);
 if (cind>0) ; cinds=[cinds cind] ; end
end
% go over all strains and pick those that contain tinds and do not contain cinds, whcih are part of a check set
okinds=[]
for i=1:vind
 if (qany)
  if ((isempty(tinds)|any(ismember(tinds,indvac(i,:)))) & (isempty(cinds)|~all(ismember(cinds,indvac(i,:)))) ) ; okinds=[okinds i]; end
 else
  if ((isempty(tinds)|all(ismember(tinds,indvac(i,:)))) & (isempty(cinds)|~any(ismember(cinds,indvac(i,:)))) ) ; okinds=[okinds i]; end
 end
end

if (1)
f=figure(1) ; hold on ;
set(f,'position',[ 200, 100, 700, 525 ]) ;
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
cptrave=mean(allca(okinds))
cptrerr=std(allca(okinds))

cptave=mean(allcta(okinds))
cpterr=std(allcta(okinds))

%l=sprintf('C_P^t %4.2f +/- %5.3f', cptrave, cptrerr')
leg={['\it C_P^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_P^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2) ]};
legend(leg, 'location', 'northwest'); legend boxoff;

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
%
 cptrave=mean(allcsa(okinds))
 cptrerr=std(allcsa(okinds))

 cptave=mean(allcsta(okinds))
 cpterr=std(allcsta(okinds))
%
 leg=[leg {['\it C_S^{train}=', num2str(cptrave,2), '+/-', num2str(cptrerr,2), '; \it C_S^{test}=', num2str(cptave,2), '+/-', num2str(cpterr,2) ]}];
 legend(leg, 'location', 'northwest'); legend boxoff;
%
 box on ;
 set(gca, 'fontsize',14)
 set(gcf, 'paperpositionmode','auto')
 print(gcf, '-dpng', [fn,'-',nptype,'-cs.png']);
 print(gcf, '-depsc2', [fn,'-',nptype,'-cs.eps']);
end

%errorbar(alle2aa(okinds), alle2ata(okinds), alle2ae(okinds),'k.')
%xlabel('E^{train}')
%ylabel('E^{test}')


end

return

%compute distances between all strains
coor=aln2coor(msamat,qgrantham);
dall=gdist(coor) ;
%reformat into matrix
nstrain=size(coor,1);
dallm=zeros(nstrain) ;

ind=1;
for i=1:nstrain
 for j=i+1:nstrain ;
  dallm(i,j)=dall(ind);
  ind=ind+1;
 end
end
dallm=dallm+dallm'; %symmetrize since the diagonal is 0

clear ichecks itests;
for i=1:vind
% test strains
 itest=indvac(i,:);
 test = alltest(itest);  % indices in alltest
 itest=gind(test); % indices of the vaccine strains in the msa/coor array
 itests(i,:) = itest ;
% check strains
 check=setdiff(alltest, test);
% table ; % need gind
 icheck=gind(check);
 ichecks(i,:)=icheck ;
% compute average distance between check and test (also, average max & average min)
 dave=0;
 dmin=0;
 dmax=0;
 for k=icheck
  davej=0;
  dminj=inf;
  dmaxj=-inf;
  for j=itest
   davej=davej+dallm(j,k) ;
   dminj=min(dminj,dallm(j,k));
   dmaxj=max(dmaxj,dallm(j,k));
  end
  dave=dave+davej/numel(itest);
  dmin=dmin+dminj;
  dmax=dmax+dmaxj;
 end
 daves(i)=dave/(numel(icheck));
 dmins(i)=dmin/(numel(icheck));
 dmaxs(i)=dmax/(numel(icheck));
end


% sort vaccines based on average correlation (or energy)
ichecks=[ichecks allcta(:) alle2aa(:)];
icheckss=sortrows(ichecks,-(ncheck+1)) ;% sort by correlation descending order

% same thing with test strains :
ntest=numel(itest);
itests=[itests allcta(:) alle2aa(:)];
itestss=sortrows(itests,-(ntest+1)) ;% sort by correlation descending order

close
showpca ;
% plot a few top models
nbest=10;
for i=1:nbest
% ii=icheckss(i,1:ncheck);
 ii=itestss(i,1:ntest);
 if q2d
  scatter( pcc(1,ii), pcc(2,ii), ms, 'rv' ) 
 else
  scatter3( pcc(1,ii), pcc(2,ii), pcc(3,ii), ms, 'rv' ) 
 end
end
% plot a few bottom models
nworst=10
for i=1:nworst
% ii=icheckss(length(ichecks)-i+1, 1:ncheck);
 ii=itestss(length(itests)-i+1, 1:ntest);
 if q2d
  scatter( pcc(1,ii), pcc(2,ii), ms, 'bs' ) 
 else
  scatter3( pcc(1,ii), pcc(2,ii), pcc(3,ii), ms, 'bs' ) 
 end
end
