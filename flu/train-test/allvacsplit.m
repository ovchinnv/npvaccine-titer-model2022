% check all possible test/train _vaccine_ (not AG) splits ; adapted from allagssplit
%
rootdir='../';
addpath(rootdir,'-END');
%
nallags=10; % number of ags
nallvacs=4 ; % number of vaccines in Cohen et al. plos 21 paper
vind = 0 ; % vaccine index
isamp = 0 ;% sample index

if (~exist('xint'))
 xint=1 ;
end
if (~exist('xp'))
 xp=1.9 ;
end
if (~exist('Diff'))
 Diff=0.49 ; % diffusion constant in regularization ; 0 to turn off
end
if (~exist('model'))
 model='dist2ave' ; % name of matlab script file that performs the fit
end
if (~exist('nptype'))
 nptype='mosaic';
end

tag='headstem';
%
qwrand=0 ; % constant initial weights
nrep=1 ; % if the weights are constant, there is no randomness so use only one replica
%
return
%
clear indvac indags vinds;
% select all possible 5-strain combinations from available strains (nalltest) ; note that order does not matter (cocktail)
for i1=1 : nallvacs
 vinds(1) = i1 ; % first vaccine
 for i2 = i1 + 1 : nallvacs
  vinds(2) = i2 ; %second vaccine
  vind=vind+1 % index of vaccination experiment
  indvac(vind,:)=vinds ;
  indags(vind,:)=1:nallags ;
% run model
% repeat nrep times :
 for irep=1:nrep
  isamp=isamp+1; % number of samples
  eval(model) ; % train
  if (isamp==1) % set (1st sample)
   allwgt=bestwgt ;
   alle2a=e2a ;
   allc=c ; % pearson correlation
   allcs=cs ; % spearman correlation
   allcpval=pval ;% p-value associated w/correlation coefficient
   allcspval=spval ;
  else
   allwgt=[allwgt;bestwgt] ; % append this sample
   alle2a=[alle2a;e2a];
   allc=[allc;c];
   allcs=[allcs;cs];
   allcpval=[allcpval pval] ;% p-value associated w/correlation coefficient
   allcspval=[allcspval spval] ;
  end
  run test; % test
  if (isamp==1) % set 1st sample
   alle2at=e2a ;
   allct=c ;
   allcst=cs ;
   allctpval=pval ;% p-value associated w/correlation coefficient
   allcstpval=spval ;
  else
   alle2at=[alle2at;e2a]; % append next sample
   allct=[allct;c];
   allcst=[allcst;cs];
   allctpval=[allctpval pval] ;% p-value associated w/correlation coefficient
   allcstpval=[allcstpval spval] ;
  end
  alltestdata(isamp,:,:) = itestsample ; % save test sample & results as well
 end % run model nrep times
%  return
 end
end

%
savename=[model,'-',enc,'-',nptype,'-vac.mat'];
save(savename, '-mat')
%
