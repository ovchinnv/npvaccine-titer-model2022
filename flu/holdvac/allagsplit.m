% check all possible test/train splits (5 strains chosen for training)
%
rootdir='../';
addpath(rootdir,'-END');
%
nallags=10 ; % number of flu strains in Cohen et al. plos 21 paper
nallvacs=4 ;
% === to exclude a vaccine panel from training :
%mkvac ;
%holdvacs={'m8'} ;% leave empty for no holdout
%holdvacs={} ;% leave empty for no holdout
%for iv=1:numel(holdvacs)
% iholdvacs(iv)=find(ismember(vacs,holdvacs(iv)));
%end
% to avoid having to call mkvac
if ~exist('iholdvacs') ; iholdvacs=[4]; end ;
itrainvacs=setdiff(1:nallvacs,iholdvacs)  % vacs used to train

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
%
if (~exist('qtestholdvac')) % whether to test on the holdout vaccine set
 qtestholdvac=0;
end
qholdvac=exist('iholdvacs') ; if qholdvac ; qholdvac=(numel(iholdvacs)>0); end ;
qtestholdvac=(qtestholdvac && qholdvac);
%
tag='headstem';
%
qwrand=0 ; % constant initial weights
nrep=1 ; % if the weights are constant, there is no randomness so use only one replica
%
return ; % without computing
%
clear indvac indags vinds;
% select all possible 5-strain combinations from available strains (nallags) ; note that order does not matter (cocktail)
for i1=1:nallags
 vinds(1) = i1 ; % first strain
 for i2 = i1 + 1 : nallags
  vinds(2) = i2 ; %second strain
  for i3 = i2 + 1 : nallags
   vinds(3) = i3 ; %third strain
   for i4 = i3 + 1 : nallags
    vinds(4) = i4 ; %fourth strain
    for i5 = i4 + 1 : nallags
     vinds(5) = i5 ; %fifth strain

    for i6 = i5 + 1 : nallags
     vinds(6) = i6 ; %sixth strain

    for i7 = i6 + 1 : nallags
     vinds(7) = i7 ; %seventh strain


     vind=vind+1 % count the vaccines -- should be "10 choose 5" at the end
     indags(vind,:)=vinds ;
     indvac(vind,:)=itrainvacs ;
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
   allerrt=err(:)';
   allctpval=pval ;% p-value associated w/correlation coefficient
   allcstpval=spval ;
  else
   alle2at=[alle2at;e2a]; % append next sample
   allct=[allct;c];
   allcst=[allcst;cs];
   allerrt=[allerrt;err(:)'];
   allctpval=[allctpval pval] ;% p-value associated w/correlation coefficient
   allcstpval=[allcstpval spval] ;
  end
  alltestdata(isamp,:,:) = itestsample ; % save test sample & results as well
 end % run model
%  return
    end
    end
    end
   end
  end
 end
end

%
%savename=[model,'-',date,'-',enc,'.mat'];
if (qholdvac)
 if (qtestholdvac)
  savename=[model,'-',enc,'-',nptype,'-ags-htest.mat'];
 else
  savename=[model,'-',enc,'-',nptype,'-ags-htrain.mat'];
 end
else
 savename=[model,'-',enc,'-',nptype,'-ags.mat'];
end
%
save(savename, '-mat')

