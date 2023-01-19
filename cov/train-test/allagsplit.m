% check all possible train/test splits with 5 strains chosen for training
%
rootdir='../';
addpath(rootdir);

nrep=1 ;
nallags=9 ;
nallvacs=4 ;
vind = 0 ; % vaccine index
isamp = 0 ;% sample index

if (~exist('xint'))
 xint=1 ;
end
if (~exist('xp'))
 xp=2 ;
end
if (~exist('Diff'))
 Diff=0.3 ; % diffusion constant in regularization ; 0 to turn off
end
if (~exist('model'))
 model='dist2ave' ; % name of matlab script file that performs the fit
end
%
qwrand=0;
wamp=0.1 ;
%
return
%
clear indvac indags vinds;
%
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
     vind=vind+1 % count the vaccines
     indags(vind,:)=vinds ;
     indvac(vind,:)=1:nallvacs ;
% run model
% repeat nrep times :
 for irep=1:nrep
  isamp=isamp+1; % number of samples
  eval(model) ; % train
  if (isamp==1)
   allwgt1=bestwgt ;
   alle2a=e2a ;
   allc=c ;
   allcs=cs ;
   allcpval=pval ;% p-value associated w/correlation coefficient
   allcspval=spval ;
  else
   allwgt1=[allwgt1;bestwgt] ; %append to the end
   alle2a=[alle2a;e2a];
   allc=[allc;c];
   allcs=[allcs;cs];
   allcpval=[allcpval pval] ;% p-value associated w/correlation coefficient
   allcspval=[allcspval spval] ;
  end
  run test; % test
  if (isamp==1)
   alle2at=e2a ;
   allct=c ;
   allcst=cs ;
   allctpval=pval ;% p-value associated w/correlation coefficient
   allcstpval=spval ;
  else
   alle2at=[alle2at;e2a];
   allct=[allct;c];
   allcst=[allcst;cs];
   allctpval=[allctpval pval] ;% p-value associated w/correlation coefficient
   allcstpval=[allcstpval spval] ;
  end
  alltestdata(isamp,:,:) = itestsample ; % save test sample & results as well
 end % run model
%return
%
    end
   end
  end
 end
end
%
savename=[model,'-',enc,'-ags.mat'];
save(savename, '-mat')

