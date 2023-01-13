%
if (~exist('modver')) % sticky
 modver=1 ; % versions of f, the similarity function
end
%
mkcoor; % make numerical coordinates from MSA
mkvac ;%
%
% whether to normalize by exp. error :
qnorm=1;
wmin=0.0; % minimum allowed weight (e.g. negative weights are unphysical)
%
if (~exist('xint'))
 xint=1. ;
end
if (~exist('xp'))
% xp=2. ; % exponent of (d^2) in model ; (i.e., corresponds to half of the exponent of d)
 xp=2.5 ; % note that the landscape is quite flat for the range 2-4 ; so this exp produces a very similar, slightly better fit
end
if (~exist('Diff'))
 Diff=0.3 ; % diffusion constant in regularization ; 0 to turn off
end
%
format long
% initial weight assignment:
% this has an important effect
if(~exist('qwrand')) ; qwrand=0 ; end
if (qwrand)
 wgt=0.5*rand(1,nres) ;% % random
else
 if (~exist('wamp')) ; wamp=0.1 ; end
 wgt=wamp*ones(1,nres) ;% uniform
end
%
dwgt=zeros(1,nres) ;
%
ndim=size(coor,2)/nres; % components per residue
nvac=length(vacs) ; % number of "vaccines"
% compute average strain in all vaccine
clear vcoor;
vcoor(j1s,:)=mean(coor(is2,:) ,1);
vcoor(j4a,:)=mean(coor(im4a,:),1);
vcoor(j4b,:)=mean(coor(im4b,:),1);
vcoor(j8ab,:)=mean(coor(im8,:),1);

maxiter=3000 ; % iterations :
sdstep=0.04 ; % steepest descent step coefficient
%
iter=1;
besq=inf;
bcorr=-inf;
%
% set error normalization
if qnorm == 1
% scale : 
 scale=0;
 iscale=0;
 for itrain=1:size(itrainsample,1) % train samples
  ia=itrainsample(itrain,1); % antigen index
  iv=itrainsample(itrain,2); % vaccine index
  scale = scale + iggemat( ia, iv ) ;
  iscale = iscale + 1 ;
 end
 enmat = iggemat / scale * iscale ;
else
 enmat = ones(size(iggmat));
end
oenmat=1./enmat;
%
while 1 %do
wgt2=wgt.^2; % squared weights
%
ind=0;
dwgt(:)=0;

for itrain=1:size(itrainsample,1) % this is a marix with two columns, ag index in first column, vaccine index in second
  ia=itrainsample(itrain,1); % antigen
  iv=itrainsample(itrain,2); % vaccine
  ind=ind+1;
%
  dcoor=reshape(vcoor(iv,:)-coor(ia,:), ndim, []);
  ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor

  d2 = sum( wgt2 .* ndcoor2 );  % squared distance

  if (modver==1)
   iggmod(ind) = 1./(xint+d2^xp) ;  % model igg signal
  elseif (modver==2)
   iggmod(ind) = 1./(xint+d2)^xp;  % model igg signal
  end

  iggexp1(ind) = iggmat(ia,iv) ;
  oenorm(ind) = oenmat(ia,iv) ;

  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error

% contribution to gradient of error wrt weight
% need to omit singular points at d2=0;
  if (modver==1)
   if (d2>0)
    dwgt = dwgt + 2 * err(ind) * oenorm(ind) * (-iggmod(ind)^2)*xp*d2^(xp-1) * 2 .* wgt .* ndcoor2 ;
   end
  elseif(modver==2)
   if (d2>0)
    dwgt = dwgt + 2 * err(ind) * oenorm(ind) * (-xp*iggmod(ind)^(1+1/xp)) * 2 .* wgt .* ndcoor2 ;
   end
  end
%return
end

ibeg=1; % start at this row (>1 to omit an igg sample)
% correlation
corrs(iter)=corr( iggmod(ibeg:end)', iggexp1(ibeg:end)' );
e2s(iter)=sum(err(ibeg:end).^2);

if ( bcorr <= corrs(iter) ) % save parameters tha correspond to the best fit
 besqs=e2s(iter);
 bcorr=corrs(iter);
 iggmodb=iggmod;
 bestwgt=wgt;
 bestiter=iter;
end

% diffusion term as regularization
d2wgt=[0 diff(diff(wgt)) 0]; % zeros on the boundaries

wgt = max(wmin, wgt - sdstep * dwgt + Diff * d2wgt);
iter=iter+1;

if iter>maxiter;break;end;end%until iter>maxiter
%wgt1=wgt ;% save weights

%iggmod2t = iggmod ; % take final model values
iggmod2t = iggmodb ; % take best model values

c=corr(iggmod2t(ibeg:end)', iggexp1(ibeg:end)')
if exist('OCTAVE_VERSION')
 cs=spearman(iggmod2t(ibeg:end)', iggexp1(ibeg:end)')
else % matlab
 cs=corr(iggmod2t(ibeg:end)', iggexp1(ibeg:end)', 'type', 'spearman')
end
err2=(iggmod2t(:) - iggexp1(:)).^2;
e2=sum(err2(ibeg:end))
e2a=mean(err2(ibeg:end))

igg=[iggmod2t(:) iggexp1(:)] ;

%clear gind1 gxpind1 % to avoid save error in octave 
save(['dist2ave_',enc,'.mat'], '-mat');

% write out weights
%
fwgt=fopen(['dist2ave_',enc,'_wgt.dat'],'w');
for i=1:nres
   fprintf(fwgt, '%d %f\n', i, bestwgt(i));
end
%
fclose(fwgt);
%
