mkcoor;
mkvac ;
qnorm=1; % whether to normalize by exp. error :
wmin=0.0; % minimum weight
%
model='avedist';
%
if (~exist('xint'))
 xint=0.8 ;
end
if (~exist('xp'))
 xp=3.8 ;
end
if (~exist('Diff'))
 Diff=0.3 ; % diffusion constant in regularization ; 0 to turn off
end
%
format long
%
% assign initial weights
% this has an important effect
if(~exist('qwrand')) ; qwrand=0 ; end
if (qwrand)
 wgt=0.5*rand(1,nres) ;% % random
else
 if (~exist('wamp')) ; wamp=0.1 ; end
 wgt=wamp*ones(1,nres) ;% uniform
end
dwgt=zeros(1,nres) ;
dwgt_this=zeros(1,nres) ;

ndim=size(coor,2)/nres; % components per residue
nvac=length(vacs) ; % number of "vaccines"

clear vcoor;
vaccines={ is2, im4a, im4b, im8 };
nvac=numel(vaccines);

% precompute ndcoor2 for faster iterations : (NOT MUCH FASTER !)
ndcoor2=zeros(nres,nstrain,nstrain);
for i=1:nstrain % all test strains
 for j=i+1:nstrain % all strains
  dcoor=reshape(coor(i,:)-coor(j,:), ndim, []);
  ndcoor2(:,i,j)=sum(dcoor.^2,1); % squared norm
  ndcoor2(:,j,i)=ndcoor2(:,i,j);
 end
end

% iterations :
maxiter=3000 ;
sdstep=0.04 ; % steepest descent step coefficient
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
% enorm = ones(size(iggexp));
 enmat = ones(size(iggmat));
end
%oenorm=1./enorm;
oenmat=1./enmat;
%
while 1 %do
wgt2=wgt.^2; % squared weights
%
ind=0;
dwgt(:)=0;
for itrain=1:size(itrainsample,1) % train samples
  ia=itrainsample(itrain,1); % antigen index
  iv=itrainsample(itrain,2); % vaccine index
%
  dave=0 ;% init average distance
  dwgt_this(:)=0;
  for k=vaccines{iv}
%   dcoor=reshape(coor(k,:)-coor(ia,:), ndim, []);
%   ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor
   d = sqrt ( sum( wgt2 .* ndcoor2(:,k,ia)' ) ) ;  % distance to this strain

   dave = dave + d ;

% contribution to gradient of error wrt weight
% to omit points at d=0; (depending on model could be singular or very large)
   if (d>0)
    dwgt_this = dwgt_this + wgt .* ndcoor2(:,k,ia)' / d ;
   end
  end % all strains in the vaccine

  nstr=numel(vaccines{iv}); % normalization

  dave = dave / nstr ; % mean squared distance between the train strain and all vaccine strains
% compute model value :
  ind=ind+1;
  iggmod(ind) = 1./(xint+dave^xp);  % model igg signal
% experimental value :
  iggexp1(ind) = iggmat(ia,iv) ;
  oenorm(ind) = oenmat(ia,iv) ;

  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error for this strain & vaccine

% gradients :
  dwgt_this = dwgt_this / nstr ;

  if (dave>0)
   dwgt = dwgt + 2 * err(ind) * (-iggmod(ind)^2)*xp*dave^(xp-1) * dwgt_this ; % add contribution from ths train strain
  end
end % over training samples

ibeg=1; % start at this row ( ibeg>1 to omit data for fitting )
% correlation
corrs(iter)=corr( iggmod(ibeg:end)', iggexp1(ibeg:end)' );
e2s(iter)=sum(err(ibeg:end).^2);

if ( bcorr <= corrs(iter) )
 besq=e2s(iter);
 bcorr=corrs(iter);
 iggmodb=iggmod;
 bestwgt=wgt;
 bestiter=iter;
end
%
% diffusion term as reg:
d2wgt=[0 diff(diff(wgt)) 0]; % zeros on the boundaries
wgt = max(wmin, wgt - sdstep * dwgt + Diff * d2wgt);
%
iter=iter+1;

if iter>maxiter;break;end;end%until iter>maxiter

iggmodt = iggmodb ; % best model

c=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
if exist('OCTAVE_VERSION')
 cs=spearman(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
else % matlab
 cs=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)', 'type', 'spearman')
end
err2=(iggmodt(:) - iggexp1(:)).^2;
e2=sum(err2(ibeg:end))
e2a=mean(err2(ibeg:end))

save(['avedist_',enc,'.mat'], '-mat');
%
% write out weights
%
fwgt=fopen(['avedist_',enc,'_wgt.dat'],'w');
for i=1:nres
   fprintf(fwgt, '%d %f\n', i, bestwgt(i));
end
%
fclose(fwgt);

