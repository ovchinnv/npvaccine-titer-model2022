% in this version, using not distance to average strain, but average strain distance
%
model='avedist';
%
if (~exist('tag'))
 tag='headstem';
end

if (~exist('nptype'))
 nptype='mosaic';
end
%
mkcoor ;
mkvac ;
qnorm=1 ;
wmin=0 ; % minimum allowed weight

if (~exist('xint'))
 xint=1 ;
end
if (~exist('xp'))
 xp=3.0 ;
end
if (~exist('Diff'))
 Diff=0.49 ; % diffusion constant in regularization ; 0 to turn off
end
format long ;
xp
xint
Diff

% this has an important effect
%wgt=1*ones(1,nres) ;% uniform
%wgt=0.25*ones(1,nres) ;% uniform
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
% compute average strain in all vaccine
vaccines={ im1, im2, im4, im8 };
nvac=numel(vaccines);

% iterations :
maxiter=3000 ;
sdstep=0.04 ; % steepest descent step coefficient
iter=1;
besq=inf;
bcorr=-1;
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
% precompute ndcoor2 for faster iterations : (NOT MUCH FASTER !)
ndcoor2=zeros(nres,nstrain,nstrain);
for i=1:nstrain % all train strains
 for j=i+1:nstrain % all strains
  dcoor=reshape(coor(i,:)-coor(j,:), ndim, []);
  ndcoor2(:,i,j)=sum(dcoor.^2,1); % squared norm
  ndcoor2(:,j,i)=ndcoor2(:,i,j);
 end
end

while 1 %do

wgt2=wgt.^2; % squared weights
% to repeat weights for all coor components :
%wgts=reshape(ones(ndim,1)*wgt,1,[]) ;
% to weight coordinates :
% wcoor=repmat(wgts,nstrain,1).*coor ;
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

   d = sqrt( sum( wgt2 .* ndcoor2(:,k,ia)' ) );  % distance to this strain

   dave = dave + d ;

% contribution to gradient of error wrt weight
% to omit points for which d2=0; (depending on model could be singular or very large)
   if (d>0)
    dwgt_this = dwgt_this +  wgt .* ndcoor2(:,k,ia)' / d ;
   end
  end % all strains in the vaccine

  nstr=numel(vaccines{iv}); % normalization

  dave = dave / nstr ; % mean squared distance between the train strain and all vaccine strains
% compute model value :
  ind=ind+1;
  iggmod(ind) = 1./(xint+dave^xp);  % model igg signal

  iggexp1(ind) = iggmat(ia,iv) ;
  oenorm(ind) = oenmat(ia,iv) ;

  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error

% gradients :
  dwgt_this = dwgt_this / nstr ;

  if (dave>0)
   dwgt = dwgt + 2 * err(ind) * (-iggmod(ind)^2)*xp*dave^(xp-1) * dwgt_this ; % add contribution from ths train strain
  end
end

ibeg=1; % start at this row
% correlation
corrs(iter)=corr(iggmod(ibeg:end)', iggexp1(ibeg:end)');
e2s(iter)=sum(err(ibeg:end).^2);

%if ( besq >= e2s(iter) )
if ( bcorr <= corrs(iter) )
 besq=e2s(iter);
 bcorr=corrs(iter);
 iggmodb=iggmod;
 bestwgt=wgt;
 bestiter=iter;
end

% diffusion term as reg
d2wgt=[0 diff(diff(wgt)) 0]; % zeros on the boundaries
wgt = max(wmin, wgt - sdstep * dwgt + Diff * d2wgt);

iter=iter+1;
if iter>maxiter;break;end;end;%until iter>maxiter

iggmodt = iggmodb ; % output best model results
%iggmodt = iggmod ; % output last model results
%
c=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)') % Pearson correlation between exp and model
if exist('OCTAVE_VERSION')
 cs=spearman(iggmodt(ibeg:end)', iggexp1(ibeg:end)') % Spearman correlation
else % matlab
 cs=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)', 'type', 'spearman')
end
err2=(iggmodt(:) - iggexp1(:)).^2; % error for each titer
e2=sum(err2(ibeg:end)) % total squared error
e2a=mean(err2(ibeg:end)) % mean squared error

save(['avedist_',tag,'_',enc,'.mat'], '-mat');

% write out weights
%
fwgt=fopen(['avedist_',tag,'_',enc,'_wgt.dat'],'w');
for i=1:nres
   fprintf(fwgt, '%d %f\n', i, wgt(i));
end
%
fclose(fwgt);
