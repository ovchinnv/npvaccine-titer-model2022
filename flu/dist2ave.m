%
if (~exist('modver')) % sticky
 modver=1 ; % versions of f, the similarity function
end
%
if (~exist('tag'))
 tag='headstem';
end

if (~exist('nptype'))
 nptype='mosaic';
end

mkcoor ;
mkvac ;

qnorm=1; % decide whether to normalize by exp. error :
wmin=0.0; % minimum allowed residue weight (e.g. to disallow negative weights)
%
if (~exist('xint'))
 xint=1 ;
end
if (~exist('xp'))
 xp=2 ;
end
if (~exist('Diff'))
 Diff=0.45 ; % diffusion constant in regularization ; 0 to turn off
end
format long ;
%
% assign initial weights (this has an important effect)
if(~exist('qwrand')) ; qwrand=0 ; end
if (qwrand)
 wgt=0.5*rand(1,nres) ;% % random
else
 if (~exist('wamp')) ; wamp=0.1 ; end
 wgt=wamp*ones(1,nres) ;% uniform
end
dwgt=zeros(1,nres) ;

ndim=size(coor,2)/nres; % components per residue
nvac=length(vacs) ; % number of "vaccines"
% compute average strain in all vaccine
clear vcoor;
vcoor(jm1,:)=mean(coor(im1,:),1);
vcoor(jm2,:)=mean(coor(im2,:),1);
vcoor(jm4,:)=mean(coor(im4,:),1);
vcoor(jm8,:)=mean(coor(im8,:),1);

if (~exist('maxiter'))
 maxiter=3000 ; % number of iterations :
end
if (~exist('sdstep'))
 sdstep=0.03 ; % steepest descent step coefficient
end

iter=1;
besq=inf;
bcorr=-inf;
% set error normalization
if qnorm == 1
% scale :
 scale=0;
 iscale=0;
 for i=1:length(train) % train strains
  for j=1:nvac % all vaccines, for now
   scale = scale + iggemat( itrain(i), j ) ;
   iscale = iscale + 1 ;
  end
 end
 enmat = iggemat / scale * iscale ;
else
 enmat = ones(size(iggmat));
end
%oenorm=1./enorm;
oenmat=1./enmat;
%
%
while 1 %do
wgt2=wgt.^2; % squared weights
ind=0; % vaccine index; needs to match the order of vaccines in the experimental dataset of Cohen21 PLoS
dwgt(:)=0;
for i=1:length(train) % all train strains
 for j=1:nvac % all vaccines
  ind=ind+1;

  dcoor=reshape(vcoor(j,:)-coor(itrain(i),:), ndim, []);
  ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor

  d2 = sum( wgt2 .* ndcoor2 );  % squared distance

  if (modver==1)
   iggmod(ind) = 1./(xint+d2^xp) ;  % model igg signal
  elseif (modver==2)
   iggmod(ind) = 1./(xint+d2)^xp;  % model igg signal
  end

  iggexp1(ind) = iggmat(itrain(i),j) ;
  oenorm(ind) = oenmat(itrain(i),j) ;

  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error

% contribution to gradient of error wrt weight
% first coded as a close translation of derived formula
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
end

ibeg=1; % start at this row
% correlation
corrs(iter)=corr(iggmod(ibeg:end)', iggexp1(ibeg:end)');
e2s(iter)=sum(err(ibeg:end).^2);

%if ( besq >= e2s(iter) ) % optimal squared error
if ( bcorr <= corrs(iter) ) % optimal correlation coefficient
 besq=e2s(iter);
 bcorr=corrs(iter);
 iggmodb=iggmod;
 bestwgt=wgt;
 bestiter=iter;
end

% include a diffusion term as ad hoc regularization
d2wgt=[0 diff(diff(wgt)) 0]; % zeros on the boundaries
wgt = max(wmin, wgt - sdstep * dwgt + Diff * d2wgt); % weights optimization step
iter=iter+1;

if iter>maxiter;break;end;end ;%until iter>maxiter
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

igg=[iggmodt(:) iggexp1(:)] ;

save(['dist2ave_',tag,'_',enc,'.mat'], '-mat');

% write out weights
%
fwgt=fopen(['dist2ave_',tag,'_',enc,'_wgt.dat'],'w');
for i=1:nres
   fprintf(fwgt, '%d %f\n', i, wgt(i));
end
%
fclose(fwgt);
