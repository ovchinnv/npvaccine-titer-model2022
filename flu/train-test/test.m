% check model fit
%
if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot') ;
end
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

% reload model fit
%load([model,'_',tag,'_',enc,'.mat']); % note that this could overwrite some variables, like isamp from allvacs
%
close all;
lims=[-0.25 1.75];
lims=[-0.25 1.25];

cols=[ 255 92 103 ; 0 255 0 ; 255 165 0 ; 0 0 250 ]/255;
%
check=allstrains ;
% take those not present in fitting set
check=setdiff(allstrains, train);
%
icheck=getind(check);
ncheck=length(check);
%
cols=repmat(cols, ncheck, 1);

if qnorm == 1
% scale : 
 scale=0;
 iscale=0;
 for i=1:length(check) % test strains
  for j=1:nvac % all vaccines, for now
   scale = scale + iggemat( icheck(i), j ) ;
   iscale = iscale + 1 ;
  end
 end
 enmat = iggemat / scale * iscale ;
else
 enmat = ones(size(iggmat));
end
oenmat=1./enmat;
%
clear iggexp1 iggexpe1 iggmod
wgt2=bestwgt.^2; % squared weights
%
ind=0;
for i=1:length(check) % all test strains
 for j=1:nvac % all vaccines
%
  if (strcmp(model,'dist2ave')) % distance to average model
   dcoor=reshape(vcoor(j,:)-coor(icheck(i),:), ndim, []);
   ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor
%
   d2 = sum( wgt2 .* ndcoor2 );  % squared distance
%
   ind=ind+1;
%
   if (modver==1)
    iggmod(ind) = 1./(xint+d2^xp) ;  % model igg signal
   elseif (modver==2)
    iggmod(ind) = 1./(xint+d2)^xp;  % model igg signal
   end
%
  elseif (strcmp(model, 'avedist')) % average distance model
%
   dave=0 ;% average distance
   for k=vaccines{j}

    dcoor=reshape(coor(k,:)-coor(icheck(i),:), ndim, []);
    ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor

    d = sqrt( sum( wgt2 .* ndcoor2 ) );  % distance to this strain
    dave = dave + d ;
   end % for k in vaccine
   nstr=numel(vaccines{j}); % normalization
   dave = dave / nstr ; % mean squared distance between the train strain and all vaccine strains
%
   ind=ind+1;
% model value :
   iggmod(ind) = 1./(xint+dave^xp);  % model igg signal
  end % if
%
  iggexp1(ind) = iggmat(icheck(i),j) ;
  iggexpe1(ind) = iggemat(icheck(i),j) ;
  oenorm(ind) = oenmat(icheck(i),j) ;
  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error
 end
end

ibeg=1; % start at this row
iggmodt=iggmod ;
c=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
if exist('OCTAVE_VERSION')
 cs=spearman(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
else % matlab
 cs=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)', 'type', 'spearman')
end
err2=(iggmodt(:) - iggexp1(:)).^2;
e2=sum(err2(ibeg:end))
e2a=mean(err2(ibeg:end)) % mean squared error

