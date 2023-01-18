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
% create test sample :
% take all those not present in fitting set
iall_mat=ones(nallags,numel(vacs));
itrainsample_mat=sparse(itrainsample(:,1), itrainsample(:,2), 1, nallags, numel(vacs));
itestsample_mat=iall_mat-itrainsample_mat ;
[ias,ivs]=find(itestsample_mat);
itestsample=[ias,ivs];
%
%
if (isempty(itestsample))
 itestsample=itrainsample ; % use all test strains
end
agtest=unique(itestsample(:,1)) ;
vactest=unique(itestsample(:,2)) ;
%
cols=repmat(cols, numel(agtest), 1);

if qnorm == 1
% scale : 
 scale=0;
 iscale=0;
 for itest=1:size(itestsample,1) % test strains
  ia=itestsample(itest,1); % antigen
  iv=itestsample(itest,2); % vaccine
  scale = scale + iggemat( ia, iv ) ;
  iscale = iscale + 1 ;
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
for itest=1:size(itestsample,1) % this is a marix with two columns, ag index in first column, vaccine index in second
  ia=itestsample(itest,1); % antigen
  iv=itestsample(itest,2); % vaccine
%
  if (strcmp(model,'dist2ave')) % distance to average model
   dcoor=reshape(vcoor(iv,:)-coor(ia,:), ndim, []);
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
   for k=vaccines{iv}

    dcoor=reshape(coor(k,:)-coor(ia,:), ndim, []);
    ndcoor2=sum(dcoor.^2,1); % squared norm of dcoor

    d = sqrt( sum( wgt2 .* ndcoor2 ) );  % distance to this strain
    dave = dave + d ;
   end % for k in vaccine
   nstr=numel(vaccines{iv}); % normalization
   dave = dave / nstr; % mean squared distance between the train strain and all vaccine strains
%
   ind=ind+1;
% model value :
   iggmod(ind) = 1./(xint+dave^xp);  % model igg signal
  end % if
%
  iggexp1(ind) = iggmat(ia,iv) ;
  iggexpe1(ind) = iggemat(ia,iv) ;
  oenorm(ind) = oenmat(ia,iv) ;
  err(ind) = ( iggmod(ind) - iggexp1(ind) ) * oenorm(ind) ; % model error
% also save testsample matrix, with results added :
  itestsample(itest,3)=iggmod(ind) ;
  itestsample(itest,4)=iggexp1(ind) ;
end

ibeg=1; % start at this row
iggmodt=iggmod ;
if exist('OCTAVE_VERSION')
 c=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
 cs=spearman(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
else % matlab
 [c,pval]=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)')
 [cs,spval]=corr(iggmodt(ibeg:end)', iggexp1(ibeg:end)', 'type', 'spearman')
end
err2=(iggmodt(:) - iggexp1(:)).^2;
e2=sum(err2(ibeg:end))
e2a=mean(err2(ibeg:end)) % mean squared error

