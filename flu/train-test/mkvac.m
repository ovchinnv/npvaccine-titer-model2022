if ~exist('rootdir')
 rootdir='.';
end
% load array of virus strains corresponding to the order of exp data of cohen et al 21. PLoS 
strains ;

% define vaccination cocktails per cohen 2021:

m1={'c09'}; % 1-antigen np
im1=getind(m1); % compute strain index

m2={'ai68', 'c09'}; % 2-antingen np
im2=getind(m2);

m4={'ai68', 'c09', 'v04', 'sh13'}; % 4-antigen np
im4=getind(m4);

m8={'ai68', 'c09', 'v04', 'sh13', 'j57', 'hk99', 'jx13', 'hb09'}; % 8-antigen np
im8=getind(m8);

allags=Virus;
nallags = length(allags) ;

% (1) to select a fraction of antigens for model fitting/training :
agtrainf=0.5 ;
nagtrain = ceil ( agtrainf * nallags ) ; % number to take a a test vs. check strains
% to pick randomly from Virus cell array :
%iagtrain = randperm(nallags,nagtrain) ; fprintf('random\n');
%
% OR
% (2) to use all strains :
%agtrain=allags ;
% OR
% (3) custom set :
iagtrain=indags(vind,:)  % from allagsplit
%
agtrain = allags(iagtrain) ; % recompute names from indices
%
% compute indices of the test strains in the exp. data array
% and/or can select certain vaccines for a train/test split :
vacs = { 'm1', 'm2', 'm4', 'm8' } ; % available vaccines
jm1=1 ; jm2=2 ; jm4=3 ; jm8=4 ; % vaccien indices

%vactrain=vacs; % take all vaccines for training ( but possibly a subset of ags, as above)
%for iv=1:numel(vactrain)
% ivactrain(iv)=find(ismember(vacs,vactrain(iv)));
%end
% or, could specify directly the vaccine indices in ivactrain
%
ivactrain=indvac(vind,:)

vactrain=vacs(ivactrain);
% build the training data matrix :
ind=0;
for ia=iagtrain(:)'
 for iv=ivactrain(:)'
  ind=ind+1;
  itrainsample(ind,1)=ia;
  itrainsample(ind,2)=iv;
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% read experimental vaccine titer data :
noxplot=1 ; % do not plot exp data
%run ../exp/showf.m; % read data from plos one cohen paper
eval(['run ',rootdir,'/exp/showf.m']);
%%%% organize cohen exp data in a matrix for convenient random access :
for ind=1:length(iggname)
 split=strsplit(char(iggname(ind)),'_') ;% split(1) has the virus ; split(2) has the vaccine
 ivir=getind(split(1));
% split(1)
 ivac= find(ismember(vacs,split(2)));
 iggmat(ivir, ivac) = iggexp(ind);
 iggemat(ivir, ivac) = iggexpe(ind);
end

%
[nstrain,nres]=size(msamat) ;
