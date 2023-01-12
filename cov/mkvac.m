% assume that coordinates & msa are loaded

strains ;

% define vaccination cocktails per cohen 2021 science :

m4a={'sars-2', 'shc014', 'ratg13', 'rs4081'};
im4a=getind(m4a);

m4b={'wiv1', 'rf1', 'rmyn02', 'pang17'};
im4b=getind(m4b);

m8=[m4a m4b];
im8=getind(m8);

s2={'sars-2'};
is2=getind(s2);

% strains are in the order of appearance in the exp data file :
allags={'sars-2', 'ratg13', 'sars', 'wiv1', 'rs4081', 'shc014', 'yun11', 'bm4831', 'btky72'} ;
nallags = length(allags) ;
%
% (1) to select a fraction of antigens for model fitting/training :
agtrainf=0.5 ;
nagtrain = ceil ( agtrainf * nallags ) ; % number to take a a test vs. check strains
% to pick randomly from Virus cell array :
%iagtrain = randperm(nallags,nagtrain) ; fprintf('random\n');
%agtrain = allags(iagtrain) ;
%
% OR
% (2) to use all strains :
agtrain=allags ;
%
% compute indices of the test strains in the exp. data array
iagtrain=getind(agtrain)
%
% and/or can select certain vaccines for a train/test split :
vacs = { '1s', '4a', '4b', '8ab' } ; % vaccine names : homotypic, mosaic-4a, mosaic-4b, mosaic-8ab
%j1s=find(ismember(vacs,'1s')); % could in principle have different order ; here will be 1-4
%j4a=find(ismember(vacs,'4a'));
%j4b=find(ismember(vacs,'4b'));
%j8ab=find(ismember(vacs,'8ab'));
j1s=1 ; j4a=2 ; j4b=3 ; j8ab=4 ; % vaccine indices
%
vactrain=vacs; % take all vaccines for training ( but possibly a subset of ags, as above)
for iv=1:numel(vactrain)
 ivactrain(iv)=find(ismember(vacs,vactrain(iv)));
end
% or, could specify directly the vaccine indices in ivactrain
%
ivactrain
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
noxplot=1 ; % read exp data but do not plot
% read exp data from Cohen 2021 Science
eval(['run ',rootdir,'/exp/showc.m']);
%eval(['run ',rootdir,'/exp/show2.m']);
%%%% organize cohen exp data in a matrix for random access (to be used in train/test version split)
for ind=1:length(iggname)
 split=strsplit(char(iggname(ind)),'_') ;% split(1) has the virus ; split(2) has the vaccine
 ivir=getind(split(1));
 ivac= find(ismember(vacs,split(2)));
 iggmat(ivir, ivac) = iggexp(ind);
 iggemat(ivir, ivac) = iggexpe(ind);
end
% NOTE : not every matrix entry is populated because the vaccination iggs do not exist for some strains
%
%%%%
[nstrain,nres]=size(msamat) ;

