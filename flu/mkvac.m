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

allstrains=Virus;
nallstrains = length(allstrains) ;

% (1) to select a fraction for model fitting/training :
trainf=0.5 ;
ntrain = ceil ( trainf * nallstrains ) ; % number to take a a test vs. check strains
% pick randomly from Virus cell array :
%itrain = randperm(nallstrains,ntrain) ; fprintf('random\n');
%train = allstraint(itrain) ;
%
% OR
% (2) to use all strains :
train=allstrains ;

% compute indices of the test strains in the exp. data array
itrain=getind(train)

%%% read exp data :
noxplot=1 ; % do not plot exp data (yet)
%run ../exp/showf.m; % read data from plos one cohen paper
eval(['run ',rootdir,'/exp/showf.m']);
%%%% organize cohen exp data in a matrix for random access :
vacs = { 'm1', 'm2', 'm4', 'm8' } ;
%
jm1=find(ismember(vacs,'m1'));
jm2=find(ismember(vacs,'m2'));
jm4=find(ismember(vacs,'m4'));
jm8=find(ismember(vacs,'m8'));
%
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
