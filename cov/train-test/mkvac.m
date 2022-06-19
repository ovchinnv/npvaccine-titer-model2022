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
allstrains={'sars-2', 'ratg13', 'sars', 'wiv1', 'rs4081', 'shc014', 'yun11', 'bm4831', 'btky72'} ;
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
%train=allstrains ;

itrain=indvac(vind,:)  % from allvacs
train = allstrains(itrain)

% compute indices of the test strains in the exp. data array
itrain=getind(train)

%%% read exp data %%
noxplot=1 ; % read exp data but do not plot
%run exp/show.m
%run ./exp/show.m ;% data from Cohen 2021 Science
eval(['run ',rootdir,'/exp/showc.m']);
%%%% organize cohen exp data in a matrix for random access (to be used in train/test version split)
vacs = { '1s', '4a', '4b', '8ab' } ; % vaccine names : homotypic, mosaic-4a, mosaic-4b, mosaic-8ab
%
j1s=find(ismember(vacs,'1s')); % could in principle have different order ; here will be 1-4
j4a=find(ismember(vacs,'4a'));
j4b=find(ismember(vacs,'4b'));
j8ab=find(ismember(vacs,'8ab'));
%
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
%
