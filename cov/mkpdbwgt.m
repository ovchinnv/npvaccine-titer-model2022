% visualize model residue weights on a CoV-RBD structure
%
if ~exist('rootdir')
 rootdir='.';
end
addpath([rootdir,'/../']); %
%
pdbname='6vxx'; % structure on which to visualize weights
qinst=1 ; % whether to use instantaneous weights (vs averaged)
%
if (~exist('model'))
 model='dist2ave';
end
if (~exist('enc'))
 enc='gr';
end
%
pdbfile=[pdbname,'.pdb']; %
%
if (~exist('qpdb')) ; qpdb=0 ; end
%
if (~qpdb) % read pdb file
% global molecule;
 disp(['==>Reading pdb structure from file ',pdbfile,' ...']);
 molecule=pdbread(pdbfile);
 qpdb=1; %
end

pdb=molecule.Model.Atom;

% manipulate PDB format and write out pdb files for charmm
anum= [pdb.AtomSerNo]';
aname={pdb.AtomName}' ;
rname={pdb.resName}' ;
segid={pdb.segID}' ;
resid=[pdb.resSeq]' ;
insertion=char({pdb.iCode});
occu =[pdb.occupancy]' ;
bet  =[pdb.tempFactor]' ;
xpdb =[pdb.X]';
ypdb =[pdb.Y]';
zpdb =[pdb.Z]';
chainid=[pdb.chainID];
element=[pdb.element];
natom=length(pdb);
%
A=ismember(chainid,'A');
typeCA=ismember(aname,'CA');
%
% prepare relevant sequence :
seqinds =A(:) & typeCA(:) ;
run ntaa ;
pdbseq1 = cellfun ( aa1,rname(seqinds) ); % grab resname for seqinds & convert to 1-letter code
pdbseq1 = pdbseq1(:)' ; % make sure a single row
%

% load weights :
if (qinst) % instantaneous
 wnam=[model, '_', enc];
 d=load([wnam, '.mat']);
 wgt=d.wgt;
else
 fits=load([wnam,'.mat']); % many fits to be averaged :
 wgt=mean(allwgt,1);
end
%
seq0=d.msamat(1,:) ;

assert(length(wgt)==length(seq0)) ; % make sure we have the correct file pair
%
% to exclude dashes in seq0
ind0=find(seq0~='-');
seq0i=seq0(ind0);
wgti=wgt(ind0);
%
%[score,aln]=nwalign(seq0i,pdbseq1,'glocal',0);
[score,aln]=nwalign(seq0i,pdbseq1,'glocal',1,'gapopen',5);

% loop over pdb seq (position 2) and transfer wgt to structure
%
ind1=0;
ind3=0;
for i=1:length(aln)
 if (aln(1,i)~='-')
  ind1=ind1+1;
  w=wgti(ind1);
 else
%  w=-1;
  w=0;
 end
%
 if (aln(3,i)~='-')
  ind3=ind3+1;
  wgt3(ind3)=w;
 end
end % for
%
% now, wgt3 should be matched to seqinds
%
%erase all beta :
bet(:)=0;
% set weights
wscale=100 ; % improve resolution since the PDB format only has 2 decimal points
bet(seqinds)=wscale*wgt3;
%
% write pdb :
%
mol.Model.Atom=pdb ;% new molecule without headers
pdbout(mol, [wnam,'.pdb'], xpdb, ypdb, zpdb, occu, bet, find( A(:) & resid(:) >= 321 & resid(:) <= 592)) ;% hardwire rbd limits

