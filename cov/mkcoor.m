% read fasta sequences for spike, isolate RBD, align, convert to numerical coordinates using an encoding
if ~exist('rootdir')
 rootdir='.';
end
%
addpath([rootdir,'/../']); % make sure aln2coor accessible
%
if ~exist('qgrantham')
 qgrantham=1; % grantham by default
end
%
if qgrantham==1
 enc='gr';
else
 enc='at';
end
%
if ~exist('qwritealn') ; qwritealn=0 ; end % whether to write alignment figures
if ~exist('qmsa') ; qmsa=0 ; end
% RBD sequence limits based on 6vxx
sbeg='qptes';
send='pcsf';
%
if (qmsa==1)
% first add the sequence of the 6VXX pdb (to isolate the RBD, as below)
f=fopen([rootdir,'/seqs/6vxx.fasta'],'r');
s=textscan(f,'%s','delimiter','\n');
s=s{1};
fclose(f);
%
seqs(1).Header = '>6VXX|PDB' ;
seq=s(2:end); % cell
seqs(1).Sequence = [seq{:}] ; %write out cell contents and concatenate

% read other sequences
[~,files]=system(['ls ',rootdir','/seqs/[A-Z]*.fasta']);
files=sort(strread(files,'%s')); % make sure to sort for compatibility between octave & matlab
%
i=1;
for file=files(:)'
 i=i+1;
 f=fopen(file{:},'r');
 s=textscan(f,'%s','delimiter','\n');
 fclose(f);
 s=s{1};
 seqs(i).Header = char(s(1));
 seqs(i).Sequence = char(s(2));
end

% to write complete fasta file (in case I want to align in a different program):
% fout=fopen('sars2-cov.fasta','w');
% seqwrite(fout,seqs,'FST');
%
msa=multialign(seqs); % not clear what/whether there are optimal parameters
save([rootdir,'/msa-cov.mat'],'msa', 'seqs') % save alignment
%
% write alignment
fid=fopen([rootdir,'/msa-cov.slx'],'w');
seqwrite(fid,msa,'slx')
fclose(fid);

qmsa=0;
return ; % quit after generating MSA

else % qmsa (whether to align)
 load([rootdir,'/msa-cov.mat']);
end
% assign strain names
for i=1:length(msa)
 h=msa(i).Header;
 ib=2;
 ie=find(h=='|',1)-1 ;
 if (isempty(ie)) ; ie=length(h); end
 msa(i).Name=h(ib:ie);
end

% extract alignment matrix:
msamat=char( {msa.Sequence}' ) ; % msa computed above or loaded from mat alignment file
s=msamat(1,:);
% but first read residue numbers that correspond to this strain (this -- s2vi -- version only)
fid=fopen([rootdir,'/seqs/6vxx.resid']);
d=textscan(fid,'%d %c') ;
fclose(fid);
resid=d{1};
s1=char(d{2})' ;
assert(1==strcmp(s1, strrep(s,'-',''))) % compare two sources , they should be the same
% now "transfer" the resids into the msa alignment (i.e. skipping gaps)
ires=0 ;
for imsa=1:numel(s)
 if( s(imsa)=='-' )
  rmsa(imsa)=-1;
 else
  ires=ires+1;
  assert(s(imsa)==s1(ires)) ;% make sure seqs align
  rmsa(imsa)=resid(ires);
 end
end

% cut between RBD sequence limits:
ibeg=strfind(s,upper(sbeg));
iend=strfind(s,upper(send));
% truncate msa to these limits; also exclude 1st strain (which was the PDB alignment template)
msamat=msamat(2:end,ibeg : iend + length(send) - 1 ) ;
% also truncate resids (include the search subsequences):
rmsa=rmsa(ibeg : iend + length(send) - 1 ) ;
% strain names :
names={msa(2:end).Name};
% to save truncated :
fid=fopen([rootdir,'/msa-rbd.slx'],'w');
for j=1:size(msamat,1)
 fprintf(fid,'%-29s %s\n', names{j}(1:min(29,length(names{j}))), msamat(j,:)) ; %left justify
end
fclose(fid);
%
% exclude missing residues :
pmissing = 0.9;
pmiss=mean((msamat=='-'),1);
imiss=find(pmiss>pmissing);
msamat(:,imiss)=''; % delete those residue positions from alignment matrix
%showalignment(msamat);

% convert to coordinates
coor=aln2coor(msamat,qgrantham);

if (qwritealn)
% write RBD alignment
% note that the first strain in the MSA is the PDB vequence 6vxx ; we first remove it from the MSA below
%
msa=msa(2:end);
alname='-cov-rbd';
run strains;
for i=1:numel(msa)
 msa(i).Header=Virus{i} ;
 msa(i).Sequence = msamat(i,:) ;
end
multialignwrite(['msa',alname,'.clu'],msa, 'header', 'Multiple Sequence Alignment of Coronavirus RBDs');
% also compute alignment score matrix
ascores=zeros(numel(msa));
for i=1:numel(msa)
 for j=i:numel(msa)
  [ascores(i,j),aln]=nwalign(strrep(msa(i).Sequence,'-',''), strrep(msa(j).Sequence,'-',''), 'GLOCAL', 0) ;
 end
% ascores(i,:)=ascores(i,:)/ascores(i,i)*100;
end
% normalize
for i=1:numel(msa)
 for j=i+1:numel(msa) % skip diag first
  ascores(i,j)=ascores(i,j)/sqrt(ascores(i,i)*ascores(j,j))*100;
 end
 ascores(i,i)=100; % diag
end

%pcolor(1-ascores) ; shading faceted ; colormap hot ;
mycolor([1:size(ascores,1)],[1:size(ascores,2)], 1-ascores) ; shading flat ; colormap hot ; box on ;
set(gca, 'xticklabel', Virus, 'yticklabel', Virus, 'tickdir','out', 'ticklength',[0,0],'xtick',[1:numel(msa)],'ytick',[1:numel(msa)],'XtickLabelRotation',90)
% write in similarity numbers :
for i=1:numel(msa)
 for j=i:numel(msa)
  sim=sprintf('%2.0f%',ascores(i,j));
  text(j-0.44,i,['\bf',sim], 'color','w', 'fontsize',12)
 end
end
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2',['msa',alname,'.eps']);
print(gcf,'-dpng',['msa',alname,'.png']);

end
