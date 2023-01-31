% make sure sqlite3.mex can be found:
% download at https://github.com/rmartinjak/mex-sqlite3, Copyright (c) by 2014 Robin Martinjak
addpath('../mex-sqlite3','../')
%
% make sure our custom flu_strains database can be found
%
db='flu_strains.db'
if ~exist('tab')
 tab='' ; % table in the strains database (90, 95, 97-99, or empty for full db)
end
%
if ~exist('qwritealn') ; qwritealn=0 ; end % whether to write alignment figures
%== align with a sequence to isolate stem & head  ==
cmd='select sequence from flu_strains where strain_id = ''ACT22035.1'''; % this seq should be ID to CA/09
ref=sqlite3(db,cmd);
seq0=ref.sequence ;
%
% specify sequence intervals and labels
%
default_flag = 0 ; % how to mark residues that are not part of the specification below
%
domains={ % just to be clear : the interior at the height of the head (i.e. HA2 helices) still considered part of the stem)
% head :
 'EDKH', 'KCPK', 1 ;
% stem :
 'MKAILV', 'VNLL',  2 ;
 'YVKS', 'KLEST', 2
% tm :
 'RIYQILA', -1 , 3
}
%
ndom=size(domains,1) ;
% create a flag aray length-matched to seq
flag = default_flag * ones( size (seq0 ) );
% populate flag array
for idom=1:ndom
 aa=domains(idom, 1) ; % beg
 if all( cell2mat(aa) == -1 )
  ibeg=1;
 else
  ibeg=strfind(seq0,aa) ;
 end
%
 bb=domains(idom, 2) ; % end
 if all( cell2mat(bb) == -1 )
  iend=length(seq0);
 else
  iend=strfind(seq0,bb) + length(char(bb))-1 ;
 end
 mark=cell2mat(domains(idom,3));
 [ibeg iend mark]
% now mark, but check for prior marks (which indicates overlap)
 if (any(flag(ibeg:iend)~=default_flag))
  warning domain overlap detected
 end
 flag(ibeg:iend)=mark;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now grab strains corresponding to the data in Cohen et al 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names={'California/04/2009', 'nam/1203/2004', 'apan/305/1957' , 'WF10', 'Aichi/%2/1968', 'Shanghai/1/2013', 'IPB13/2013', ...
    'HuBei/06/2009', 'Peru/033/2010', '2576/1979'  }

table=['flu_strains',num2str(tab)];

% use pairwise alignment to align domains rather than whole seqs

%alnflags=[2]; % aas with these flags will be used for alignment
%alname='stem' ; % domain name for file output

%alnflags=[2 3];
%alname='stem-tm' ;

%alnflags=[1];
%alname='head' ;

%alnflags=[0 1 2 3];
%alname='full' ;

alnflags=[1 2]; % head and stem
alname='headstem';


okflags=ismember(flag, alnflags)  %binary array which marks selected residues
okflags=diff([0 okflags 0]) % find bounds ; note that aray is 1 larger now -- that's deliberate to ID the intervals : see below
ibegs=find(okflags==1);
iends=find(okflags==-1)-1;

istr=0;
for name=names
 istr=istr+1;
 str=['''%',char(name),'%'''] ;
 cmd= ['select * from ',table,' where name LIKE ',str,' and length(sequence)>550 and Nunk=0 and sequence not LIKE ''%J%''']
% cmd= ['select * from ',table,' where name LIKE ',str];
 strains=sqlite3(db, cmd); % in case of multiple matches, can take the first
 full(istr).Header = strains(1).strain_id ;
 full(istr).Name = strains(1).name ;
 seq=strains(1).sequence ;
 full(istr).Sequence = seq ;

% align to template seq0 :
 [score,aln]=nwalign(seq0, seq, 'GLOCAL', 0) ;
% transfer flag array onto alignment (this is simple, just omit hyphens):
 ind=0;
 for i=1:length(aln)
  if(aln(1,i)~='-')
   ind=ind+1;
   ialn(i)=ind; %index map
  else
   ialn(i)=-999 ;
  end
 end

% build the sequence for MSAs
 assert(all(size(ibegs)==size(iends)))
 tmpseq='';
 for i=1:length(ibegs)
  ibeg=ibegs(i) ; iend=iends(i) ;
  tmpseq=[tmpseq,aln(3,find(ialn==ibeg):find(ialn==iend))];
 end
%
 seqs(istr).Sequence = strrep(tmpseq,'-',''); % create sequence without deleted parts
 seqs(istr).Name = full(istr).Name;
 seqs(istr).Header = full(istr).Header;
end
%
% align and save
msa=multialign(seqs); % note that am not using any custom alignment params here
save(['msa',alname,'.mat'],'msa') ;

% save in text format
if (qwritealn)
% set header to abbreviated names, for simplicity
clear strains; run strains;
for i=1:numel(msa)
 msa(i).Header=Virus{i} ;
end
multialignwrite(['msa',alname,'.clu'],msa, 'header', 'Multiple Sequence Alignment of Influenza Hemagglutinin Ectodomains');
% also compute alignment score matrix
ascores=zeros(numel(seqs));
for i=1:numel(seqs)
 for j=i:numel(seqs)
  [ascores(i,j),aln]=nwalign(seqs(i), seqs(j), 'GLOCAL', 0) ;
 end
% ascores(i,:)=ascores(i,:)/ascores(i,i)*100;
end
% normalize
for i=1:numel(seqs)
 for j=i+1:numel(seqs) % skip diag first
  ascores(i,j)=ascores(i,j)/sqrt(ascores(i,i)*ascores(j,j))*100;
 end
 ascores(i,i)=100; % diag
end
%
%pcolor(1-ascores) ; shading faceted ; colormap hot ;
mycolor([1:size(ascores,1)],[1:size(ascores,2)], 1-ascores) ; shading flat ; colormap hot ;
set(gca, 'xticklabel', Virus, 'yticklabel', Virus, 'tickdir','out', 'ticklength',[0,0],'xtick',[1:10])
% write in similarity numbers :
for i=1:numel(seqs)
 for j=i:numel(seqs)
  sim=sprintf('%2.0f%',ascores(i,j));
  text(j-0.33,i,['\bf',sim], 'color','w', 'fontsize',12)
 end
end
set(gcf,'paperpositionmode','auto')
print(gcf,'-depsc2',['msa',alname,'.eps']);
print(gcf,'-dpng',['msa',alname,'.png']);

end