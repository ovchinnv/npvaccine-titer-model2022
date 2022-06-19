% lookup table
global Virus inds
% Alex Cohen's data:
Virus={'SARS-2', 'RaTG13', 'SHC014', 'Rs4081', 'Pang17', 'RmYN02', 'Rf1', 'WIV1', 'SARS', 'Yun11', 'BM4831', 'BtKY72'};
Acc={'MN985325.1', 'QHR63300', 'KC881005', 'KY417143', 'QIA48632', 'EPI_ISL_412977', 'DQ412042', 'KF367457', 'AAP13441.1', 'JX993988', 'NC_014470', 'KY352407'};

%
% in matlab data searches do not include added paths ; therefore, need to cd
wd=cd;
cd('~/cov/cohen21');
%
[~,files]=system('ls ../[A-Z]*.fasta');
files=sort(strread(files,'%s'));
cd(wd)

%loop over all Viruses, and find its index in the coordinate/alignment array

j=0; % accession/Virus index

for acc=Acc
j=j+1;

i=1; 
for file=files(:)'
 fname=char(file);fname=fname(4:end);
 ie=find(fname=='.',1);
 nm=fname(1:ie-1);
 nm2=char(names(i)); % note : names correspond directly to msa(2:end)
 ie=find(nm2=='.',1);
 nm2=nm2(1:ie-1);
 if ~(isempty(strfind(char(acc),nm))&isempty(strfind(char(acc),nm2)))
  inds(j)=i; % virus at position j has its sequence on the msa/coor array at position i
 end
 i=i+1;
end
end

