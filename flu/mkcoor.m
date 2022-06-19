if ~exist('rootdir')
 rootdir='.';
end
%
addpath([rootdir,'/../']); % make sure aln2coor accessible
% convert aa letters in MSA to numerical coordinates
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
if (~exist('tag'))
 tag='headstem';
end

% load alignment matrix:
msaname=['msa',tag,'.mat'];
load(msaname);
msamat=char( {msa.Sequence}' ) ;
% exclude missing residues :
pmissing = 0.9; % fraction of strains in which the residue in the alignment does not exist ( '-' ) ;
pmiss=mean((msamat=='-'),1);
imiss=find(pmiss>pmissing);
msamat(:,imiss)=''; % delete those residue positions from alignment matrix
%showalignment(msamat);

% convert to numerical coordinates
coor=aln2coor(msamat,qgrantham);

return
