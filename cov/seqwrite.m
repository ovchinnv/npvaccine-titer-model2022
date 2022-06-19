%write sequence to a file; would use this to write an aa/nt sequence
function[]=seqwrite(fid,seqs,fmt)
if nargin<3
  seqfmt='fasta';
else
  seqfmt=fmt;
end

if (strcmp(upper(seqfmt),'FASTA')|strcmp(upper(seqfmt),'FST'))
 for i=1:numel(seqs)
  header=seqs(i).Header;
  if(header(1)=='>')
   fprintf(fid,'%s\n',strtrim(header));
  else
   fprintf(fid,'>%s\n',strtrim(header));
  end
  fprintf(fid,'%s\n',strtrim(seqs(i).Sequence));
 end
elseif ( strcmp(upper(seqfmt),'SELEX') | strcmp(upper(seqfmt),'SLX') )
 for i=1:numel(seqs)
  header=seqs(i).Header;
  if(header(1)=='>')
   fprintf(fid,'%-30s',strtrim(header(2:min(30,length(header))))); % remove '>" assuming it got here in error
  else
   fprintf(fid,'30%s',strtrim(header(1:min(29,length(header)))));
  end
  fprintf(fid,'%s\n',strtrim(seqs(i).Sequence));
 end
else
 error (['format ',seqfmt,' not recognized']);
end
