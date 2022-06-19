function oinds=getind(x);
global Virus
% return index of the strain "x" in the Virus array
% x can be vector
% need Virus to be defined
for i=1:numel(x) ;
 oinds(i)=find(ismember(upper(Virus(:)), upper(x(i))));
end
