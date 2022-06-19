% find index corresponding to a strain name abbrv 
function oinds=getind(x);
global Virus
global inds
for i=1:numel(x) ;
 oinds(i)=inds(find(ismember(upper(Virus(:)), upper(x(i)))));
end
