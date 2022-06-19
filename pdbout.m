% output x,y,z coordinates in PDB format
%
function pdbout(molecule,pdbname,x,y,z,occupancy,temp,inds);
%
 tempModel=molecule;
 if (exist('inds','var'))
  disp(['Writing a subset of atoms : ',num2str(inds(:)')]);
  tempModel.Model.Atom=tempModel.Model.Atom(inds);
  xx=x(inds);
  yy=y(inds);
  zz=z(inds);
  if (exist('occupancy','var') && ~isempty(occupancy) ) ; occu=occupancy(inds) ; end
  if (exist('temp','var') && ~isempty(temp) ) ; tfact=temp(inds) ; end
 else
  xx=x;
  yy=y;
  zz=z;
  if (exist('occupancy','var') && ~isempty(occupancy) ) ; occu=occupancy ; end
  if (exist('temp','var') && ~isempty(temp) ) ; tfact=temp ; end
 end
%
 n=length(xx);
%
 for i=1:n
  tempModel.Model.Atom(i).X=xx(i);
  tempModel.Model.Atom(i).Y=yy(i);
  tempModel.Model.Atom(i).Z=zz(i);
 end
 if (exist('occupancy','var') && ~isempty(occupancy))
  for i=1:n
   tempModel.Model.Atom(i).occupancy=occu(i);
  end
 end
%
 if (exist('temp','var') && ~isempty(temp))
  for i=1:n
   tempModel.Model.Atom(i).tempFactor=tfact(i);
  end
 end
%
 pdbwrite(pdbname,tempModel)
 fid=fopen(pdbname,'A'); fprintf(fid,'END\n'); fclose(fid); % append END to the file because CHARMM needs at least one line after the last coordinate
%
end
