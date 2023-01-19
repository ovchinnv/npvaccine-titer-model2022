function[]=mycolor(x,y,matrix);
%plots contours with data in the _CENTER_ of the cell
[ny,nx]=size(matrix);
my=length(y);
mx=length(x);
if (my~=ny|mx~=nx); 'plot error: matrix dimensions must agree'
return
end
%matrixn=matrix;
for j=1:ny;matrix(j,nx+1)=matrix(j,nx);end;
for i=1:nx;matrix(ny+1,i)=matrix(ny,i);end;
for i=2:mx;xn(i)=0.5*(x(i)+x(i-1));end
for j=2:my;yn(j)=0.5*(y(j)+y(j-1));end

%matrix

xn(1)=1.5*x(1)-0.5*x(2);
yn(1)=1.5*y(1)-0.5*x(2);
%mx
%my
xn(mx+1)=1.5*x(mx)-0.5*x(mx-1);
yn(my+1)=1.5*y(my)-0.5*y(my-1);
%yn


figure(gcf);
pcolor(xn,yn,matrix);
return;
