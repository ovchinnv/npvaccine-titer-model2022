% plot MSE slice as a function of parameters
%
if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot')
end

% first, load fit file
name = 'dist2ave_headstem_gr' ; escale=2 ; % in the model 1 code we are optimizing eps/2 (eps as defined in the paper, so times 2 here for vis.)
load([name,'.mat'])

close all ; figure(1) ;  hold on ; box on ;

xlabel(opt(2).name);
ylabel(opt(1).name);

ylabel('\it\alpha', 'fontsize',15)
xlabel('\it\beta', 'fontsize',15)

% rr could be a 3d array if we are also scanning diff (i.e. depends on what params optimized)
% note that we can plot rr, epc, css, as fit quality measures ;
check=rr;
%
%
ndim=numel(size(check));
if (ndim==3)
 islice=5; % slice corresponding to a particular Diff constant
 Diff=opt(3).range(islice)
 slice=check(:,:,islice);
else
 slice=check; % assume 2d
end

pcolor(opt(2).range * escale, opt(1).range, slice);shading interp; colorbar ; view(2) ; title('MSE')
% plot as a 3D mesh
%mesh(opt(2).range, opt(1).range, slice); view(3) ; zlabel('MSE');

axis tight;
c=[0 0.25];
c=[0.0 0.1];
c=[0.01 0.11];
caxis(c);
%zlim(c);

% plot minimum curve
[a,inds]=min(slice) ;
% this comes out too jagged
%scatter3(opt(2).range, opt(1).range(inds),a, 'k.')
% smooth :
xx=opt(2).range * escale;
yy=opt(1).range(inds);
np=100 ;
xxs=linspace(xx(1),xx(end),np);
yys=interp1(xx,yy,xxs, 'cubic')
plot(xxs,yys,'w-', 'linewidth', 2)

text(-1,2.1, 'A', 'fontsize',20)

set(gca, 'fontsize', 17)
set(gca, 'tickdir','out'); box on ;

set(gcf, 'paperpositionmode','auto')
%print(gcf, '-depsc2', '-zbuffer', [name,'-D',num2str(Diff),'.eps']) ;%
print(gcf,  '-depsc2', [name,'-D',num2str(Diff),'.eps']) ;%
%print(gcf, '-dpng', [name,'-D',num2str(Diff),'.png']) ;
print(gcf,'-dtiff', [name,'-D',num2str(Diff),'.tif']) ; % good but uncompressed
