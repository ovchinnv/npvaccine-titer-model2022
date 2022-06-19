if exist('OCTAVE_VERSION')
 graphics_toolkit('gnuplot')
end

% plot minimum MSE curve vs exponent
% load data file
name = 'dist2ave_headstem_gr' ; escale=2 ; % in the model 1 code we are optimizing eps/2 (eps as defined in the paper, so times 2 here for vis.)

load([name,'.mat'])

xlabel(opt(2).name);
ylabel(opt(1).name);

check=rr;
%check=cps;
%check=css;
%
styles={'' 'rx','go','bv','m*','ks','-'};
leg={};
figure(2); clf ; hold on ; box on;
for islice=2:6
%
ndim=numel(size(check));
if (ndim==3)
 Diff=opt(3).range(islice)
 leg=[leg {['D=',num2str(Diff)]}];
 slice=check(:,:,islice);
else
 slice=check; % assume 2d
end

[a,inds]=min(slice) ;
%plot(opt(2).range*escale, a, [char(styles(islice)),'.-'])
lw=2; ms=10;
plot(opt(2).range*escale, a, [char(styles(islice)),'-'], 'linewidth', lw, 'markersize', ms)
%xlabel(opt(2).name);
xlabel('\it\beta','fontsize',15);
ylabel('\it MSE') ;
%ylim(c); hold on ;
ylim([0 0.1])
end

l=legend(leg, 'location', 'northeast'); legend boxoff ;
set(l, 'fontsize',14)

text(-0.75*escale, 0.1, 'B', 'fontsize',20)

set(gca, 'fontsize', 18)
set(gca, 'tickdir','out'); box on ;

set(gcf, 'paperpositionmode','auto')
%print(gcf, '-dtiff', [name,'-1D','.tif']) ; % uncompressed
print(gcf, '-depsc2', [name,'-1D','.eps']) ;
