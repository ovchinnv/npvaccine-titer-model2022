rootdir='../';
addpath(rootdir,'-END');

% whether to initialize (otherwise, continue from prev. state) :
if (~exist('qinit'))
 qinit=1;
end
%
qwrand=0; % whether residue weights are initialized randomly
wamp=0.1; % scaling amplitude for initial weights
qmode='scan'; % scan
nrep = 1 ; % number of repetitions to compute error (i.e., if stochastic)
%
if (~exist('model'))
 model='dist2ave' ; % name of matlab script file that performs the fit
end
%
if (qinit)
 clear rave;
 clear opt;
 maxoptiter=100;
 i=0; % opt var index
% create optimization structure
% specify structure
 i=i+1;
 opt(i).name='xint'; % must match parameter name in error evaluator
 opt(i).ival=1 ; % initial value
 opt(i).dh=0.005 ;% FD step in computing derivative
 opt(i).step=0.02 ;% evolution step in gradient descent
 opt(i).maxchange=0.1 ;% maximum change between iterations
 opt(i).minval=0.1 ; % minimum allowed value (-inf to disable for mini)
 opt(i).maxval=2 ; % maximum allowed value (inf to disable for mini)
%
 i=i+1;
 opt(i).name='xp';
 opt(i).ival=1 ;
 opt(i).dh=0.005 ;
 opt(i).step=0.02 ;
 opt(i).maxchange=0.1 ;
 opt(i).minval=0.1 ;
 opt(i).maxval=4 ;
%
 i=i+1;
 opt(i).name='Diff';
 opt(i).ival='0.3';
 opt(i).maxchange=0.1 ;
 opt(i).minval=0. ;
 opt(i).maxval=0.5 ;
%% generic inits below
 for j=1:i
  opt(j).val=opt(j).ival ;
 end
 npar=length(opt); % number of parameters to optimize
 qmini=strcmp(upper(qmode),'MINI') ;
 qscan=strcmp(upper(qmode),'SCAN') ;
 if (qmini)
  xpar=zeros(0,npar);
  rr=zeros(1,0);
  cps=zeros(1,0);
  css=zeros(1,0);
 elseif (qscan)
  shape='';
  for k=1:npar
   opt(k).range=[opt(k).minval : opt(k).maxchange : opt(k).maxval];
   shape=[shape, num2str( floor ( ( opt(k).maxval - opt(k).minval ) / opt(k).maxchange ) + 1 )];
   if (k<npar) ; shape=[shape,',']; end
  end
  rr=eval(['zeros(',shape,')']);
  cps=eval(['zeros(',shape,')']);
  css=eval(['zeros(',shape,')']);
 else
  error(['Unknown mode : ',qmode]);
 end
%
 qinit=0;
end
%
fmt='%.16f'; % format for num2str

if (qmini) % minimization :

for optiter=1:maxoptiter
% set initial values to make sure that all parems are defined:
 for ipar=1:npar
  eval( [opt(ipar).name,'=',sprintf(fmt,opt(ipar).val),';']);
 end
% compute error derivative by FD
 for ipar=1:npar
  dd=[-opt(ipar).dh +opt(ipar).dh];
%
  for k=1:2 % compute diff
   eval( [opt(ipar).name,'=',sprintf(fmt, opt(ipar).val + dd(k) ),';']);
% custom:
   for irepe=1:nrepe
    eval(model) % run model
    rerrs(k, irepe)=e2a ; % error
    corrps(k, irepe)=c  ; % pearson
    corrss(k, irepe)=cs ; % spearman
   end
   rerr=mean(rerrs,2);
   corrp=mean(corrps,2);
   corrs=mean(corrss,2);
% custom
  end
  drerr(ipar)= (rerr(2)-rerr(1))/(dd(2)-dd(1))
  rave(ipar) = 0.5*(rerr(1)+rerr(2))
 end % over all params
%
% compute new parameter values
%
 for ipar=1:npar
  dh = min( opt(ipar).maxchange, opt(ipar).step * abs(drerr(ipar))); % evolution size
  opt(ipar).val = opt(ipar).val - dh*sign(drerr(ipar));
  opt(ipar).val = min(opt(ipar).maxval, max( opt(ipar).val, opt(ipar).minval )); % make sure to keep inside a prescribed interval
 end

 xpar=[xpar ; [opt(:).val] ];

 rr=[rr mean(rave)]; % keep track of error
 cps=[cps mean(corrp)]; % avg. pearson
 css=[css mean(corrs)]; % avg. spearman

 fprintf('%s %d\n', ' === iteration',optiter);
 fprintf( ['%s ', repmat('%17.12f', [1,npar]), '\n'], 'New parameters: ',[opt(:).val]);

end % iterations
%
savename=[model,'_',date,tag,'_',enc,'.mat'];

elseif (qscan)
 % build multidimensional loop
 pre='' ;scaninds=''; post='';
 for ipar=1:npar
  loop=sprintf('for i%d=1:%d ; %s=opt(%d).range(i%d), ;', ipar, length(opt(ipar).range), opt(ipar).name, ipar, ipar );
  pre=[pre , loop];
  scanind=sprintf('i%d',ipar);
  scaninds=[scaninds,scanind]; if (ipar<npar) ; scaninds=[scaninds,',']; end
  post=[post,'end;'];
 end
 body=['errs=[];crps=[];crss=[];',scaninds,', ; for irep=1:nrep ; ',model,'; errs=[errs e2a]; crps=[crps c]; crss=[crss cs]; end ;',... % continuation
  'rr(',scaninds,')=mean(errs) ;cps(',scaninds,')=mean(crps) ;css(',scaninds,')=mean(crss) ; '];
 script=[pre,'', body, '', post] ;
% shall we ?
 disp('Running script:');
 disp(script)
 eval(script) ;% run it

 savename=[model,'_',tag,'_',enc,'.mat'];
end

save(savename, '-mat')
