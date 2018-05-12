function varargout = rsa_tutorial(what, varargin)
% Script for performing representational similarity analysis on single
% finger emg data


baseDir     = '/Users/naveed/Dropbox (Personal)/Presentations/2017/rsa_tutorial';

emgChan     = {'APB','FPB','2Lum','3Lum','4Lum','ADM',...
               '1DO','2DO','3DO','4DO',...
               'WEr','WEu','WFr','WFu'};
           
switch(what)
    case 'visualize_data'
        D = load(fullfile(baseDir,'sf_emg.mat'));
        
        % visualize data across all channels
        plt.trace(1:14,D.emg,'split',D.targetForce);
        plt.set('xtick',1:14);
        
    case 'estimate_distances'
        D = load(fullfile(baseDir,'sf_emg.mat'));
        
        T=[];
        S=[];
        for s=unique(D.SN)'                 % loop over all subj
            for f=unique(D.targetForce)'    % loop over all forces
                
                Di = getrow(D,D.SN==s & D.targetForce==f);
                
                % setup variables of interest
                Y           = Di.emg;
                design      = indicatorMatrix('identity',Di.digit);
                partition   = Di.BN;
                
                % get estimate of averaged within class scatter matrix
                [~,~,Sw]=distance_mahalanobis(Y',Di.digit);
                pwY     = Y*Sw^(-1/2);
                
                % get contrast matrix
                C       = indicatorMatrix('allpairs',[1:5]);
                
                % estimate crossvalidated, mahalanobis distance
                d2      = distance_ldc(pwY,design,C,partition);
                
                Ti.SN           = s;
                Ti.targetForce  = f;
                Ti.d2           = d2;
                T               = addstruct(T,Ti);
            end;
        end;
        varargout={T};
    
    case 'visualize_distances'
        T = rsa_tutorial('estimate_distances');
        
        plt.subplot(1,3,1);
        plt.line(T.targetForce,mean(T.d2,2));
        plt.subplot(1,3,[2 3]);
        plt.trace(1:10,T.d2,'split',T.targetForce);
        
    case 'visualize_rdm'
        T = rsa_tutorial('estimate_distances');
        
        rdm = pivottablerow(T.targetForce,T.d2,'mean(x,1)');
        
        for i=1:3
            r = squareform(rdm(i,:));
            plt.subplot(1,3,i);
            plt.image(r);
            axis square;
        end;
        
    case 'visualize_mds'
        T = rsa_tutorial('estimate_distances');
        
        pltno = 1;
        for tf=unique(T.targetForce)'
            Ti      = getrow(T,T.targetForce == tf);
            loading = rsa_tutorial('make_mds',Ti.d2);
            
            plt.subplot(1,3,pltno);
            plt.xy(loading.vec1,loading.vec2,loading.digit,'split',loading.digit,'leglocation','northeast');
            pltno = pltno + 1;
        end;        
        
    case 'make_mds'
        mdsdist     = varargin{1};
        E=[];
                
        % Determine the template from group distances
        [y,eigvals] = cmdscale(mean(mdsdist));
        Aa          = [y zeros(size(y,1),4-size(y,2))];
        
        for n=1:size(mdsdist,1)
            [y,eigvals] = cmdscale(mdsdist(n,:));
            Y=[y zeros(size(y,1),4-size(y,2))];
            
            [d(n),Bp,transform]=procrustes(Aa,Y,'Scaling',false);
            A(:,:,n)=Bp;
            
            D.vec1=A(:,1,n);
            D.vec2=A(:,2,n);
            D.vec3=A(:,3,n);
            D.digit=[1:5]';
            E=addstruct(E,D);
        end;
        
        varargout = {E};

    case 'modelcomparison'
        corrtype = 'spearman';
        T = rsa_tutorial('estimate_distances');
        
        % get emperical model
        E = getrow(T,T.targetForce == 0.75);
        
        % get model H0 (muscle activation for lowest force level)
        h0 = pivottablerow(T.targetForce,T.d2,'mean(x,1)','subset',T.targetForce==0.25);
    
        % get model H1 (muscle activation for middle force level)
        h1 = pivottablerow(T.targetForce,T.d2,'mean(x,1)','subset',T.targetForce==0.5);
    
        D = [];
        for i=unique(E.SN)'
            r_h0 = corr(h0',E.d2(i,:)','type',corrtype);
            r_h1 = corr(h1',E.d2(i,:)','type',corrtype);
            
            Di.SN   = i;
            Di.r_h0 = r_h0;
            Di.r_h1 = r_h1;
            D       = addstruct(D,Di);
        end;
        
        % calculate statistics between the two different models
        ttest(fisherz(D.r_h0),fisherz(D.r_h1),2,'paired');
        
        % calculate noise ceilings and plot data
        [~,c] = crossval_correlation(E.d2,'type',corrtype);
        
        plt.bar([],[D.r_h0, D.r_h1]);
        plt.drawline(c','dir','horz','style',style.custom('black'));
        
    
end;