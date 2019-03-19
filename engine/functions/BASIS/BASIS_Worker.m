%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: local workers at each cpu that calculate likelihood and prior
%           values in log-scale, perform MCMC chains
% created by Stephen Wu, 2015 April 28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Notes
% Note 1: prior annealing not available 

% Last update: 2016 June 14, by Daniel Waelchli
%   2016-06-14: (Daniel)
%       1. reviewed old version, restructured worker (pMALA partly implemented) 
%   2016-04-25: (Daniel)
%       1. included position dependent methods
%   2016-03-18: (Daniel)
%       1. included auto tuning of scaling param eps
%   2016-02-12: (Daniel)
%       1. renamed local_para -> sys_para
%       2. implemented EIG check, adapting EVals x EVecs to extended bounds
%       3. included rescaling, can be turned on/off in sys_para
%   2016-02-08: (Daniel)
%       1. N_accept and N_reject do not include burn-in
%   2015-06-09:
%       1. add handling of burn-in period for each MCMC chain [line 73]
%   2015-06-09:
%       1. check work.N_job = 0 [line 25]
%       2. add flexibility on controling burn-in stages [line 71]
%   2015-10-24:
%       1. support len(N_burn)!=Gen_burn [line 73]

function [out_worker] = BASIS_Worker( ind_worker, work, sys_para, local_info)

% Initialize empty output
out_worker = struct([]);
if isempty(work) || (work.N_job == 0)
    return
else
    out_worker(work.N_job).theta = [];
end

% Message from worker
if sys_para.TF_msg
    fprintf('Worker #%d: My assigned work is %d chains, tot len = %d\n',...
        ind_worker,work.N_job,sum(work.len));
end

if local_info.gen == 1
    % Calculate prior for all sample at once
    tmp = PDF_Pri(work.theta,sys_para.pri);
    % Calculate likelihood and store outputs
    for i = 1:work.N_job
        out_worker(i).theta = work.theta(i,:);
        out_worker(i).Ns = 1;
        % Note: out_lik store as a cell here for consistency
        out_worker(i).out_gradient = zeros(1,sys_para.N_dim);
        [out_worker(i).lik,out_worker(i).out_lik{1}] = ...
            feval(sys_para.lik.name,work.theta(i,:),sys_para.lik.para);
        out_worker(i).out_sig{1} = zeros(sys_para.N_dim);
        out_worker(i).out_gamma = zeros(1,sys_para.N_dim);
        out_worker(i).pri = tmp(i);
        out_worker(i).eps = sys_para.eps;
        % Initialize statistics
        out_worker(i).y = zeros(1,sys_para.N_dim);
        out_worker(i).correction_sigma = 0;
        out_worker(i).theta_inbds = 1;
        out_worker(i).accept = 1;
        out_worker(i).bds_flg = zeros(1,sys_para.N_dim);
        out_worker(i).posdef = 1;
        out_worker(i).burnin_reps = 0;
    end
else
    
    % Take known likelihood values and do MCMC for each job assigned
    for i = 1:work.N_job
        
        % initialize the MCMC chain
        thetao = work.theta(i,:);
        gradiento = calculate_Gradient(work.out_lik{i}, sys_para, local_info);
        [SIG,c_flg,bds_flg,~] = ...
            calculate_Sigma(work.out_lik{i}, work.theta(i,:), sys_para, local_info);
        GAMMA = calculate_Gamma(work.out_lik{i}, SIG, c_flg, sys_para, local_info);
        lnfo_lik = work.lik(i);
        lnfo_pri = work.pri(i);
        o_check = work.out_lik{i}.check;
        outo_lik = work.out_lik{i};
        eps = work.eps(i);
        reps = 0;
        
        % initialize output from worker (delete later if Ns = 0 at the end)
        out_worker(i).Ns = 0;
        out_worker(i).theta = thetao;
        out_worker(i).lik = lnfo_lik;
        out_worker(i).pri = lnfo_pri;
        out_worker(i).out_lik{1} = outo_lik; % store the cell
        out_worker(i).out_gradient = gradiento';
        out_worker(i).out_sig{1} = SIG; % store the cell
        out_worker(i).out_gamma = GAMMA;
        out_worker(i).correction_sigma = c_flg>0;
        out_worker(i).posdef = outo_lik.posdef;
        out_worker(i).eps = eps;
        out_worker(i).burnin_reps = reps;
        out_worker(i).accept = nan; %placeholder
        out_worker(i).bds_flg = bds_flg;
     
        % Auto tune burn-in period
        if sys_para.auto_rescale  && local_info.gen <= sys_para.Gen_rescale
            
            acr = 1;
            
            while reps < sys_para.max_rep && ...
                    (acr<sys_para.tacr_bds(1) || acr>sys_para.tacr_bds(2))
                
                reps = reps + 1;
                n_ac = 0;
                
                for j = 1:sys_para.N_rescale(...
                        min(local_info.gen,length(sys_para.N_burn)))
                    
                    % Propose a new sample based on current theta
                    [TF_inBDs,thetac] = ...
                        Rnd_Prop( thetao, eps, gradiento, SIG, GAMMA, o_check, sys_para);
      
                    if TF_inBDs
                        % Calc. likelihood and prior for proposed sample
                        [lnfc_lik,outc_lik] = ...
                            feval(sys_para.lik.name,thetac,sys_para.lik.para);
                        lnfc_pri = PDF_Pri(thetac,sys_para.pri);
                        c_check = outc_lik.check;
                        
                        % Calc. ratio of acceptance
                        r = accept_ratio(lnfo_lik, lnfo_pri, lnfc_lik, lnfc_pri, eps, thetao, thetac, gradiento, gradientc, SIG, GAMMA, outc_lik, o_check, c_check, sys_para, local_info);
                                                           
                        r = min(1,r);
                        % Accept | Reject
                        state = find(mnrnd(1,[r,1-r]));
                        
                        % Accept or reject proposed sample
                        if state == 1
                            % count acceptance
                            n_ac = n_ac + 1;
                            
                            % update
                            thetao = thetac;
                            gradiento = gradientc;
                            lnfo_lik = lnfc_lik;     
                            lnfo_pri = lnfc_pri;
                            outo_lik = outc_lik;
                            o_check = c_check;
                            
                            [SIG,c_flg,bds_flg,~] = ...
                                calculate_Sigma(outc_lik, thetac, sys_para, local_info);
                            GAMMA = calculate_Gamma(outc_lik, SIG, c_flg, sys_para, local_info);
                            
                            % Update Statistics
                            if ~sys_para.discard_rescaling
                                out_worker(i).theta = thetao;
                                out_worker(i).lik = lnfo_lik;
                                out_worker(i).pri = lnfo_pri;
                                out_worker(i).out_lik{1} = outo_lik; % store the cell
                                out_worker(i).out_gradient = gradiento';
                                out_worker(i).out_sig{1} = SIG; % store the cell
                                out_worker(i).out_gamma = GAMMA;
                                out_worker(i).correction_sigma = c_flg>0;
                                out_worker(i).posdef = outo_lik.posdef;
                            end
                        end
                    end
                end
                % tune scaling parameter
                acr = n_ac/sys_para.N_burn(min(local_info.gen,length(sys_para.N_burn)));
                scale = norminv(sys_para.tacr*0.5)/norminv(acr*0.5);
                eps = eps*min(max(scale,0.1),100);
            end
            % update scaling parameter
            out_worker(i).eps = eps;
            out_worker(i).burnin_reps = reps;
            
        end
        
        % Standard burn-in period
        if sys_para.burnin
            
            if sys_para.discard_rescaling
                % load initial position if adaptive scaling discarded
                thetao = work.theta(i,:);
                [SIG,c_flg,bds_flg,~] = ...
                    calculate_Sigma(work.out_lik{i},work.theta(i,:),sys_para,local_info);
                GAMMA = calculate_Gamma(work.out_lik{i}, SIG, c_flg, sys_para, local_info);
                lnfo_lik = work.lik(i);
                lnfo_pri = work.pri(i);
            end
            
            for j = 1:sys_para.N_burn(...
                    min(local_info.gen,length(sys_para.N_burn)))
                
                [TF_inBDs,thetac] = ...
                       Rnd_Prop( thetao, eps, gradiento, SIG, GAMMA, o_check, sys_para);
                    
                if TF_inBDs
                    % Calc. likelihood and prior for proposed sample
                    [lnfc_lik,outc_lik] = ...
                        feval(sys_para.lik.name,thetac,sys_para.lik.para);
                    c_check = outc_lik.check;
                    lnfc_pri = PDF_Pri(thetac,sys_para.pri);
                    gradientc = calculate_Gradient(outc_lik, sys_para, local_info);
                    
                    % Calc. ratio of acceptance
                    r = accept_ratio(lnfo_lik, lnfo_pri, lnfc_lik, lnfc_pri, eps, thetao, thetac, gradiento, gradientc, SIG, GAMMA, outc_lik, o_check, c_check, sys_para, local_info);
                     
                    r = min(1,r);
                    % Accept | Reject
                    state = find(mnrnd(1,[r,1-r]));
                    
                    % Accept or reject proposed sample
                    if state == 1
                        % update
                        thetao = thetac;
                        gradiento = gradientc;
                        lnfo_lik = lnfc_lik;
                        lnfo_pri = lnfc_pri;
                        outo_lik = outc_lik;
                        o_check = c_check;
                        % update Sigma
                        [SIG,c_flg,bds_flg,~] = ...
                            calculate_Sigma(outc_lik, thetac, sys_para, local_info);
                        GAMMA = calculate_Gamma(outc_lik, SIG, c_flg, sys_para, local_info);

                        % Update Statistics
                        out_worker(i).accept = 1;
                        out_worker(i).theta = thetao;
                        out_worker(i).lik = lnfo_lik;
                        out_worker(i).pri = lnfo_pri;
                        out_worker(i).out_lik{1} = outo_lik; % store the cell
                        out_worker(i).out_gradient = gradiento';
                        out_worker(i).out_sig{1} = SIG; % store the cell
                        out_worker(i).out_gamma = GAMMA;
                        out_worker(i).correction_sigma = c_flg>0;
                        out_worker(i).posdef = outo_lik.posdef;
                    else
                        out_worker(i).accept = 0;
                    end
                else
                    out_worker(i).accept = 0;
                end
            end
            
        end
        
        % Actual MCMC period
        for j = 1:work.len(i)
            
            [TF_inBDs,thetac] = ...
                Rnd_Prop( thetao, eps, gradiento, SIG, GAMMA, o_check, sys_para);
            
            if TF_inBDs
                % Calc. likelihood and prior for proposed sample
                [lnfc_lik,outc_lik] = ...
                    feval(sys_para.lik.name,thetac,sys_para.lik.para);
                lnfc_pri = PDF_Pri(thetac,sys_para.pri);
                c_check = outc_lik.check;
                gradientc = calculate_Gradient(outc_lik, sys_para, local_info);
                
                % Calc. ratio of acceptance
                r = accept_ratio(lnfo_lik, lnfo_pri, lnfc_lik, lnfc_pri, eps, thetao, thetac, gradiento, gradientc, SIG, GAMMA, outc_lik, o_check, c_check, sys_para, local_info);
                
                r = min(1,r);
                state = find(mnrnd(1,[r,1-r]));
                
                % Accept or reject proposed sample
                if state == 1
                    [SIG,c_flg,bds_flg,~] = ...
                        calculate_Sigma(outc_lik, thetac, sys_para, local_info);
                    GAMMA = calculate_Gamma(outc_lik, SIG, c_flg, sys_para, local_info);
                    
                    % add new theta to the record & outc_lik
                    out_worker(i).theta(end+1,:) = thetac;
                    out_worker(i).Ns(end+1) = 1;
                    out_worker(i).lik(end+1) = lnfc_lik;
                    out_worker(i).pri(end+1) = lnfc_pri;
                    out_worker(i).out_lik{end+1} = outc_lik;
                    out_worker(i).out_gradient(end+1,:) = gradiento';
                    out_worker(i).out_sig{end+1} = SIG;
                    out_worker(i).out_gamma(end+1,:) = GAMMA;
                    
                    % update
                    thetao = thetac;
                    gradiento = gradientc;
                    lnfo_lik = lnfc_lik;
                    lnfo_pri = lnfc_pri;
                    o_check = c_check;
                    
                    % update statistics
                    out_worker(i).correction_sigma(end+1) = c_flg>0;
                    out_worker(i).posdef(end+1) = outc_lik.posdef;
                    out_worker(i).eps(end+1) = eps;
                    out_worker(i).burnin_reps(end+1) = reps;
                    out_worker(i).accept(end+1) = 1;
                    out_worker(i).bds_flg(end+1,:) = bds_flg;
                else
                    % update statistics
                    out_worker(i).accept(end) = 0;
                    out_worker(i).bds_flg(end,:) = bds_flg;
                    % increase current theta count by one when rejected
                    out_worker(i).Ns(end) = out_worker(i).Ns(end) + 1;
                end
            else
                % update statistics
                out_worker(i).accept(end) = 0;
                out_worker(i).bds_flg(end,:) = bds_flg;
                % increase current theta count by one when out of bounds
                out_worker(i).Ns(end) = out_worker(i).Ns(end) + 1;
            end
        end
        
        % delete first sample if the count is 0 (from initialization)
        if out_worker(i).Ns(1) == 0
            out_worker(i).theta(1,:) = [];
            out_worker(i).Ns(1) = [];
            out_worker(i).lik(1) = [];
            out_worker(i).pri(1) = [];
            out_worker(i).out_lik(1) = [];
            out_worker(i).accept(1) = [];
            out_worker(i).bds_flg(1,:) = [];
            out_worker(i).posdef(1) = [];
            out_worker(i).eps(1) = [];
            out_worker(i).burnin_reps(1) = [];
            out_worker(i).correction_sigma(1) = [];
            out_worker(i).out_gradient(1,:) = [];
            out_worker(i).out_sig(1) = [];
            out_worker(i).out_gamma(1,:) = [];
        end
    end
end

end

%% Function to calculate Gradient
function [GRADIENT] = calculate_Gradient(out_lik, local_para, local_info)

if out_lik.check == 1
    % ode failed, gradient not available
    GRADIENT = nan(local_para.N_dim,1);
else
    switch local_para.method
        
        case {'smMALA','spmMALA','pMALA'}
            GRADIENT = local_info.p(local_info.gen)*out_lik.gradient(:);

        otherwise
            GRADIENT = nan(local_para.N_dim,1);
            
    end
end

end

%% Function to calculate Sigma
function [SIGMA, do_correction, bds_flg, err_flg] = calculate_Sigma(out_lik,theta,local_para,local_info)

bds_flg = zeros(1,local_para.N_dim);

if out_lik.check ~= 0
    % error in loglike
    % sigma not available 
    SIGMA = local_info.cov_s;
    do_correction = 1;
else
    switch local_para.cov_check.method


        case {'BASIS','NONE'}
            do_correction = 1;
            
        case 'EIG'
            posdef = out_lik.posdef;
            V = out_lik.eig.V;
            D = out_lik.eig.D;
            Dp = D/local_info.p(local_info.gen);

            do_correction = do_eigendirection_check( theta, Dp, V, local_para);
            do_correction = do_correction || ~posdef;

        otherwise
            error('Unknown covariance check method. Choose NONE or EIG');

    end

    if do_correction
        % do the correction in the covariance matrix
        switch local_para.cov_check.method

            case {'BASIS','NONE'}
                SIGMA = local_info.cov_s ;

            case 'EIG'
                Dp_new = Dp;
                if ~posdef
                    E = sort(eig(local_info.cov_s));
                    Dp_new(Dp_new<=0) = E(1:sum(Dp_new<=0));
                    Dp_new(Dp_new<0) = 0;
                end
                if do_eigendirection_check( theta, Dp_new, V, local_para)
                    [Dp_new,bds_flg] = adapt_evals_to_bounds(theta, Dp_new, V, local_para);
                end
                do_eigendirection_sanity_check( theta, Dp_new, V, local_para);
                SIGMA =  V * diag(Dp_new) * V';     
        end
    else
        % no corrections
        % accept proposal and scale it
        SIGMA = out_lik.inv_G/local_info.p(local_info.gen);
    end
end

% force symmetry
SIGMA = (SIGMA + SIGMA.') / 2;

% sanity check
[~,err_flg] = cholcov(SIGMA);
if err_flg ~= 0
    % correction failed (SIGMA not symmetric posdef)
    % use sample covariance matrix as SIGMA
    SIGMA = local_info.cov_s;
    warning('cov correction failed, setting SIGMA = local_info.cov_s');
end

end

%% Function to calculate gamma in position dependent manifold MALA (pMALA)
function [GAMMA, err_flg] = calculate_Gamma(out_lik, SIG, do_correction, local_para, local_info)

Nth = local_para.N_dim;
err_flg = 0;
if (out_lik.check > 0)
    GAMMA = nan(1,Nth);
    err_flg = 1;
else
    
    switch local_para.method
        
        case 'pMALA'
            Dp = (out_lik.eig.D/local_info.p(local_info.gen))';
            [Vc,Dc] = eig(SIG);
            [I,p] = compareCol(Vc,out_lik.eig.V);
            if p == 1
                Dc = diag(Dc)'*I;
                DC = repmat(Dc,local_para.N_dim,1);
                DP = repmat(Dp,local_para.N_dim,1);
                J = DC-DC'./(DP-DP');
                J(logical(eye(local_para.N_dim))) = Dc./Dp;
                GAMMA = zeros(1,Nth);
                %optimize?
                for i = 1:Nth
                    sum_term = zeros(Nth,1);
                    for j = 1:Nth
                        G_hat = out_lik.eig.V*(J.*(out_lik.eig.V'*out_lik.dGK(:,(j-1)*Nth+1:j*Nth)*out_lik.eig.V))*out_lik.eig.V';
                        sum_term = sum_term+G_hat*SIG(:,j);
                    end
                    GAMMA(i) = -0.5*SIG(:,i)'*sum_term;
                end
                GAMMA = local_info.p(local_info.gen)*GAMMA; % multiply by p because out_lik.dGK not scaled yet

            else
                warning('column comparison failed, eigenvectors not found (gamma = zero)');
                GAMMA = zeros(1,Nth);
            end
            
        case 'spmMALA'
            GAMMA = zeros(1,Nth);
            if do_correction == 0
                for j = 1:Nth
                    %optimize?
                    sum_term = SIG*out_lik.dGK(:,(j-1)*Nth+1:j*Nth)*SIG;
                    GAMMA = GAMMA + sum_term(:,j)';
                end
                GAMMA = local_info.p(local_info.gen)*GAMMA;
            else
                GAMMA = zeros(1,Nth);
            end
            
        otherwise
            GAMMA = nan(1,Nth);
            
    end
        
end

if any(isnan(GAMMA(:))) || any(isinf(GAMMA(:)))
    GAMMA = nan(1,Nth);
    err_flg = 1;
end

%disp(GAMMA);

end

%% Function to calculate prior pdf values in log scale
% Note 1: only pass sys_para.pri into para
% Note 2: theta is matrix with size - #samples x #dimension
% Note 3: ln_f(k) = -inf if prior = 0 for the sample (k)
function [ln_f] = PDF_Pri(theta,para)
if para.TF_custom
    ln_f = feval(para.name_custom_pdf,theta,para.para_custom_pdf);
else
    ln_f = zeros(size(theta,1),1);
    for i = 1:size(theta,2)
        switch length(para.para{i})
            case 1
                ln_f = ln_f + ...
                    log(pdf(para.name{i},theta(:,i),para.para{i}));
            case 2
                ln_f = ln_f + log(pdf(para.name{i},theta(:,i),...
                    para.para{i}(1),para.para{i}(2)));
            case 3
                ln_f = ln_f + log(pdf(para.name{i},theta(:,i),...
                    para.para{i}(1),para.para{i}(2),para.para{i}(3)));
        end
    end
end
end

%% Function to propose next sample based on input theta
function [TF_inBDs,thetac] = Rnd_Prop( theta, eps, gradient, SIG, GAMMA, o_check, sys_para)
% propose sample with BASIS
bds = sys_para.hard_bds;
theta = theta(:);
GAMMA = GAMMA(:);

switch sys_para.method
    case 'BASIS'
        thetac = mvnrnd( theta, eps*SIG );
        
    case {'smMALA','SN'}
        if o_check == 0 || o_check == 2 % smMALA or gBASIS (SIG = cov_s)
            thetac = mvnrnd( theta + 0.5*eps*SIG*gradient, eps*SIG );
        elseif o_check == 1 % BASIS
            thetac = mvnrnd( theta, eps*SIG );
        else
            error('invalid state: o_check');
        end

    case 'pMALA'
        thetac = mvnrnd( theta + 0.5*eps*SIG*gradient + eps*GAMMA, eps*SIG );
        
    otherwise
        error('Unknown method. Choose BASIS or smMALA or pMALA or SN')
        
end

% check hard bounds (assume theta is row vector)
% bds(1,:) = row vector of min; bds(2,:) = row vector of max
if any(thetac < bds(1,:)) || any(thetac > bds(2,:))
    TF_inBDs = false;
else
    TF_inBDs = true;
end
end


%% Functions to calculate acceptance ratio
% Note 1: not annealing prior
function r = accept_ratio(lnfo_lik,  lnfo_pri, lnfc_lik, lnfc_pri, eps, thetao, thetac, gradiento, gradientc, SIG, GAMMA, outc_lik, o_check, c_check, sys_para, local_info)


switch sys_para.method
    case 'BASIS'
        r = exp(local_info.p(local_info.gen)*...
            (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri));
        
    case {'smMALA','SN'}
        if o_check == 0 || o_check == 2
            if c_check == 0 || c_check == 2 % smMALA or gBASIS (SIG = cov_s)
                
                tmpc = thetac(:) - thetao(:) - 0.5*eps*SIG*gradiento(:);
                tmpo = thetao(:) - thetac(:) - 0.5*eps*SIG*gradientc(:);
                
                qc = -(tmpc'*(eps*SIG\tmpc))/2;
                qo = -(tmpo'*(eps*SIG\tmpo))/2;
                
                r = exp(local_info.p(local_info.gen)*...
                    (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri) + ...
                    qo - qc);
                
            elseif c_check == 1 % reject
                r = 0.0;
            else
                error('invalid state: c_check');
            end
            
        elseif o_check == 1 % BASIS
            
            r = exp(local_info.p(local_info.gen)*...
                (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri));
        else
            error('invalid state: o_check');
        end
        
    case 'pMALA'
        
        if o_check == 0
            if c_check == 0 || c_check == 3
                [SIGc,c_flg,~] = ...
                    calculate_Sigma(outc_lik, thetac, sys_para, local_info);
                
                if c_flg == 0 % pMALA
                    
                    tmpc = thetac(:) - thetao(:) - 0.5*eps*SIG*gradiento(:) - eps*GAMMA(:);
                    tmpo = thetao(:) - thetac(:) - 0.5*eps*SIGc*gradientc(:) - eps*GAMMA(:);
                    
                    qc = -(tmpc'*(eps*SIG\tmpc))/2;
                    qo = -(tmpo'*(eps*SIGc\tmpo))/2;
                    
                    r = exp(local_info.p(local_info.gen)*...
                        (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri) + ...
                        qo - qc);
                    
                else % gBASIS
                    
                    tmpc = thetac(:) - thetao(:) - 0.5*eps*SIG*gradiento(:);
                    tmpo = thetao(:) - thetac(:) - 0.5*eps*SIG*gradientc(:);
                    
                    qc = -(tmpc'*(eps*SIG\tmpc))/2;
                    qo = -(tmpo'*(eps*SIG\tmpo))/2;
                    
                    r = exp(local_info.p(local_info.gen)*...
                        (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri) + ...
                        qo - qc);
                    
                end
                
            elseif c_check == 1 % reject
                
                r = 0.0;
                
            elseif c_chek == 2 % smMALA or gBASIS (SIG = cov_s)
                
                tmpc = thetac(:) - thetao(:) - 0.5*eps*SIG*gradiento(:);
                tmpo = thetao(:) - thetac(:) - 0.5*eps*SIG*gradientc(:);
                
                qc = -(tmpc'*(eps*SIG\tmpc))/2;
                qo = -(tmpo'*(eps*SIG\tmpo))/2;
                
                r = exp(local_info.p(local_info.gen)*...
                    (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri) + ...
                    qo - qc);
                
            else
                error('invalid state: c_check');
            end
            
        elseif o_check == 1 % BASIS
            
            r = exp(local_info.p(local_info.gen)*...
                (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri));
            
        elseif o_chek == 2
            if c_check == 0 || c_check == 2 || c_check == 3 % gBASIS (SIG = cov_s)
                
                tmpc = thetac(:) - thetao(:) - 0.5*eps*SIG*gradiento(:);
                tmpo = thetao(:) - thetac(:) - 0.5*eps*SIG*gradientc(:);
                
                qc = -(tmpc'*(eps*SIG\tmpc))/2;
                qo = -(tmpo'*(eps*SIG\tmpo))/2;
                
                r = exp(local_info.p(local_info.gen)*...
                    (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri) + ...
                    qo - qc);
                
            elseif c_check == 1 % reject
                
                r = 0.0;
                
            else
                error('invalid state: c_check');
            end
            
        elseif o_check == 3
            if c_check == 0 || c_check == 2 || c_check == 3 % smMALA or gBASIS (SIG = cov_s)
                
                tmpc = thetac(:) - thetao(:) - 0.5*eps*SIG*gradiento(:);
                tmpo = thetao(:) - thetac(:) - 0.5*eps*SIG*gradientc(:);
                
                qc = -(tmpc'*(eps*SIG\tmpc))/2;
                qo = -(tmpo'*(eps*SIG\tmpo))/2;
                
                r = exp(local_info.p(local_info.gen)*...
                    (lnfc_lik - lnfo_lik) + (lnfc_pri - lnfo_pri) + ...
                    qo - qc);
                
            elseif c_check == 1 % reject
                
                r = 0.0;
                
            else
                error('invalid state: c_check');
            end
            
        else
            error('invalid state: o_check');
        end
        
    otherwise
        error('Unknown method. Choose BASIS or smMALA or pMALA or SN')
        
end

end


%% Helper Functions
function do_correction = do_eigendirection_check( x, Dp, V, local_para)

x = x(:);
df = length(x);
chi = chi2inv(local_para.conf,df);

do_correction = 0;

for l = 1:local_para.N_dim
    
    sc = sqrt(Dp(l)/chi);
    
    pp = x + sc*V(:,l);
    if ~in_the_box(pp,local_para.hard_bds)
        do_correction = 1;
        break;
    end
    
    pm = x - sc*V(:,l);
    if ~in_the_box(pm,local_para.hard_bds)
        do_correction = 1;
        break;
    end
    
end


end

function err_flg = do_eigendirection_sanity_check( x, Dp, V, local_para)

eps = 1e-3;
x = x(:);
df = length(x);
chi = chi2inv(local_para.conf,df);

bds = local_para.hard_bds;
if local_para.use_extended_bds
    width = bds(2,:)-bds(1,:);
    bds(2,:) = bds(2,:)+width.*local_para.extend_bds;
    bds(1,:) = bds(1,:)-width.*local_para.extend_bds;
end

err_flg = 0;

for l = 1:local_para.N_dim
    
    sc = sqrt(Dp(l)/chi);
    
    pp = x + sc*V(:,l);
    if any(pp'<(bds(1,:)-eps) | pp'>(bds(2,:)+eps));
        warning('adaption failed');
        err_flg = 1;
        
        break;
    end
    
    pm = x - sc*V(:,l);
    if any(pm'<(bds(1,:)-eps) | pm'>(bds(2,:)+eps) )
        warning('adaption failed');
        err_flg = 1;
        break;
    end
    
end


end

function [Dp,bds_flag] = adapt_evals_to_bounds( x, Dp, V, local_para)

x = x(:);
df = length(x);
chi2 = chi2inv(local_para.conf,df);

bds_flag = zeros(1,local_para.N_dim);

bds = local_para.hard_bds;
if local_para.use_extended_bds
    width = bds(2,:)-bds(1,:);
    bds(2,:) = bds(2,:)+width.*local_para.extend_bds;
    bds(1,:) = bds(1,:)-width.*local_para.extend_bds;
end

for l = 1:local_para.N_dim
    
    sc = sqrt(Dp(l)*chi2);
    
    pp = x + sc*V(:,l);
    flag = pp' <= bds(1,:);
    bds_flag = bds_flag | flag;
    if any(flag)>0 %exceeding lower bds in pos dir
        c = abs(1./V(flag,l) .* (bds(1,flag)'-x(flag)));
        sc = min(c(:));
    end
    
    pp = x + sc*V(:,l);
    flag = pp' >= bds(2,:);
    bds_flag = bds_flag | flag;
    if any(flag)>0 %exceeding lower bds in neg dir
        c = abs(1./V(flag,l) .* (bds(2,flag)'-x(flag)));
        sc = min(c(:));
    end
    
    
    
    pm = x - sc*V(:,l);
    flag = pm' <= bds(1,:);
    bds_flag = bds_flag | flag;
    if any(flag)>0 %exceeding lower bds in neg dir
        c = abs(1./V(flag,l) .* (bds(1,flag)'-x(flag)));
        sc = min(c(:));
    end
    
    pm = x - sc*V(:,l);
    flag = pm' >= bds(2,:);
    bds_flag = bds_flag | flag;
    if any(flag)>0 %exceeding upper bds in neg dir
        c = abs(1./V(flag,l) .* (bds(2,flag)'-x(flag)));
        sc = min(c(:));
    end
    
    Dp(l) = sc^2/chi2;
end

end


function flag = in_the_box( x, Box )

x = x(:)';

flag = ~any( x<Box(1,:) | x>Box(2,:) );

end


function [I,p] = compareCol(V,W)
I = zeros(size(V));
p = 0;
for i = 1:size(V,2)
    for j = 1:size(V,2)
        if isalmost(V(:,i),W(:,j)) || isalmost(V(:,i),-W(:,j))
            I(i,j) = 1;
            break;
        end
    end
end

if all(sum(I,1) == ones(1,size(V,1))) && all(sum(I,2) == ones(size(V,2),1))
    p = 1;
end

end

function b = isalmost(x,y)
tol = 1e-6;
b = norm(x-y,2) < tol;
end


