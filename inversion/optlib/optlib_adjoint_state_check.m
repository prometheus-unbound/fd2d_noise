function [dcheck, dcheck_struct] = optlib_adjoint_state_check(m,df,hpmin,hpmax,step,usr_par)
    
% check if the calculation of a adjoint state is working correctly
%

[j, ~, ~, adj_state] = eval_objective_and_gradient(m, optlib_generate_random_string(8), usr_par);

j
norm_adjst = sum(sum(sum( adj_state )))
norm_df = sum(sum(sum( df )))
djdf = sum(sum(sum( adj_state .* df )))


dcheck = zeros(hpmax-hpmin+1,6);
it=0;
jh_vec=zeros((hpmax-hpmin)/step+1, 1);
for hp=hpmin:step:hpmax
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, (hpmax-hpmin)/step+1 );
    dfh = 10^hp * df;
    usr_par.debug.df = dfh;
    [jh] = eval_objective(m,optlib_generate_random_string(8), usr_par);
    jh_vec(it) = jh;
    
    djdfh = (jh-j) / (10^hp);
    dcheck(it,:) = [10^hp, djdf, djdfh, abs(djdf - djdfh), abs(djdf - djdfh) / abs(djdf), djdf / djdfh];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).djdm_LHS = djdf;
    dcheck_struct(it).djdmh_RHS = djdfh;
    dcheck_struct(it).absdif_LR  = abs(djdf - djdfh);
    dcheck_struct(it).reldif_LR  = abs(djdf - djdfh) / abs(djdf);
    dcheck_struct(it).ratio_LR   = djdf / djdfh;
    
    
end

h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
title('check adjoint state - relative error')
j
jh_vec

end