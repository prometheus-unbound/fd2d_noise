function [dcheck, dcheck_struct] = optlib_adjoint_stf_check(u,u0,du,t,src,rec,hpmin,hpmax,step,usr_par)
    
% check if the calculation of a adjoint source time function is working correctly
%

[j, adj_stf(1,:)] = make_adjoint_sources_inversion( u, u0, t, usr_par.veldis, usr_par.measurement, src, rec );

j
norm_adj_stf = norm(adj_stf)
norm_du = norm(du)
djdu = fliplr(adj_stf) * du'


dcheck = zeros(hpmax-hpmin+1,6);
it=0;
jh_vec=zeros((hpmax-hpmin)/step+1, 1);
for hp=hpmin:step:hpmax
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, (hpmax-hpmin)/step+1 );
    duh = u + 10^hp * du;
    [jh] = make_adjoint_sources_inversion( duh, u0, t, usr_par.veldis, usr_par.measurement, src, rec );
    jh_vec(it) = jh;
    
    djduh = (jh-j) / (10^hp);
    dcheck(it,:) = [10^hp, djdu, djduh, abs(djdu - djduh), abs(djdu - djduh) / abs(djdu), djdu / djduh];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).djdm_LHS = djdu;
    dcheck_struct(it).djdmh_RHS = djduh;
    dcheck_struct(it).absdif_LR  = abs(djdu - djduh);
    dcheck_struct(it).reldif_LR  = abs(djdu - djduh) / abs(djdu);
    dcheck_struct(it).ratio_LR   = djdu / djduh;
    
    
end

h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
title('check adjoint stf - relative error')
j
jh_vec

end