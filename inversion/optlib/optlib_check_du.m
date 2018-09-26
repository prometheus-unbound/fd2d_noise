function [dcheck, dcheck_struct] = optlib_check_du( m, dm, hpmin, hpmax, step, usr_par )

[~, du, u] = eval_hessian_vector_product( m, dm, optlib_generate_random_string(8), usr_par);

norm_du = norm(du)
norm_dm = norm(dm)

dcheck = zeros(hpmax-hpmin+1,6);
it=0;

for hp=hpmin:step:hpmax
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, (hpmax-hpmin)/step+1 );
    mh = m + 10^hp * dm;
    
    [~, ~, uh] = eval_hessian_vector_product( mh, 0, optlib_generate_random_string(8), usr_par);
    uh_vec(it,:) = (uh-u) / 10^hp;
    
    figure(1)
    clf
    hold on
    plot( du, 'k' )
    plot( (uh-u) / 10^hp , 'r' )
    drawnow
    
    norm_duh = norm( (uh-u) / (10^hp) );
    dcheck(it,:) = [10^hp, norm_du, norm_duh, abs(norm_du - norm_duh), abs(norm_du - norm_duh) / abs(norm_du), norm_du / norm_duh];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).dc_LHS = norm_du;
    dcheck_struct(it).dch_RHS = norm_duh;
    dcheck_struct(it).absdif_LR  = abs(norm_du - norm_duh);
    dcheck_struct(it).reldif_LR  = abs(norm_du - norm_duh) / abs(norm_du);
    dcheck_struct(it).ratio_LR   = norm_du / norm_duh;
    
end


h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
title('check derivative - relative error')
uh_vec;

end