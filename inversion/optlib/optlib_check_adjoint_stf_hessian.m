function [dcheck, dcheck_struct] = optlib_check_adjoint_stf_hessian(u_ref,u_0,du1,du2,t,src,rec,hpmin,hpmax,step,usr_par)
    
[~, d_adjstf_du1] = make_adjoint_sources( u_ref, u_0, du1, t, usr_par.veldis, usr_par.measurement, src, rec, '2nd' );
[~, adjstf_ref] = make_adjoint_sources( u_ref, u_0, 0, t, usr_par.veldis, usr_par.measurement, src, rec, '1st' );

norm_dadjstf_du1 = norm(d_adjstf_du1)
norm_du1 = norm(du1)
norm_du2 = norm(du2)
d_adjstf_du1du2 = fliplr(d_adjstf_du1) * du2'


dcheck = zeros(hpmax-hpmin+1,6);
it=0;
for hp=hpmin:step:hpmax
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, (hpmax-hpmin)/step+1 );
    duh = u_ref + 10^hp * du1;
    [~, adjstf_h] = make_adjoint_sources( duh, u_0, 0, t, usr_par.veldis, usr_par.measurement, src, rec, '1st' );
    
    d_adjstf_h_du2 = ( fliplr( adjstf_h - adjstf_ref ) / 10^hp ) * du2';
    adjstf_h_vec(it,:) = adjstf_h;
    
    figure(1)
    clf
    hold on
    plot(d_adjstf_du1,'b*')
    plot((adjstf_h-adjstf_ref)/10^hp,'r+')
    
%     figure(2)
%     plot((adjstf_h-adjstf_ref)/10^hp - d_adjstf_du1)
    
    dcheck(it,:) = [10^hp, d_adjstf_du1du2, d_adjstf_h_du2, abs(d_adjstf_du1du2 - d_adjstf_h_du2), abs(d_adjstf_du1du2 - d_adjstf_h_du2) / abs(d_adjstf_du1du2), d_adjstf_du1du2 / d_adjstf_h_du2];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).d_adjstf_du1du2_LHS = d_adjstf_du1du2;
    dcheck_struct(it).d_adjstf_h_du2_RHS = d_adjstf_h_du2;
    dcheck_struct(it).absdif_LR  = abs(d_adjstf_du1du2 - d_adjstf_h_du2);
    dcheck_struct(it).reldif_LR  = abs(d_adjstf_du1du2 - d_adjstf_h_du2) / abs(d_adjstf_du1du2);
    dcheck_struct(it).ratio_LR   = d_adjstf_du1du2 / d_adjstf_h_du2;
    
    
end

h1 = figure;
set(gca,'FontSize',12);
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
title('adstf hessian check - relative error')

end