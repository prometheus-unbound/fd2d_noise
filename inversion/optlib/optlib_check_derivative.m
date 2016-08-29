function [dcheck, dcheck_struct] = optlib_check_derivative(m,dm,hpmin,hpmax,step,usr_par)
    
% check if the calculation of a derivative is working correctly
%

[j, g] = eval_objective_and_gradient(m, optlib_generate_random_string(8), usr_par);

nx = usr_par.config.nx; nz = usr_par.config.nz; 
[Lx,Lz] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
cm = cbrewer('div','RdBu',120,'PCHIP');
angle = [0 90];

tmp2 = reshape(g,nx,nz,2);

figure(222)
subplot(1,2,1)
cla
mesh( X, Z, tmp2(:,:,1)' )
hold on
plot3( usr_par.network.array(:,1), usr_par.network.array(:,2),  0*usr_par.network.array(:,2) + max(max(abs( tmp2(:,:,1) ))), 'x' )
colormap(cm)
colorbar
caxis([ -max(max(abs( tmp2(:,:,1) )))  max(max(abs( tmp2(:,:,1) )))]);
axis square
view(angle)
% zlim([-2 2]*1e-5)
drawnow

subplot(1,2,2)
cla
mesh( X, Z, tmp2(:,:,end)' )
hold on
colormap(cm)
colorbar
caxis([ -max(max(abs( tmp2(:,:,end) )))  max(max(abs( tmp2(:,:,end) )))]);
axis square
view(angle)
% zlim([-5 5]*1e-5)
drawnow

j
normg = norm(g)
norm_dm = norm(dm)
djdm = g' * dm


dcheck = zeros(hpmax-hpmin+1,6);
it=0;

for hp=hpmin:step:hpmax
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, (hpmax-hpmin)/step+1 );
    mh = m + 10^hp * dm;
    
    [jh] = eval_objective(mh, optlib_generate_random_string(8), usr_par);
    jh_vec(it) = jh;
    
    djdmh = (jh-j) / (10^hp);
    dcheck(it,:) = [10^hp, djdm, djdmh, abs(djdm - djdmh), abs(djdm - djdmh) / abs(djdm), djdm / djdmh];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).djdm_LHS = djdm;
    dcheck_struct(it).djdmh_RHS = djdmh;
    dcheck_struct(it).absdif_LR  = abs(djdm - djdmh);
    dcheck_struct(it).reldif_LR  = abs(djdm - djdmh) / abs(djdm);
    dcheck_struct(it).ratio_LR   = djdm / djdmh;
    
    
end


h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
title('check derivative - relative error')
jh_vec

end