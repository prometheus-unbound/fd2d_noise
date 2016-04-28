function [dcheck, dcheck_struct] = optlib_check_hessian( m, dm1, dm2, hpmin, hpmax, step, usr_par )
    
nx = usr_par.config.nx; nz = usr_par.config.nz; 

[Hdm1] = eval_hessian_vector_product( m, dm1, optlib_generate_random_string(8), usr_par);
[~, g] = eval_objective_and_gradient( m, optlib_generate_random_string(8), usr_par );

Hdm1dm2 = Hdm1' * dm2;
normg = norm(g)

angle = [0 90];

fig666 = figure(666);
set(fig666,'units','normalized','position',[.1 .5 0.5 0.4])
clf
[Lx,Lz] = input_parameters();
[X,Z] = define_computational_domain(Lx,Lz,nx,nz);
uiui = reshape( Hdm1, nx, nz, 2 );

subplot(1,2,1)
mesh( X, Z, uiui(:,:,1)' )
hold on
plot3( usr_par.network.array(:,1), usr_par.network.array(:,2),  0*usr_par.network.array(:,2) + max(max(abs( uiui(:,:,1) ))), 'x' )
cm = cbrewer('div','RdBu',120,'PCHIP');
colormap(cm);
caxis([ -max(max(abs( uiui(:,:,1) )))  max(max(abs( uiui(:,:,1) )))]);
axis square
colorbar
view(angle)
% zlim([-2 2]*1e-5)
drawnow

subplot(1,2,2)
mesh( X, Z, uiui(:,:,end)' )
hold on
cm = cbrewer('div','RdBu',120,'PCHIP');
colormap(cm);
caxis([ -max(max(abs( uiui(:,:,end) )))  max(max(abs( uiui(:,:,end) )))]);
axis square
colorbar
view(angle)
% zlim([-5 5]*1e-5)
drawnow


norm_dm1 = norm(dm1)
norm_dm2 = norm(dm2)

dcheck = zeros(hpmax-hpmin+1,6);
it=0;


for hp=hpmin:step:hpmax
    
    it=it+1;
    fprintf('current power of 10: %5d \niteration:           %2d/%2d\n', hp, it, (hpmax-hpmin)/step+1 );

    mh = m + 10^hp * dm1;
    
    [~, gh] = eval_objective_and_gradient(mh, optlib_generate_random_string(8), usr_par);    
    
    gh_vec(it,:) = gh; 
    tmp = ( gh - g ) / (10^hp);
    dgdhdm = tmp' * dm2;
    
    
    fig999 = figure(999);
    set(fig999,'units','normalized','position',[.1 .1 0.5 0.4])
    tmp2 = reshape(tmp,nx,nz,[]);
    
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
%     zlim([-2 2]*1e-5)
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
%     zlim([-5 5]*1e-5)
    drawnow
    
    
    dcheck(it,:) = [10^hp, Hdm1dm2, dgdhdm, abs(Hdm1dm2 - dgdhdm), abs(Hdm1dm2 - dgdhdm) / abs(Hdm1dm2), Hdm1dm2 / dgdhdm];
    dcheck_struct(it).powerof10 = 10^hp;
    dcheck_struct(it).Hdm1dm2_LHS = Hdm1dm2;
    dcheck_struct(it).dgdhdm_RHS = dgdhdm;
    dcheck_struct(it).absdif_LR  = abs(Hdm1dm2 - dgdhdm);
    dcheck_struct(it).reldif_LR  = abs(Hdm1dm2 - dgdhdm) / abs(Hdm1dm2);
    dcheck_struct(it).ratio_LR   = Hdm1dm2 / dgdhdm;

    
end


h1 = figure;
loglog(dcheck(:,1),dcheck(:,5))
xlabel('h','Interpreter','Latex')
ylabel('error','Interpreter','Latex')
title('Hessian vector product - relative error')

end
