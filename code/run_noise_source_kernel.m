function [X,Z,K_s] = run_noise_source_kernel(simulation_mode,i_ref)

%==========================================================================
% run simulation to compute sensitivity kernel for noise power-spectral
% density distribution
%
% input:
%--------
% simulation_mode: 'noise_source_kernel'
% i_ref: number of reference station
%
% output:
%--------
% X, Z: coordinate axes
% K_s: sensitivity kernel
%
%==========================================================================

cm = cbrewer('div','RdBu',100,'PCHIP');


%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------
[Lx,Lz,nx,nz,dt,nt,order,model_type] = input_parameters();
[X,Z,x,z,dx,dz] = define_computational_domain(Lx,Lz,nx,nz);
[mu,rho] = define_material_parameters(nx,nz,model_type); 
output_specs


%- time axis --------------------------------------------------------------    
t=-(nt-1)*dt:dt:(nt-1)*dt;
nt=length(t);


%- read adjoint source locations ------------------------------------------
fid=fopen([adjoint_source_path 'source_locations_' num2str(i_ref)],'r');

adsrc = zeros(1,1);
k = 1;
while (feof(fid)==0)
    adsrc(k,1) = fscanf(fid,'%g',1);
    adsrc(k,2) = fscanf(fid,'%g',1);
    fgetl(fid);
    k = k+1;
end

fclose(fid);


%- read adjoint source time functions -------------------------------------
ns = size(adsrc,1);
stf = zeros(ns,nt);

for n=1:ns
    fid = fopen([adjoint_source_path '/src_' num2str(i_ref) '_' num2str(n)],'r');
    stf(n,1:1:nt) = fscanf(fid,'%g',nt);
    fclose(fid);
end
    

%- compute indices for adjoint source locations ---------------------------    
adsrc_id = zeros(ns,2);
for i=1:ns
    adsrc_id(i,1) = min( find( min(abs(x-adsrc(i,1))) == abs(x-adsrc(i,1)) ) );
    adsrc_id(i,2) = min( find( min(abs(z-adsrc(i,2))) == abs(z-adsrc(i,2)) ) );
end


%- initialise interferometry ----------------------------------------------       
[~, n_sample, w_sample] = input_interferometry();
G_1_vel = zeros(nx,nz,length(t));
           
ifft_coeff = zeros(length(t),n_sample) + 1i*zeros(length(t),n_sample);
for k = 1:n_sample
    ifft_coeff(:,k) = 1/sqrt(2*pi) * exp( 1i*w_sample(k)*t' ) * dw;
end


%- dynamic fields and absorbing boundary field ----------------------------
v = zeros(nx,nz);
sxy = zeros(nx-1,nz);
szy = zeros(nx,nz-1);


%- initialise absorbing boundary taper a la Cerjan ------------------------
[absbound] = init_absbound();


%==========================================================================
% iterate
%==========================================================================

if( strcmp(make_plots,'yes') )
    figure;
    set(gca,'FontSize',20);
end


for n=1:length(t)
    
    %- compute divergence of current stress tensor ------------------------    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    
    
    %- add point sources --------------------------------------------------    
    for i=1:ns
        DS(adsrc_id(i,1),adsrc_id(i,2)) = DS(adsrc_id(i,1),adsrc_id(i,2)) + stf(i,n);
    end
    
    
    %- update velocity field ----------------------------------------------    
    v = v + dt*DS./rho;
    
    
    %- apply absorbing boundary taper -------------------------------------    
    v = v .* absbound;
    
    
    %- compute derivatives of current velocity and update stress tensor ---    
    sxy = sxy + dt*mu(1:nx-1,:) .* dx_v(v,dx,dz,nx,nz,order);
    szy = szy + dt*mu(:,1:nz-1) .* dz_v(v,dx,dz,nx,nz,order);
     
    
    %- save adjoint state -------------------------------------------------
    G_1_vel(:,:,n) = v;
    
    
    %- plot velocity field ------------------------------------------------
    if (strcmp(make_plots,'yes'))
        plot_velocity_field;
    end

end


%==========================================================================
% compute noise source kernels
%==========================================================================

G_1_dis = flip( cumsum(G_1_vel,3), 3) * dt;


%- load Fourier transformed Greens function from forward simulation
load(['../output/interferometry/G_2_' num2str(i_ref) '.mat']);


K_s = zeros(nx,nz);
for n = 1:length(t)
    
    M_tn = zeros(nx,nz) + 1i*zeros(nx,nz);
    
    if( mod(n,5)==0 && t(n)<0.0 )
        for k = 1:n_sample
            M_tn(:,:) = M_tn(:,:) + spectrum(k) * conj(G_2(:,:,k)) * ifft_coeff(n,k);
        end
        
        K_s = K_s + real( M_tn .* G_1_dis(:,:,n) );
    end
    
end


%==========================================================================
% output 
%==========================================================================

%- store the movie if wanted ----------------------------------------------
if strcmp(make_movie,'yes')
    writerObj=VideoWriter(movie_file,'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end

