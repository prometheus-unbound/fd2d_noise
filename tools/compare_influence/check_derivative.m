
clear all
close all

addpath(genpath('../../'))
[~,~,~,~,dt,nt] = input_parameters();
t = -(nt-1)*dt:dt:(nt-1)*dt;


folder_1 = '~/Desktop/runs/inversion/data/';
load([folder_1 'array_16_ref.mat'])
n_rec = size(array,1)-1;
u_h = load([folder_1 'data_16_ref_uniform_homog_structure.mat']);
u_1 = load([folder_1 'data_16_ref_uniform_structure_1.mat']);
 

% folder_2 = '~/Desktop/runs/inversion_basis_fct/data/';
% load([folder_2 'array_16_ref.mat'])
% n_rec = size(array,1)-1;
% u1_h_0 = load([folder_2 'data_16_ref_0_uniform1_homogeneous.mat');
% u2_h_0 = load([folder_2 data_16_ref_0_uniform2_homogeneous.mat');

% u1b_h_0 = load([folder_2 data_16_ref_0_uniform_1gaussian_homogeneous.mat');
% u2b_h_0 = load([folder_2 data/data_16_ref_0_uniform_2gaussian_homogeneous.mat');


% measurement = 'log_amplitude_ratio';
% measurement = 'amplitude_difference';
measurement = 'waveform_difference';
% measurement = 'cc_time_shift';


veldis = 'vel';
% veldis = 'dis';

i_ref = 2;
rec_id = 11;


if( strcmp(measurement,'cc_time_shift') )
    tau = 1:7;
    i_zero = find( t==0 );
else
    du = 2 * rand(1,length(t)) - 0.5;
    tau = logspace(-10,-2,100);
end



% each reference station will act as a source once
src = ref_stat(i_ref,:);
rec = array( find(~ismember(array,src,'rows') ) , :);


% choose one receiver
indices = (i_ref-1)*n_rec + 1 : i_ref*n_rec;
indices = indices(rec_id);
rec = rec(rec_id,:);

first(1,:) = u_h.c_data( indices , : );
second(1,:) = u_1.c_data( indices , : );


% calculate misfit for u_0 and u
[misfit, adjoint(1,:)] = misfits( first, second, t, veldis, measurement, src, rec );


% calculate derivative with adjoint source time function and du
if( ~strcmp(measurement,'cc_time_shift') )
    addu = dot( fliplr(adjoint), du(1,:) );
end


% shift only acausal branch 
if( strcmp(measurement,'cc_time_shift') )
    du = repmat( [ t>=0 ] .* first, length(tau), 1 );
    
    figure
    plot(t,first);
    hold on
end


% calculate misfit for u_0 and perturbed u
misfit_tau_du = zeros(1,length(tau));
diff_quotient = zeros(1,length(tau));

for j = 1:length(tau)
    
    
    if( strcmp(measurement,'cc_time_shift') )
        
        du( j, tau(j):i_zero ) = first( 1, 1:i_zero-tau(j)+1);
        du( j, 1:tau(j)-1) = 0.0;
        
        plot(t,du(j,:));
        
        [misfit_tau_du(j), ~] = misfits( du(j,:), second, t, veldis, measurement, src, rec );
        
        % approximate derivative with finite difference approximation
        diff_quotient(j) = ( misfit_tau_du(j) - misfit );
        
        % calculate derivative with adjoint source time function and du
        addu(j) = dot( fliplr(adjoint), du(j,:) - first );
        
    else
        
        [misfit_tau_du(j), ~] = misfits( first + tau(j) * du, second, t, veldis, measurement, src, rec );
        
        % approximate derivative with finite difference approximation
        diff_quotient(j) = ( misfit_tau_du(j) - misfit ) / tau(j);
        
    end
    
    
end


% summary
table(:,1) = tau;
table(:,2) = addu;
table(:,3) = diff_quotient;


% plot absolute difference 
figure
loglog( tau, abs( (diff_quotient-addu) / addu ) * 100 , 'b' )


rmpath(genpath('../../'))
run ../startup.m

