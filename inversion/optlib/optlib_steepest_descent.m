function [flag, mfinal, usr_par] = optlib_steepest_descent(m0, options, usr_par)
    

    options = optlib_init_default_parameters_steepest_descent(options);

    if (options.verbose)
        options
    end

    % constant 0<del<1/2 for Armijo condition
    delta = 0.001;
    % constant del<theta<1 for Wolfe condition
    theta = 0.6;

    sigma = options.init_step_length;
    m = m0;
    [j, g] = eval_objective_and_gradient(m, 'yo', usr_par);

    normg0 = norm(g);
    normg = normg0;

    it = 0;
    fid = fopen([fd2d_path() filesep 'inversion' filesep options.output_file], 'a+');

    fprintf(fid, 'it=%d   j=%e   ||g||=%e  \n', it, j, normg);
    fprintf('it=%d   j=%e   ||g||=%e  \n', it, j, normg);

    % main loop
    while (normg > options.tolerance * normg0 && it < options.max_iterations)

        it = it + 1;
        sigma0 = sigma;
        s = - g;
        stg = s' * g;

        if (options.stepsize == 0)
            % choose sigma by Armijo stepsize rule starting with previous
            % stepsize sig0. If sig0 is acceptable try sig=2^k*sig0.
            [sigma, jn] = optlib_armijo(m, s, stg, j, delta, sigma0, usr_par);
        else
            % choose sigma by Powell-Wolfe stepsize rule starting with previous
            % stepsize sig0.
            [sigma, jn, gn] = optlib_stepsize_wolfe(m, s, stg, j, delta, theta, sigma0, usr_par);
        end
        m = m + sigma * s;
        j = jn;
        if (options.stepsize == 0)
            [gn] = eval_grad_objective(m, 'yo', usr_par);
        end
        g = gn;
        normg = norm(g);

        % save(sprintf([fd2d_path(), filesep, 'inversion', filesep, 'models', filesep, 'model_%i.mat'], it), 'm', 'gn', 'jn')
        fprintf(fid, 'it=%3.d   f=%e   ||g||=%e   sig=%5.3f\n', it, j, normg, sigma);
        fprintf('it=%3.d   f=%e   ||g||=%e   sig=%5.3f\n', it, j, normg, sigma);
    end

    if (normg <= options.tolerance * normg0)
        fprintf(fid, 'Successful termination with ||g||<%e*||g0||:\n', options.tolerance);
        fprintf('Successful termination with ||g||<%e*||g0||:\n', options.tolerance);
        flag = 0;
    else
        fprintf(fid, 'Maximum number of iterations reached.\n');
        fprintf('Maximum number of iterations reached.\n');
        flag = 1;
    end
    
    mfinal = m;
    fclose(fid);

end
