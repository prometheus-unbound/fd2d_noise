function name = filename( choice, n_ref )
    
    if( strcmp(choice,'array') )
        name = [ fd2d_path() 'output' filesep choice '_' num2str(n_ref) '_ref.mat' ];
    else
        [~,~,~,~,~,~,~,model_type,source_type] = input_parameters();
        name = [ fd2d_path() 'output' filesep choice '_' num2str(n_ref) '_ref_model_' num2str(model_type) '_source_' source_type '.mat' ];
    end

end