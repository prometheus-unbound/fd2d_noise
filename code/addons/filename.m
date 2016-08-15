function name = filename( choice, number )
    
    if( strcmp(choice,'array') )
        name = [ fd2d_path() 'output' filesep 'array_nref-' num2str(number) '.mat' ];
    elseif( strcmp(choice,'G_fft') )
        [~,~,~,~,~,~,~,model_type] = input_parameters();
        name = [ fd2d_path() 'output' filesep 'G_fft_iref-' num2str(number) '_model-' num2str(model_type) '.mat' ];
    else
        [~,~,~,~,~,~,~,model_type,source_type] = input_parameters();
        name = [ fd2d_path() 'output' filesep choice '_nref-' num2str(number) '_model-' num2str(model_type) '_source-' source_type '.mat' ];
    end

end