function process_example_data(base_dir, example_name, reference_rates, export, ylims)
    % Construct the file path
    data_path = fullfile(base_dir,  example_name, 'data', 'output.dat');
    
    % Open the file for reading
    fid = fopen(data_path, 'r');
    if fid == -1
        warning('Cannot open file: %s. Skipping this dataset.', data_path);
        return;
    end
    
    % Read the header lines
    filename_base = strtrim(string(fgetl(fid)));
    tau_line = fgetl(fid);
    tau = str2num(tau_line)'; 
    
    % Read the remaining numerical data
    raw_data = textscan(fid, '%f');
    fclose(fid);
    
    % Reshape the 1D array into a matrix
    num_cols = 2 + length(tau);
    data = reshape(raw_data{1}, num_cols, [])';
    
    % Parse the columns
    h = data(:, 1);
    n_err = data(:, 2); 
    error_mat = data(:, 3:end);
    
    % Calculate convergence rates using your provided function
    rates = get_rates(h, error_mat, floor(length(h)/2));
    n_rate = get_rates(h, n_err, floor(length(h)/2));
    
    % Print the normal error convergence rate to the console
    fprintf('Convergence rate for n_err (%s): %.4f\n', example_name, n_rate);
    
    % Construct final filename and generate the plot
    filename = example_name + ".pdf";
    plot_conv(tau, error_mat, h, rates, filename, reference_rates, export, ylims);

    % plot_n_conv(n_err, n_rate, h, 1)
end