using DelimitedFiles
using SparseArrays
using LinearAlgebra

const COND_FIGURE_DIR = normpath(joinpath(@__DIR__, "..", "figures"))
mkpath(COND_FIGURE_DIR)

function condest(A::SparseMatrixCSC)
    n = size(A, 1)
    
    # Power iteration for the dominant eigenvalue of A
    v = rand(n)
    v ./= norm(v, 2)
    λ_max = 0.0
    for _ in 1:100
        w = A * v
        λ_max = norm(w, 2)
        v = w ./ λ_max
    end
    
    # Inverse power iteration for the dominant eigenvalue of A^{-1}
    F = lu(A)
    v = rand(n)
    v ./= norm(v, 2)
    λ_inv_max = 0.0
    for _ in 1:100
        w = F \ v
        λ_inv_max = norm(w, 2)
        v = w ./ λ_inv_max
    end
    
    return λ_max * λ_inv_max
end

function process_condition_number(base_dir::String, example_name::String, do_export::Bool)
    data_dir = joinpath(base_dir, example_name, "data")
    output_file = joinpath(data_dir, "output.dat")
    
    println("==================================================")
    println("Analysing Example: $example_name")
    println("==================================================")
    
    if !isfile(output_file)
        @warn "Data file not found for $example_name. Skipping."
        return
    end
    
    # Parse the output.dat header to extract taus
    lines = readlines(output_file)
    taus = parse.(Float64, split(lines[2]))
    stab_cases = length(taus)
    
    # Parse the remaining numerical data
    raw_data = Float64[]
    for line in lines[3:end]
        append!(raw_data, parse.(Float64, split(line)))
    end
    
    num_cols = 2 + stab_cases
    num_rows = div(length(raw_data), num_cols)
    data = transpose(reshape(raw_data, num_cols, num_rows))
    
    h_all = data[:, 1]
    n_export = length(h_all) # Processing all available rows
    h_vals = h_all[1:n_export]
    
    # Pre-allocate condition number matrix
    cond_nums = zeros(n_export, stab_cases)
    completed_refs = 0
    
    # Process Matrices
    for ref in 1:n_export
        println("  Refinement level $ref/$n_export")
        level_successful = true
        
        for stab in 1:stab_cases
            mat_filename = joinpath(data_dir, "$(example_name)_ref$(ref)_stab$(stab)_mat.dat")
            
            if !isfile(mat_filename)
                level_successful = false
                continue
            end
            
            # Load the COO text data
            B = readdlm(mat_filename)
            if isempty(B)
                level_successful = false
                continue
            end
            
            # Convert to a sparse matrix (Julia's equivalent of spconvert)
            I = Int.(B[:, 1])
            J = Int.(B[:, 2])
            V = Float64.(B[:, 3])
            A_global = sparse(I, J, V)
            
            # Direct condition estimation (Simplified per your request)
            try
                cond_nums[ref, stab] = condest(A_global)
            catch e
                @warn "Estimation failed for stab=$stab (Likely singular matrix). $e"
                cond_nums[ref, stab] = NaN
            end
        end
        
        if !level_successful
            break 
        end
        completed_refs = ref
    end
    
    # Calculate Rates and Plot
    if completed_refs > 0
        cond_nums_cropped = cond_nums[1:completed_refs, :]
        h_vals_cropped = h_vals[1:completed_refs]
        
        skip_idx = div(length(h_vals_cropped), 2)
        
        if completed_refs > skip_idx + 1
            rates = get_rates(h_vals_cropped, cond_nums_cropped, skip_idx)
        else
            rates = fill(NaN, stab_cases)
        end
        
        fig_filename = joinpath(COND_FIGURE_DIR, "$(example_name)_cond.pdf")
        reference_slopes = [-2.0] # O(h^-2) growth line
        
        plot_cond(taus, cond_nums_cropped, h_vals_cropped, rates, fig_filename, reference_slopes, do_export)
    else
        @warn "No valid refinement levels completed for $example_name. Skipping plot."
    end
end


function plot_cond(taus, A_mat, h, rates, filename::String, refs, do_export::Bool)
    markers = [:circle, :cross, :x, :square, :utriangle, :star5]
    
    taus_str = map(taus) do t
        t == 0 ? "0" : "10^{$(round(Int, log10(t)))}"
    end
    
    p = plot(xaxis = :log, yaxis = :log, 
             legend = :topright, 
             xlabel = L"h", 
             ylabel = L"\text{Approximate Condition Number } \kappa_2(A)",
             framestyle = :box,
             extra_kwargs = Dict(:subplot => Dict("ylabel shift" => "-3mm", "xlabel shift" => "-3mm"))
    )
    
    for k in 1:length(taus)
        current_marker = markers[mod1(k, length(markers))]
        m_size = current_marker in [:cross, :x, :plus] ? 8 : 5
        
        box_width = "2.0cm"
        leg_label = L"\makebox[%$box_width][l]{$\tau = %$(taus_str[k])$} \text{rate} \approx %$(round(rates[k], digits=2))"
        
        plot!(p, h, A_mat[:, k], 
              marker = current_marker,
              markersize = m_size,         
              label = leg_label)
    end
    
    if !isempty(refs)
        first_row = A_mat[1, :]
        base_val = exp(sum(log.(first_row)) / length(first_row))
        
        for ref in refs
            ref_y = base_val .* [(h[end]/h[1])^ref, 1.0]
            plot!(p, [h[end], h[1]], ref_y, 
                  color = :black, 
                  linestyle = :dash, 
                  label = "")
        end
    end
    
    if do_export
        mkpath(dirname(filename))
        savefig(p, filename)
    end
    
    return p
end