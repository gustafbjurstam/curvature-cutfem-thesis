using Plots
using Plots.Measures 
using LaTeXStrings
using LinearAlgebra

pgfplotsx() 

push!(Plots.PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{amsmath}")

const FIGURE_DIR = normpath(joinpath(@__DIR__, "..", "figures"))
mkpath(FIGURE_DIR)

matlab_colors = [
    "#0072BD", # Blue
    "#D95319", # Orange
    "#EDB120", # Yellow
    "#8E35A1", # Purple (Lightened from #7E2F8E)
    "#77AC30", # Green
    "#4DBEEE", # Cyan
    "#A2142F"  # Red
]

default(
    palette = matlab_colors,
    
    markercolor       = :auto, 
    markerstrokecolor = :auto, 
    markerstrokewidth = 1.5,    
    
    guidefontsize    = 18, # Axis labels
    titlefontsize    = 18, 
    tickfontsize     = 15, # Axis tick numbers
    legendfontsize   = 14, 
    legendfonthalign = :left, 
    
    linewidth  = 2.0,     
    framestyle = :box,
    grid       = false,
    gridalpha  = 0.3,
    minorticks = true,

    tick_direction = :in,
    
    size = (600, 450)
)

function get_rates(h::AbstractVector, error_mat::AbstractMatrix, skip::Int)
    h_fit = h[skip+1:end]
    err_fit = error_mat[skip+1:end, :]
    
    A = hcat(log.(h_fit), ones(length(h_fit)))
    
    rates = A \ log.(err_fit)
    
    return rates[1, :]
end

# Multiple dispatch: an overloaded version specifically for 1D error vectors
function get_rates(h::AbstractVector, err_vec::AbstractVector, skip::Int)
    h_fit = h[skip+1:end]
    err_fit = err_vec[skip+1:end]
    
    A = hcat(log.(h_fit), ones(length(h_fit)))
    rates = A \ log.(err_fit)
    
    return rates[1]
end

function plot_conv(taus, A_mat, h, rates, filename::String, refs, do_export::Bool, y_limits)
    markers = [:circle, :cross, :x, :square, :utriangle, :star5]
    
    taus_str = map(taus) do t
        t == 0 ? "0" : "10^{$(round(Int, log10(t)))}"
    end
    
    y_label_str = L"\lVert H_h - H^e \rVert_{L^2(\Gamma_h)}"
    
    p = plot(xaxis = :log, yaxis = :log, 
             legend = :bottomright, 
             xlabel = L"h", 
             ylabel = y_label_str,
             framestyle = :box,
             extra_kwargs = Dict(:subplot => Dict("ylabel shift" => "-3mm", "xlabel shift" => "-3mm"))
            )
    
    for k in 1:length(taus)
        current_marker = markers[mod1(k, length(markers))]
        m_size = current_marker in [:cross, :x, :plus] ? 8 : 5
        
        box_width = "2.0cm"
        
        leg_label = L"\makebox[%$box_width][l]{$\tau = %$(taus_str[k])$} $\text{rate} \approx %$(round(rates[k], digits=2))$"
        
        plot!(p, h, A_mat[:, k], 
              marker = current_marker,
              markersize = m_size,         
              label = leg_label)
    end
    
    if !isempty(refs)
        col_maxes = maximum(A_mat, dims=1)
        base_val = exp(sum(log.(col_maxes)) / length(col_maxes))
        
        for ref in refs
            ref_y = base_val .* [(h[end]/h[1])^ref, 1.0]
            plot!(p, [h[end], h[1]], ref_y, 
                  color = :black, 
                  linestyle = :dash, 
                  label = "")
        end
    end
    
    if !isnothing(y_limits)
        ylims!(p, y_limits)
    end
    
    if do_export
        mkpath(dirname(filename))
        savefig(p, filename)
    end
    
    return p
end

function process_example_data(base_dir::String, example_name::String, reference_rates, do_export::Bool, y_limits=nothing)
    data_path = joinpath(base_dir, example_name, "data", "output.dat")
    
    if !isfile(data_path)
        @warn "Cannot open file: $data_path. Skipping this dataset."
        return
    end
    
    lines = readlines(data_path)
    filename_base = strip(lines[1])
    tau = parse.(Float64, split(lines[2]))
    
    raw_data = Float64[]
    for line in lines[3:end]
        append!(raw_data, parse.(Float64, split(line)))
    end
    
    num_cols = 2 + length(tau)
    num_rows = div(length(raw_data), num_cols)
    data = transpose(reshape(raw_data, num_cols, num_rows))
    
    h = data[:, 1]
    n_err = data[:, 2]
    error_mat = data[:, 3:end]
    
    skip = floor(Int, length(h) / 2)
    rates = get_rates(h, error_mat, skip)
    n_rate = get_rates(h, n_err, skip)
    
    println("Convergence rate for n_err ($example_name): $(round(n_rate, digits=4))")
    
    filename = joinpath(FIGURE_DIR, "$example_name.pdf")
    plot_conv(tau, error_mat, h, rates, filename, reference_rates, do_export, y_limits)
end