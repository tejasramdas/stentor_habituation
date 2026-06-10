using FillArrays, ReverseDiff, Turing, Random, Mooncake, DataFrames, Colors, GLMakie

param_lims = Dict(:a_1_0 => (0, 30), :a_2_0 => (0, 30), :a_1_stdev => (0, 5), :a_2_stdev => (0, 5), :b_1_0 => (0, 60), :c_1_0 => (-0.1, 1.1), :d_1_0 => (-0.1, 1.1), :b_2_0 => (0, 60), :c_2_0 => (-0.1, 1.1), :d_2_0 => (-0.1, 1.1), :b_1_conc => (0, 30), :b_2_conc => (0, 30), :c_1_conc => (-0.1, 1.1), :c_2_conc => (-0.1, 1.1), :d_1_conc => (-0.1, 1.1), :d_2_conc => (-0.1, 1.1))
param_labels = Dict(:a_1_0 => "Hill coeff (T1)", :a_2_0 => "Hill coeff (T2)", :a_1_stdev => "Hill coeff spread (T1)", :a_2_stdev => "Hill coeff spread (T2)", :b_1_0 => "N50 (trial 1)", :c_1_0 => "Initial response (trial 1)", :d_1_0 => "Asymptote (trial 1)", :b_2_0 => "N50 (trial 2)", :c_2_0 => "Initial response (trial 2)", :d_2_0 => "Asymptote (trial 2)", :b_1_conc => "N50 conc (T1)", :b_2_conc => "N50 conc (T2)", :c_1_conc => "T1 initial conc", :c_2_conc => "T2 initial conc", :d_1_conc => "T1 asymptote conc", :d_2_conc => "T2 asymptote conc", :b_1 => "N50 (trial 1)", :c_1 => "Initial response (trial 1)", :d_1 => "Asymptote (trial 1)", :b_2 => "N50 (trial 2)", :c_2 => "Initial response (trial 2)", :d_2 => "Asymptote (trial 2)")
# =============================================================================
# Math Functions and Models
# =============================================================================
begin
    """Hill function with parameters a (slope), b (N50), c (max), d (min)."""
    function hill(x::Int, a::Float64=1.0, b::Float64=30.0, c::Float64=1.0, d::Float64=0.0)
        return (c - d) / (1 + ((x - 0.99) / b)^a) + d
    end
    function hill(x::Float64, a::Float64=1.0, b::Float64=30.0, c::Float64=1.0, d::Float64=0.0)
        return (c - d) / (1 + ((x - 0.99) / b)^a) + d
    end

    """Analytical derivative of the Hill function with respect to x."""
    function hill_derivative(x, a=1, b=30, c=1, d=0)
        u = (x - 0.99) / b
        z = 1 + u^a
        return -(c - d) * (a / b) * u^(a - 1) / (z^2)
    end

    """Analytical derivative of the Hill function with respect to hill(x)."""
    function hill_derivative_y(y, a=1, b=30, c=1, d=0)
        return -(a / (b * (c - d))) * ((y - d)^((a + 1) / a)) * ((c - y)^((a - 1) / a))
    end

end

# =============================================================================
# Turing models
# =============================================================================
begin
    @model function single_trial(num_stim::Int, num_cells::Int)
        a_1_0 ~ truncated(Normal(0, 4), -4, 4)
        a_1_stdev ~ Gamma(1, 1)

        #=b_1_0 ~ Gamma(2, 10)=#
        b_1_0 ~ truncated(Gamma(2, 10), 2, 60)
        b_1_conc ~ Gamma(2, 5)

        c_1_0 ~ Beta(1, 1)
        c_1_conc ~ Gamma(5, 2)

        d_1_0 ~ Beta(1, 1)
        d_1_conc ~ Gamma(5, 2)

        a_1_raw ~ filldist(Normal(0, 1), num_cells)
        # a_1 = exp.(2.5 .+ a_1_stdev .* a_1_raw)
        a_1 = exp.(a_1_0 .+ a_1_stdev .* a_1_raw)
        #=b_1 ~ filldist(Gamma(b_1_conc, b_1_0 / b_1_conc), num_cells)=#
        b_1 ~ filldist(truncated(Gamma(b_1_conc, b_1_0 / b_1_conc), 2, 60), num_cells)
        c_1 ~ filldist(Beta((c_1_0) * c_1_conc + 1e-10, (1 - c_1_0) * c_1_conc + 1e-10), num_cells)
        d_1 ~ filldist(Beta((d_1_0) * d_1_conc + 1e-10, (1 - d_1_0) * d_1_conc + 1e-10), num_cells)

        function cell_dist(a::Real, b::Real, c::Real, d::Real)
            cell_func(i::Int) = hill(i, a, b, c, d)
            series_dist = arraydist(map(i -> Bernoulli(cell_func(i)), 1:num_stim))
        end

        t1_contractions ~ arraydist(map(i -> cell_dist(a_1[i], b_1[i], c_1[i], d_1[i]), 1:num_cells))
        # t1_contractions ~ arraydist(map(i -> cell_dist(10.0, b_1[i], c_1[i], d_1[i]), 1:num_cells))
    end
    cond_param_single_trial(num_stim, num_cells, params) = condition(single_trial(num_stim, num_cells), (a_1_0=params[:a_1_0], a_1_stdev=params[:a_1_stdev], b_1_0=params[:b_1_0], b_1_conc=params[:b_1_conc], c_1_0=params[:c_1_0], c_1_conc=params[:c_1_conc], d_1_0=params[:d_1_0], d_1_conc=params[:d_1_conc]))
    cond_data_single_trial(responses) = condition(single_trial(size(responses, 1), size(responses, 2)), (t1_contractions = responses))

end


# =============================================================================
# Data Extraction and Parameter Processing Functions
# =============================================================================
begin
    # UNUSED: extract_responses_turing — never called from paper.jl
    """Extract response data from Turing chain for specified number of stimuli and cells."""
    function extract_responses_turing(chain, num_stim, num_cells)
        t1 = reshape(Array(group(chain, :t1_contractions)), num_stim, num_cells)
        t2 = reshape(Array(group(chain, :t2_contractions)), num_stim, num_cells)
        return vcat(t1, t2)
    end

    function compute_qt_b(b_1_0, b_1_conc; q=[0.5])
        b_1 = quantile.(truncated.(Gamma.(b_1_conc, b_1_0 ./ b_1_conc), 2, 60), q)
    end

    function compute_qt_c(c_1_0, c_1_conc; q=[0.5])
        c_1 = quantile.(Beta.(c_1_0 .* c_1_conc, (1 .- c_1_0) .* c_1_conc), q)
    end

    function compute_qt_d(d_1_0, d_1_conc; q=[0.5])
        d_1 = quantile.(Beta.(d_1_0 .* d_1_conc, (1 .- d_1_0) .* d_1_conc), q)
    end

    """Extract all population-level parameters from Turing chain as vectors (not means)."""
    function extract_params_samples_turing(chain)
        cell_params = [:b_1, :c_1, :d_1]

        samps = Dict{Symbol,Any}(
            :cell => Dict(k => Array(group(chain, k)) for k in cell_params)
        )
        samps[:pop] = Dict{Symbol,Vector{Float64}}()
        a_1_0 = vec(Array(group(chain, :a_1_0)))
        b_1_0 = vec(Array(group(chain, :b_1_0)))
        b_1_conc = vec(Array(group(chain, :b_1_conc)))
        c_1_0 = vec(Array(group(chain, :c_1_0)))
        c_1_conc = vec(Array(group(chain, :c_1_conc)))
        d_1_0 = vec(Array(group(chain, :d_1_0)))
        d_1_conc = vec(Array(group(chain, :d_1_conc)))
        samps[:pop][:a_1_0] = exp.(a_1_0)
        samps[:pop][:b_1_0] = vec(compute_qt_b(b_1_0, b_1_conc, q=[0.5]))
        samps[:pop][:c_1_0] = vec(compute_qt_c(c_1_0, c_1_conc, q=[0.5]))
        samps[:pop][:d_1_0] = vec(compute_qt_d(d_1_0, d_1_conc, q=[0.5]))
        samps[:pop][:b_1_0_low] = vec(compute_qt_b(b_1_0, b_1_conc, q=[0.025]))
        samps[:pop][:c_1_0_low] = vec(compute_qt_c(c_1_0, c_1_conc, q=[0.025]))
        samps[:pop][:d_1_0_low] = vec(compute_qt_d(d_1_0, d_1_conc, q=[0.025]))
        samps[:pop][:b_1_0_high] = vec(compute_qt_b(b_1_0, b_1_conc, q=[0.975]))
        samps[:pop][:c_1_0_high] = vec(compute_qt_c(c_1_0, c_1_conc, q=[0.975]))
        samps[:pop][:d_1_0_high] = vec(compute_qt_d(d_1_0, d_1_conc, q=[0.975]))
        a_1_true = exp.(Array(group(chain, :a_1_0)) .+ Array(group(chain, :a_1_stdev)) .* Array(group(chain, :a_1_raw)))
        samps[:cell][:a_1] = a_1_true
        return samps
    end
    """Extract and compute median parameter values from Turing chain."""
    function extract_params_qt_turing(chain; q=[0.025, 0.5, 0.975])
        # pop_params = [:a_1_0, :a_2_0, :a_1_stdev, :a_2_stdev, :b_1_0, :b_1_conc, :b_2_0, :b_2_conc, :c_1_0, :c_1_conc, :c_2_0, :c_2_conc, :d_1_0, :d_1_conc, :d_2_0, :d_2_conc]
        # cell_params = [:a_1_raw, :a_2_raw, :b_1, :c_1, :d_1, :b_2, :c_2, :d_2]
        pop_params = [:a_1_0, :b_1_0, :c_1_0, :d_1_0]
        cell_params = [:a_1, :b_1, :c_1, :d_1]
        samps = extract_params_samples_turing(chain)
        params_dict = Dict(
            :pop => Dict(v => quantile(samps[:pop][v], q) for v in pop_params),
            :cell => Dict(v => hcat([quantile(samps[:cell][v][:, i], q) for i in 1:size(samps[:cell][v], 2)]...) for v in cell_params)
        )
        return params_dict
    end

    """Extract and compute median parameter values from Turing chain."""
    function extract_params_median_turing(chain)
        # pop_params = [:a_1_0, :a_2_0, :a_1_stdev, :a_2_stdev, :b_1_0, :b_1_conc, :b_2_0, :b_2_conc, :c_1_0, :c_1_conc, :c_2_0, :c_2_conc, :d_1_0, :d_1_conc, :d_2_0, :d_2_conc]
        # cell_params = [:a_1_raw, :a_2_raw, :b_1, :c_1, :d_1, :b_2, :c_2, :d_2]
        pop_params = [:a_1_0, :b_1_0, :c_1_0, :d_1_0]
        cell_params = [:a_1, :b_1, :c_1, :d_1]
        qts = extract_params_qt_turing(chain)
        params_dict = Dict(
            :pop => Dict(v => qts[:pop][v][2] for v in pop_params),
            :cell => Dict(v => qts[:cell][v][2, :] for v in cell_params)
        )
        return params_dict
    end

    function samples_to_param_vec(samps)
        n_samples = length(first(values(samps[:pop])))
        return [Dict(
            :pop => Dict(k => v[i] for (k, v) in samps[:pop]),
            :cell => Dict(k => v[i, :] for (k, v) in samps[:cell])
        ) for i in 1:n_samples]
    end




    function get_inferred_hill(params; stim=1:60, num_trials=1)
        if num_trials == 1
            a_1 = params[:pop][:a_1_0]
        else
            a_1 = params[:pop][:a_1_0]
            a_2 = params[:pop][:a_2_0]
        end
        trial_1 = hill.(collect(stim), a_1, params[:pop][:b_1_0], params[:pop][:c_1_0], params[:pop][:d_1_0])
        return trial_1
    end

    function get_inferred_hill_sc(params; stim=1:60)
        n_cells = length(params[:cell][:b_1])
        a_1 = params[:cell][:a_1]
        mat_1 = zeros(length(collect(stim)), n_cells)
        for i in 1:n_cells
            mat_1[:, i] = hill.(collect(stim), a_1[i], params[:cell][:b_1][i], params[:cell][:c_1][i], params[:cell][:d_1][i])
        end
        return mat_1
    end


    function calculate_response_diff(params, stim1, stim2)
        curve = get_inferred_hill(params)
        return curve[stim1] .- curve[stim2]
    end

    """
    Compute the ratio of habituation trajectories between two trials.
    For each trial, calculates hill(stim_i) - hill(1) for all stimuli,
    then returns (T2 diffs) ./ (T1 diffs) per stimulus.
    """
    function trial_diff_ratio(params_t1, params_t2; stim=5, norm=true)
        t1_diff = calculate_response_diff(params_t1, 1, stim)
        t2_diff = calculate_response_diff(params_t2, 1, stim)
        if norm
            t1_diff = t1_diff ./ get_inferred_hill(params_t1, stim=stim)[1]
            t2_diff = t2_diff ./ get_inferred_hill(params_t2, stim=stim)[1]
        end
        return log.(t2_diff ./ t1_diff)
        # return t1_diff ./ t2_diff
    end

    function trial_discrete_derivative_ratio(params_t1, params_t2; step=5)
        stim = collect(1:55)
        t1_diffs = calculate_response_diff(params_t1, stim, stim .+ step)
        t2_diffs = calculate_response_diff(params_t2, stim, stim .+ step)
        return log.(t2_diffs ./ t1_diffs)
    end

    function calculate_statistic_ci(chain, stat_func::Function; prob=0.95)
        params = extract_params_samples_turing(chain)
        vals = stat_func(params)
        # Convert to array and flatten (handles Chains objects or multi-chain arrays)
        vals_flat = vec(Array(vals))
        lower_q = (1 - prob) / 2
        upper_q = 1 - lower_q
        return quantile(vals_flat, [lower_q, 0.5, upper_q])
    end

    function calculate_statistic_ci(chain1, chain2, stat_func::Function; prob=0.95)
        params1 = extract_params_samples_turing(chain1)
        params2 = extract_params_samples_turing(chain2)
        vals = stat_func(params1, params2)
        # Convert to array and flatten (handles Chains objects or multi-chain arrays)
        vals_flat = vec(Array(vals))
        lower_q = (1 - prob) / 2
        upper_q = 1 - lower_q
        return quantile(vals_flat, [lower_q, 0.5, upper_q])
    end

    function calculate_statistic_ci_vec(chain1, chain2, stat_func::Function; prob=0.95)
        params1 = samples_to_param_vec(extract_params_samples_turing(chain1))
        params2 = samples_to_param_vec(extract_params_samples_turing(chain2))
        vals = stat_func.(params1, params2)
        vals_flat = vec(Array(vals))
        lower_q = (1 - prob) / 2
        upper_q = 1 - lower_q
        return quantile(vals_flat, [lower_q, 0.5, upper_q])
    end
    function calculate_statistic_ci_vec(chain1, chain2, chain3, chain4, stat_func::Function, compare_func::Function; prob=0.95)
        params1 = samples_to_param_vec(extract_params_samples_turing(chain1))
        params2 = samples_to_param_vec(extract_params_samples_turing(chain2))
        params3 = samples_to_param_vec(extract_params_samples_turing(chain3))
        params4 = samples_to_param_vec(extract_params_samples_turing(chain4))
        vals_12 = stat_func.(params1, params2)
        vals_34 = stat_func.(params3, params4)
        vals = compare_func(vals_12, vals_34)
        # Convert to array and flatten (handles Chains objects or multi-chain arrays)
        vals_flat = vec(Array(vals))
        lower_q = (1 - prob) / 2
        upper_q = 1 - lower_q
        return quantile(vals_flat, [lower_q, 0.5, upper_q])
    end


    function calculate_statistic_ci_corr(chain1, chain2, p1_func::Function, p2_func::Function; prob=0.95)
        # params[:cell][p1] is (num_cells, num_samples)
        function corr_samples(params1, params2)
            arr1 = p1_func(params1, params2)
            arr2 = p2_func(params1, params2)
            num_samples = size(arr1, 2)
            corrs = [cor(arr1[i, :], arr2[i, :]) for i in 1:num_samples]
        end
        return calculate_statistic_ci(chain1, chain2, corr_samples; prob=prob)
    end
end




# =============================================================================
# Visualization and Plotting Functions
# =============================================================================
begin
    """Plot Hill function curves with specified parameters."""
    function plot_params(ax, a, b, c, d; color=nothing, label="")
        x = 1:0.1:60
        y = hill.(x, a, b, c, d)
        lines!(ax, x, y, label=label, color=(isnothing(color) ? :black : color))
        ylims!(ax, -0.1, 1.1)
    end

    function plot_params_derivative(a, b, c, d, ax; color=nothing, label="", discrete=false)
        if discrete
            x = 1:0.1:60.9
            y = hill.(x, a, b, c, d)
            dy = diff(y) ./ 0.1
            x = x[1:end-1]
        else
            x = 1:0.1:60
            dy = hill_derivative.(x, a, b, c, d)
        end
        lines!(ax, x, dy, label=label, color=(isnothing(color) ? :black : color))
    end

    function plot_hill_vs_derivative(a, b, c, d, ax; color=nothing, label="", discrete=false, plot_points=true)
        if discrete
            x = 1:0.1:60.9
            y = hill.(x, a, b, c, d)
            dy = diff(y) ./ 0.1
            x = x[1:end-1]
            y = y[1:end-1]
        else
            x = 1:0.1:60
            y = hill.(x, a, b, c, d)
            dy = hill_derivative.(x, a, b, c, d)
        end
        lines!(ax, y, dy, label=label, color=isnothing(color) ? :black : color)
        if discrete && plot_points
            scatter!(ax, y[1:10:end], dy[1:10:end], color=isnothing(color) ? :black : color)
        end
        # ax.xlabel = "Response probability"
        # ax.ylabel = "Learning rate"
    end

    """Plot Hill function curves for both trials with different parameters."""
    function plot_curves(ax, params; label="Trial 1", color=:black, axlabel=true, half_max=false)
        plot_params(ax, params[:pop][:a_1_0], params[:pop][:b_1_0], params[:pop][:c_1_0], params[:pop][:d_1_0], color=color, label=label)
        # Half-max dot at x = b (N50), y = (c + d) / 2
        if half_max
            b = params[:pop][:b_1_0]
            halfmax_y = (params[:pop][:c_1_0] + params[:pop][:d_1_0]) / 2
            scatter!(ax, [b + 1], [halfmax_y], color=color, markersize=10)
        end
        if axlabel
            ax.xlabel = "Stimulus number"
            ax.ylabel = "Inferred single-cell \nresponse probability"
        end
    end

    function plot_band(ax, params; color=:black, plot_x=1:0.1:60)
        x_pts = 1:0.1:60
        curves = hcat([hill.(x_pts, Ref(params[:pop][:a_1_0][i]), Ref(params[:pop][:b_1_0][i]), Ref(params[:pop][:c_1_0][i]), Ref(params[:pop][:d_1_0][i])) for i in 1:size(params[:pop][:a_1_0], 1)]...)
        lo = [quantile(curves[i, :], 0.025) for i in 1:size(curves, 1)]
        hi = [quantile(curves[i, :], 0.975) for i in 1:size(curves, 1)]
        band!(ax, plot_x, lo, hi, color=(color, 0.3))
    end

    # curves = hcat([hill.(1:0.1:60, Ref(b1[:pop][:a_1_0][i]), Ref(b2[:pop][:b_1_0][i]), Ref(b1[:pop][:c_1_0][i]), Ref(b1[:pop][:d_1_0][i])) for i in 1:size(b1[:pop][:a_1_0], 1)]...)
    # lo = [quantile(curves[i, :], 0.025) for i in 1:size(curves, 1)]
    # hi = [quantile(curves[i, :], 0.975) for i in 1:size(curves, 1)]
    # med = [quantile(curves[i, :], 0.5) for i in 1:size(curves, 1)]
    """Plot Hill function curves for both trials with different parameters using derivatives."""
    function plot_curves_derivative(ax, params; label="Trial 1", color=:black, axlabel=true)
        plot_params_derivative(params[:pop][:a_1_0], params[:pop][:b_1_0], params[:pop][:c_1_0], params[:pop][:d_1_0], ax, color=color, label=label)
        if axlabel
            ax.xlabel = "Stimulus number"
            ax.ylabel = "Learning rate"
        end
    end

    """Plot Hill function curves vs their derivatives for both trials with different parameters."""
    function plot_curves_vs_derivatives(ax, params; label="Trial 1", color=:black, axlabel=false)
        plot_hill_vs_derivative(params[:pop][:a_1_0], params[:pop][:b_1_0], params[:pop][:c_1_0], params[:pop][:d_1_0], ax, color=color, label=label, discrete=true)
        if axlabel
            ax.xlabel = "Response probability"
            ax.ylabel = "Learning rate"
        end
    end

    function plot_cell_curves(ax, chain; mode="median", color=:black)
        if mode == "median"
            mat_1 = get_inferred_hill_sc(extract_params_median_turing(chain))
        elseif mode == "mean"
            mat_1 = get_inferred_hill_sc(extract_params_mean_turing(chain))
        end
        for i in 1:size(mat_1, 2)
            lines!(ax, 1:60, mat_1[:, i], color=color, alpha=0.1)
        end
    end

    """Plot inferred sc median response curves for both trials."""
    function plot_inferred_population_sc_median(ax, chain; chain2=nothing)
        trial_1 = get_inferred_hill(extract_params_median_turing(chain))
        lines!(ax, 1:60, trial_1, color=:black, label="Trial 1")
        if !isnothing(chain2)
            trial_2 = get_inferred_hill(extract_params_median_turing(chain2))
            lines!(ax, 1:60, trial_2, color=:gray, label="Trial 2")
        end
    end

    """Plot inferred population median response curves for both trials."""
    function plot_inferred_population_mean_median(ax, chain; chain2=nothing)
        mat_1 = get_inferred_hill_sc(extract_params_median_turing(chain))
        lines!(ax, 1:60, vec(mean(mat_1, dims=2)), color=:black, label="Trial 1")
        if !isnothing(chain2)
            mat_2 = get_inferred_hill_sc(extract_params_median_turing(chain2))
            lines!(ax, 1:60, vec(mean(mat_2, dims=2)), color=:gray, label="Trial 2")
        end
    end

    """Plot inferred population mean response curves for both trials."""
    function plot_inferred_population_mean_mean(ax, chain)
        mat_1, mat_2 = get_inferred_hill_sc(extract_params_mean_turing(chain))
        lines!(ax, 1:60, vec(mean(mat_1, dims=2)), color=:black, label="Trial 1")
        lines!(ax, 1:60, vec(mean(mat_2, dims=2)), color=:gray, label="Trial 2")
    end

    """Single cell curves as histogram. Each row is response probability, each column is stimulus number."""
    function plot_inferred_population_hist(ax, chain)
        mat_1, mat_2 = get_inferred_hill(extract_params_median_turing(chain), stim=1:0.1:60.9)
        bins = collect(0.0:0.05:1)
        # Helper to get (n_stim, n_bins) histogram array from (n_stim, n_cells)
        function histmat(mat)
            hmat = hcat([fit(Histogram, mat[i, :], bins, closed=:right).weights for i in axes(mat, 1)]...)' / size(mat, 1)
        end
        # Stack trial 1 and 2 vertically for heatmap
        hmat_1 = histmat(mat_1)
        hmat_2 = histmat(mat_2)
        # (Row 1..60: trial 1, 61..120: trial 2)
        hmat_stacked = vcat(hmat_1, hmat_2)
        # Make the y axis: 1:120, x axis: bin_centers
        heatmap!(ax, 1:0.1:120.9, 0:0.05:1, hmat_stacked; colormap=:viridis)
        ax.xlabel = "Stimulus number"
        ax.ylabel = "Response probability"
        ylims!(ax, 0, 1)
        xlims!(ax, 1, 120.9)
        # ax.yticks = (collect(0:0.2:0.99) .- 0.025, collect(0:0.2:0.99))
        ax.xticks = ([1, 20, 40, 60, 80, 100, 120], ["1", "20", "40", "60", "80", "100", "120"])
    end

    """Single cell curves as heatmap. Each row is cell, each column is stimulus number."""
    function plot_inferred_sc_heatmap(ax, chain; chain2=nothing, data=nothing, plot_data=true)
        mat_1 = get_inferred_hill_sc(extract_params_median_turing(chain), stim=1:0.1:60.9)
        if !isnothing(chain2)
            mat_2 = get_inferred_hill_sc(extract_params_median_turing(chain2), stim=1:0.1:60.9)
            stacked_mat = vcat(mat_1, mat_2)
        else
            stacked_mat = mat_1
        end
        # Sort columns by overall sum (descending)
        if plot_data
            sort_inds = sortperm(vec(sum(data, dims=1)))
        else
            sort_inds = sortperm(vec(sum(stacked_mat, dims=1)))
        end
        sorted_data = data[:, sort_inds]
        rects = [Rect(i[1], i[2] - 0.5, 1, 1) for i in findall(sorted_data .== 1)]
        sorted_mat = stacked_mat[:, sort_inds]
        if !isnothing(chain2)
            heatmap!(ax, 1:0.1:120.9, 1:size(sorted_mat, 2), sorted_mat; colormap=:viridis, colorrange=(0, 1))
        else
            heatmap!(ax, 1:0.1:60.9, 1:size(sorted_mat, 2), sorted_mat; colormap=:viridis, colorrange=(0, 1))
        end
        if plot_data
            poly!(ax, rects; color=RGBA(0, 0, 0, 0.0), strokecolor=:red, strokewidth=1.5)
        end
        ax.xlabel = "Stimulus number"
        ax.xticks = ([1, 20, 40, 60, 80, 100, 120], ["1", "20", "40", "60", "80", "100", "120"])
        ax.yticks = 1:20:size(sorted_mat, 2)
        ylims!(ax, 0.5, size(sorted_mat, 2) + 0.5)
        ax.ylabel = "Cell number"
        return sort_inds
    end

    """Plot inferred single-cell heatmaps for all ISI/ITI conditions in a grid."""
    function plot_all_inferred_sc_heatmaps(chains, responses; isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5])
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="Cell number", horz_label="Stimulus number", fig_size=(1600, 1200))
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                key = get_key(isi * 60, iti * 3600)
                data = responses[key]["control_data"]
                plot_inferred_sc_heatmap(axs[isi, j], chains[key]; data=data, plot_data=true)
            end
        end
        return fig, axs
    end

    """Plot log(T2 diff / T1 diff) relative to stimulus 1 for all ISI/ITI conditions."""
    function plot_all_trial_diff_ratio(chains; isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5], num_stim=60)
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="log(T2 / T1) habituation ratio", horz_label="Stimulus number")
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                key = get_key(isi * 60, iti * 3600)
                t1_params = extract_params_median_turing(chains[key][:trial_1])
                t2_params = extract_params_median_turing(chains[key][:trial_2])
                ratio = trial_diff_ratio(t1_params, t2_params; num_stim=num_stim)
                lines!(axs[isi, j], 2:num_stim, ratio, color=:black)
                hlines!(axs[isi, j], [0.0], color=:gray, linestyle=:dash)
                vlines!(axs[isi, j], [10.0], color=:gray, linestyle=:dash)
                xlims!(axs[isi, j], 1, num_stim)
            end
        end
        return fig, axs
    end

    """Plot habituation curves (T1/T2) and trial diff ratio side by side for a single condition."""
    function plot_condition_diff_ratio(chain; num_stim=60, fig=nothing)
        if isnothing(fig)
            fig = Figure(size=(900, 400), fontsize=20)
        end
        ax1 = Axis(fig[1, 1], xlabel="Stimulus number", ylabel="Response probability",
            xgridvisible=false, ygridvisible=false)
        ax2 = Axis(fig[1, 2], xlabel="Stimulus number", ylabel="log(T2/T1) habituation ratio",
            xgridvisible=false, ygridvisible=false)

        # Habituation curves with CI bands
        t1_med = extract_params_median_turing(chain[:trial_1])
        t2_med = extract_params_median_turing(chain[:trial_2])
        t1_samps = extract_params_samples_turing(chain[:trial_1])
        t2_samps = extract_params_samples_turing(chain[:trial_2])
        plot_curves(ax1, t1_med, label="Trial 1", color=:black, axlabel=false)
        plot_curves(ax1, t2_med, label="Trial 2", color=:gray, axlabel=false)
        plot_band(ax1, t1_samps, color=:black)
        plot_band(ax1, t2_samps, color=:gray)
        ylims!(ax1, 0, 1)
        xlims!(ax1, 1, num_stim)

        # Trial diff ratio
        ratio = trial_discrete_derivative_ratio(t1_med, t2_med; step=5)
        lines!(ax2, 1:size(ratio, 1), ratio, color=:black)
        hlines!(ax2, [0.0], color=:gray, linestyle=:dash)
        xlims!(ax2, 1, num_stim)

        return fig, [ax1, ax2]
    end

    """Plot a scatter plot of two single-cell parameters against each other, using extract_params_mean."""
    function plot_sc_param_scatter(ax, param1, param2; norm_param_1=nothing, norm_param_2=nothing, plot_quantiles=false, skip_n=0)
        if norm_param_1 !== nothing
            param1 = param1 ./ norm_param_1
        end
        if norm_param_2 !== nothing
            param2 = param2 ./ norm_param_2
        end
        if plot_quantiles
            # sort param1 and param2 by param1
            sort_inds = sortperm(param1)[skip_n+1:end]
            sort_param1 = param1[sort_inds]
            sort_param2 = param2[sort_inds]
            num_per_quantile = floor(Int, length(sort_param2) / 5)
            #split into 5 quantiles and plot mean with error bars
            quantiles_medians_param1 = median.([sort_param1[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param1)])
            quantiles_medians_param2 = median.([sort_param2[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param2)])
            quantiles_means_param1 = mean.([sort_param1[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param1)])
            quantiles_means_param2 = mean.([sort_param2[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param2)])
            quantiles_low_param1 = quantile.([sort_param1[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param1)], 0.025)
            quantiles_low_param2 = quantile.([sort_param2[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param2)], 0.025)
            quantiles_high_param1 = quantile.([sort_param1[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param1)], 0.975)
            quantiles_high_param2 = quantile.([sort_param2[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param2)], 0.975)
            quantiles_err_param1 = std.([sort_param1[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param1)])
            quantiles_err_param2 = std.([sort_param2[i:i+num_per_quantile-1] for i in 1:num_per_quantile:length(sort_param2)])
            scatter!(ax, quantiles_medians_param1, quantiles_medians_param2, color=:black)
            errorbars!(ax, quantiles_medians_param1, quantiles_medians_param2, quantiles_medians_param1 .- quantiles_low_param1, quantiles_high_param1 .- quantiles_medians_param1, direction=:x, color=:black)
            errorbars!(ax, quantiles_medians_param1, quantiles_medians_param2, quantiles_medians_param2 .- quantiles_low_param2, quantiles_high_param2 .- quantiles_medians_param2, direction=:y, color=:black)
            # scatter!(ax, quantile_means_param1, quantile_means_param2, color=:black)
            # errorbars!(ax, quantile_means_param1, quantile_means_param2, quantiles_err_param1, direction=:x, color=:black)
            # errorbars!(ax, quantile_means_param1, quantile_means_param2, quantiles_err_param2, direction=:y, color=:black)
        else
            scatter!(ax, param1, param2, color=:black)# color=1:num_cells, colormap=:viridis)
        end
        # ax.title = "$(param1) vs $(param2)"
    end


    """Plot comprehensive data analysis with inferred parameters and curves."""
    function plot_param_hist(chain, param, ax; xlims=nothing, ylims=nothing, xlabel=nothing, ylabel=nothing)
        if param in [:a, :b_1_0, :c_1_0, :d_1_0, :b_m_0, :c_2_0, :d_2_0, :b_1_conc, :b_m_var, :c_1_conc, :c_2_conc, :d_1_conc, :d_2_conc]
            vals = extract_params_all_turing(chain)[:pop][param]
        else
            vals = extract_params_median_turing(chain)[:cell][param]
        end
        hist!(ax, vals, color=:black, normalization=:probability)
        ylims = Dict(:a => (0, 0.3), :b_1_0 => (0, 0.3), :c_1_0 => (0, 0.3), :d_1_0 => (0, 0.3), :b_m_0 => (0, 0.3), :c_2_0 => (0, 0.3), :d_2_0 => (0, 0.3), :b_1_conc => (0, 0.3), :b_m_var => (0, 0.3), :c_1_conc => (0, 0.3), :c_2_conc => (0, 0.3), :d_1_conc => (0, 0.3), :d_2_conc => (0, 0.3))
        xlims!(ax, xlims[param])
        ylims!(ax, ylims[param])
        ax.xlabel = xlabels[param]
        ax.ylabel = "Posterior probability"
    end

    """Make ISI/ITI axis grid."""
    function make_isi_iti_axis_grid(isi_vals, iti_vals; vert_label="", horz_label="", fig_size=(1200, 900), label_pos=:bottom)
        fig = Figure(size=fig_size)
        axs = [Axis(fig[(i-1)*3+2:(i-1)*3+4, (j-1)*3+2:(j-1)*3+4]) for i in isi_vals, j in 1:length(iti_vals)]
        for isi in isi_vals
            Label(fig[(isi-1)*3+2:(isi-1)*3+4, 1], "ISI = $(isi) min", rotation=pi / 2, tellwidth=false, tellheight=false, fontsize=24)
        end
        for (j, iti) in enumerate(iti_vals)
            Label(fig[1, (j-1)*3+2:(j-1)*3+4], "ITI = $(iti) hr", tellheight=false, tellwidth=false, fontsize=24)
        end
        [ax.xgridvisible = false for ax in axs]
        [ax.ygridvisible = false for ax in axs]
        Label(fig[2:length(isi_vals)*3+1, 0], vert_label, rotation=pi / 2, tellwidth=false, tellheight=false, fontsize=24)
        if label_pos == :bottom
            Label(fig[length(isi_vals)*3+2, 2:length(iti_vals)*3+1], horz_label, tellwidth=false, tellheight=false, fontsize=24)
        else
            Label(fig[1, 2:length(iti_vals)*3+1], horz_label, tellwidth=false, tellheight=false, fontsize=24)
        end
        return fig, axs
    end

    """Plot all single-cell curves across different ISI/ITI conditions."""
    function plot_all_sc_curves(chains; isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5])
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="Inferred single-cell\nresponse probability", horz_label="Stimulus number")
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                plot_curves(axs[isi, j], extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)]), label=false)
                # band plot of 95% CI
                params = extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)])
                samps = extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)])
                plot_band(axs[isi, j], params)
                xlims!(axs[isi, j], 1, 60)
                ylims!(axs[isi, j], 0, 1)
            end
        end
        # Add a single y label in column 0, centered vertically across all rows
        return fig, axs
    end

    function plot_all_sc_curves_split(chains; responses=nothing, isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5])
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="Inferred single-cell\nresponse probability", horz_label="Stimulus number")
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                key = get_key(isi * 60, iti * 3600)
                for (trial_key, clr) in [(:trial_1, :black), (:trial_2, :gray)]
                    meds = extract_params_median_turing(chains[key][trial_key])
                    plot_curves(axs[isi, j], meds, axlabel=false, color=clr)
                    samps = extract_params_samples_turing(chains[key][trial_key])
                    plot_band(axs[isi, j], samps, color=clr)
                end
                xlims!(axs[isi, j], 1, 60)
                ylims!(axs[isi, j], 0, 1)
            end
        end
        return fig, axs
    end

    """Plot inferred population mean (average of per-cell median Hill curves) for split single_trial chains.
    Overlays raw data population mean as blue scatter markers."""
    function plot_all_sc_population_mean_split(chains, responses; isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5])
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="Population mean\nresponse probability", horz_label="Stimulus number")
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                key = get_key(isi * 60, iti * 3600)
                data = responses[key]["control_data"]
                num_stim = size(data, 1) ÷ 2
                for (trial_key, rows, clr) in [(:trial_1, 1:num_stim, :black), (:trial_2, num_stim+1:2*num_stim, :gray)]
                    ch = chains[key][trial_key]
                    # Median population params
                    a_1_0_med = median(vec(Array(ch[:a_1_0])))
                    a_1_std_med = median(vec(Array(ch[:a_1_stdev])))
                    # Median cell-level params (median across samples for each cell)
                    a_1_raw_med = vec(median(Array(group(ch, :a_1_raw)), dims=1))
                    b_1_med = vec(median(Array(group(ch, :b_1)), dims=1))
                    c_1_med = vec(median(Array(group(ch, :c_1)), dims=1))
                    d_1_med = vec(median(Array(group(ch, :d_1)), dims=1))
                    a_cells = exp.(a_1_0_med .+ a_1_std_med .* a_1_raw_med)
                    n_cells = length(b_1_med)
                    # Per-cell Hill curves, then average → median population mean
                    cell_curves = hcat([hill.(1:num_stim, a_cells[i], b_1_med[i], c_1_med[i], d_1_med[i]) for i in 1:n_cells]...)
                    pop_mean = vec(mean(cell_curves, dims=2))
                    lines!(axs[isi, j], 1:num_stim, pop_mean, color=clr, label=(trial_key == :trial_1 ? "Inferred T1" : "Inferred T2"))
                    # 5-95% credible interval on population mean across all samples
                    a_1_0_s = vec(Array(ch[:a_1_0]))
                    a_1_std_s = vec(Array(ch[:a_1_stdev]))
                    a_1_raw_arr = Array(group(ch, :a_1_raw))  # (n_samples_total, n_cells)
                    b_1_arr = Array(group(ch, :b_1))
                    c_1_arr = Array(group(ch, :c_1))
                    d_1_arr = Array(group(ch, :d_1))
                    n_samps = length(a_1_0_s)
                    pop_mean_curves = zeros(num_stim, n_samps)
                    for s in 1:n_samps
                        a_cells_s = exp.(a_1_0_s[s] .+ a_1_std_s[s] .* a_1_raw_arr[s, :])
                        cell_curves_s = hcat([hill.(1:num_stim, a_cells_s[i], b_1_arr[s, i], c_1_arr[s, i], d_1_arr[s, i]) for i in 1:n_cells]...)
                        pop_mean_curves[:, s] = vec(mean(cell_curves_s, dims=2))
                    end
                    lo = [quantile(pop_mean_curves[i, :], 0.05) for i in 1:num_stim]
                    hi = [quantile(pop_mean_curves[i, :], 0.95) for i in 1:num_stim]
                    band!(axs[isi, j], 1:num_stim, lo, hi, color=(clr, 0.15))
                end
                # Overlay raw data population mean
                mean_t1 = vec(mean(data[1:num_stim, :], dims=2))
                mean_t2 = vec(mean(data[num_stim+1:end, :], dims=2))
                scatter!(axs[isi, j], 1:num_stim, mean_t1, color=(:dodgerblue, 0.9), markersize=5, label="Data T1")
                scatter!(axs[isi, j], 1:num_stim, mean_t2, color=(:dodgerblue, 0.4), markersize=5, label="Data T2")
                xlims!(axs[isi, j], 1, num_stim)
                ylims!(axs[isi, j], 0, 1)
            end
        end
        return fig, axs
    end

    """Plot all single-cell derivative curves across different ISI/ITI conditions."""
    function plot_all_sc_curves_derivatives(data_chains; isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5])
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="Learning rate", horz_label="Stimulus number")
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                plot_curves_derivative(axs[isi, j], extract_params_median_turing(data_chains[get_key(isi * 60, iti * 3600)]), label=false)
            end
        end
        return fig, axs
    end

    """Plot all single-cell curves vs their derivatives across different ISI/ITI conditions."""
    function plot_all_sc_curves_vs_derivatives(chains; isi_vals=[1, 2, 3], iti_vals=[1, 2, 3, 5])
        fig, axs = make_isi_iti_axis_grid(isi_vals, iti_vals, vert_label="Learning rate", horz_label="Response probability")
        for isi in isi_vals
            for (j, iti) in enumerate(iti_vals)
                plot_curves_vs_derivatives(axs[isi, j], extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1]), label=false, color=:black)
                plot_curves_vs_derivatives(axs[isi, j], extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)][:trial_2]), label=false, color=:gray)
            end
        end
        xlims!.(axs, 0, 1)
        return fig, axs
    end


    """Plot single curves for trial 1 using population-level averaged parameters, averaging across ITIs for each ISI."""
    function plot_trial1_curves_all_conditions(chains, isi_list, iti, ax; time_axis=false, leg=true)
        # Plot curves for each ISI, averaging across ITIs
        for isi in isi_list
            a_0_avg = extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][:a_1_0]
            b_0_avg = extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][:b_1_0]
            c_0_avg = extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][:c_1_0]
            d_0_avg = extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][:d_1_0]

            # Only change x values for plotting if time_axis is true
            if time_axis
                # x is time in minutes, step is # isi
                x = 0:isi:59*isi
                stim_numbers = (x ./ isi) .+ 1
                y = hill.(stim_numbers, a_0_avg, b_0_avg, c_0_avg, d_0_avg)
                lines!(ax, x, y, color=isi_colors[isi], linewidth=1.5,
                    label="$isi min ISI")
                plot_band(ax, extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1]), color=isi_colors[isi], plot_x=(isi .* (0:0.1:59)))
            else
                x = 1:0.1:60
                y = hill.(x, a_0_avg, b_0_avg, c_0_avg, d_0_avg)
                lines!(ax, x, y, color=isi_colors[isi], linewidth=1.5,
                    label="$isi min ISI")
                plot_band(ax, extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1]), color=isi_colors[isi])
            end
        end

        # Add axis labels and title
        ax.xlabel = time_axis ? "Time (min)" : "Stimulus number"
        ax.ylabel = "Response probability"
        ax.xgridvisible = false
        ax.ygridvisible = false
        # ax.title = "Trial 1 Curves - Population Level Parameters (Averaged Across ITIs)"

        # Set axis limits
        if time_axis
            max_time = maximum(isi_list) * 60
            xlims!(ax, 0, max_time)
        else
            xlims!(ax, 0, 60)
        end
        ylims!(ax, 0, 1)

        # Add legend if leg is true
        if leg
            axislegend(ax, position=:rt, fontsize=10)
        end

        return ax
    end

    function plot_param_by_isi(chains, ax; isi_list=[1, 2, 3], iti=1, param=:c_1)
        param_medians = []
        ci_upper = []
        ci_lower = []
        for isi in isi_list
            key = get_key(isi * 60, iti * 3600)
            if param == :recovery
                params_all_1 = extract_params_median_turing(chains[key][:trial_1])
                params_all_2 = extract_params_median_turing(chains[key][:trial_2])
                params_all_vals_1 = extract_params_samples_turing(chains[key][:trial_1])
                params_all_vals_2 = extract_params_samples_turing(chains[key][:trial_2])
                push!(param_medians, params_all_2[:pop][:c_1_0] ./ params_all_1[:pop][:c_1_0])
                cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> p2[:pop][:c_1_0] ./ p1[:pop][:c_1_0]; prob=0.95)
            elseif param == :hill
                push!(param_medians, log.(extract_params_median_turing(chains[key][:trial_1])[:pop][:a_1_0]))
                cis = calculate_statistic_ci(chains[key][:trial_1], (p1) -> log.(p1[:pop][:a_1_0]); prob=0.95)
            else
                push!(param_medians, extract_params_median_turing(chains[key][:trial_1])[:pop][param])
                cis = calculate_statistic_ci(chains[key][:trial_1], (p1) -> p1[:pop][param]; prob=0.95)
            end
            push!(ci_lower, cis[2] - cis[1])
            push!(ci_upper, cis[3] - cis[2])
        end
        lines!(ax, isi_list, param_medians, color=:black, linewidth=2, label="$iti min ITI")
        errorbars!(ax, isi_list, param_medians, ci_lower, ci_upper, color=:black, linewidth=2)
        return ax
    end

    function plot_param_hist_by_isi(chains, ax; isi_list=[1, 2, 3], iti=1, param=:c_1_0)
        for isi in isi_list
            if param == :recovery
                samps = extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)][:trial_2])[:pop][:c_1_0] ./ extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][:c_1_0]
                hist!(ax, samps, color=isi_colors[isi], bins=0:0.1:1, normalization=:probability, transparency=true, alpha=0.3)
            elseif param == :b_1_0
                samps = extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][:b_1_0]
                hist!(ax, samps, color=isi_colors[isi], bins=0:2:30, normalization=:probability, transparency=true, alpha=0.3)
            else
                samps = extract_params_samples_turing(chains[get_key(isi * 60, iti * 3600)][:trial_1])[:pop][param]
                hist!(ax, samps, color=isi_colors[isi], normalization=:probability, transparency=true, alpha=0.3)
            end
        end
        return ax
    end

    function PLOT_INIT_by_iti(
        chains, axs; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5],
        error_bars=true
    )
        # Create three axes: c_m_0, c_0, c_multiplier
        params_all = extract_params_all_turing(chains)
        for isi in isi_list
            init_1 = [params_all[get_key(isi * 60, iti * 3600)][:pop][:c_1_0] for iti in iti_list]
            init_2 = [params_all[get_key(isi * 60, iti * 3600)][:pop][:c_2_0] for iti in iti_list]
            recovery = init_2 ./ init_1
            lines!(axs[1], iti_list, init_1, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            lines!(axs[2], iti_list, init_2, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            lines!(axs[3], iti_list, recovery, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            # errorbars!(ax, iti_list, recovery_values_isi, sqrt.(recovery_vars_isi), color=isi_colors[isi], linewidth=2)
        end
        xlims!(axs[1], 0, 6)
        ylims!(axs[1], 0.4, 1.2)
        # axislegend(ax_c_m, position=:rt, fontsize=10)
        ax.xlabel = "Rest duration (hrs)"
        ax.ylabel = "Recovery"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end


    function plot_recovery_by_iti(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5],
        error_bars=true
    )
        # Create three axes: c_m_0, c_0, c_multiplier
        for isi in isi_list
            recovery_values_isi = Float64[]
            recovery_ci_lower = Float64[]
            recovery_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                params_all_1 = extract_params_median_turing(chains[key][:trial_1])
                params_all_2 = extract_params_median_turing(chains[key][:trial_2])
                params_all_vals_1 = extract_params_samples_turing(chains[key][:trial_1])
                params_all_vals_2 = extract_params_samples_turing(chains[key][:trial_2])
                push!(recovery_values_isi, params_all_2[:pop][:c_1_0] ./ params_all_1[:pop][:c_1_0])
                cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> p2[:pop][:c_1_0] ./ p1[:pop][:c_1_0]; prob=0.95)
                push!(recovery_ci_lower, cis[2] - cis[1])
                push!(recovery_ci_upper, cis[3] - cis[2])
                # push!(recovery_values_isi, (params_all[:pop][:c_2_0]) ./ (params_all[:pop][:c_1_0] - params_all[:pop][:d_1_0]))
                # push!(recovery_values_isi, params_all[:pop][:c_2_0])
                # recovery = params_all[:cell][:c_2] ./ params_all[:cell][:c_1]
                # push!(recovery_values_isi, median(recovery))
                # push!(recovery_vars_isi, var(recovery))
            end
            lines!(ax, iti_list .+ 0.1 * (isi - 1), recovery_values_isi, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            errorbars!(ax, iti_list .+ 0.1 * (isi - 1), recovery_values_isi, recovery_ci_lower, recovery_ci_upper, color=isi_colors[isi], linewidth=2)
        end
        hlines!(ax, [1.0], color=:gray, linestyle=:dash)
        xlims!(ax, 0, 6)
        ylims!(ax, 0, 1.5)
        axislegend(ax, position=:rb, fontsize=10)
        ax.xlabel = "Rest duration (hrs)"
        ax.ylabel = "Initial response ratio"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """Plot population trial 2 parameters b_m_0 and c_m_0 for each ITI with error bars, separate line for each ISI."""
    function plot_population_t2_params_by_iti(
        chains, isi_list, iti_list;
        fig=Figure(size=(1500, 400)),
        error_bars=true
    )
        # Create three axes: b_m_0, c_m_0, c_multiplier
        ax_b_m = Axis(fig[1, 1])
        ax_c_m = Axis(fig[1, 2])
        ax_c_mult = Axis(fig[1, 3])

        # Plot for each ISI
        for isi in isi_list
            b_m_values = [extract_params_median_turing(chains[get_key(isi * 60, iti * 3600)])[:pop][:b_m] for iti in iti_list]
            b_m_var = [extract_params_quantiles_turing(chains[get_key(isi * 60, iti * 3600)])[:pop][:b_m] for iti in iti_list]
            # Plot b_m_0 parameter with error bars
            errorbars!(ax_b_m, iti_list, b_m_values, b_m_vars,
                color=isi_colors[isi], linewidth=2)
            lines!(ax_b_m, iti_list, b_m_values, color=isi_colors[isi],
                linewidth=2, label="$isi min ISI")

            # Plot c_m_0 parameter with error bars
            if error_bars
                errorbars!(ax_c_m, iti_list, c_m_values, sqrt.(c_m_vars),
                    color=isi_colors[isi], linewidth=2)
            end
            lines!(ax_c_m, iti_list, c_m_values, color=isi_colors[isi],
                linewidth=2, label="$isi min ISI")

            # Plot c_multiplier (c_m_0 / c_0) with error bars
            if error_bars
                errorbars!(ax_c_mult, iti_list, c_mult_values, sqrt.(c_mult_vars),
                    color=isi_colors[isi], linewidth=2)
            end
            lines!(ax_c_mult, iti_list, c_mult_values, color=isi_colors[isi],
                linewidth=2, label="$isi min ISI")
        end

        # Set axis labels and titles
        ax_b_m.xlabel = "ITI (hours)"
        ax_b_m.ylabel = "N50 multiplier (potentiation)"
        ax_b_m.title = "Trial 2 N50 Multiplier by ITI/ISI"

        ax_c_m.xlabel = "ITI (hours)"
        ax_c_m.ylabel = "Trial 2 initial response"
        ax_c_m.title = "Trial 2 initial response by ITI/ISI"

        ax_c_mult.xlabel = "ITI (hours)"
        ax_c_mult.ylabel = "T2/T1 initial response (cₘ₀/c₀) → more recovery"
        ax_c_mult.title = "Recovery (cₘ₀/c₀) by ITI/ISI"

        # Set axis limits
        xlims!(ax_b_m, 0, maximum(iti_list) + 1)
        xlims!(ax_c_m, 0, maximum(iti_list) + 1)
        xlims!(ax_c_mult, 0, maximum(iti_list) + 1)

        ylims!(ax_b_m, 0, 2)
        ylims!(ax_c_m, -0.1, 1.1)
        ylims!(ax_c_mult, -0.1, 1.1)

        # Add legends
        axislegend(ax_b_m, position=:rt, fontsize=10)
        axislegend(ax_c_m, position=:rt, fontsize=10)
        axislegend(ax_c_mult, position=:rb, fontsize=10)

        # Add overall title
        title_text = "Trial 2 Population Parameters by ITI/ISI"
        Label(fig[0, :], title_text, tellwidth=false, fontsize=16)

        return fig, [ax_b_m, ax_c_m, ax_c_mult]
    end


    function plot_b_m_by_iti(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5],
        error_bars=true
    )
        # Create three axes: c_m_0, c_0, c_multiplier
        for isi in isi_list
            b_m_values_isi = Float64[]
            b_m_ci_lower = Float64[]
            b_m_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> log.(p1[:pop][:b_1_0] ./ p2[:pop][:b_1_0]); prob=0.95)
                # cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> (p2[:pop][:b_1_0] ./ p1[:pop][:b_1_0]); prob=0.95)
                push!(b_m_values_isi, cis[2])
                push!(b_m_ci_lower, cis[2] - cis[1])
                push!(b_m_ci_upper, cis[3] - cis[2])
            end
            lines!(ax, iti_list .+ 0.1 * (isi - 1), b_m_values_isi, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            if error_bars
                errorbars!(ax, iti_list .+ 0.1 * (isi - 1), b_m_values_isi, b_m_ci_lower, b_m_ci_upper, color=isi_colors[isi], linewidth=2)
            end
        end
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        xlims!(ax, 0, 6)
        ylims!(ax, 0, 1.5)
        # axislegend(ax_c_m, position=:rt, fontsize=10)
        ax.xlabel = "Rest duration (hrs)"
        ax.ylabel = "log(N50 ratio)"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """Plot the trial diff ratio (log T2/T1 decrease) by ITI for each ISI condition."""
    function plot_trial_diff_ratio_by_iti(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5],
        stim=5, error_bars=true
    )
        for isi in isi_list
            ratio_values = Float64[]
            ratio_ci_lower = Float64[]
            ratio_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                cis = calculate_statistic_ci_vec(
                    chains[key][:trial_1], chains[key][:trial_2],
                    (p1, p2) -> trial_diff_ratio(p1, p2; stim=stim);
                    prob=0.95
                )
                push!(ratio_values, cis[2])
                push!(ratio_ci_lower, cis[2] - cis[1])
                push!(ratio_ci_upper, cis[3] - cis[2])
            end
            lines!(ax, iti_list .+ 0.1 * (isi - 1), ratio_values, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            if error_bars
                errorbars!(ax, iti_list .+ 0.1 * (isi - 1), ratio_values, ratio_ci_lower, ratio_ci_upper, color=isi_colors[isi], linewidth=2)
            end
        end
        xlims!(ax, 0, 6)
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        ax.xlabel = "Rest duration (hrs)"
        ax.ylabel = "Initial rate ratio"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """Plot the difference in hill_derivative_y (learning rate at a given response level)
    between T2 and T1 by ITI for each ISI condition."""
    function plot_derivative_y_diff_by_iti(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5],
        y_val=0:0.01:1, error_bars=true
    )
        for isi in isi_list
            diff_values = Float64[]
            diff_ci_lower = Float64[]
            diff_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                cis = calculate_statistic_ci_vec(
                    chains[key][:trial_1], chains[key][:trial_2],
                    (p1, p2) -> sum(abs.(hill_derivative_y.(max(p2[:pop][:d_1_0], p1[:pop][:d_1_0]):0.01:min(p2[:pop][:c_1_0], p1[:pop][:c_1_0]), p2[:pop][:a_1_0], p2[:pop][:b_1_0], p2[:pop][:c_1_0], p2[:pop][:d_1_0]) .-
                                         hill_derivative_y.(max(p2[:pop][:d_1_0], p1[:pop][:d_1_0]):0.01:min(p2[:pop][:c_1_0], p1[:pop][:c_1_0]), p1[:pop][:a_1_0], p1[:pop][:b_1_0], p1[:pop][:c_1_0], p1[:pop][:d_1_0])));
                    prob=0.95
                )
                push!(diff_values, cis[2])
                push!(diff_ci_lower, cis[2] - cis[1])
                push!(diff_ci_upper, cis[3] - cis[2])
            end
            lines!(ax, iti_list .+ 0.1 * (isi - 1), diff_values, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            if error_bars
                errorbars!(ax, iti_list .+ 0.1 * (isi - 1), diff_values, diff_ci_lower, diff_ci_upper, color=isi_colors[isi], linewidth=2)
            end
            hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        end
        xlims!(ax, 0, 6)
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        ax.xlabel = "Rest duration (hrs)"
        ax.ylabel = "Phase plot diff"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    function plot_recovery_vs_b_m_0(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5], plot_ci=false
    )
        for isi in isi_list
            recovery_values = Float64[]
            recovery_ci_lower = Float64[]
            recovery_ci_upper = Float64[]
            b_m_values = Float64[]
            b_m_ci_lower = Float64[]
            b_m_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                rec_cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> log.(p1[:pop][:c_1_0] ./ p2[:pop][:c_1_0]); prob=0.95)
                push!(recovery_values, rec_cis[2])
                push!(recovery_ci_lower, rec_cis[2] - rec_cis[1])
                push!(recovery_ci_upper, rec_cis[3] - rec_cis[2])

                # b_m_0 and its CI
                b_m_cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> log.(p1[:pop][:b_1_0] ./ p2[:pop][:b_1_0]); prob=0.95)
                push!(b_m_values, b_m_cis[2])
                push!(b_m_ci_lower, b_m_cis[2] - b_m_cis[1])
                push!(b_m_ci_upper, b_m_cis[3] - b_m_cis[2])
            end
            if plot_ci
                # Plot crosshairs (error bars in both directions)
                errorbars!(ax, recovery_values, b_m_values, recovery_ci_lower, recovery_ci_upper,
                    direction=:x, color=isi_colors[isi], linewidth=1.5)
                errorbars!(ax, recovery_values, b_m_values, b_m_ci_lower, b_m_ci_upper,
                    direction=:y, color=isi_colors[isi], linewidth=1.5)
            end
            scatter!(ax, recovery_values, b_m_values, color=isi_colors[isi], markersize=10, label="$isi min ISI")
        end
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        vlines!(ax, [0.0], color=:gray, linestyle=:dash)
        # ablines!(ax, 0, 1, color=:gray, linestyle=:dot)
        ax.xlabel = "log(Initial response ratio)"
        ax.ylabel = "log(N50 ratio)"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """Plot recovery (c₂₀/c₁₀) vs trial diff ratio for each ISI across ITIs."""
    function plot_recovery_vs_trial_diff_ratio(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5], stim=5, plot_ci=false
    )
        for isi in isi_list
            recovery_values = Float64[]
            recovery_ci_lower = Float64[]
            recovery_ci_upper = Float64[]
            ratio_values = Float64[]
            ratio_ci_lower = Float64[]
            ratio_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                # Recovery and its CI
                rec_cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> log.(p1[:pop][:c_1_0] ./ p2[:pop][:c_1_0]); prob=0.95)
                push!(recovery_values, rec_cis[2])
                push!(recovery_ci_lower, rec_cis[2] - rec_cis[1])
                push!(recovery_ci_upper, rec_cis[3] - rec_cis[2])

                # Trial diff ratio and its CI
                ratio_cis = calculate_statistic_ci_vec(
                    chains[key][:trial_1], chains[key][:trial_2],
                    (p1, p2) -> trial_diff_ratio(p1, p2; stim=stim);
                    prob=0.95
                )
                push!(ratio_values, ratio_cis[2])
                push!(ratio_ci_lower, ratio_cis[2] - ratio_cis[1])
                push!(ratio_ci_upper, ratio_cis[3] - ratio_cis[2])
            end
            if plot_ci
                errorbars!(ax, recovery_values, ratio_values, recovery_ci_lower, recovery_ci_upper,
                    direction=:x, color=isi_colors[isi], linewidth=1.5)
                errorbars!(ax, recovery_values, ratio_values, ratio_ci_lower, ratio_ci_upper,
                    direction=:y, color=isi_colors[isi], linewidth=1.5)
            end
            scatter!(ax, recovery_values, ratio_values, color=isi_colors[isi], markersize=10, label="$isi min ISI")
        end
        vlines!(ax, [0.0], color=:gray, linestyle=:dash)
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        # ablines!(ax, 0, 1, color=:gray, linestyle=:dot)
        ax.xlabel = "log(Recovery)"
        ax.ylabel = "Initial rate ratio"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """Plot recovery (c₂₀/c₁₀) vs derivative diff (sum |dy₂ - dy₁|) for each ISI across ITIs."""
    function plot_recovery_vs_derivative_diff(
        chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5], plot_ci=false
    )
        for isi in isi_list
            recovery_values = Float64[]
            recovery_ci_lower = Float64[]
            recovery_ci_upper = Float64[]
            diff_values = Float64[]
            diff_ci_lower = Float64[]
            diff_ci_upper = Float64[]
            for iti in iti_list
                key = get_key(isi * 60, iti * 3600)
                # Recovery and its CI
                rec_cis = calculate_statistic_ci(chains[key][:trial_1], chains[key][:trial_2], (p1, p2) -> log.(p1[:pop][:c_1_0] ./ p2[:pop][:c_1_0]); prob=0.95)
                push!(recovery_values, rec_cis[2])
                push!(recovery_ci_lower, rec_cis[2] - rec_cis[1])
                push!(recovery_ci_upper, rec_cis[3] - rec_cis[2])

                # Derivative diff and its CI
                diff_cis = calculate_statistic_ci_vec(
                    chains[key][:trial_1], chains[key][:trial_2],
                    (p1, p2) -> sum(abs.(hill_derivative_y.(max(p2[:pop][:d_1_0], p1[:pop][:d_1_0]):0.01:min(p2[:pop][:c_1_0], p1[:pop][:c_1_0]), p2[:pop][:a_1_0], p2[:pop][:b_1_0], p2[:pop][:c_1_0], p2[:pop][:d_1_0]) .-
                                         hill_derivative_y.(max(p2[:pop][:d_1_0], p1[:pop][:d_1_0]):0.01:min(p2[:pop][:c_1_0], p1[:pop][:c_1_0]), p1[:pop][:a_1_0], p1[:pop][:b_1_0], p1[:pop][:c_1_0], p1[:pop][:d_1_0])));
                    prob=0.95
                )
                push!(diff_values, diff_cis[2])
                push!(diff_ci_lower, diff_cis[2] - diff_cis[1])
                push!(diff_ci_upper, diff_cis[3] - diff_cis[2])
            end
            if plot_ci
                errorbars!(ax, recovery_values, diff_values, recovery_ci_lower, recovery_ci_upper,
                    direction=:x, color=isi_colors[isi], linewidth=1.5)
                errorbars!(ax, recovery_values, diff_values, diff_ci_lower, diff_ci_upper,
                    direction=:y, color=isi_colors[isi], linewidth=1.5)
            end
            scatter!(ax, recovery_values, diff_values, color=isi_colors[isi], markersize=10, label="$isi min ISI")
        end
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
        vlines!(ax, [0.0], color=:gray, linestyle=:dash)
        ax.xlabel = "log(Recovery)"
        ax.ylabel = "Phase plot diff"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """
    Compute normalized phase curves (learning rate vs response probability) for a given ISI/ITI.
    Returns y_vals, dy1_norm (trial 1), dy2_norm (trial 2).
    """
    function get_normalized_phase_curves(chains, isi::Int, iti::Int; y_range=0.0:0.01:1.0)
        key = get_key(isi * 60, iti * 3600)
        params = extract_params_median_turing(chains[key])

        a = params[:pop][:a]
        b_1 = params[:pop][:b_1_0]
        c_1 = params[:pop][:c_1_0]
        d_1 = params[:pop][:d_1_0]
        b_m = params[:pop][:b_m_0]
        c_2 = params[:pop][:c_2_0]
        d_2 = params[:pop][:d_2_0]

        # Define y range within valid bounds for both trials
        min_y = max(minimum(y_range), maximum([d_1, d_2]) + 0.01)
        max_y = min(maximum(y_range), minimum([c_1, c_2]) - 0.01)
        y_vals = collect(min_y:0.01:max_y)

        # Compute derivatives for both trials
        dy1 = abs.(hill_derivative_y.(y_vals, a, b_1, c_1, d_1))
        dy2 = abs.(hill_derivative_y.(y_vals, a, b_1 * b_m, c_2, d_2))

        # Normalize by the max of both curves
        max_val = max(maximum(dy1), maximum(dy2))
        dy1_norm = dy1 ./ max_val
        dy2_norm = dy2 ./ max_val

        return y_vals, dy1_norm, dy2_norm
    end

    """
    Plot normalized phase difference between trial 1 and trial 2.
    Uses get_normalized_phase_curves and plots the difference (trial2 - trial1).
    """
    function plot_normalized_phase_diff(chains, ax, isis::Vector{Int}, itis::Vector{Int}; y_range=0.0:0.01:1.0)
        for isi in isis
            for iti in itis
                y_vals, dy1_norm, dy2_norm = get_normalized_phase_curves(chains, isi, iti; y_range=y_range)
                diff = dy2_norm .- dy1_norm
                lines!(ax, y_vals, diff, color=isi_colors[isi], alpha=iti_alphas[iti], linewidth=2, label="$iti hr ITI")
            end
        end
        ax.xlabel = "Response probability"
        ax.ylabel = "Learning rate difference"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

    """
    Plot the peak of the normalized phase difference for each ISI across all ITIs.
    """
    function plot_peak_phase_diff_by_iti(chains, ax; isi_list=[1, 2, 3], iti_list=[1, 2, 3, 5], y_range=0.0:0.01:1.0)
        for isi in isi_list
            peak_values = Float64[]
            for iti in iti_list
                y_vals, dy1_norm, dy2_norm = get_normalized_phase_curves(chains, isi, iti; y_range=y_range)
                diff = dy2_norm .- dy1_norm
                push!(peak_values, maximum(diff))
            end
            lines!(ax, iti_list, peak_values, color=isi_colors[isi], linewidth=2, label="$isi min ISI")
            # scatter!(ax, iti_list, peak_values, color=isi_colors[isi], markersize=10)
        end
        ax.xlabel = "Rest duration (hrs)"
        ax.ylabel = "Peak learning rate difference"
        ax.xgridvisible = false
        ax.ygridvisible = false
    end

end

# =============================================================================
# Debugging
# =============================================================================
begin
    function print_inferred_params(inferred_chains)
        param_df = DataFrame(ISI=Int[], ITI=Int[], a=Float64[], b_1_0=Float64[], c_1_0=Float64[], d_1_0=Float64[], b_m_0=Float64[], c_2_0=Float64[], d_2_0=Float64[], b_1_conc=Float64[], b_m_var=Float64[], c_1_conc=Float64[], c_2_conc=Float64[], d_1_conc=Float64[], d_2_conc=Float64[], recovery=Float64[])
        for isi in 1:3
            println("ISI: ", isi * 60)
            for iti in [1, 2, 3, 5]
                println("ITI: ", iti * 3600)
                key = get_key(isi * 60, iti * 3600)
                params_mean = extract_params_mean_turing(inferred_chains[key])[:pop]
                params_cell = extract_params_mean_turing(inferred_chains[key])[:cell]
                push!(param_df, [isi * 60, iti * 3600, params_mean[:a], params_mean[:b_1_0], params_mean[:c_1_0], params_mean[:d_1_0], params_mean[:b_m_0], params_mean[:c_2_0], params_mean[:d_2_0], params_mean[:b_1_conc], params_mean[:b_m_var], params_mean[:c_1_conc], params_mean[:c_2_conc], params_mean[:d_1_conc], params_mean[:d_2_conc], params_mean[:c_2_0] ./ params_mean[:c_1_0]])
            end
        end
        println(param_df)
        return param_df
    end


    """
    Interactive parameter explorer with ISI/ITI dropdowns, population histograms, 
    scatter plots of single-cell means, and a slider for single-cell posterior histograms.
    """
    function chain_visualization(chains)
        f = Figure(size=(1200, 900))
        all_params = extract_params_all_turing(chains)

        # =========================================
        # Controls
        # =========================================
        isi_opts = [1, 2, 3]
        iti_opts = [1, 2, 3, 5]

        # Menus
        Label(f[1, 1], "ISI:", halign=:right, tellwidth=false)
        isi_menu = Menu(f[1, 2], options=zip(["1 min", "2 min", "3 min"], isi_opts), default="1 min", tellwidth=false)
        Label(f[1, 3], "ITI:", halign=:right, tellwidth=false)
        iti_menu = Menu(f[1, 4], options=zip(["1 hr", "2 hr", "3 hr", "5 hr"], iti_opts), default="1 hr", tellwidth=false)

        # Slider for Cell Index
        Label(f[5, 1], "Cell Index:", halign=:right, tellwidth=false)
        cell_slider = Slider(f[5, :], range=1:100, startvalue=1, tellwidth=false)
        cell_label = Label(f[5, 7], @lift("Cell: $(lpad($(cell_slider.value), 2, "0"))"))

        # =========================================
        # Layout & Axes
        # =========================================

        # Row 2: Population Parameter Histograms
        pop_params = [:a, :b_1_0, :c_1_0, :d_1_0, :b_m_0, :c_2_0]
        xlims = [(0, 20), (0, 30), (0, 1), (0, 1), (0, 2), (0, 1)]
        pop_labels = ["a (Hill coeff)", "b_1_0 (N50)", "c_1_0 (Init)", "d_1_0 (Asym)", "b_m_0 (N50 ratio)", "c_2_0 (T2 Init)"]
        pop_axes = [Axis(f[2, i], title=lab, xlabel=lab) for (i, lab) in enumerate(pop_labels)]
        [xlims!(pop_axes[i], xlims[i]) for i in 1:length(pop_params)]

        # Row 3: Scatter Plots (Single Cell Means)
        # b_1 vs b_m, c_1 vs c_2, d_1 vs d_2
        scatter_pairs = [(:b_1, :b_m), (:c_1, :c_2), (:d_1, :d_2)]
        scatter_labels = ["b_1 vs b_m", "c_1 vs c_2", "d_1 vs d_2"]
        scatter_axes = [Axis(f[3, i*2-1:i*2]) for (i, lab) in enumerate(scatter_labels)]

        # Row 4: Single Cell Parameter Histograms (Posterior for selected cell)
        sc_params = [:b_1, :c_1, :d_1, :b_m, :c_2]
        sc_labels = ["b_1 (N50)", "c_1 (Init)", "d_1 (Asym)", "b_m (Ratio)", "c_2 (T2 Init)"]
        sc_axes = [Axis(f[4, i+1], title=lab, xlabel=lab) for (i, lab) in enumerate(sc_labels)]

        [hist!(pop_axes[i], @lift(all_params[:samples][get_key($(isi_menu.selection) * 60, $(iti_menu.selection) * 3600)][:pop][pop_params[i]]), color=:black) for i in 1:length(pop_params)]

        [scatter!(scatter_axes[i], @lift(all_params[:median][get_key($(isi_menu.selection) * 60, $(iti_menu.selection) * 3600)][:cell][scatter_pairs[i][1]]), @lift(all_params[:median][get_key($(isi_menu.selection) * 60, $(iti_menu.selection) * 3600)][:cell][scatter_pairs[i][2]]), color=:black) for i in 1:length(scatter_pairs)]

        [hist!(sc_axes[i], @lift(all_params[:samples][get_key($(isi_menu.selection) * 60, $(iti_menu.selection) * 3600)][:cell][sc_params[i]][$(cell_slider.value), :]), color=:black) for i in 1:length(sc_params)]
        display(f)
    end

    """
    Interactive visualization of trial 1 single-cell priors.
    Shows 4 parameters (a, b, c, d) × 3 columns (population location prior, population spread prior,
    reactive single-cell prior). Sliders control population parameter values; single-cell PDFs update live.

    Arguments:
    - `model`: :potentiation or :potentiation_new (default)
    - `fig_size`: figure size, default (1800, 1000)
    - `n_points`: PDF evaluation points, default 500
    """
    function plot_priors(; model::Symbol=:potentiation_new, fig_size=(1800, 1000), n_points::Int=500)
        # Mode helpers
        _mean_gamma(α, θ) = α * θ
        _mode_beta(α, β) = (α > 1 && β > 1) ? (α - 1) / (α + β - 2) : 0.5
        _median_lognormal(μ, σ) = exp(μ)

        if model == :potentiation_new
            configs = [
                (label="a₁", loc_name="a_1_0", spread_name="a_1_stdev",
                    loc_prior=LogNormal(0, 4), spread_prior=Gamma(1, 1),
                    loc_init=_median_lognormal(0, 4), spread_init=1.0,
                    loc_range=0.01:0.01:20.0, spread_range=0.1:0.1:15.0,
                    pop_loc_x=(0.0, 20.0), pop_spread_x=(0.0, 5.0),
                    sc_x=(0.0, 20.0),
                    sc_dist=(l, s) -> LogNormal(log(max(l, 1e-6)), max(s, 1e-6))),
                (label="b₁", loc_name="b_1_0", spread_name="b_1_conc",
                    loc_prior=Gamma(2, 10), spread_prior=Gamma(2, 5),
                    loc_init=_mean_gamma(2, 10), spread_init=_mean_gamma(2, 5),
                    loc_range=0.1:0.1:60.0, spread_range=0.1:0.1:30.0,
                    pop_loc_x=(0.0, 60.0), pop_spread_x=(0.0, 30.0),
                    sc_x=(0.0, 40.0),
                    sc_dist=(l, s) -> Gamma(max(s, 0.1), l / max(s, 0.1))),
                (label="c₁", loc_name="c_1_0", spread_name="c_1_conc",
                    loc_prior=Beta(7, 3), spread_prior=Gamma(5, 2),
                    loc_init=_mode_beta(7, 3), spread_init=_mean_gamma(5, 2),
                    loc_range=0.01:0.01:0.99, spread_range=0.1:0.1:30.0,
                    pop_loc_x=(0.0, 1.0), pop_spread_x=(0.0, 30.0),
                    sc_x=(0.0, 1.0),
                    sc_dist=(l, s) -> Beta(max(l * s, 0.1), max((1 - l) * s, 0.1))),
                (label="d₁", loc_name="d_1_0", spread_name="d_1_conc",
                    loc_prior=Beta(3, 7), spread_prior=Gamma(5, 2),
                    loc_init=_mode_beta(3, 7), spread_init=_mean_gamma(5, 2),
                    loc_range=0.01:0.01:0.99, spread_range=0.1:0.1:30.0,
                    pop_loc_x=(0.0, 1.0), pop_spread_x=(0.0, 30.0),
                    sc_x=(0.0, 1.0),
                    sc_dist=(l, s) -> Beta(max(l * s, 0.1), max((1 - l) * s, 0.1))),
            ]
        elseif model == :potentiation
            configs = [
                (label="b₁", loc_name="b_1_0", spread_name="b_1_conc",
                    loc_prior=Gamma(2, 10), spread_prior=Gamma(2, 5),
                    loc_init=_mean_gamma(2, 10), spread_init=_mean_gamma(2, 5),
                    loc_range=0.1:0.1:60.0, spread_range=0.1:0.1:30.0,
                    pop_loc_x=(0.0, 60.0), pop_spread_x=(0.0, 30.0),
                    sc_x=(0.0, 40.0),
                    sc_dist=(l, s) -> Gamma(max(s + 1, 0.1), l / max(s, 0.1))),
                (label="c₁", loc_name="c_1_0", spread_name="c_1_conc",
                    loc_prior=Beta(7, 3), spread_prior=Gamma(10, 2),
                    loc_init=_mode_beta(7, 3), spread_init=_mode_gamma(10, 2),
                    loc_range=0.01:0.01:0.99, spread_range=0.1:0.1:50.0,
                    pop_loc_x=(0.0, 1.0), pop_spread_x=(0.0, 50.0),
                    sc_x=(0.0, 1.0),
                    sc_dist=(l, s) -> Beta(max((l + 1e-6) * s, 0.1), max((1 - l + 1e-6) * s, 0.1))),
                (label="d₁", loc_name="d_1_0", spread_name="d_1_conc",
                    loc_prior=Beta(3, 7), spread_prior=Gamma(10, 2),
                    loc_init=_mode_beta(3, 7), spread_init=_mode_gamma(10, 2),
                    loc_range=0.01:0.01:0.99, spread_range=0.1:0.1:50.0,
                    pop_loc_x=(0.0, 1.0), pop_spread_x=(0.0, 50.0),
                    sc_x=(0.0, 1.0),
                    sc_dist=(l, s) -> Beta(max((l + 1e-6) * s, 0.1), max((1 - l + 1e-6) * s, 0.1))),
            ]
        else
            error("Unknown model: $model. Use :potentiation or :potentiation_new.")
        end

        n_params = length(configs)

        # Build slider specs (alternating location / spread for each parameter)
        slider_specs = NamedTuple[]
        for cfg in configs
            push!(slider_specs, (label=cfg.loc_name, range=cfg.loc_range, startvalue=cfg.loc_init))
            push!(slider_specs, (label=cfg.spread_name, range=cfg.spread_range, startvalue=cfg.spread_init))
        end

        fig = Figure(size=fig_size)
        sg = SliderGrid(fig[1:n_params, 4], slider_specs...)

        for (i, cfg) in enumerate(configs)
            loc_obs = sg.sliders[2i-1].value
            spread_obs = sg.sliders[2i].value

            # Col 1: Population location prior + slider marker
            ax1 = Axis(fig[i, 1], title="$(cfg.loc_name)", ylabel="Density")
            ax1.xgridvisible = false
            ax1.ygridvisible = false
            xs1 = collect(range(cfg.pop_loc_x..., length=n_points))
            ys1 = cdf.(cfg.loc_prior, xs1)
            lines!(ax1, xs1, ys1, color=:black, linewidth=2)
            band!(ax1, xs1, zeros(n_points), ys1, color=(:black, 0.15))
            vlines!(ax1, loc_obs, color=:red, linewidth=2, linestyle=:dash)
            let cfg = cfg
                cdf_label1 = @lift(string("CDF: ", round(cdf(cfg.loc_prior, $loc_obs), digits=2)))
                text!(ax1, cfg.pop_loc_x[2], 0.0, text=cdf_label1, color=:red, fontsize=14, align=(:right, :bottom))
            end

            # Col 2: Population spread prior + slider marker
            ax2 = Axis(fig[i, 2], title="$(cfg.spread_name)", ylabel="Density")
            ax2.xgridvisible = false
            ax2.ygridvisible = false
            xs2 = collect(range(cfg.pop_spread_x..., length=n_points))
            ys2 = cdf.(cfg.spread_prior, xs2)
            lines!(ax2, xs2, ys2, color=:black, linewidth=2)
            band!(ax2, xs2, zeros(n_points), ys2, color=(:black, 0.15))
            vlines!(ax2, spread_obs, color=:red, linewidth=2, linestyle=:dash)
            let cfg = cfg
                cdf_label2 = @lift(string("CDF: ", round(cdf(cfg.spread_prior, $spread_obs), digits=2)))
                text!(ax2, cfg.pop_spread_x[2], 0.0, text=cdf_label2, color=:red, fontsize=14, align=(:right, :bottom))
            end

            # Col 3: Single-cell prior (reactive to sliders)
            ax3 = Axis(fig[i, 3], title="$(cfg.label) (single-cell)", ylabel="Density")
            ax3.xgridvisible = false
            ax3.ygridvisible = false
            xs3 = collect(range(cfg.sc_x..., length=n_points))
            let cfg = cfg, xs3 = xs3
                sc_ys = @lift(cdf.(cfg.sc_dist($loc_obs, $spread_obs), xs3))
                lines!(ax3, xs3, sc_ys, color=:dodgerblue, linewidth=2)
                band!(ax3, xs3, zeros(n_points), sc_ys, color=(:dodgerblue, 0.15))
            end
            ylims!(ax1, 0, 1)
            ylims!(ax2, 0, 1)
            ylims!(ax3, 0, 1)
        end

        Label(fig[0, 1:3], "Trial 1 Priors — $(model)", fontsize=20, tellwidth=false)
        return fig
    end
end
