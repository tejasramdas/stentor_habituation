using JSON, HDF5, Random, Dates, StatsBase, CairoMakie, JLD2, Dates #GLMakie 


# =======================
# Helper functions
# =======================
function empty_figure(; x=(0, 60), y=(0, 1.1), size=(1000, 800))
    fig = Figure(size=size)
    ax = Axis(fig[1, 1])
    xlims!(ax, x...)
    ylims!(ax, y...)
    ax.xgridvisible = false
    ax.ygridvisible = false
    return fig, ax
end

function clean_axis!(ax; xlims=(0, 60), ylims=(0, 1.1))
    empty!(ax)
    ax.xgridvisible = false
    ax.ygridvisible = false
    xlims!(ax, xlims...)
    ylims!(ax, ylims...)
    return ax
end

# =======================
# Data loading functions
# =======================

function get_key(isi, iti)
    return "hab_ISI$(isi)_ITI$(iti)"
end
cond_keys = [get_key(isi * 60, iti * 3600) for isi in 1:3 for iti in [1, 2, 3, 5]]
function get_annot(fold)
    try
        f = read(h5open("$fold_redirect/data/$fold/annotated/$(fold)_contractions.h5", "r")["manual"])
        return f
    catch
        run(`datalad get -d $fold_redirect $fold_redirect/data/$fold/annotated/$(fold)_contractions.h5`)
        f = read(h5open("$fold_redirect/data/$fold/annotated/$(fold)_contractions.h5", "r")["manual"])
        return f
    end
end
function get_locs(fold)
    run(`datalad get -d $fold_redirect $fold_redirect/data/$fold/tiled/$(fold)_tiled_data.h5`)
    f_locs = read(h5open("$fold_redirect/data/$fold/tiled/$(fold)_tiled_data.h5", "r")["cell_locs"])
    return f_locs
end
function get_info(fold)
    return JSON.parsefile("$fold_redirect/data/$fold/trial_1/dat.json")
end
function get_all(folds)
    return [get_annot(f) for f in folds], [get_info(f) for f in folds], folds
end
function get_all_with_locs(folds)
    return [get_annot(f) for f in folds], [get_locs(f) for f in folds], [get_info(f) for f in folds], folds
end
function process_data(all, info, folds, locs)
    locs = [map(x -> (x[1], x[2]), locs[i]) for i in 1:size(locs, 1)]
    all_filt = [all[i][:, findall(x -> !isnan(x), dropdims(sum(all[i], dims=1), dims=1))] for i in 1:size(all, 1)]
    locs_filt = [locs[i][findall(x -> !isnan(x), dropdims(sum(all[i], dims=1), dims=1))] for i in 1:size(locs, 1)]
    nums = [size(a, 2) for a in all]
    stages = map(i -> i["stage_number"], info)
    nums_filt = [size(a, 2) for a in all_filt]
    comb = hcat(all...)
    comb_locs_filt = vcat(locs_filt...)
    filt_nan = findall(x -> !isnan(x), dropdims(sum(comb, dims=1), dims=1))
    comb_filt = comb[:, filt_nan]
    sort_inds = sortperm(vec(sum(comb_filt[1:60, :], dims=1)))
    comb_sort_t1 = comb_filt[:, sort_inds]
    sort_inds_all = sortperm(vec(sum(comb_filt[1:120, :], dims=1)))
    comb_sort_all = comb_filt[:, sort_inds_all]
    rand_inds = shuffle(1:size(filt_nan, 1))
    comb_shuf = comb_filt[:, rand_inds]
    cell_indx = vcat([zeros(nums[i]) .+ i for i in 1:size(all, 1)]...)[filt_nan]
    return Dict(["annots" => all, "locs_filt" => locs_filt, "folders" => folds, "annots_filt" => annots_filt, "nums_filt" => nums_filt, "nums" => nums, "comb" => comb, "comb_locs_filt" => comb_locs_filt, "comb_filt" => comb_filt, "sort_inds" => sort_inds, "comb_sort_t1" => comb_sort_t1, "rand_inds" => rand_inds, "comb_shuf" => comb_shuf, "cell_indx" => cell_indx, "sort_inds_all" => sort_inds_all, "comb_sort_all" => comb_sort_all, "stages" => stages])
end
function process_data_simple(annots, infos, folds)
    annots_filt = [annots[i][:, findall(x -> !isnan(x), dropdims(sum(annots[i], dims=1), dims=1))] for i in eachindex(annots)]
    nums = [size(a, 2) for a in annots]
    stages = map(i -> i["stage_number"], infos)
    nums_filt = [size(a, 2) for a in annots_filt]
    comb = hcat(annots...)
    comb_filt = hcat(annots_filt...)
    sort_inds_all = sortperm(vec(sum(comb_filt[1:120, :], dims=1)))
    comb_sort_all = comb_filt[:, sort_inds_all]
    return Dict(["annots" => annots, "folders" => folds, "annots_filt" => annots_filt, "nums_filt" => nums_filt, "nums" => nums, "comb" => comb, "comb_filt" => comb_filt, "comb_sort_all" => comb_sort_all, "sort_inds_all" => sort_inds_all, "stages" => stages])
end
function get_cond(cond, expts_data)
    process_data_simple(get_all((expts_data[cond]))...)
end

function get_control_data(proc_data; num_responders=70, num_nonresponders=30, seed=42, control=true)
    println("Getting control data with $num_responders responders and $num_nonresponders nonresponders")
    for (i, j) in proc_data
        if control
            responders = findall(j["comb_filt"][1, :] .== 1)
            nonresponders = findall(j["comb_filt"][1, :] .== 0)
            num_responders_this = num_responders
            num_nonresponders_this = num_nonresponders
            println("$i: $(length(responders)) responders, $(length(nonresponders)) nonresponders")
            if size(j["comb_filt"], 2) < num_responders + num_nonresponders
                println("$i: $(size(j["comb_filt"], 2)) cells < $num_responders + $num_nonresponders")
            end
            if length(responders) < num_responders
                println("$i: $(length(responders)) responders < $num_responders. Using $(num_responders - length(responders)) extra nonresponders")
                num_responders_this = length(responders)
                num_nonresponders_this = num_nonresponders + num_responders - length(responders)
            end
            if length(nonresponders) < num_nonresponders
                println("$i: $(length(nonresponders)) nonresponders < $num_nonresponders. Using $(num_nonresponders - length(nonresponders)) extra responders")
                num_nonresponders_this = length(nonresponders)
                num_responders_this = num_responders + num_nonresponders - length(nonresponders)
            end
            Random.seed!(seed)
            random_responders = sort(sample(responders, num_responders_this, replace=false))
            random_nonresponders = sort(sample(nonresponders, num_nonresponders_this, replace=false))
            proc_data[i]["control_data"] = hcat(j["comb_filt"][:, random_responders], j["comb_filt"][:, random_nonresponders])
            proc_data[i]["responders_inds"] = random_responders
            proc_data[i]["nonresponders_inds"] = random_nonresponders
        end
    end
    return proc_data
end

function get_run_info()
    run_info = Dict{String,Any}()
    open("folder_info.txt", "w") do file
        num_folders = length(readdir("$fold_redirect/data"))
        for (j, i) in enumerate(readdir("$fold_redirect/data"))
            this_run_info = Dict{String,Any}()
            try
                this_run_info["date"] = Date(i, "yyyy_mm_dd_HH_MM_SS")
            catch
                this_run_info["error"] = "Not habituation"
                run_info[i] = this_run_info
                continue
            end
            println(j, " of ", num_folders)
            # println(file, "------")
            try
                info = get_info(i)
                this_run_info["volts"] = info["stimulus_voltage"]
                this_run_info["stage"] = info["stage_number"]
                this_run_info["delay"] = info["delay"]
                this_run_info["isi"] = info["isi"]
                this_run_info["iti"] = info["iti"]
                this_run_info["trial_length"] = info["trial_length"]
            catch
                this_run_info["error"] = "No info"
            end
            this_run_info["date"] = Date(i, "yyyy_mm_dd_HH_MM_SS")
            try
                contractions = get_annot(i)
                this_run_info["contractions"] = contractions
                this_run_info["initial_response"] = mean(filter(!isnan, contractions[1, :]))
                this_run_info["num_cells"] = size(contractions, 2)
                this_run_info["num_cells_nan"] = sum(isnan.(vec(mean(contractions, dims=1))))
            catch
                this_run_info["error"] = "No contraction data"
                run_info[i] = this_run_info
                continue
            end
            # println(file, "Volts: $(info["stimulus_voltage"])")
            # println(file, "Stage: $(info["stage_number"])")
            # println(file, "delay: $(info["delay"])")
            # println(file, "isi: $(info["isi"])")
            # println(file, "iti: $(info["iti"])")
            # println(file, "trial length: $(info["trial_length"])")
            # println(file, "--------------------------------")
            this_run_info["error"] = ""
            run_info[i] = this_run_info
        end
    end
    return run_info
end
# =======================
# Data check plotting
# =======================

function plot_folders_by_dates(proc_data)
    f = Figure(size=(1000, 1000))
    ax = Axis(f[1, 1])
    ax2 = Axis(f[2, 1])
    ax3 = Axis(f[3, 1])
    ax4 = Axis(f[4, 1])
    # ax4 = Axis(f[:,2])
    # Extract all folder names and conditions
    rows = []
    dates = []
    ylabels = []
    init_resps = []
    t2_resps = []
    isis = []
    itis = []
    nums = []
    for isi in 1:3
        for iti in [1, 2, 3, 5]
            num_folds = length(proc_data[get_key(isi * 60, iti * 3600)]["folders"])
            push!(itis, fill(iti, num_folds))
            push!(isis, fill(isi, num_folds))
            push!(rows, fit(Histogram, Dates.value.(Date.(proc_data[get_key(isi * 60, iti * 3600)]["folders"], "yyyy_mm_dd_HH_MM_SS") .- Date("2024-08-01")), 0:450).weights)
            push!(init_resps, map(i -> mean(i[1, :]), proc_data[get_key(isi * 60, iti * 3600)]["annots_filt"]))
            push!(t2_resps, map(i -> mean(i[61, :]), proc_data[get_key(isi * 60, iti * 3600)]["annots_filt"]))
            push!(ylabels, "$(lpad(isi*60, 4)) - $(lpad(iti*3600, 5))")
            push!(dates, Dates.value.(Date.(proc_data[get_key(isi * 60, iti * 3600)]["folders"], "yyyy_mm_dd_HH_MM_SS") .- Date("2024-08-01")))
            push!(nums, proc_data[get_key(isi * 60, iti * 3600)]["nums_filt"])
            scatter!(ax2, Date("2024-08-01") .+ Day.(dates[end]), init_resps[end], color=[:red, :black, :blue][isi])#,marker=[:circle,:triangle,:square,:star,:diamond][iti])
            scatter!(ax3, Date("2024-08-01") .+ Day.(dates[end]), t2_resps[end] ./ init_resps[end], color=[:red, :black, :blue][isi])
            scatter!(ax4, Date("2024-08-01") .+ Day.(vcat(dates...)), vcat(nums...), color=[:red, :black, :blue][isi])
            println(ylabels[end])
        end
    end
    heatmap!(ax, hcat(rows...))
    ax.yticks = (1:length(ylabels), ylabels)
    ax.xticks = (0:100:450, string.(Date("2024-08-01") .+ Day.(0:100:450)))
    # scatter!(ax4, vcat(nums...), vcat(init_resps...), color=:blue)
    println(size(vcat(isis...)))
    return f
end

function plot_condition_analysis(proc_data, cond_keys)
    f = Figure(size=(1200, 800))
    ax = Axis(f[1, 1], xlabel="Condition (ISI_ITI)", ylabel="Initial Response",
        title="Initial Response vs Number of Cells by Condition")


    # Store data for each condition
    all_x = Float64[]
    all_y = Float64[]
    all_sizes = Float64[]
    all_labels = String[]
    all_dates = Float64[]

    for (idx, key) in enumerate(cond_keys)
        data = proc_data[key]
        folders = data["folders"]
        all_filt = data["all_filt"]

        # Extract ISI from key for color coding
        isi_val = parse(Int, split(split(key, "ISI")[2], "_")[1]) ÷ 60

        today = Dates.today()
        for (folder_idx, folder) in enumerate(folders)
            if folder_idx <= length(all_filt)
                try
                    folder_date = Date(folder, "yyyy_mm_dd_HH_MM_SS")
                    days_since = Dates.value(today - folder_date)
                    push!(all_dates, days_since)
                catch
                    push!(folder_dates, 0)  # Handle folders that don't match date format
                end
                folder_data = all_filt[folder_idx]
                initial_response = mean(folder_data[1, :])
                num_cells = size(folder_data, 2)

                push!(all_x, idx)  # X position based on condition index
                push!(all_y, initial_response)
                push!(all_sizes, num_cells)
                push!(all_labels, key)
            end
        end
    end

    # Create scatter plot with dates as color axis
    scatter!(ax, all_x, all_y, markersize=all_sizes, color=all_dates, alpha=0.6, label="Data points")

    ax.xticks = (1:length(cond_keys), [replace(key, "hab_ISI" => "ISI", "_ITI" => "_ITI") for key in cond_keys])
    ax.xticklabelrotation = 45
    return f
end


function proc_data_stats(proc_data, cond_keys)
    # Create figure and layout
    f = Figure(size=(1200, 1200))

    # Create dropdown menus for ISI and ITI
    isi_options = [1, 2, 3]  # Available ISI values in minutes
    iti_options = [1, 2, 3, 5]  # Available ITI values in hours

    isi_menu = Menu(f[1, 1], options=zip(["1 min", "2 min", "3 min"], isi_options), default="1 min")
    iti_menu = Menu(f[1, 2], options=zip(["1 hr", "2 hr", "3 hr", "5 hr"], iti_options), default="1 hr")

    # Create axis for the stacked bar chart
    ax = Axis(f[2:5, 1:2], xlabel="Experiment", ylabel="Count", title="Data Statistics")

    # Create label for displaying counts
    count_label = Label(f[6, 1:2], "", tellwidth=false, tellheight=false)

    ax2 = Axis(f[2:3, 3:4], xlabel="Condition", ylabel="Count")

    stats = Dict{String,Any}()
    for k in cond_keys
        stats[k] = Dict{String,Any}()
        stats[k]["R"] = sum(map(i -> sum(i[1, :]), proc_data[k]["annots_filt"]))
        stats[k]["N"] = sum(proc_data[k]["nums_filt"]) - stats[k]["R"]
        stats[k]["P"] = round(stats[k]["R"] / (stats[k]["R"] + stats[k]["N"]), digits=2)
        stats[k]["X"] = sum(proc_data[k]["nums"]) - sum(proc_data[k]["nums_filt"])
    end

    barplot!(ax2, 1:length(cond_keys), [stats[k]["R"] for k in cond_keys])
    barplot!(ax2, 1:length(cond_keys), [stats[k]["N"] for k in cond_keys], offset=[stats[k]["R"] for k in cond_keys])
    barplot!(ax2, 1:length(cond_keys), [stats[k]["X"] for k in cond_keys], offset=[stats[k]["R"] for k in cond_keys] .+ [stats[k]["N"] for k in cond_keys])
    ax2.xticks = (1:length(cond_keys), [replace(key, "hab_ISI" => "ISI", "_ITI" => "_ITI") for key in cond_keys])
    ax2.xticklabelrotation = π / 2

    ax3 = Axis(f[4:5, 3:4], xlabel="Condition", ylabel="Proportion (P)", title="Proportion of Responders")
    lines!(ax3, 1:length(cond_keys), [stats[k]["P"] for k in cond_keys])
    ylims!(ax3, 0, 1)
    ax3.xticks = (1:length(cond_keys), [replace(key, "hab_ISI" => "ISI", "_ITI" => "_ITI") for key in cond_keys])
    ax3.xticklabelrotation = π / 2


    # Function to update the plot based on selected condition
    function update_plot(key)
        empty!(ax)

        data = proc_data[key]
        nums = data["nums"]
        nums_filt = data["nums_filt"]
        responder_counts = map(i -> sum(i[1, :]), data["annots_filt"])

        n_expts = length(nums)
        x_pos = 1:n_expts

        # Calculate the three categories for stacked bars
        nan_cells = nums .- nums_filt  # nums - nums_filt
        filtered_non_responders = nums_filt .- responder_counts  # nums_filt - responders

        # Create stacked bar chart
        barplot!(ax, x_pos, responder_counts, color=:blue, label="Responders")
        barplot!(ax, x_pos, filtered_non_responders, offset=responder_counts, color=:orange, label="Filtered non-responders")
        barplot!(ax, x_pos, nan_cells, offset=responder_counts .+ filtered_non_responders, color=:lightgray, label="Nan cells")

        # Add legend
        axislegend(ax, position=:rt)

        # Update title with current condition
        ax.title = "Data Statistics - $key"

        # Update count label
        count_label.text = "R: $(sum(responder_counts)) \t N: $(sum(filtered_non_responders)) \t P: $(round(sum(responder_counts)/((sum(responder_counts)+sum(filtered_non_responders))), digits=2)) \t X: $(sum(nan_cells))"
    end

    # Initial plot
    update_plot(get_key(isi_menu.selection[] * 60, iti_menu.selection[] * 3600))

    # Connect dropdown to update function
    on(isi_menu.selection) do selected_isi
        selected_iti = iti_menu.selection[]
        update_plot(get_key(selected_isi * 60, selected_iti * 3600))
    end
    on(iti_menu.selection) do selected_iti
        selected_isi = isi_menu.selection[]
        update_plot(get_key(selected_isi * 60, selected_iti * 3600))
    end

    return f
end


function show_control_heatmaps(processed_data, cond_keys; init_prob=0.75)
    # Create figure and layout
    f = Figure(size=(900, 800))

    # Create separate dropdown menus for ISI and ITI
    isi_options = [1, 2, 3]
    iti_options = [1, 2, 3, 5]
    isi_menu = Menu(f[1, 5], options=zip(string.(isi_options), isi_options), default=1, tellheight=false, tellwidth=false)
    iti_menu = Menu(f[2, 5], options=zip(string.(iti_options), iti_options), default=1, tellheight=false, tellwidth=false)

    # Create axis for the heatmap
    ax = Axis(f[1, 1:4], xlabel="Stimulus number", ylabel="Cell number",
        title="Control Data Heatmap")

    # Create axis for the mean plot
    ax_mean = Axis(f[2, 1:4], xlabel="Stimulus number", ylabel="Mean response",
        title="Mean Population Response")

    # Create axis to show recovery by folders alls_filt
    ax_recovery = Axis(f[3, 1:4], xlabel="Stimulus number", ylabel="Recovery",
        title="Recovery by folders alls_filt")
    ylims!(ax_recovery, 0, 1.2)
    ax_recovery.xticklabelrotation = pi / 4

    # Helper to get the current key from menu selections
    function get_current_key()
        isi = isi_menu.selection[]
        iti = iti_menu.selection[]
        return get_key(isi * 60, iti * 3600)
    end

    # Function to update the heatmap based on selected condition
    function update_heatmap()
        empty!(ax)
        empty!(ax_mean)
        key = get_current_key()
        data = processed_data[key]["control_data"]
        # Sort data by sum along axis 1 (rows)
        sort_indices = sortperm(vec(sum(data, dims=1)))
        sorted_data = data[:, sort_indices]
        heatmap!(ax, sorted_data)
        ax.title = "Control Data Heatmap - $key"

        # Plot mean across axis 2 for trial 1 (1-60) and trial 2 (61-120)
        trial1_mean = vec(mean(sorted_data[1:60, :], dims=2))
        trial2_mean = vec(mean(sorted_data[61:120, :], dims=2))

        scatter!(ax_mean, 1:60, trial1_mean, color=:black, label="Trial 1")
        scatter!(ax_mean, 1:60, trial2_mean, color=:gray, label="Trial 2")
        axislegend(ax_mean, position=:rt)
        ax_mean.title = "Mean Population Response - $key"
        ylims!(ax_mean, -0.1, 1.1)

        empty!(ax_recovery)
        folders = processed_data[key]["folders"]
        folder_data = processed_data[key]["annots_filt"]
        t2_inits = [mean(i[61, :]) for i in folder_data]
        recovery_vals = [mean(i[61, :]) / mean(i[1, :]) for i in folder_data]
        thresh = quantile.(Binomial.(size.(folder_data, 2), init_prob), 0.01) ./ size.(folder_data, 2)
        thresh_upper = quantile.(Binomial.(size.(folder_data, 2), init_prob), 0.99) ./ size.(folder_data, 2)
        inits = [mean(i[1, :]) for i in folder_data]
        scatter!(ax_recovery, (1:length(folders)) .+ 0.2, t2_inits, color=:black)
        scatter!(ax_recovery, 1:length(folders), inits, color=:red)
        scatter!(ax_recovery, 1:length(folders), thresh, color=:blue)
        # scatter!(ax_recovery, 1:length(folders), thresh_upper, color=:yellow)
        ax_recovery.xticks = (1:length(folders), folders .* " (" .* string.(1:length(folders)) .* ", " .* string.(size.(folder_data, 2)) .* ")")
        ax_recovery.title = "Recovery - $key"
    end

    # Initial heatmap
    update_heatmap()

    # Connect dropdowns to update function
    on(isi_menu.selection) do _
        update_heatmap()
    end
    on(iti_menu.selection) do _
        update_heatmap()
    end

    return f
end

function show_all_filtered_responses(proc_data)

    # Create the figure
    fig = Figure(size=(1000, 1000))

    # Create menus for ISI and ITI selection
    isi_options = [1, 2, 3]  # Available ISI values in minutes
    iti_options = [1, 2, 3, 5]  # Available ITI values in hours

    isi_menu = Menu(fig[1, 1], options=zip(["1 min", "2 min", "3 min"], isi_options), default="1 min")
    iti_menu = Menu(fig[1, 2], options=zip(["1 hr", "2 hr", "3 hr", "5 hr"], iti_options), default="1 hr")

    Label(fig[1, 1], "ISI:", tellwidth=false, halign=:right)
    Label(fig[1, 2], "ITI:", tellwidth=false, halign=:right)

    # Function to update the plot
    function update_plot(selected_isi, selected_iti)
        # Clear existing axes
        delete!.(fig.content[3:end])

        # Get new data
        new_key = get_key(selected_isi * 60, selected_iti * 3600)
        new_data = proc_data[new_key]
        new_all_filt_data = new_data["annots_filt"]
        new_num_filt_values = new_data["nums_filt"]

        # Recalculate grid dimensions
        new_n_plots = length(new_num_filt_values)
        new_n_cols = min(4, new_n_plots)
        new_n_rows = ceil(Int, new_n_plots / new_n_cols)

        # Create new axes
        for (idx, num_cells) in enumerate(new_num_filt_values)
            row = ceil(Int, idx / new_n_cols) + 1  # +1 to account for menu row
            col = ((idx - 1) % new_n_cols) + 1

            ax = Axis(fig[row, col],
                title="$(num_cells) cells\n$(new_data["folders"][idx])",
                xlabel="Stimulus number",
                ylabel="Mean response")

            # Take the corresponding data from all_filt_data
            selected_data = new_all_filt_data[idx]

            # Calculate means for trial 1 (stimuli 1-60) and trial 2 (stimuli 61-120)
            trial1_mean = vec(mean(selected_data[1:60, :], dims=2))
            trial2_mean = vec(mean(selected_data[61:120, :], dims=2))

            # Plot trial 1 and trial 2
            lines!(ax, 1:60, trial1_mean, color=:black, linewidth=2, label="Trial 1")
            lines!(ax, 1:60, trial2_mean, color=:gray, linewidth=2, label="Trial 2")

            # Set axis properties
            xlims!(ax, 0, 60)
            ylims!(ax, 0, 1)

            # Add legend only to first plot to avoid clutter
            if idx == 1
                axislegend(ax, position=:rt)
            end
        end

        # Update title
        Label(fig[0, :], "Filtered Response Visualization - ISI: $(selected_isi)min, ITI: $(selected_iti)hr", fontsize=18, tellwidth=false)
    end

    # Set up menu callbacks
    on(isi_menu.selection) do selected_isi
        selected_iti = iti_menu.selection[]
        update_plot(selected_isi, selected_iti)
    end

    on(iti_menu.selection) do selected_iti
        selected_isi = isi_menu.selection[]
        update_plot(selected_isi, selected_iti)
    end

    # Initial plot
    update_plot(1, 1)

    return fig
end


function plot_mean_responses_after_date_by_condition(proc_data, date_str; min_cells::Int=0, max_cells::Int=typemax(Int))
    # Convert input date string to Date
    input_date = Date(date_str, "yyyy_mm_dd")
    isi_vals = [1, 2, 3]
    iti_vals = [1, 2, 3, 5]
    n_isi = length(isi_vals)
    n_iti = length(iti_vals)
    n_conditions = n_isi * n_iti

    # Create the figure and a single axis
    f = Figure(size=(max(900, 60 * n_conditions), 450))
    ax = Axis(f[1, 1], xlabel="Condition", ylabel="Mean Response", title="Mean Initial Response by Condition")

    # Add a dropdown menu for stimulus selection (first or 61st)
    stim_options = Dict("First stimulus (1)" => 1, "61st stimulus" => 61, "Recovery (61/1)" => 3)
    stim_menu = Menu(f[0, 1], options=collect(keys(stim_options)), width=200, tellwidth=false)

    # Function to update the plot based on selected stimulus
    function update_plot(selected_stim_label)
        selected_stim = stim_options[selected_stim_label]
        all_x = Float64[]
        all_y = Float64[]
        all_sizes = Int[]
        all_colors = Symbol[]
        all_labels = String[]
        weighted_means = Float64[]
        weighted_x = Float64[]
        weighted_colors = Symbol[]
        xtick_labels = String[]
        xtick_positions = Float64[]
        all_cell_means = Float64[]

        # Assign a unique color for each condition (by ISI/ITI)
        condition_colors = [:blue, :green, :orange, :purple, :red, :cyan, :magenta, :brown, :pink, :teal, :olive, :gold]
        cond_idx = 1

        for (i_isi, isi) in enumerate(isi_vals)
            for (j_iti, iti) in enumerate(iti_vals)
                key = get_key(isi * 60, iti * 3600)
                if !haskey(proc_data, key)
                    cond_idx += 1
                    continue
                end
                folders = proc_data[key]["folders"]
                all_filt = proc_data[key]["annots_filt"]
                nums_filt = proc_data[key]["nums_filt"]
                means = Float64[]
                nums = Int[]
                folder_labels = String[]
                for (idx, folder) in enumerate(folders)
                    # Extract date part from folder name
                    folder_date_str = join(split(folder, "_")[1:3], "_")
                    folder_date = Date(folder_date_str, "yyyy_mm_dd")
                    n_cells = nums_filt[idx]
                    mat = all_filt[idx]
                    stim_idx = selected_stim
                    # Defensive: check matrix size
                    if stim_idx > size(mat, 1)
                        continue
                    end
                    if stim_idx == 3
                        mean_val = mean(mat[61, :]) / mean(mat[1, :])
                    else
                        mean_val = mean(mat[stim_idx, :])
                    end
                    push!(means, mean_val)
                    push!(nums, n_cells)
                    push!(folder_labels, folder)
                    # Collect all cell means for global mean (all points, regardless of cell count)
                    if stim_idx == 3
                        append!(all_cell_means, vec(mat[61, :]) / vec(mat[1, :]))
                    else
                        append!(all_cell_means, vec(mat[stim_idx, :]))
                    end
                end

                # For the condition mean, only use points that pass the filter
                weighted_mean = sum(means .* nums) / sum(nums)

                # X value for this condition
                xval = cond_idx
                # Add jitter to x positions for each folder
                jitter = 0.18
                xvals = xval .+ (rand(length(means)) .- 0.5) * jitter

                append!(all_x, xvals)
                append!(all_y, means)
                append!(all_sizes, nums)
                append!(all_colors, fill(condition_colors[mod1(cond_idx, length(condition_colors))], length(means)))
                append!(all_labels, fill("ISI: $(isi) min, ITI: $(iti) hr", length(means)))

                push!(weighted_means, weighted_mean)
                push!(weighted_x, xval)
                push!(weighted_colors, condition_colors[mod1(cond_idx, length(condition_colors))])

                push!(xtick_labels, "ISI: $(isi)\nITI: $(iti)")
                push!(xtick_positions, xval)

                cond_idx += 1
            end
        end

        # Compute global mean across all cells across all conditions
        global_mean = isempty(all_cell_means) ? NaN : mean(all_cell_means)

        # Clear the axis before plotting
        empty!(ax)

        # Plot individual folder means with jitter
        scatter!(ax, all_x, all_y; markersize=all_sizes, color=all_colors, label="Folders")

        # Plot the weighted mean as a single larger dot at the center of the column
        scatter!(ax, weighted_x, weighted_means; color=:black, markersize=20, marker=:circle, label="Weighted Mean")

        # Plot a horizontal line for the mean across all cells across all conditions
        if !isnan(global_mean)
            hlines!(ax, global_mean; color=:red, linestyle=:dash, linewidth=2, label="Global Mean")
        end

        # Set x-axis ticks to the condition labels
        ax.xticks = (xtick_positions, xtick_labels)
        ax.yticks = 0:0.1:1
        # Update y-label and title
        if selected_stim == 1
            ax.ylabel = "Mean Response (first stimulus)"
            ax.title = "Mean Initial Response by Condition"
        elseif selected_stim == 61
            ax.ylabel = "Mean Response (stimulus 61)"
            ax.title = "Mean Response at Stimulus 61 by Condition"
        else
            ax.ylabel = "Mean Response (stimulus $selected_stim)"
            ax.title = "Mean Response Recovery by Condition"
        end
        axislegend(ax, position=:rb)
    end

    # Initial plot
    update_plot(stim_menu.selection[])

    # Connect menu to update function
    on(stim_menu.selection) do selected_stim_label
        update_plot(selected_stim_label)
    end

    f
end
