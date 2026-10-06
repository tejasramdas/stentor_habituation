fold_redirect = "/Users/tejasramdas/Desktop/repos/cell_learning/stentor_habituation"
fig_folder = "$fold_redirect/analysis/main_analysis/figs/paper_figs/"
# pkg_folder = "$fold_redirect/analysis/"

using Pkg
Pkg.activate(".")

println("Loading packages...")

using Revise, JLD2, FileIO

# using GLMakie

isi_colors = Dict(1 => :blue, 2 => :green, 3 => :orange)
iti_alphas = Dict(1 => 0.9, 2 => 0.7, 3 => 0.5, 5 => 0.3)


includet("$fold_redirect/analysis/main_analysis/loader.jl")
includet("$fold_redirect/analysis/main_analysis/plotting.jl")
includet("$fold_redirect/analysis/main_analysis/turing.jl")

println("Loading data...")

# DATA_FOLDER = "processed_data"
# CHAINS_FOLDER = "chains"

inferred_chains = load("$fold_redirect/analysis/main_analysis/processed_data/inferred_chains.jld2", "chains")
processed_data = load("$fold_redirect/analysis/main_analysis/processed_data/collated.jld2", "data")

data_1_isi_1_iti = processed_data[get_key(60, 3600)]
responses = processed_data[get_key(60, 3600)]["control_data"]
chain = inferred_chains[get_key(60, 3600)]
pop_mean_demo = vec(mean(responses, dims=2))

function make_fig1()
    # Figure 1: Intro to the problem and recap of prior work
    fig1 = Figure(size=(1200, 400), fontsize=24)
    schematic_file = load(fig_folder * "schematic.pdf")
    stentor_contraction_file = load(fig_folder * "stentor_contraction.pdf")
    ax1_1 = Axis(fig1[1:6, 1:8])
    ax1_2 = Axis(fig1[1:3, 9:14])
    # ax1_3 = Axis(fig1[4:6, 9:14])
    # ax1_4 = Axis(fig1[1:6, 15:20])
    # Label(fig1[1:2, 1], "Schematic of Stentor contraction", fontsize = 16,tellwidth=false, tellheight=false)
    image!(ax1_1, rotr90(schematic_file))
    image!(ax1_2, rotr90(stentor_contraction_file))
    _, ax1_3 = make_protocol_heatmap(data_1_isi_1_iti["control_data"]; f=fig1[4:6, 9:14])
    hidexdecorations!(ax1_3, label=false)
    hideydecorations!(ax1_3, label=false)
    _, ax1_4 = make_t2_plot(data_1_isi_1_iti["control_data"]; f=fig1[1:6, 15:20])
    [ax.ygridvisible = false for ax in [ax1_1, ax1_2, ax1_3, ax1_4]]
    [ax.xgridvisible = false for ax in [ax1_1, ax1_2, ax1_3, ax1_4]]
    hidexdecorations!.([ax1_1, ax1_2])
    hideydecorations!.([ax1_1, ax1_2])

    # Align the axes with the image by setting alignmode
    [ax.alignmode = Outside() for ax in [ax1_1, ax1_2, ax1_3, ax1_4]]
    rowgap!(fig1.layout, 20)
    colgap!(fig1.layout, 20)
    ylims!(ax1_4, 0, 1)

    Label(fig1[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig1[1, 9, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig1[4, 9, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig1[1, 15, TopLeft()], "D", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig1
end

function make_fig2(chain=chain, responses=responses)
    # Figure 2 proves that our methods works well
    fig2 = Figure(size=(1300, 800), fontsize=24)
    ax2_1 = Axis(fig2[1, 1])
    ax2_2 = Axis(fig2[1, 2:3])
    ax2_3 = Axis(fig2[2, 1])
    ax2_4 = Axis(fig2[2, 2])
    ax2_5 = Axis(fig2[2, 3])
    make_hill_schematic(ax2_1)
    plot_inferred_sc_heatmap(ax2_2, chain[:trial_1], data=responses, plot_data=true, chain2=chain[:trial_2])
    # plot_inferred_population_hist(ax2_2, chain)
    plot_cell_curves(ax2_3, chain[:trial_1])
    plot_inferred_population_sc_median(ax2_3, chain[:trial_1])
    make_population_mean_plot(ax2_4, responses; label_lines=false)
    plot_inferred_population_mean_median(ax2_4, chain[:trial_1], chain2=chain[:trial_2])
    plot_curves(ax2_5, extract_params_median_turing(chain[:trial_1]), label="Trial 1", color=:black, axlabel=true)
    plot_curves(ax2_5, extract_params_median_turing(chain[:trial_2]), label="Trial 2", color=:gray, axlabel=false)
    plot_band(ax2_5, extract_params_samples_turing(chain[:trial_1]))
    plot_band(ax2_5, extract_params_samples_turing(chain[:trial_2]), color=:gray)
    ylims!(ax2_3, 0, 1.0)
    ylims!(ax2_4, 0, 1.0)
    ylims!(ax2_5, 0, 1.0)
    xlims!(ax2_3, 1, 60.0)
    ax2_3.xlabel = "Stimulus number (= mins)"
    ax2_5.xlabel = "Stimulus number (= mins)"
    ax2_3.ylabel = "Response probability"
    xlims!(ax2_4, 1, 60.0)
    xlims!(ax2_5, 1, 60.0)
    [ax.xgridvisible = false for ax in [ax2_1, ax2_2, ax2_3, ax2_4, ax2_5]]
    [ax.ygridvisible = false for ax in [ax2_1, ax2_2, ax2_3, ax2_4, ax2_5]]
    hidexdecorations!.([ax2_1, ax2_2, ax2_3, ax2_4, ax2_5], ticklabels=false, label=false, ticks=true)
    hideydecorations!.([ax2_1, ax2_2, ax2_3, ax2_4, ax2_5], ticklabels=false, label=false, ticks=true)
    axislegend(ax2_4)
    axislegend(ax2_5)
    [ax.alignmode = Outside() for ax in [ax2_1, ax2_2, ax2_3, ax2_4, ax2_5]]

    Label(fig2[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig2[1, 2, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig2[2, 1, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig2[2, 2, TopLeft()], "D", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig2[2, 3, TopLeft()], "E", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig2
end

function make_fig3()
    # Figure 3 expands to show potentiation across all conditions. Also discuss frequency dependency and decay of potentiation effect with larger ITI.
    fig3 = Figure(size=(1200, 600), fontsize=24)
    pop_params_1 = extract_params_samples_turing(chain[:trial_1])
    sc_params_1 = extract_params_median_turing(chain[:trial_1])
    pop_params_2 = extract_params_samples_turing(chain[:trial_2])
    sc_params_2 = extract_params_median_turing(chain[:trial_2])
    axs3_1 = [Makie.Axis(fig3[1, i], xgridvisible=false, ygridvisible=false) for i in 1:4]
    hist!(axs3_1[1], sc_params_1[:cell][:a_1], color=:black, normalization=:probability, bins=0:3:30, transparency=true, alpha=0.8, label="Trial 1")
    hist!(axs3_1[1], sc_params_2[:cell][:a_1], color=:gray, normalization=:probability, bins=0:3:30, transparency=true, alpha=0.8, label="Trial 2")
    hist!(axs3_1[2], sc_params_1[:cell][:b_1], color=:black, normalization=:probability, bins=0:3:30, transparency=true, alpha=0.8)
    hist!(axs3_1[2], sc_params_2[:cell][:b_1], color=:gray, normalization=:probability, bins=0:3:30, transparency=true, alpha=0.8)
    hist!(axs3_1[3], sc_params_2[:cell][:c_1] ./ sc_params_1[:cell][:c_1], color=:black, normalization=:probability, bins=0:0.2:2)
    hist!(axs3_1[4], sc_params_2[:cell][:b_1] ./ sc_params_1[:cell][:b_1], color=:black, normalization=:probability, bins=0:0.2:2)
    ylims!.(axs3_1, 0, 0.7)
    xlims!(axs3_1[2], 0, 30)
    xlims!(axs3_1[3], 0, 1.5)
    xlims!(axs3_1[4], 0, 2)
    [axs3_1[i].xlabel = label for (i, label) in enumerate(["Hill coefficient", "N50", "Initial response ratio", "N50 ratio"])]
    axs3_1[1].ylabel = "Proportion of cells"
    axislegend(axs3_1[1], position=:rt, fontsize=12)
    # axs3_1[2].ylabel = "Proportion of cells"

    axs3_2 = [Axis(fig3[2, i], xgridvisible=false, ygridvisible=false) for i in 1:4]
    plot_sc_param_scatter(axs3_2[1], sc_params_1[:cell][:c_1], sc_params_2[:cell][:b_1], plot_quantiles=true)
    xlims!(axs3_2[1], 0, 1.1)
    ylims!(axs3_2[1], 0, 30)
    axs3_2[1].xlabel = "Trial 1 Initial response"
    axs3_2[1].ylabel = "Trial 2 N50"
    plot_sc_param_scatter(axs3_2[2], sc_params_1[:cell][:c_1], sc_params_2[:cell][:c_1], norm_param_2=sc_params_1[:cell][:c_1], plot_quantiles=true)
    xlims!(axs3_2[2], 0.0, 1.1)
    ylims!(axs3_2[2], 0, 1.1)
    axs3_2[2].xlabel = "Trial 1 Initial response"
    axs3_2[2].ylabel = "Initial response ratio"
    plot_sc_param_scatter(axs3_2[3], sc_params_1[:cell][:b_1], sc_params_2[:cell][:b_1], norm_param_2=sc_params_1[:cell][:b_1], plot_quantiles=true)
    xlims!(axs3_2[3], 0, 30)
    ylims!(axs3_2[3], 0, 2)
    axs3_2[3].xlabel = "Trial 1 N50"
    axs3_2[3].ylabel = "N50 ratio"
    plot_sc_param_scatter(axs3_2[4], sc_params_2[:cell][:c_1], sc_params_2[:cell][:b_1], norm_param_1=sc_params_1[:cell][:c_1], norm_param_2=sc_params_1[:cell][:b_1], plot_quantiles=true)
    axs3_2[4].ylabel = "N50 ratio"
    axs3_2[4].xlabel = "Recovery"
    xlims!(axs3_2[4], 0.0, 1.1)
    [ax.alignmode = Outside() for ax in [axs3_1; axs3_2]]

    Label(fig3[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[1, 2, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[1, 3, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[1, 4, TopLeft()], "D", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[2, 1, TopLeft()], "E", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[2, 2, TopLeft()], "F", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[2, 3, TopLeft()], "G", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig3[2, 4, TopLeft()], "H", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig3
end

function make_fig4()
    fig4 = Figure(size=(1200, 800), fontsize=24)
    ax4_1 = [Makie.Axis(fig4[1, i]) for i in 1:3]
    ax4_2 = [Makie.Axis(fig4[2, i]) for i in 1:3]
    [ax.xgridvisible = false for ax in [ax4_1; ax4_2]]
    [ax.ygridvisible = false for ax in [ax4_1; ax4_2]]
    [hidexdecorations!.(ax, ticks=true, ticklabels=false, label=false) for ax in [ax4_1, ax4_2]]
    [hideydecorations!.(ax, ticks=true, ticklabels=false, label=false) for ax in [ax4_1, ax4_2]]
    _ = make_freq_plot(processed_data, ax4_1[1], time=false)
    # _ = plot_trial1_curves_all_conditions(inferred_chains, [1, 2, 3], [1, 2, 3, 5], ax4_2; time_axis=false)
    # _ = plot_trial1_curves_all_conditions(inferred_chains, [1, 2, 3], [1, 2, 3, 5], ax4_3; time_axis=true, leg=false)
    _ = plot_trial1_curves_all_conditions(inferred_chains, [1, 2, 3], 1, ax4_1[2]; time_axis=false, leg=false)
    _ = plot_trial1_curves_all_conditions(inferred_chains, [1, 2, 3], 1, ax4_1[3]; time_axis=true, leg=false)
    # _ = plot_recovery_by_iti(inferred_chains, ax4_2[1])
    # _ = plot_param_hist_by_isi(inferred_chains, ax4_2[1], param=:a_1_0)
    # _ = plot_param_hist_by_isi(inferred_chains, ax4_2[2], param=:b_1_0)
    # _ = plot_param_hist_by_isi(inferred_chains, ax4_2[3], param=:recovery)
    _ = plot_param_by_isi(inferred_chains, ax4_2[1], param=:hill)
    _ = plot_param_by_isi(inferred_chains, ax4_2[2], param=:b_1_0)
    _ = plot_param_by_isi(inferred_chains, ax4_2[3], param=:recovery)
    # _ = plot_trial1_curves_all_conditions_single_cell(inferred_chains, [1,2,3], [1,2,3,5];fig=fig4[3,1], time_axis=false, top_n_cells=30)
    # _ = plot_trial1_curves_all_conditions_single_cell(inferred_chains, [1,2,3], [1,2,3,5];fig=fig4[3,2], time_axis=true, top_n_cells=30)

    xlims!(ax4_1[1], 1, 60)
    xlims!(ax4_1[2], 1, 60)
    xlims!(ax4_1[3], 1, 180)
    ax4_1[3].xticks = [60, 120, 180]


    xlims!(ax4_2[1], 0, 4)
    xlims!(ax4_2[2], 0, 4)
    xlims!(ax4_2[3], 0, 4)
    ylims!(ax4_2[1], 0, 4)
    ylims!(ax4_2[2], 0, 30)
    ylims!(ax4_2[3], 0, 1.1)
    [ax.xticks = ([1, 2, 3], ["1", "2", "3"]) for ax in ax4_2]
    [ax4_2[i].xlabel = "ISI (mins)" for i in [1, 2, 3]]
    ax4_2[1].ylabel = "Log Hill coefficient"
    ax4_2[2].ylabel = "N50"
    ax4_2[3].ylabel = "Initial response ratio"

    [ax.alignmode = Outside() for ax in ax4_1]
    [ax.alignmode = Outside() for ax in ax4_2]

    Label(fig4[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig4[1, 2, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig4[1, 3, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig4[2, 1, TopLeft()], "D", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig4[2, 2, TopLeft()], "E", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig4[2, 3, TopLeft()], "F", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig4
end

function make_fig7()
    fig6 = Figure(size=(1500, 450), fontsize=24)
    # ax6_1 = [Makie.Axis(fig6[1, 1]) for i in 1:3]
    ax6_1 = [Makie.Axis(fig6[1, 1]), Makie.Axis(fig6[1, 2]), Makie.Axis(fig6[1, 3])]
    # make_derivative_example_plot(ax6_1, ax6_2, ax6_3)
    [ax.xgridvisible = false for ax in ax6_1]
    [ax.ygridvisible = false for ax in ax6_1]
    # hidexdecorations!.([ax6_1], ticks=true, ticklabels=false, label=false)
    plot_curves(ax6_1[1], extract_params_median_turing(chain[:trial_1]), half_max=true)
    plot_curves(ax6_1[1], extract_params_median_turing(chain[:trial_2]), label="Trial 2", color=:gray, half_max=true)
    axislegend(ax6_1[1], position=:rt, fontsize=10)
    # Vertical dotted lines from x-axis up to half-max scatter point
    p1 = extract_params_median_turing(chain[:trial_1])[:pop]
    p2 = extract_params_median_turing(chain[:trial_2])[:pop]
    b1 = p1[:b_1_0]
    b2 = p2[:b_1_0]
    halfmax_y1 = (p1[:c_1_0] + p1[:d_1_0]) / 2
    halfmax_y2 = (p2[:c_1_0] + p2[:d_1_0]) / 2
    linesegments!(ax6_1[1], [Point2f(b1 + 1, 0), Point2f(b1 + 1, halfmax_y1)], color=:black, linestyle=:dot)
    linesegments!(ax6_1[1], [Point2f(b2 + 1, 0), Point2f(b2 + 1, halfmax_y2)], color=:gray, linestyle=:dot)
    # Double-headed arrow between N50 values
    # arrow_y = min(halfmax_y1, halfmax_y2) - 0.05
    arrow_y = 0.1
    arrows2d!(ax6_1[1], [b1 + b2] / 2 .+ 1, [arrow_y * 0.5], [b2 - b1] / 3, [0.0], color=:red, shaftwidth=5, tipwidth=20)
    arrows2d!(ax6_1[1], [b1 + b2] / 2 .+ 1, [arrow_y * 0.5], [b1 - b2] / 3, [0.0], color=:red, shaftwidth=5, tipwidth=20)

    plot_curves_derivative(ax6_1[2], extract_params_median_turing(chain[:trial_1]))
    plot_curves_derivative(ax6_1[2], extract_params_median_turing(chain[:trial_2]), label="Trial 2", color=:gray)
    # Vertical dotted line at x=10
    vlines!(ax6_1[2], [5], color=:black, linestyle=:dot, alpha=0.5)
    # Sloping lines from x=1 to x=10 for each derivative curve
    dy1_at_1 = hill_derivative(1, p1[:a_1_0], p1[:b_1_0], p1[:c_1_0], p1[:d_1_0])
    dy1_at_5 = hill_derivative(5, p1[:a_1_0], p1[:b_1_0], p1[:c_1_0], p1[:d_1_0])
    dy2_at_1 = hill_derivative(1, p2[:a_1_0], p2[:b_1_0], p2[:c_1_0], p2[:d_1_0])
    dy2_at_5 = hill_derivative(5, p2[:a_1_0], p2[:b_1_0], p2[:c_1_0], p2[:d_1_0])
    linesegments!(ax6_1[2], [Point2f(1, dy1_at_1), Point2f(5, dy1_at_5)], color=:black, linestyle=:dash)
    linesegments!(ax6_1[2], [Point2f(1, dy2_at_1), Point2f(5, dy2_at_5)], color=:gray, linestyle=:dash)
    # Dot at midpoint of each line segment
    scatter!(ax6_1[2], [3], [(dy1_at_1 + dy1_at_5) / 2], color=:black, markersize=10, marker=:diamond)
    scatter!(ax6_1[2], [3], [(dy2_at_1 + dy2_at_5) / 2], color=:gray, markersize=10, marker=:diamond)
    plot_curves_vs_derivatives(ax6_1[3], extract_params_median_turing(chain[:trial_1]))
    plot_curves_vs_derivatives(ax6_1[3], extract_params_median_turing(chain[:trial_2]), label="Trial 2", color=:gray)
    # Red arrow at (0.6, -0.12) pointing at -45 degrees
    arrow_len = 0.04
    arrows2d!(ax6_1[3], [0.6], [-0.12], [arrow_len * cos(-π / 4)], [arrow_len * sin(-π / 4)], color=:red, shaftwidth=3, tipwidth=10)
    arrows2d!(ax6_1[3], [0.6], [-0.12], [-arrow_len * cos(-π / 4)], [-arrow_len * sin(-π / 4)], color=:red, shaftwidth=3, tipwidth=10)
    ylims!(ax6_1[1], 0, 1)
    xlims!(ax6_1[3], 0, 1)
    xlims!.(ax6_1[1:2], 1, 60)
    # xlims!.([ax6_1[4], ax6_1[5]], 1, 60)
    # [ax.xlabel = "Stimulus number" for ax in [ax6_5]]
    # [ax.ylabel = "Response" for ax in [ax6_1, ax6_4]]
    # [ax.ylabel = "Learning rate" for ax in [ax6_2, ax6_3, ax6_5, ax6_6]]
    # axislegend(ax6_1)
    [ax.alignmode = Outside() for ax in ax6_1]

    Label(fig6[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig6[1, 2, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig6[1, 3, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig6
end

function make_fig8(; potentiation_isi_list=[1, 2, 3])
    # Figure 7 shows inferred parameters for all conditions with error bars
    fig7 = Figure(size=(1500, 500), fontsize=24)
    ax7_1 = Makie.Axis(fig7[1:2, 1:2])
    ax7_2 = [Makie.Axis(fig7[1, i]) for i in 3:5]
    ax7_3 = [Makie.Axis(fig7[2, i]) for i in 3:5]
    ax7_1.xgridvisible = false
    ax7_1.ygridvisible = false
    [ax.xgridvisible = false for ax in ax7_2]
    [ax.ygridvisible = false for ax in ax7_3]
    hidexdecorations!(ax7_1, ticks=true, ticklabels=false, label=false)
    hidexdecorations!.(ax7_2, ticks=true, ticklabels=false, label=false)
    hidexdecorations!.(ax7_3, ticks=true, ticklabels=false, label=false)
    hideydecorations!(ax7_1, ticks=true, ticklabels=false, label=false)
    hideydecorations!.(ax7_2, ticks=true, ticklabels=false, label=false)
    hideydecorations!.(ax7_3, ticks=true, ticklabels=false, label=false)
    ax7_1.alignmode = Outside()
    [ax.alignmode = Outside() for ax in ax7_2]
    [ax.alignmode = Outside() for ax in ax7_3]
    plot_recovery_by_iti(inferred_chains, ax7_1, isi_list=[1, 2, 3])
    plot_b_m_by_iti(inferred_chains, ax7_2[1], isi_list=potentiation_isi_list)
    plot_trial_diff_ratio_by_iti(inferred_chains, ax7_2[2], isi_list=potentiation_isi_list)
    plot_derivative_y_diff_by_iti(inferred_chains, ax7_2[3], isi_list=potentiation_isi_list)
    plot_recovery_vs_b_m_0(inferred_chains, ax7_3[1], isi_list=[1, 2, 3], plot_ci=true)
    plot_recovery_vs_trial_diff_ratio(inferred_chains, ax7_3[2], isi_list=[1, 2, 3], plot_ci=true)
    plot_recovery_vs_derivative_diff(inferred_chains, ax7_3[3], isi_list=[1, 2, 3], plot_ci=true)
    ylims!(ax7_1, -0.0, 1.2)
    ylims!.([ax7_2[1], ax7_3[1]], -0.5, 1)
    ylims!.([ax7_2[2], ax7_3[2]], -5, 15)
    ylims!.([ax7_2[3], ax7_3[3]], -1, 15)
    xlims!.(ax7_3, -0.2, 1.6)
    # [ax.xlabel = " " for ax in ax7_2[[1, 3]]]
    # [ax.xlabel = " " for ax in ax7_3[[1, 3]]]
    # lines!(ax7_3, 0:0.01:1, 0:0.01:1, color=:black, linestyle=:dash)
    # lines!(ax7_3, [0.0, 1], [1, 1], color=:black, linestyle=:dash)
    # lines!(ax7_1[3], [1, 1], [0.0, 1], color=:black, linestyle=:dash)
    # axislegend(ax7_1[1], position=:rb, fontsize=10)
    # axislegend(ax7_3, position=:lb, labelsize=12)

    Label(fig7[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig7[1, 3, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig7[1, 4, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig7[1, 5, TopLeft()], "D", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig7[2, 3, TopLeft()], "E", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig7[2, 4, TopLeft()], "F", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig7[2, 5, TopLeft()], "G", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig7
end

function make_fig_s1()
    fig_s1 = Figure(size=(1600, 520), fontsize=22)
    ax1 = Axis(fig_s1[1, 1], xlabel="Stimulus number", ylabel="Response probability", title="Habituation curve")
    ax2 = Axis(fig_s1[1, 2], xlabel="Stimulus number", ylabel="Learning rate (dp/ds)", title="Learning-rate curve")
    ax3 = Axis(fig_s1[1, 3], xlabel="Response probability", ylabel="Learning rate (dp/ds)", title="Phase portrait")
    colgap!(fig_s1.layout, 35)

    xs = 1:0.1:60
    marker_xs = 1:60
    K = 20.0
    c = 0.9
    d = 0.1

    examples = [
        (n=5.0, color=colorant"#0072B2", label="Switch-like (n > 1)"),
        (n=0.7, color=colorant"#D55E00", label="Gradual (n < 1)"),
    ]

    for ex in examples
        y = hill.(xs, ex.n, K, c, d)
        dy = hill_derivative.(xs, ex.n, K, c, d)
        marker_y = hill.(marker_xs, ex.n, K, c, d)
        marker_dy = hill_derivative.(marker_xs, ex.n, K, c, d)

        lines!(ax1, xs, y, color=ex.color, linewidth=1.5, label=ex.label)
        scatter!(ax1, marker_xs, marker_y, color=ex.color, markersize=8)

        lines!(ax2, xs, dy, color=ex.color, linewidth=1.5)
        scatter!(ax2, marker_xs, marker_dy, color=ex.color, markersize=8)

        lines!(ax3, y, dy, color=ex.color, linewidth=1.5)
        scatter!(ax3, marker_y, marker_dy, color=ex.color, markersize=8)
    end

    axislegend(ax1, position=:rt, fontsize=10)
    xlims!(ax1, 1, 60)
    xlims!(ax2, 1, 60)
    xlims!(ax3, 0, 1)
    ylims!(ax1, 0, 1)
    ylims!(ax2, -0.06, 0.005)
    ylims!(ax3, -0.06, 0.005)
    [ax.xgridvisible = false for ax in [ax1, ax2, ax3]]
    [ax.ygridvisible = false for ax in [ax1, ax2, ax3]]
    [ax.alignmode = Outside() for ax in [ax1, ax2, ax3]]

    Label(fig_s1[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig_s1[1, 2, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig_s1[1, 3, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig_s1
end

function make_fig_s2()
    fig_s2 = Figure(size=(1600, 520), fontsize=22)
    ax1 = Axis(fig_s2[1, 1], xlabel="Stimulus number", ylabel="Response probability", title="Habituation curve")
    ax2 = Axis(fig_s2[1, 2], xlabel="Stimulus number", ylabel="Learning rate (dp/ds)", title="Learning-rate curve")
    ax3 = Axis(fig_s2[1, 3], xlabel="Response probability", ylabel="Learning rate (dp/ds)", title="Phase portrait")
    colgap!(fig_s2.layout, 35)

    xs = 1:0.1:60
    marker_xs = 1:60
    shift = 10.0

    trial_1 = (n=5.0, K=24.0, c=0.9, d=0.1, color=:black, label="Trial 1")
    trial_2 = (n=5.2, K=24.0, c=0.88, d=0.12, color=:gray, label="Trial 2")

    y1 = hill.(xs, trial_1.n, trial_1.K, trial_1.c, trial_1.d)
    dy1 = hill_derivative.(xs, trial_1.n, trial_1.K, trial_1.c, trial_1.d)
    marker_y1 = hill.(marker_xs, trial_1.n, trial_1.K, trial_1.c, trial_1.d)
    marker_dy1 = hill_derivative.(marker_xs, trial_1.n, trial_1.K, trial_1.c, trial_1.d)

    shifted_xs = xs .+ shift
    shifted_marker_xs = marker_xs .+ shift
    y2 = hill.(shifted_xs, trial_2.n, trial_2.K, trial_2.c, trial_2.d)
    dy2 = hill_derivative.(shifted_xs, trial_2.n, trial_2.K, trial_2.c, trial_2.d)
    marker_y2 = hill.(shifted_marker_xs, trial_2.n, trial_2.K, trial_2.c, trial_2.d)
    marker_dy2 = hill_derivative.(shifted_marker_xs, trial_2.n, trial_2.K, trial_2.c, trial_2.d)

    lines!(ax1, xs, y1, color=trial_1.color, linewidth=1.5, label=trial_1.label)
    scatter!(ax1, marker_xs, marker_y1, color=trial_1.color, markersize=8)
    lines!(ax1, xs, y2, color=trial_2.color, linewidth=1.5, label=trial_2.label)
    scatter!(ax1, marker_xs, marker_y2, color=trial_2.color, markersize=8)

    lines!(ax2, xs, dy1, color=trial_1.color, linewidth=1.5)
    scatter!(ax2, marker_xs, marker_dy1, color=trial_1.color, markersize=8)
    lines!(ax2, xs, dy2, color=trial_2.color, linewidth=1.5)
    scatter!(ax2, marker_xs, marker_dy2, color=trial_2.color, markersize=8)

    lines!(ax3, y1, dy1, color=trial_1.color, linewidth=1.5)
    scatter!(ax3, marker_y1, marker_dy1, color=trial_1.color, markersize=8)
    lines!(ax3, y2, dy2, color=trial_2.color, linewidth=1.5)
    scatter!(ax3, marker_y2, marker_dy2, color=trial_2.color, markersize=8)

    axislegend(ax1, position=:rt, fontsize=10)
    xlims!(ax1, 1, 60)
    xlims!(ax2, 1, 60)
    xlims!(ax3, 0, 1)
    ylims!(ax1, 0, 1)
    ylims!(ax2, -0.06, 0.005)
    ylims!(ax3, -0.06, 0.005)
    [ax.xgridvisible = false for ax in [ax1, ax2, ax3]]
    [ax.ygridvisible = false for ax in [ax1, ax2, ax3]]
    [ax.alignmode = Outside() for ax in [ax1, ax2, ax3]]

    Label(fig_s2[1, 1, TopLeft()], "A", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig_s2[1, 2, TopLeft()], "B", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    Label(fig_s2[1, 3, TopLeft()], "C", fontsize=24, padding=(5, 5, 5, 5), font="Arial Bold")
    return fig_s2
end

fig1 = make_fig1()
save(fig_folder * "/svg/fig1.svg", fig1; backend=CairoMakie)
save(fig_folder * "/pdf/fig1.pdf", fig1; backend=CairoMakie)
println("Saved figure 1...")

fig2 = make_fig2()
save(fig_folder * "/svg/fig2.svg", fig2; backend=CairoMakie)
save(fig_folder * "/pdf/fig2.pdf", fig2; backend=CairoMakie)
println("Saved figure 2...")

fig3 = make_fig3()
save(fig_folder * "/svg/fig3.svg", fig3; backend=CairoMakie)
save(fig_folder * "/pdf/fig3.pdf", fig3; backend=CairoMakie)
println("Saved figure 3...")

fig4 = make_fig4()
save(fig_folder * "/svg/fig4.svg", fig4; backend=CairoMakie)
save(fig_folder * "/pdf/fig4.pdf", fig4; backend=CairoMakie)
println("Saved figure 4...")

fig5, _ = plot_all_sc_curves_split(inferred_chains)
display(fig5)
save(fig_folder * "/svg/fig5.svg", fig5; backend=CairoMakie)
save(fig_folder * "/pdf/fig5.pdf", fig5; backend=CairoMakie)
println("Saved figure 5...")

# samps = extract_params_samples_turing(inferred_chains[get_key(60, 3600)][:trial_1])
# curves = hcat([hill.(1:0.1:60, Ref(samps[:pop][:a_1_0][i]), Ref(samps[:pop][:b_1_0][i]), Ref(samps[:pop][:c_1_0][i]), Ref(samps[:pop][:d_1_0][i])) for i in 1:size(samps[:pop][:a_1_0], 1)]...)
# lo = [quantile(curves[i, :], 0.05) for i in axes(curves, 1)]
fig6, _ = plot_all_sc_curves_vs_derivatives(inferred_chains; isi_vals=[1, 2, 3])
display(fig6)
save(fig_folder * "/svg/fig6.svg", fig6; backend=CairoMakie)
save(fig_folder * "/pdf/fig6.pdf", fig6; backend=CairoMakie)
println("Saved figure 6...")

fig7 = make_fig7()
save(fig_folder * "/svg/fig7.svg", fig7; backend=CairoMakie)
save(fig_folder * "/pdf/fig7.pdf", fig7; backend=CairoMakie)
println("Saved figure 7...")

fig8 = make_fig8()
save(fig_folder * "/svg/fig8.svg", fig8; backend=CairoMakie)
save(fig_folder * "/pdf/fig8.pdf", fig8; backend=CairoMakie)
println("Saved figure 8...")

fig8_v2 = make_fig8(potentiation_isi_list=[1])
save(fig_folder * "/svg/fig8_v2.svg", fig8_v2; backend=CairoMakie)
save(fig_folder * "/pdf/fig8_v2.pdf", fig8_v2; backend=CairoMakie)
println("Saved figure 8 v2...")

fig8_v3 = make_fig8(potentiation_isi_list=[2])
save(fig_folder * "/svg/fig8_v3.svg", fig8_v3; backend=CairoMakie)
save(fig_folder * "/pdf/fig8_v3.pdf", fig8_v3; backend=CairoMakie)
println("Saved figure 8 v3...")

fig8_v4 = make_fig8(potentiation_isi_list=[3])
save(fig_folder * "/svg/fig8_v4.svg", fig8_v4; backend=CairoMakie)
save(fig_folder * "/pdf/fig8_v4.pdf", fig8_v4; backend=CairoMakie)
println("Saved figure 8 v4...")


fig_s1 = make_fig_s1()
save(fig_folder * "/svg/fig_s1.svg", fig_s1; backend=CairoMakie)
save(fig_folder * "/pdf/fig_s1.pdf", fig_s1; backend=CairoMakie)
println("Saved figure s1...")

fig_s2 = make_fig_s2()
save(fig_folder * "/svg/fig_s2.svg", fig_s2; backend=CairoMakie)
save(fig_folder * "/pdf/fig_s2.pdf", fig_s2; backend=CairoMakie)
println("Saved figure s2...")

fig_s2, axs_s2 = plot_all_sc_curves_vs_derivatives(inferred_chains; isi_vals=[1, 2, 3])
display(fig_s2)

#############
# Reported statistics
#############

println("Population a")
println(round.(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:a_1_0]; prob=0.95), digits=1))

println("Population b median")
println(round.(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:b_1_0]; prob=0.95), digits=2))
println("Population b CI")
println(round.(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:b_1_0_low]; prob=0.95), digits=2))
println(round.(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:b_1_0_high]; prob=0.95), digits=2))

print("Cell halfmax range")
params = extract_params_median_turing(chain[:trial_1])
println(round.(quantile(params[:cell][:b_1], [0.025, 0.5, 0.975]), digits=2))

println("Population c")
println(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:c_1_0]; prob=0.95))
println("Population cell CI c")
println(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:c_1_0_low]; prob=0.95))
println(calculate_statistic_ci(chain[:trial_1], c -> c[:pop][:c_1_0_high]; prob=0.95))

println("Distribution of single cell b")
println(calculate_statistic_ci(chain[:trial_1], c -> [median(c[:cell][:b_1][i, :]) for i in 1:size(c[:cell][:b_1], 1)]; prob=0.95))
println(calculate_statistic_ci(chain[:trial_1], c -> [quantile(c[:cell][:b_1][i, :], 0.025) for i in 1:size(c[:cell][:b_1], 1)]; prob=0.95))
println(calculate_statistic_ci(chain[:trial_1], c -> [quantile(c[:cell][:b_1][i, :], 0.975) for i in 1:size(c[:cell][:b_1], 1)]; prob=0.95))

print("Population recovery")
println(calculate_statistic_ci(chain[:trial_1], chain[:trial_2], (c1, c2) -> c2[:pop][:c_1_0] ./ c1[:pop][:c_1_0]; prob=0.95))

print("Population potentiation")
println(calculate_statistic_ci(chain[:trial_1], chain[:trial_2], (c1, c2) -> c2[:pop][:b_1_0] ./ c1[:pop][:b_1_0]; prob=0.95))

println("N50 correlation with N50 ratio across cells")
calculate_statistic_ci_corr(chain[:trial_1], chain[:trial_2], (c1, c2) -> c1[:cell][:b_1], (c1, c2) -> (c2[:cell][:b_1] ./ c1[:cell][:b_1]); prob=0.95)

println("Initial response correlation with N50 across cells")
calculate_statistic_ci_corr(chain[:trial_1], chain[:trial_2], (c1, c2) -> c1[:cell][:c_1], (c1, c2) -> (c1[:cell][:b_1]); prob=0.95)

println("Initial response correlation with recovery ratio across cells")
calculate_statistic_ci_corr(chain[:trial_1], chain[:trial_2], (c1, c2) -> c1[:cell][:c_1], (c1, c2) -> (c2[:cell][:c_1] ./ c1[:cell][:c_1]); prob=0.95)

println("Recovery ratio correlation with potentiation across cells")
calculate_statistic_ci_corr(chain[:trial_1], chain[:trial_2], (c1, c2) -> (c2[:cell][:c_1] ./ c1[:cell][:c_1]), (c1, c2) -> c2[:cell][:b_1] ./ c1[:cell][:b_1]; prob=0.95)

median_a = zeros(3, 5)
low_a = zeros(3, 5)
high_a = zeros(3, 5)

median_b = zeros(3, 5)
low_b = zeros(3, 5)
high_b = zeros(3, 5)

median_c = zeros(3, 5)
low_c = zeros(3, 5)
high_c = zeros(3, 5)

median_d = zeros(3, 5)
low_d = zeros(3, 5)
high_d = zeros(3, 5)

median_recovery_ratio = zeros(3, 5)
low_recovery_ratio = zeros(3, 5)
high_recovery_ratio = zeros(3, 5)

median_potentiation = zeros(3, 5)
low_potentiation = zeros(3, 5)
high_potentiation = zeros(3, 5)

for isi in [1, 2, 3]
    for iti in [1, 2, 3, 5]
        vals = calculate_statistic_ci(inferred_chains[get_key(isi * 60, iti * 3600)][:trial_1], c -> c[:pop][:a_1_0]; prob=0.95)
        median_a[isi, iti] = vals[2]
        low_a[isi, iti] = vals[1]
        high_a[isi, iti] = vals[3]

        vals = calculate_statistic_ci(inferred_chains[get_key(isi * 60, iti * 3600)][:trial_1], c -> c[:pop][:b_1_0]; prob=0.95)
        median_b[isi, iti] = vals[2]
        low_b[isi, iti] = vals[1]
        high_b[isi, iti] = vals[3]

        vals = calculate_statistic_ci(inferred_chains[get_key(isi * 60, iti * 3600)][:trial_1], c -> c[:pop][:c_1_0]; prob=0.95)
        median_c[isi, iti] = vals[2]
        low_c[isi, iti] = vals[1]
        high_c[isi, iti] = vals[3]

        vals = calculate_statistic_ci(inferred_chains[get_key(isi * 60, iti * 3600)][:trial_1], c -> c[:pop][:d_1_0]; prob=0.95)
        median_d[isi, iti] = vals[2]
        low_d[isi, iti] = vals[1]
        high_d[isi, iti] = vals[3]

        vals = calculate_statistic_ci(inferred_chains[get_key(isi * 60, iti * 3600)][:trial_1], inferred_chains[get_key(isi * 60, iti * 3600)][:trial_2], (c1, c2) -> c2[:pop][:c_1_0] ./ c1[:pop][:c_1_0]; prob=0.95)
        median_recovery_ratio[isi, iti] = vals[2]
        low_recovery_ratio[isi, iti] = vals[1]
        high_recovery_ratio[isi, iti] = vals[3]

        vals = calculate_statistic_ci(inferred_chains[get_key(isi * 60, iti * 3600)][:trial_1], inferred_chains[get_key(isi * 60, iti * 3600)][:trial_2], (c1, c2) -> c2[:pop][:b_1_0] ./ c1[:pop][:b_1_0]; prob=0.95)
        median_potentiation[isi, iti] = vals[2]
        low_potentiation[isi, iti] = vals[1]
        high_potentiation[isi, iti] = vals[3]
    end
end

println("\n=== Hill coefficients (a) by ISI×ITI ===")
for isi in [1, 2, 3]
    println("ISI $(isi) min: ", round.(median_a[isi, [1, 2, 3]], digits=1), " (5hr: $(round(median_a[isi, 5], digits=1)))")
end

println("\n=== N50 (b) by ISI, 1hr ITI ===")
for isi in [1, 2, 3]
    println("ISI $(isi): median=$(round(median_b[isi, 1], digits=1)), 95% CI $(round(low_b[isi, 1], digits=1))–$(round(high_b[isi, 1], digits=1))")
end

println("\n=== Initial response (c) by ISI, 1hr ITI ===")
for isi in [1, 2, 3]
    println("ISI $(isi): median=$(round(median_c[isi, 1], digits=2)), 95% CI $(round(low_c[isi, 1], digits=2))–$(round(high_c[isi, 1], digits=2))")
end

println("\n=== Recovery ratio by ISI, 1hr ITI ===")
for isi in [1, 2, 3]
    println("ISI $(isi): median=$(round(median_recovery_ratio[isi, 1], digits=2)), 95% CI $(round(low_recovery_ratio[isi, 1], digits=2))–$(round(high_recovery_ratio[isi, 1], digits=2))")
end

println("\n=== Potentiation (N50 ratio) by ISI, 1hr ITI ===")
for isi in [1, 2, 3]
    println("ISI $(isi): median=$(round(median_potentiation[isi, 1], digits=2)), 95% CI $(round(low_potentiation[isi, 1], digits=2))–$(round(high_potentiation[isi, 1], digits=2))")
end

println("\n=== Cross-ISI recovery comparison ===")
iti = 1
println(calculate_statistic_ci_vec(inferred_chains[get_key(1 * 60, iti * 3600)][:trial_1], inferred_chains[get_key(1 * 60, iti * 3600)][:trial_2],
    inferred_chains[get_key(3 * 60, iti * 3600)][:trial_1], inferred_chains[get_key(3 * 60, iti * 3600)][:trial_2],
    (c1, c2) -> c2[:pop][:c_1_0] ./ c1[:pop][:c_1_0], (v1, v2) -> v2 ./ v1; prob=0.95))

iti = 5
println(calculate_statistic_ci_vec(inferred_chains[get_key(1 * 60, iti * 3600)][:trial_1], inferred_chains[get_key(1 * 60, iti * 3600)][:trial_2],
    inferred_chains[get_key(3 * 60, iti * 3600)][:trial_1], inferred_chains[get_key(3 * 60, iti * 3600)][:trial_2],
    (c1, c2) -> c2[:pop][:c_1_0] ./ c1[:pop][:c_1_0], (v1, v2) -> v2 ./ v1; prob=0.95))
