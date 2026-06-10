fold_redirect = "."

fig_folder = "$fold_redirect/analysis/main_analysis/figs/paper_figs/"

using Pkg
Pkg.activate(fold_redirect)
using Revise

using GLMakie
# ============================================
# Data loading
# ============================================
includet("$fold_redirect/analysis/main_analysis/loader.jl")

### Process raw data, remove NaN, control for initial response
# =============================================================
function get_processed_data(folder_json="$fold_redirect/analysis/experiment_folders/habituation_folders.json"; save_data=false, seed=42)
    expts_data = JSON.parsefile(folder_json)
    proc_data = Dict{String,Any}()
    for i in cond_keys
        proc_data[i] = get_cond(i, expts_data)
    end
    giant = hcat([proc_data[k]["comb_filt"] for k in cond_keys]...)
    println("Aggregate initial responses: ", round(mean(giant[1, :]), digits=3), " with $(size(giant,2)) cells")
    total_cells = 100
    num_responders = round(Int, 100 * mean(giant[1, :]))
    processed_data = get_control_data(proc_data; seed=seed, num_responders=num_responders, num_nonresponders=total_cells - num_responders)
    if save_data
        folder = Dates.format(now(), "yyyy_mm_dd")
        if !isdir("$fold_redirect/analysis/main_analysis/processed_data/$folder")
            mkdir("$fold_redirect/analysis/main_analysis/processed_data/$folder")
            save("$fold_redirect/analysis/main_analysis/processed_data/$folder/processed_data.jld2", Dict("data" => processed_data))
        else
            println("Folder $folder already exists")
        end
    end
    return processed_data
end

# processed_data = get_processed_data("$fold_redirect/analysis/experiment_folders/habituation_folders_new.json", save_data=false, seed=20)
# processed_data_filt = get_processed_data("$fold_redirect/analysis/experiment_folders/habituation_folders_filt.json", save_data=false, seed=1)

FOLDER = "processed_data"

processed_data = load("$fold_redirect/analysis/main_analysis/$FOLDER/collated.jld2")["data"]


# ============================================
# Probabilistic modeling
# ============================================
includet("$fold_redirect/analysis/main_analysis/turing.jl")

# Inference

function get_init_params(processed_data, key)
    init_c_1 = mean(processed_data[key]["control_data"][1, :])
    init_d_1 = mean(processed_data[key]["control_data"][60, :])
    init_c_2 = mean(processed_data[key]["control_data"][61, :])
    init_d_2 = mean(processed_data[key]["control_data"][120, :])
    b_1_thresh = init_d_1 + (init_c_1 - init_d_1) / 2
    b_2_thresh = init_d_2 + (init_c_2 - init_d_2) / 2
    init_b_1 = findfirst(x -> x < b_1_thresh, vec(mean(processed_data[key]["control_data"][1:60, :], dims=2)))
    init_b_2 = findfirst(x -> x < b_2_thresh, vec(mean(processed_data[key]["control_data"][61:120, :], dims=2)))
    init_from_param_1 = (a_1_0=2, b_1_0=init_b_1, c_1_0=init_c_1, d_1_0=init_d_1)
    init_from_param_2 = (a_1_0=2, b_1_0=init_b_2, c_1_0=init_c_2, d_1_0=init_d_2)
    return init_from_param_1, init_from_param_2
end


# Inference
function do_inference_single_trial(processed_data; save_chains=false, folder=nothing, num_samples=5000, nadapts=5000, num_chains=4)
    if folder === nothing && save_chains
        println("No folder specified, not saving chains")
    end
    inferred_chains = Dict{String,Any}()
    for isi in 1:3
        for iti in [1, 2, 3, 5]
            println("Inferring ISI: ", isi * 60, " ITI: ", iti * 3600)
            cond_chains = Dict{Symbol,Any}()
            init_param_1, init_param_2 = get_init_params(processed_data, get_key(isi * 60, iti * 3600))
            init_from_param_1 = InitFromParams(init_param_1)
            init_from_param_2 = InitFromParams(init_param_2)
            cond_chains[:trial_1] = sample(cond_data_single_trial(processed_data[get_key(isi * 60, iti * 3600)]["control_data"][1:60, :]), NUTS(; adtype=AutoMooncake()), MCMCThreads(), num_samples, num_chains, nadapts=nadapts, progress=true, initial_params=repeat([init_from_param_1], num_chains))
            cond_chains[:trial_2] = sample(cond_data_single_trial(processed_data[get_key(isi * 60, iti * 3600)]["control_data"][61:120, :]), NUTS(; adtype=AutoMooncake()), MCMCThreads(), num_samples, num_chains, nadapts=nadapts, progress=true, initial_params=repeat([init_from_param_2], num_chains))
            inferred_chains[get_key(isi * 60, iti * 3600)] = cond_chains
        end
    end
    if save_chains
        if folder === nothing
            println("No folder specified, not saving chains")
        else
            if !isdir("$fold_redirect/analysis/main_analysis/processed_data/$folder")
                mkdir("$fold_redirect/analysis/main_analysis/processed_data/$folder")
            end
            if !isfile("$fold_redirect/analysis/main_analysis/processed_data/$folder/inferred_chains.jld2")
                save("$fold_redirect/analysis/main_analysis/processed_data/$folder/inferred_chains.jld2", "chains", inferred_chains)
            else
                println("Chains already exist for $folder")
            end
        end
    end
    return inferred_chains
end

inferred_chains_single_trial = do_inference_single_trial(processed_data; num_samples=2500, nadapts=2500, num_chains=4)

# chain_folder = "chains"

# save("$fold_redirect/analysis/main_analysis/processed_data/inferred_chains.jld2", "chains", inferred_chains_single_trial)

inferred_chains = load("$fold_redirect/analysis/main_analysis/processed_data/inferred_chains.jld2")["chains"]

summarystats(inferred_chains[get_key(60, 3600)][:trial_1])[:a_1_0, :ess_bulk]
describe(inferred_chains[get_key(120, 10800)][:trial_1])

