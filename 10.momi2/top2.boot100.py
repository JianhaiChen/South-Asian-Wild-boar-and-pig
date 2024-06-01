import momi
import autograd.numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import logging
import sys
import pandas as pd

logging.basicConfig(level=logging.INFO, filename="migration_models_bootstrap.log")

# Load the SFS data
sfs = momi.Sfs.load("sfstop2.gz")

# Define yticks for the plot
yticks = [2e5, 4e5, 6e5, 1e6, 2e6]
pop_x_positions = {"sus": 0, "SAS": 1, "ISEAWB": 2, "mseaw": 3}

# Function to plot and save demographic models
def plot_demography(model, filename, bootstrap_results=None):
    fig = momi.DemographyPlot(model=model,
                              pop_x_positions=pop_x_positions,
                              figsize=(7, 10),
                              major_yticks=yticks,
                              linthreshy=5e4, pulse_color_bounds=(0, .25),
                              draw=False)
    if bootstrap_results:
        for params in bootstrap_results:
            fig.add_bootstrap(params, alpha=0.3)
    fig.draw()
    fig.draw_N_legend(loc="upper right")
    plt.savefig(f"scebog.pakibangledesh.{filename}")

# 1. No migration model
no_migration_model = momi.DemographicModel(N_e=1e5)
no_migration_model.set_data(sfs)
no_migration_model.add_size_param("n_sus")
no_migration_model.add_size_param("n_SAS")
no_migration_model.add_size_param("n_mseaw")
no_migration_model.add_size_param("n_ISEAWB")
no_migration_model.add_leaf("sus", N="n_sus")
no_migration_model.add_leaf("SAS", N="n_SAS")
no_migration_model.add_leaf("mseaw", N="n_mseaw")
no_migration_model.add_leaf("ISEAWB", N="n_ISEAWB")
no_migration_model.move_lineages("SAS", "sus", t=1360000)
no_migration_model.move_lineages("mseaw", "SAS", t=685000)
no_migration_model.move_lineages("ISEAWB", "mseaw", t=360000)
no_migration_model.optimize(method="TNC")
no_migration_log_lik = no_migration_model.log_likelihood()
print("No migration model log likelihood:", no_migration_log_lik)
plot_demography(no_migration_model, "no_migration_modeltop2.pdf")

# 2. Migration model from sus to ISEAWB
migration_sus_to_ISEAWB_model = no_migration_model.copy()
migration_sus_to_ISEAWB_model.add_time_param("t_migration", upper=350000)
migration_sus_to_ISEAWB_model.add_pulse_param("mfrac_sus_ISEAWB", lower=0, upper=0.3)
migration_sus_to_ISEAWB_model.move_lineages("ISEAWB", "sus", t="t_migration", p="mfrac_sus_ISEAWB")
migration_sus_to_ISEAWB_model.optimize(method="TNC")
migration_sus_to_ISEAWB_log_lik = migration_sus_to_ISEAWB_model.log_likelihood()
print("Migration model (sus to ISEAWB) log likelihood:", migration_sus_to_ISEAWB_log_lik)
plot_demography(migration_sus_to_ISEAWB_model, "migration_sus_to_ISEAWB_modeltop2.pdf")

# Bootstrap 100 runs for Migration model (sus to ISEAWB)
n_bootstraps = 100
bootstrap_results = []
for i in range(n_bootstraps):
    print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps}")
    resampled_sfs = sfs.resample()
    migration_sus_to_ISEAWB_model.set_data(resampled_sfs)
    migration_sus_to_ISEAWB_model.set_params(randomize=True)
    migration_sus_to_ISEAWB_model.optimize(method="TNC")
    bootstrap_results.append(migration_sus_to_ISEAWB_model.get_params())

# Save the plot with bootstrap results
plot_demography(migration_sus_to_ISEAWB_model, "migration_sus_to_ISEAWB_model_bootstraptop2.pdf", bootstrap_results=bootstrap_results)

# Calculate median and 95% interval for each parameter
bootstrap_df = pd.DataFrame(bootstrap_results)
median_values = bootstrap_df.median()
ci_lower = bootstrap_df.quantile(0.025)
ci_upper = bootstrap_df.quantile(0.975)

# Save the results to log
logging.info("Bootstrap Results for Migration model (sus to ISEAWB)")
logging.info(f"Median Values:\n{median_values}")
logging.info(f"95% Confidence Intervals:\nLower:\n{ci_lower}\nUpper:\n{ci_upper}")

# Print number of parameters (degrees of freedom) and parameter estimates for each model
def print_params(model, model_name):
    params = model.get_params()
    print(f"{model_name} parameters ({len(params)}):")
    for param, value in params.items():
        print(f"  {param}: {value}")

print("No migration model log likelihood:", no_migration_log_lik)
print_params(no_migration_model, "No migration model")

print("Migration model (sus to ISEAWB) log likelihood:", migration_sus_to_ISEAWB_log_lik)
print_params(migration_sus_to_ISEAWB_model, "Migration model (sus to ISEAWB)")

print("Bootstrap median values for Migration model (sus to ISEAWB):")
print(median_values)
print("95% Confidence Intervals for Migration model (sus to ISEAWB):")
print(f"Lower:\n{ci_lower}")
print(f"Upper:\n{ci_upper}")