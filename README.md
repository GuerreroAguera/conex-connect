# conex-connect
Conex-Connect is a statistical extreme-value approach to study brain connectivity based on multi-channel EEG data.

## Organizing the data frame. 

To fully apply the Connex-Connect method, the data frame must have the EEG channels arranged in columns; column names are the EEG channels' names. Also, one must provide the following additional columns: 
Phase: 0 for pre-seizure phase; 1 for post-seizure phase.
Lead: Time lag to which you want to fit the model (e.g., 0, 10, 25, 50, ... seconds).
Time: time in seconds of the EEG records.

## Marginal Fit

To fit the generalized Pareto marginal model to each EEG channel, one must use the script "fit_marginal.R."

Provide the data frame to the object df_eeg following the instructions in Organizing the data frame. 

The fit setup is given in the section General Inputs of the script. The script is self-explicative; all input objects have a description with its aim.

Pay attention to the folder where to save the outputs (out_dir). This folder, with its files, is needed during the fit of the dependence structure model. **Do not change the pattern of the folder and file names**.

The objects nep and lambdas work as tunning parameters for the goodness-of-fit. If the diagnostic plots reveal a poor fit, one might need to change these objects' values.

## Conditional Extremal Dependence Model Fit

After fitting the marginal models, both for the pre- and post-seizure onset phases, the script "fit_ht.R" fits the conditional extremal dependence model.

The model setup is given in the section General Inputs of the script. The script is self-explicative; all input objects have a description with its aim.

Pay attention to the folder from where the marginal outputs (marginal_fit_dir) are loaded. Also, pay attention to the folder where to save the outputs (ht_fit_dir). **Do not change the pattern of the folder and file names*.

The objects nep and lbda_rng work as tunning parameters for the goodness-of-fit. If the diagnostic plots reveal a poor fit, one might need to change these objects' values.

Note that the number of bootstrap samples (boot_n) must be the same for both the marginal fits and conditional extremal dependence model.

## Minimum Working Example 

To run a complete analysis based on the provided data set (mwe_dataset.csv), run the following scripts in the given order: 

1. fit_marginal_lead_0_pre.R;

2. fit_marginal_lead_0_post.R;

3. fit_ht.R.
