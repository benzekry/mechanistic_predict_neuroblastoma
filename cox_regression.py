import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.stats.api as sms
import os
import lifelines
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation
from shutil import copyfile
import warnings
import pdb
from packaging import version

float_format   = "%.3g"
pd.options.display.float_format = '{:,.3g}'.format
# Adapt the names of columns that will be used in cph.summary depending on lifelines version
if version.parse(lifelines.__version__) >= version.parse("0.22.0"):
    col_names_ci = ['coef lower 95%', 'coef upper 95%']
else:
    col_names_ci = ['lower 0.95', 'upper 0.95']
col_names_export_cox = ['exp(coef)', 'p'] + col_names_ci
def cox_regression(
    features_input,   # pd.DataFrame or <csv file> (index has to be first column)
    outcomes_input,   # pd.DataFrame or <csv file> index has to be first column
    duration_col,     # <str> name of column with time-to-event
    event_col,        # <str> name of column with censoring variable
    k_cross=10,       # number of folds in cross-validation
    n_repeats_cv=100,
    output_dir='cox_regression',
    signif_threshold=0.05,
    penalizer=0
    ):
    '''
    Implements Cox regression of each individual features and all together.
    '''
    output = dict()
    #-----------------------------------------------------------------------------
    # Load data
    #-----------------------------------------------------------------------------
    if isinstance(features_input, str):
        df_features_input = pd.read_csv(features_input, index_col=0)
    elif isinstance(features_input, pd.DataFrame):
        df_features = features_input
    df_features = df_features.astype('float64') # require for cross validation
    if isinstance(outcomes_input, str):
        df_outcomes = pd.read_csv(outcomes_input, index_col=0)
    elif isinstance(outcomes_input, pd.DataFrame):
        df_outcomes = outcomes_input
    #-----------------------------------------------------------------------------
    # Fill nan values of features with median
    #-----------------------------------------------------------------------------
    for column in df_features.columns:
        if sum(df_features[column].apply(np.isnan)) > 0:
            warnings.warn('Nan values found in column %s. Filled with median' % column)
        df_features[column] = df_features[column].fillna(value=df_features[column].median());
    df = pd.concat([df_features, df_outcomes[[duration_col, event_col]]], axis=1);
    #-----------------------------------------------------------------------------
    # Drop nan values of outcomes
    #-----------------------------------------------------------------------------
    if (sum(df[duration_col].isna() > 0)) | (sum(df[event_col].isna() > 0)):
        idxs_na_duration = df.index[df[duration_col].isna()]
        idxs_na_event    = df.index[df[event_col].isna()]
        idxs_na          = idxs_na_duration.append(idxs_na_event).unique()
        warnings.warn('There were %d missing values in outcomes (considering entries that had feature). They are removed from analysis' % len(idxs_na))
        df = df.dropna()
    #-----------------------------------------------------------------------------
    # Create output dir if does not exist
    #-----------------------------------------------------------------------------
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    #-----------------------------------------------------------------------------
    # Loop over individual features
    #-----------------------------------------------------------------------------
    print('-----------------------------------------------------------------------')
    print('Individual features')
    print('-----------------------------------------------------------------------')
    # file            = open(output_dir+'/cox_summary.tex', 'w')
    p_values        = dict()
    features_signif = []
    df_indiv        = pd.DataFrame(columns=col_names_export_cox)
    for feature in df_features.columns:
        cph = CoxPHFitter(penalizer=penalizer)
        cph.fit(df[[feature, duration_col, event_col]], duration_col=duration_col, event_col=event_col)
        print(feature)
        print('p = %.3g' % cph.summary['p'].values[0])
        p_values[feature] = cph.summary['p'].values[0]
        if p_values[feature] < signif_threshold:
            features_signif.append(feature)
        df_indiv = pd.concat([df_indiv, cph.summary[col_names_export_cox]])
        # LaTeX export
        # cph.summary.to_latex(output_dir+'/cox_summary_'+feature+'.tex', float_format='%.3g')
        # file.write('\subsection{\\texorpdfstring{\\nolinkurl{'+feature+'}}{'+feature+'}} \n')
        # file.write('\input{cox_summary_'+feature+'.tex} \n')
    # file.close()
    export_cox(df_indiv, output_file=os.path.join(output_dir, "cox_summary_univariate.tex"))
    #-----------------------------------------------------------------------------
    # All features
    #-----------------------------------------------------------------------------
    print('-----------------------------------------------------------------------')
    print('All features multivariate')
    print('-----------------------------------------------------------------------')
    cph_all = CoxPHFitter(penalizer=penalizer)
    cph_all.fit(df, duration_col=duration_col, event_col=event_col)
    # cph_all.print_summary()
    cph_all.summary.to_csv(output_dir+'/cox_summary_all_features.csv', float_format='%.3g')
    df_out_reduced = export_cox(cph_all.summary, os.path.join(output_dir, "cox_summary_all_features.tex"))
    print(df_out_reduced)
    output['cph_all'] = cph_all
    #-----------------------------------------------------------------------------
    ## Plot of HR
    #-----------------------------------------------------------------------------
    cph_all.plot();
    plt.xlabel("$\log(HR)$ (95% CI)")
    plt.savefig(os.path.join(output_dir, 'hazard_ratios_all_features.pdf'), bbox_inches='tight')
    plt.show()
    #-----------------------------------------------------------------------------
    ## Cross-validation
    #-----------------------------------------------------------------------------
    scores = []
    for replicate in range(n_repeats_cv):
        np.random.seed(replicate)
        cph_all_cv = CoxPHFitter(penalizer=penalizer)
        scores_loc = k_fold_cross_validation(cph_all_cv,
                                             df,
                                             duration_col=duration_col,
                                             event_col=event_col,
                                             scoring_method='concordance_index',
                                             k=k_cross)
        scores.append(scores_loc)
    scores = np.array(scores).flatten()
    print("")
    l,u = sms.DescrStatsW(scores).tconfint_mean() # confidence interval
    print("Mean (95%% C.I) c-index in %d fold cross-validation: %.3g (%.3g - %.3g)" % (k_cross, np.mean(scores), l, u))
    file_score = open(output_dir+'/mean_k_fold_score_all.tex', 'w')
    file_score.write('%.3g' % np.mean(scores))
    file_score.close()
    output['scores_k_fold_all_multivariate'] = scores
    #-----------------------------------------------------------------------------
    # Only features identified as significant in univariate analysis
    #-----------------------------------------------------------------------------
    print('-----------------------------------------------------------------------')
    print('Only significant features in univariate (p < %.3g): ' % signif_threshold)
    print(features_signif)
    print('-----------------------------------------------------------------------')
    cph_signif = CoxPHFitter(penalizer=penalizer)
    cph_signif.fit(
            df[features_signif+[duration_col, event_col]],
            duration_col=duration_col,
            event_col=event_col)
    print('Only significant features')
    # cph_signif.print_summary()
    df_out_reduced = export_cox(cph_signif.summary, output_file=os.path.join(output_dir, "cox_summary_features_signif.tex"))
    print(df_out_reduced)
    output['cph_signif_univariate'] = cph_signif
    #-----------------------------------------------------------------------------
    ## Plot of HRs
    #-----------------------------------------------------------------------------
    cph_signif.plot();
    plt.xlabel("$\log(HR)$ (95% CI)")
    plt.savefig(os.path.join(output_dir, 'hazard_ratios_only_signif.pdf'), bbox_inches='tight')
    plt.show()
    #-----------------------------------------------------------------------------
    ## Cross-validation
    #-----------------------------------------------------------------------------
    scores = []
    for replicate in range(n_repeats_cv):
        np.random.seed(replicate)
        cph_signif_cv = CoxPHFitter(penalizer=penalizer)
        scores_loc = k_fold_cross_validation(cph_signif_cv,
                                             df[features_signif+[duration_col, event_col]],
                                             duration_col=duration_col,
                                             event_col=event_col,
                                             scoring_method='concordance_index',
                                             k=k_cross)
        scores.append(scores_loc)
    scores = np.array(scores).flatten()
    print("")
    l,u = sms.DescrStatsW(scores).tconfint_mean() # confidence interval
    print("Mean (95%% C.I) c-index in %d fold cross-validation: %.3g (%.3g - %.3g)" % (k_cross, np.mean(scores), l, u))
    file_score = open(output_dir+'/mean_k_fold_score_signif.tex', 'w')
    file_score.write('%.3g' % np.mean(scores))
    file_score.close()
    output["features_signif"]                 = features_signif
    output["scores_k_fold_signif_univariate"] = scores
    #-----------------------------------------------------------------------------
    # Only features identified as significant in multivariate analysis
    #-----------------------------------------------------------------------------
    print('-----------------------------------------------------------------------')
    print('Only significant features in multivariate (p < %.3g): ' % signif_threshold)
    p_values_multivariate        = cph_all.summary['p']
    features_signif_multivariate = cph_all.summary.index[p_values_multivariate < signif_threshold].tolist()
    print(features_signif_multivariate)
    print('-----------------------------------------------------------------------')
    cph_signif_mv = CoxPHFitter(penalizer=penalizer)
    cph_signif_mv.fit(
            df[features_signif_multivariate+[duration_col, event_col]],
            duration_col=duration_col,
            event_col=event_col)
    print('Only significant features')
    # cph_signif_mv.print_summary()
    df_out_reduced = export_cox(cph_signif_mv.summary, output_file=os.path.join(output_dir, "cox_summary_features_signif_multivariate.tex"))
    print(df_out_reduced)
    output['cph_signif_multivariate'] = cph_signif_mv
    #-----------------------------------------------------------------------------
    ## Plot of HRs
    #-----------------------------------------------------------------------------
    cph_signif_mv.plot();
    plt.xlabel("$\log(HR)$ (95% CI)")
    plt.savefig(os.path.join(output_dir, 'hazard_ratios_only_signif_multivariate.pdf'), bbox_inches='tight')
    plt.show()
    #-----------------------------------------------------------------------------
    ## Cross-validation
    #-----------------------------------------------------------------------------
    scores = []
    for replicate in range(n_repeats_cv):
        np.random.seed(replicate)
        cph_signif_mv_cv = CoxPHFitter(penalizer=penalizer)
        scores_loc = k_fold_cross_validation(cph_signif_mv_cv,
                                             df[features_signif_multivariate+[duration_col, event_col]],
                                             duration_col=duration_col,
                                             event_col=event_col,
                                             scoring_method='concordance_index',
                                             k=k_cross)
        scores.append(scores_loc)
    scores = np.array(scores).flatten()
    print("")
    l,u = sms.DescrStatsW(scores).tconfint_mean() # confidence interval
    print("Mean (95%% C.I) c-index in %d fold cross-validation: %.3g (%.3g - %.3g)" % (k_cross, np.mean(scores), l, u))
    file_score = open(output_dir+'/mean_k_fold_score_signif_multivariate.tex', 'w')
    file_score.write('%.3g' % np.mean(scores))
    file_score.close()
    output["features_signif_multivariate"]      = features_signif_multivariate
    output["scores_k_fold_signif_multivariate"] = scores
    #-----------------------------------------------------------------------------
    # Compile output TeX file
    #-----------------------------------------------------------------------------
    if os.path.exists("/Users/benzekry/work/code/stat_tools_python/cox_results.tex"):
        current_dir = os.getcwd()
        os.chdir(output_dir)
        copyfile("/Users/benzekry/work/code/stat_tools_python/cox_results.tex", "results.tex")
        os.system("pdflatex results.tex")
        os.system("pdflatex results.tex")
        os.system('rm *.aux *.log *.out *.toc')
        os.chdir(current_dir)
    return output
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
# Utils
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
def export_cox(df, output_file):
    df_out_reduced = df.copy()[col_names_export_cox]
    df_out_reduced[col_names_ci] = np.exp(df_out_reduced[col_names_ci])
    df_out_reduced = df_out_reduced[col_names_export_cox]
    df_out_reduced = df_out_reduced.rename(columns={"exp(coef)": "Hazard ratio",
                                                    'coef lower 95%': 'coef lower 95\%',
                                                    'coef upper 95%': 'coef upper 95\%'
                                                    })
    df_out_reduced.to_latex(output_file, float_format=float_format, escape=False)
    return df_out_reduced
