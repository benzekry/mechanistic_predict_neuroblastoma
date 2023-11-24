import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pdb
import os
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

def log_rank_threshold(df,                          # dataframe
             variable,                              # string
             tte_name='PFS',
             threshold=None,
             event_name='progression',
             labels=["Low", "High"],
             output_file="survival",
             ylabel="Survival (%)",
             ci_show=True
             ):
    '''
    Dichotomized Kaplan-Meier plots for <variable>  from dataframe <df> at value <threshold>
    '''
    if threshold is None:
        threshold = df[variable].median()
    group_high_bool = df[variable] >= threshold
    group_low_bool = df[variable] < threshold
    # pdb.set_trace()
    if (sum(group_high_bool)==0) or  (sum(group_low_bool)==0):
        print("Warning. One of the two groups is empty. Analysis canceled")
        return
    fig, ax = plt.subplots();
    km = KaplanMeierFitter()
    km.fit(df.loc[group_high_bool, tte_name], df.loc[group_high_bool, event_name], label=labels[1])
    km.plot(ax=ax, ci_show=ci_show, show_censors=True)
    km.fit(df.loc[group_low_bool, tte_name], df.loc[group_low_bool, event_name], label=labels[0])
    km.plot(ax=ax, ci_show=ci_show, show_censors=True)
    # Cosmetics
    plt.ylim(0,1);
    ax.set_yticklabels([0, 20, 40, 60, 80, 100]);
    plt.ylabel(ylabel);
    plt.xlabel("Months");
    legend = ax.get_legend();
    plt.setp(legend.get_texts(), fontsize=18);
    legend.get_frame().set_linewidth(0.0);
    results = logrank_test(
        df.loc[group_high_bool, tte_name],
        df.loc[group_low_bool, tte_name],
        df.loc[group_high_bool, event_name],
        df.loc[group_low_bool, event_name]
        )
    # print('Logrank test, p = %.3g' % results.p_value)
    plt.text(0.05, 0.1, s="p = {:.3g}".format(results.p_value), transform=ax.transAxes, fontsize=16)
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))
    plt.savefig(output_file+".pdf", bbox_inches='tight')
    plt.show();

def log_rank(df,                          # dataframe
             group1,                      # string
             group2,                      # string
             tte_name='PFS',
             event_name='progression',
             labels=["Low", "High"],
             output_file="survival",
             ylabel="Survival (%)",
             xlabel="Months"
              ):
    '''
    Kaplan-Meier plots and logrank test between groups <group1> and <group2> of dataframe <df>
    '''
    group1_bool = (df[group1]==1)
    group2_bool = (df[group2]==1)
    fig, ax = plt.subplots();
    km1 = KaplanMeierFitter()
    km1.fit(df.loc[group1_bool, tte_name], df.loc[group1_bool, event_name], label=labels[0])
    km1.plot(ax=ax, ci_show=True, show_censors=True)
    km2 = KaplanMeierFitter()
    km2.fit(df.loc[group2_bool, tte_name], df.loc[group2_bool, event_name], label=labels[1])
    km2.plot(ax=ax, ci_show=True, show_censors=True)
    # Cosmetics
    plt.ylim(0,1);
    ax.set_yticklabels([0, 20, 40, 60, 80, 100]);
    plt.ylabel(ylabel);
    plt.xlabel(xlabel);
    legend = ax.get_legend();
    plt.setp(legend.get_texts(), fontsize=18);
    legend.get_frame().set_linewidth(0.0);
    plt.savefig(output_file+".pdf", bbox_inches='tight')
    results = logrank_test(
        df.loc[group1_bool, tte_name],
        df.loc[group2_bool, tte_name],
        df.loc[group1_bool, event_name],
        df.loc[group2_bool, event_name]
        )
    print('Logrank test, p = %.3g' % results.p_value)
    print("Median "+tte_name+" "+group1+" = %.3g months" % km1.median_)
    print("Median "+tte_name+" "+group2+" = %.3g months" % km2.median_)
    plt.show();
