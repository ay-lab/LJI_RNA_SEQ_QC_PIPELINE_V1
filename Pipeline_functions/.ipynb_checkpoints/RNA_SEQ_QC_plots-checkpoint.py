####################### Plot functions for RNA Seq pipeline - Paired end #######################

# Version: N/A (02/04/2020)

# Niu Du (ndu [at] lji.org)
# La Jolla Institute for Immunology (LJI)
# La Jolla, CA USA

# User environments
# python 3.X 
# seaborn 0.9.0
# plotly 4.4.1
# scanpy 1.4.1

import json
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotly.offline as py
from sklearn.decomposition import PCA,TruncatedSVD
import numpy as np
from scipy.stats import norm
import plotly.graph_objs as go
from RNA_SEQ_Func import read_bamcoverage,platt,fit_platt

def QC_plot(vis_parameters,dict_threshold):
    
    try:
        df_QC_report = pd.read_csv('QC_report.csv').set_index('Unnamed: 0')

        QC_category = list(set(df_QC_report['recommendation']))
        
        colors = sns.color_palette("Set1")


        # Making dot plots for QC parameters and save to png
        
        if (len(QC_category) > 1):
            if ('1.Good' in QC_category): 
                category_length = list(df_QC_report.groupby('recommendation').count().mean(axis = 1)[1:].values)
                category_length = [sum(category_length)/5] + category_length
            else:
                category_length = list(df_QC_report.groupby('recommendation').count().mean(axis = 1).values)
                
                
            f,ax = plt.subplots(1,len(QC_category),figsize = (20,3),gridspec_kw={'width_ratios': category_length},sharey = True)
            for i,category in enumerate(sorted(QC_category)):
                for c,parameter in enumerate(vis_parameters):
                    ax[i].plot(df_QC_report[df_QC_report['recommendation'] == category][parameter],'o',label = parameter,color = colors[c])
                    ax[i].hlines(y = dict_threshold[parameter],xmin = -0.5,xmax = len(df_QC_report[df_QC_report['recommendation'] == category])-0.5,linestyles = 'dashed',colors = colors[c],label = dict_threshold[parameter])

                xlabels = df_QC_report[df_QC_report['recommendation'] == category].index
                ax[i].set_title(category)
            
            
                if ('1.Good' in QC_category):     
                    if i%2 == 0:   
                        ax[i].set_title(category)
                    if i%2 == 1:   
                        ax[i].set_title(category+'\n')
                    if i != 0:
                        ax[i].set_xticklabels(xlabels,rotation = 90)
                    else:
                        ax[0].set_xticks('');
                        ax[0].set_xticklabels('');   
                else:
                    ax[i].set_xticklabels(xlabels,rotation = 90)
                    if i%2 == 0:   
                        ax[i].set_title(category)
                    if i%2 == 1:   
                        ax[i].set_title(category+'\n')
                        
            plt.legend()
            plt.subplots_adjust(wspace=0.0, hspace=0)
            plt.legend(bbox_to_anchor=(1., 1.))
            plt.savefig(f'QC_plots/{vis_parameters}.png',format='png', bbox_inches='tight')
        
        elif len(QC_category) == 1:
            f,ax = plt.subplots(figsize = (20,3))
            for c,parameter in enumerate(vis_parameters):
                ax.plot(df_QC_report[parameter],'o',label = parameter,color = colors[c])
                ax.hlines(y = dict_threshold[parameter],xmin = -0.5,xmax = len(df_QC_report)-0.5,linestyles = 'dashed',colors = colors[c],label = dict_threshold[parameter])
            xlabels = df_QC_report.index
            ax.set_title(QC_category[0])
            ax.set_xticklabels(xlabels,rotation = 90)
            plt.legend()
            plt.subplots_adjust(wspace=0.0, hspace=0)
            plt.legend(bbox_to_anchor=(1., 1.))
            plt.savefig(f'QC_plots/{vis_parameters}.png',format='png', bbox_inches='tight')
            
    except FileNotFoundError:
        print('Please generate QC report')
                   
            
        
def bam_plot():
    try:
        df_QC_report = pd.read_csv('QC_report.csv',index_col=0)
        
        # Making dot plots for bam coverage QC and save to png
        QC_category = list(set(df_QC_report['recommendation']))
        df_bamcoverage = pd.DataFrame(columns = [x.split('/')[-1] for x in glob('qualimap/*')],index = range(100))
        for col in df_bamcoverage.columns:
            df_bamcoverage[col] = read_bamcoverage(col) 
        for category in QC_category:
            sample_list = list(df_QC_report[df_QC_report['recommendation'] == category].index)
            df_bamcoverage[sample_list].plot(figsize = (8,5))
            plt.title(category)
            plt.legend(bbox_to_anchor=(1., 1.),ncol=max(1,int(len(sample_list)/12)))
            plt.savefig(f'QC_plots/bam_coverage_{category}.png',format='png', bbox_inches='tight')
            
    except FileNotFoundError:
        print('Please generate QC report')

def RNA_QC_spearman():
    '''Using spearman correlation to detect outliers'''
    try:
        df_corr = pd.read_csv('counts/TPM_counts.csv',index_col = 0).corr(method = 'spearman')
        sns.clustermap(df_corr,vmin = 0, vmax = 1,figsize = (int(len(df_corr)/3),int(len(df_corr)/3)))
        plt.savefig(f'QC_plots/Spearman_correlation.png',format='png', bbox_inches='tight')
    except FileNotFoundError:
        print('Please generate TPM table')
        
        
def soft_threshold(dict_conf):
    '''Detects minimal soft threshold for total STAR counts'''
    try:
        gene_recovery_perc = dict_conf['QC_threshold']['gene_recovery_perc']
        df_QC_report = pd.read_csv('QC_report.csv',index_col = 0)
        para = fit_platt([df_QC_report['final_STAR_counts'].values,df_QC_report['Total_genes'].values])
        fit_X = np.arange(0,max(df_QC_report['final_STAR_counts'].values),max(df_QC_report['final_STAR_counts'].values)/100)
        Y_max = max(platt(para,fit_X,0))
        Y_threshold = Y_max*gene_recovery_perc # Y percentage threshold based on input
        X_threshold = -np.log(1-gene_recovery_perc)*Y_max/para[1] # Calculated X threshold based on Y threshold value
        
        data_init_filter = df_QC_report[df_QC_report['final_STAR_counts']>X_threshold]['Total_genes']
        para_norm = norm.fit(data_init_filter)
        Y_lower = norm.ppf(0.05,para_norm[0],para_norm[1]) # New Y lower limit based on normal distribution
        
        plt.subplots(figsize = (12,6))
        plt.plot(df_QC_report['final_STAR_counts'].values,df_QC_report['Total_genes'].values,'o')
        X_fit = np.arange(0,max(df_QC_report['final_STAR_counts'].values),max(df_QC_report['final_STAR_counts'].values)/100)
        plt.plot(X_fit,platt(para,X_fit,0),'-')
        plt.axvline(x = max(X_threshold,dict_conf['QC_threshold']['final_STAR_counts']),linestyle = '--',color = 'r')
        plt.axhline(y = Y_lower,linestyle = '--',color = 'r')
        perc_threshold ="{:.1%}".format(gene_recovery_perc)
        
        if X_threshold > int(dict_conf['QC_threshold']['final_STAR_counts']):
            text = f'{int(X_threshold):,}'
            plt.text(X_threshold*1.1, Y_max/2, f'<- {gene_recovery_perc} gene recovery threshold = {text}')
        else:
            X_threshold = int(dict_conf['QC_threshold']['final_STAR_counts'])
            text = f'{int(X_threshold):,}'
            plt.text(int(X_threshold*1.1), Y_max/2, f'<- minimal set threshold = {text}')
            
        plt.text(X_threshold*1.5, Y_max/1.3, f'^ minimal gene threshold of {int(Y_lower)} at 0.95 confidence interval',color = 'r')
        plt.xlabel('STAR counts')
        plt.ylabel('Number of genes')
        plt.savefig(f'QC_plots/STAR_minimal_counts_soft_threshold.png',format='png', bbox_inches='tight')
    except FileNotFoundError:
        print('Please generate TPM table')
    
        
def Plot_3D(df_meta_input = None,n_comps = 10,dim = 3):
    try:
        df_QC_report = pd.read_csv('QC_report.csv',index_col=0)
        
        # Import metadata; if not provided only QC will be used
        try:
            df_meta_plot = pd.read_csv(df_meta_input,index_col=0)
            df_meta_plot = df_meta_plot.select_dtypes('object').merge(df_QC_report[['recommendation','Outlier']],left_index = True,right_index = True)
            df_meta_plot = df_meta_plot.astype(str)
            
        except FileNotFoundError:
            print('No metadata file found')
            df_meta_plot = df_QC_report[['recommendation','Outlier']]
            
        df_plot = pd.read_csv('./counts/TPM_counts.csv',index_col = 0)[df_meta_plot.index]
        
        # Compute PCA results for plotiing
        data_array = df_plot.apply(lambda x: np.log2(x+1)).values.T
        pca = PCA(n_components=3)
        PCA_result = pca.fit_transform(data_array)
        Ratio = pca.explained_variance_ratio_
        mu = PCA_result.mean(axis=0)
        PCA_result = PCA_result - mu
        # eigenvectors, eigenvalues, V = np.linalg.svd(PCA_result.T, full_matrices=False)
        # projected_data = np.dot(PCA_result, eigenvectors)
        # Scalar = projected_data.std(axis=0).mean()
        X = [xx[0] for xx in PCA_result]
        Y = [xx[1] for xx in PCA_result]
        Z = [xx[2] for xx in PCA_result]


        df_PCA = pd.DataFrame({'sample':df_meta_plot.index.values,'X':X,'Y':Y,'Z':Z}).set_index('sample')
        group_factors = df_meta_plot.columns

        scatter_data = []
        Scatter_data_length = []
        if dim == 3:
            for select_factor in group_factors:
                trace = []
                for N in sorted(set(df_meta_plot[select_factor])):
                    trace.append(go.Scatter3d(
                    x = df_PCA[df_meta_plot[select_factor] == N]['X'].values,
                    y = df_PCA[df_meta_plot[select_factor] == N]['Y'].values,
                    z = df_PCA[df_meta_plot[select_factor] == N]['Z'].values,
                    mode='markers',
                    name = N,
                    text=df_meta_plot[df_meta_plot[select_factor] == N].index.values,
                    hoverinfo='text',
                    marker=dict(
                        size=4,
                        line=dict(
                            width=0.5
                        ),
                        opacity=0.8
                    )
                ))

                scatter_data+=trace
                Scatter_data_length.append(len(trace))

        if dim == 2:
            for select_factor in group_factors:
                trace = []
                for N in sorted(set(df_meta_plot[select_factor])):
                    trace.append(go.Scatter(
                    x = df_PCA[df_meta_plot[select_factor] == N]['X'].values,
                    y = df_PCA[df_meta_plot[select_factor] == N]['Y'].values,
                    mode='markers',
                    name = N,
                    text=df_meta_plot[df_meta_plot[select_factor] == N].index.values,
                    hoverinfo='text',
                    marker=dict(
                        size=4,
                        line=dict(
                            width=0.5
                        ),
                        opacity=0.8
                    )
                ))

                scatter_data+=trace
                Scatter_data_length.append(len(trace))

        Vis_table = []
        for i in np.arange(len(Scatter_data_length)):
            Vis_list = []
            for j,elements in enumerate(Scatter_data_length):
                if j == i:
                    Vis_list.append([True]*elements)
                else:
                    Vis_list.append([False]*elements)
            Vis_table.append([item for sublist in Vis_list for item in sublist])
        dict_Vis = dict(zip(group_factors,Vis_table))

        updatemenus = [dict(active=0,
                 buttons=list([   
                    dict(label = Group_factor,
                         method = 'update',
                         args = [ {'visible': dict_Vis[Group_factor]+[True]},{'title': str(Group_factor)}]) for Group_factor in group_factors
                                     ])
                     )
            ]

        layout = go.Layout(
                legend=dict(x=-0.15, y=0),
                scene = dict(
                xaxis=dict(
                title='PC1  '+'%.2f'%Ratio[0],
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                    )
                ),

                yaxis=dict(
                title='PC2  '+'%.2f'%Ratio[1],
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                    )
                ),

                zaxis=dict(
                title='PC3  '+'%.2f'%Ratio[2],
                titlefont=dict(
                    family='Courier New, monospace',
                    size=18,
                    color='#7f7f7f'
                    )
                ),),

                updatemenus=updatemenus
                )
        fig = go.Figure(data=scatter_data,layout = layout)
        py.plot(fig,filename='check_QC_PCA.html', auto_open=False)
        
    except FileNotFoundError:
        print('Please generate QC report')