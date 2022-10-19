# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 15:58:59 2021

@author: Shiwei Liu @ Harvard University

scripts for making figures for compiled mitotic measurements
"""

import pandas as pd
import numpy as np
import os, sys, glob, pathlib
import matplotlib.pyplot as plt
import math
import yaml, tqdm


# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 15:58:59 2021

@author: Shiwei Liu @ Harvard University

scripts for making figures for compiled mitotic measurements
"""

import pandas as pd
import numpy as np
import os, sys, glob, pathlib
import matplotlib.pyplot as plt
import math
import yaml, tqdm





# function for the boxplot
def boxplot_for_NE_assembly (df_plot, 
                             _exp_name = None, 
                             _NE_region = None,
                            _marker=None, 
                            _cell_line=None,
                             #_time_interval=1, # min
                             _show_datapoint=True,
                             _xrange = [],
                             _less_xticks=2,
                             _high_ylimit = 5,
                             _low_ylimit = 0.75,
                             save_fig = False,
                             save_folder = None,
                            ):
    # common params for figure and style
    fig, ax = plt.subplots(figsize=(3,2), dpi=300)
    ax.tick_params(direction='out', length=0.5, width=0.5, colors='k',
               grid_color='k', grid_alpha=0.5, size=2,labelsize=7.5)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    medianprops = {'color': 'red', 'linewidth': 1.2,'linestyle':'-'}
    meanprops = {'color': 'red', 'linewidth': 0,'linestyle':'-'}
    boxprops = {'color': 'black', 'linestyle': '-','linewidth': 0.5}
    whiskerprops = {'color': 'black', 'linewidth': 0.5,'linestyle': '-'}
    capprops = {'color': 'black', 'linewidth': 0.5,'linestyle': '-'}
    flierprops = {'color': 'black', 'marker': 'x'}
    _other_kwargs = {'whis':(10,90),'showfliers':False, 'showbox':True,'showmeans':True,'showcaps':True,'meanline':True,
                'grid':False}

    # retrieve only the measurement columns
    num_cols = [_col for _col in df_plot.columns if type(_col) is not str ]
    df_plot = df_plot [num_cols]           
    
    if len(_xrange)>0:
        _xrange_start,_xrange_end=_xrange[0],_xrange[1]
        df_plot=df_plot.iloc[:,_xrange_start:_xrange_end]

    # boxplot for the data df
    ax =  df_plot.boxplot(
           medianprops=medianprops,
           meanprops =meanprops,
           boxprops=boxprops,
           whiskerprops=whiskerprops,
           capprops=capprops,
           flierprops=flierprops, **_other_kwargs)

    # adjust the xyticks
    # use 3 for the core-noncore ratio;
    if np.max(df_plot.max())>1:
        _ylimit = _high_ylimit
    # use 1 for the other measurements;
    else:
        _ylimit = _low_ylimit
    ax.set_ylim(-0.1,_ylimit)
    # set tick before set ticklabel
    ax.set_yticks(np.arange(0,_ylimit,step=0.25))
    #ax.set_xticks(np.arange(round(df_plot.columns[0]-1),len(df_plot.columns)))
    #ax.set_xticks(np.arange(round(df_plot.columns[0]-1),round(df_plot.columns[-1]+1)))
    plt.xticks(rotation=-45,size=6)
    
    #ax.set_xticklabels(labels=np.arange(round(df_plot.columns[0]-1),len(df_plot.columns)),rotation=-45,size=6.5)
    _yticklabels=np.arange(0,_ylimit,step=0.25,dtype='float32')
    ax.set_yticklabels([str(round(float(label), 2)) for label in _yticklabels],size=6)
    # if reduce the xtick number
    for label in ax.xaxis.get_ticklabels()[1::_less_xticks]:
        label.set_visible(False)
    # reduce ytick for core-noncore ratio:
    if _ylimit == _high_ylimit:
        for label in ax.yaxis.get_ticklabels()[1::_less_xticks]:
            label.set_visible(False)
        
    
        
    # if show scattered points
    if _show_datapoint:
        for _col_index, _col in enumerate(df_plot.columns):
            y = df_plot[_col].dropna()
            # Add some random "jitter" to the x-axis
            x = np.zeros([len(y),])
            x[:]=_col_index+1
            ax.plot(x, y, 'k.', alpha=0.5, markersize=1)
    
    # set title 
    if _exp_name is not None and type(_exp_name) is str:
        _title = _exp_name.replace('_',' ')
        if _NE_region is not None and type(_NE_region) is str:
            _title = _title + ' ' +  _NE_region.replace('_',' ')
            if type(_marker) is str:
                _title = _marker+ ' ' +_title
                if type(_cell_line) is str:
                    _title = _cell_line+ ' ' +_title
                
    else:
        _title = 'analysis'
    
    _sample_size =0
    for _row_df in df_plot.iloc():
        if _row_df.isnull().sum() != len(_row_df):
            _sample_size+=1
        
    _title =_title+f' (N={_sample_size})'
    ax.set_title(f'{_title}', fontdict={'fontsize':9})

    
    # set legend
    ax.set_ylabel('Normalized accumulation level',fontdict={'fontsize':7.5})
    ax.set_xlabel('Normalized timepoint (min)',fontdict={'fontsize':7.5})
    
    # save figures (optional)
    if save_fig and os.path.exists(save_folder):
        fig_folder = os.path.join(save_folder,'Figures')
        if not os.path.exists(fig_folder):
            os.mkdir(fig_folder)   
        fig_name = _title.replace(' ', '_')
        #fig_name = '_'.join(fig_name.split('_')[:-1])
        fig.savefig(os.path.join(fig_folder,f'{fig_name}.pdf'),bbox_inches='tight')
        print (f'-- Saving analysis Figures for {fig_name}.')
        
            
    return ax