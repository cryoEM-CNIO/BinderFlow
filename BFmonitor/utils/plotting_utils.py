import pandas as pd
import os
import re
import time 
import numpy as np
import glob 
import argparse
import subprocess
import plotly.express as px
import plotly.graph_objects as go
from Bio import SeqIO

# Function to update scatter plot, histograms, and row count

def update_scatter_plot(directory, merged_df,filtered_df, x_variable, y_variable, input_pdb_path):
    '''
    updates the scatter plot (and plot it indeed)

    Input:

    directory --> Directory where all the data we want to plot is stored. In this case the parent directory where the output folder and Scoring_Stats.csv are located

    filtered_df --> Filtered df with only those hits that fulfill all the conditions imposed

    x_variable --> Variable that is plotted in the x axis

    y_variable --> Variable that is plotted in the y axis

    input_pdb_path --> Path where the input is stored (for scaffold or PD comparisons)

    Output:

    A beautiful and dynamic scatter plot, where the gray points are those that not fullfil the metrics (except no design can be considered a hit, then all are colored). The size of the points indicate if they are the input (bigger) or new desigs (smaller)
    '''
    try:
        input_pdb_path=glob.glob(f'{directory}/{input_pdb_path}')[0]
    except:
        input_pdb_path=None
    if not merged_df.empty:
        for description in merged_df['description']:
            if not filtered_df.empty:
                if (filtered_df['description']==description).any():
                    continue
                else:
                    merged_df.loc[merged_df['description']==description, "length"] = np.nan

        merged_df['size'] = merged_df.apply(
        lambda row: float(1) if (
            row['original_design'] == input_pdb_path
            ) else float(0.5),
        axis=1
        )

        row_count_text = f"{len(merged_df)}"
        scatter_plot = px.scatter(merged_df,
                                  y=y_variable,
                                  x=x_variable,
                                  color='length',
                                  marginal_x = "violin",
                                  marginal_y = "violin",
                                  render_mode = 'webgl',
                                  hover_data=['description','plddt_binder','pae_interaction','length'],
                                  size='size'
                                )
        if x_variable == 'pae_interaction' or x_variable =='CUTRE' or x_variable == 'dG' or x_variable == 'binder_surf_hyd':
            scatter_plot.update_xaxes(autorange="reversed")        
        
        scatter_plot.update_traces(marker=dict(opacity=0.7))
        return scatter_plot, row_count_text
    else:
        return None, "No data available."

def radar_plot(designs_list, merged_df, pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    '''
    Input:

    designs_list--> List of the two designs to compare. Original design is being used for the posterior filtering 

    merged_df --> Dataframe with all the stats. This can be computed inside the function, but i doubt the legibility of the function will change

    Thresholds --> All thresholds for the parameters. This are needed for the normalization of the variables and the comparisons with the thresholds

    Output:

    A beautiful radar plot whose thresholds can be changed. The values are normalized respect to the thresholds (plotly.express does not yet support an assymetric radar chart). If the value is above the metric, it is fulfilled (and the further from the line, the better is that metric).
    '''

    if len(designs_list) == 2:
        if designs_list[0] == designs_list[1]:
            designs_list=[designs_list[0]]  
    
    selected_variables=['pae_interaction', 'CUTRE', 'plddt_binder', 'dSASA', 'Shape_complementarity', 'interface_hbonds', 'interface_unsat_hbonds', 'binder_surf_hyd', 'description']
    #Add normalization respect to the parameters used in the sliders
    normalized_df=merged_df.drop([column for column in merged_df.columns if column not in selected_variables], axis=1)

    #####################This is a defect for development, remove for deployment###############
    try:
        normalized_df.drop('SCORE:', axis=1)
    except KeyError:
        normalized_df=normalized_df
    ###############################################################################################
    normalized_df['pae_interaction']=max(pae_interaction_thres)/normalized_df['pae_interaction']
    normalized_df['CUTRE']=max(CUTRE_thres)/normalized_df['CUTRE']
    normalized_df['plddt_binder']=normalized_df['plddt_binder']/min(plddt_binder_thres)
    normalized_df['dSASA']=normalized_df['dSASA']/min(dsasa_thres)
    normalized_df['Shape_complementarity']=normalized_df['Shape_complementarity']/min(shape_complementarity_thres)
    normalized_df['interface_hbonds']=normalized_df['interface_hbonds']/min(interface_hbond_thres)
    normalized_df['interface_unsat_hbonds']=max(interface_unsat_hbond_thres)/normalized_df['interface_unsat_hbonds']
    normalized_df['binder_surf_hyd']=max(binder_surf_hyd_thres)/normalized_df['binder_surf_hyd']

    #If any value == inf change it by 2

    normalized_df=normalized_df[normalized_df['description'].isin(designs_list)]
    normalized_df.replace([np.inf, -np.inf], 2, inplace=True)
    normalized_df.replace([np.NaN], 1, inplace=True)

    # Use pd.Categorical to set the order of the 'original_design' column based on designs_list
    normalized_df['description'] = pd.Categorical(
            normalized_df['description'],
            categories=designs_list,
            ordered=True
        )

    # Sort the DataFrame by the ordered 'original_design' column
    normalized_df = normalized_df.sort_values('description')
    
    #All the text variables
    text_variables=['description']
    
    #Change to the format the radar plot needs 
    reformatted_dict={
        'values1':[],
        'variables':[],
        'values2':[]
    }
    
    for variable in normalized_df.columns:
        if variable not in text_variables:
            reformatted_dict['variables']   .append(variable)
            if len(normalized_df.index) > 0:    
                reformatted_dict['values1'].append(normalized_df[variable].iloc[0])
            else:
                reformatted_dict['values1'].append(np.nan)
            if len(normalized_df.index)==2:
                reformatted_dict['values2'].append(normalized_df[variable].iloc[1])
            else:
                reformatted_dict['values2'].append(np.nan)
        else:
            continue

    reformatted_df=pd.DataFrame(reformatted_dict)
    
    AU=compute_area(normalized_df, designs_list)

    #Plot!
    radar_figure=go.Figure()
    #First design
    if len(designs_list) > 0:
        radar_figure.add_trace(
                go.Scatterpolar(
                r=list(reformatted_df['values1'])+ [reformatted_df['values1'][0]],
                theta=list(reformatted_df['variables'])+[reformatted_df['variables'][0]],
                name=f'{AU["design"][0].replace("_substituted_dldesign","").replace("_af2pred","")} --> {AU["AreaUnits"][0]} AU',
                fillcolor = 'rgba(134,48,113,0.5)', fill = 'toself',
                line=dict(color='#863071'))
        )
    #Second design
    if len(designs_list)==2:
        radar_figure.add_trace(
            go.Scatterpolar(
                r=list(reformatted_df['values2'])+ [reformatted_df['values2'][0]],
                theta=list(reformatted_df['variables'])+[reformatted_df['variables'][0]],
                name=f'{AU["design"][1].replace("_substituted_dldesign","").replace("_af2pred","")} --> {AU["AreaUnits"][1]} AU',
                fillcolor = 'rgba(237,176,129,0.5)', fill ='toself',
                line=dict(color='#edb081'))
        )
    #Threhsolds 
    radar_figure.add_trace(
        go.Scatterpolar(
            r=[1]*len(reformatted_df['variables'])+ [1],
            theta=list(reformatted_df['variables'])+[reformatted_df['variables'][0]],
            name='Reference Line',
            line=dict(color='black', width=2, dash='dot',),
        )
    )

    return radar_figure


def compute_area(normalized_df, designs_list):
    '''
    Function to compute the area of the radar plot, which can be used as metric for the binder selection

    Input:

    normalized_df --> DF with the values of each variable normalized by the thresholds values

    designs_list --> List with the designs selected 

    Output:

    AU --> Dictionary of the area units of the radar plot of the different designs selected
    '''
    AU={
        'design':designs_list,
        'AreaUnits':[]
    }

    normalized_df=normalized_df.drop('description', axis=1)

    r_values=[]
    
    for _,row in normalized_df.iterrows():
        area=0
        r_values=row.values.tolist() 
        n = len(r_values)  # Number of vertices
        theta = 2 * np.pi / n  # Angle between consecutive vertices
        sin_theta = np.sin(theta)
    
        for i in range(n):
            r_i = r_values[i]
            r_next = r_values[(i + 1) % n]  # Wrap around to the first vertex
            area += r_i * r_next
    
        area *= 0.5 * sin_theta
        AU['AreaUnits'].append(round(area,3))
    
    return AU

def update_designs_list(designs_list, design_to_plot):
    """
    Updates the designs_list based on the provided design_to_plot data.
    
    Input:
        designs_list (list): A list to store selected designs. It should ideally have up to 2 elements.
        design_to_plot (dict): A dictionary containing the design data to be added, 
                               specifically in 'points[0]["customdata"][0]'.
    Output:
        None: The function updates the designs_list in place.
    """
    try:
        # Extract the design data from the input dictionary
        new_design = design_to_plot["points"][0]["customdata"][0]
    except (KeyError, IndexError, TypeError) as e:
        print(f"Error extracting design data: {e}")
        print("No valid design selected.")
        return

    # Update the designs_list based on its length
    if len(designs_list) < 2:
        # Add the new design if the list has fewer than 2 items
        designs_list.append(new_design)
    else:
        # Replace both designs in the list
        designs_list[1] = designs_list[0]
        designs_list[0] = new_design

def scatter_plot_AF3(working_dir, xaxis, yaxis):
    '''
    Function to plot data from AF3 and compare with AF2

    Input:
    
    - working_dir: Parent directory of the campaign
    - xaxis: Variable of the xaxis
    - yaxis: Variable of the yaxis

    Output:

    - Scatterplot 
    '''

    data_df = pd.read_csv(os.path.join(working_dir, 'Scoring_Stats_AF3.csv'))

    scatter_plot = px.scatter(data_df,
                            y=yaxis,
                            x=xaxis,
                            color='length',
                            marginal_x = "violin",
                            marginal_y = "violin",
                            render_mode = 'webgl',
                            hover_data=['description','plddt_binder','pae_interaction','length']
                            )
    if xaxis == 'pae_interaction':
        scatter_plot.update_xaxes(autorange="reversed")

        scatter_plot.update_traces(marker=dict(opacity=0.7))
        return scatter_plot
    else:
        return scatter_plot