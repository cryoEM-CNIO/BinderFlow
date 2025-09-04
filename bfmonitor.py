#!/usr/bin/env python3

'''
POSSIBLE PROBLEMS

- If you don't use the FastRelax cycle to improve the sequences, the names of the files change. In this case, you should remove the cycle_1 from the csv of Rosetta and CUTRE
This could have an easy fix

'''

'''
things to do:

- Creating files similar to johnbercow for synthesis
- Implementing johnbercow for our use in this webapp
- Pointing out hotspot residue in the PDB viewer from molstar
- Creating campaigns from the webapp
'''

import os
import pandas as pd
import dash
import dash_molstar as dashmolstar
from dash_molstar.utils import molstar_helper
from dash import dash_table
import dash_bio as dashbio
from dash_molstar.utils.representations import Representation
from dash import html
import dash_bootstrap_components as dbc
from dash import Dash, dcc, Input, Output, State, callback
from dash.exceptions import PreventUpdate
import dash_bio.utils.ngl_parser as ngl_parser
import plotly.express as px
import subprocess
import re
from Bio import PDB
import argparse
import socket
import numpy as np
import time
import os.path
from BFmonitor.utils.hits_utils import *
from BFmonitor.utils.generic_utils import *
from BFmonitor.utils.plotting_utils import *




# Parsing port number and host
parser = argparse.ArgumentParser()
parser.add_argument("--port", "-p", help = "choose port to run the webapp")
parser.add_argument("--host", "-host", help = "choose host to run the webapp")
parser.add_argument("--debug", "-d", help = "launch app in debug mode")

args, unknown = parser.parse_known_args()

# Set localhost and 8051 as host and port by default
if not args.port: port_number = 8050
else: port_number = args.port
if not args.host: hostname = socket.gethostname()
else: hostname = args.host 
if not args.debug: debug_mode = False
else: debug_mode = True

# Messsssssinesssss
working_dir = os.getcwd()
output_dir = os.path.join(working_dir, 'output')
directories_list=get_working_directories(working_dir)
if not directories_list:
    directories_list=[working_dir]
designs_list=[]

script_dir = os.path.dirname(os.path.abspath(__file__))
#Initial lists for extraction dropdowns

initial_organisms=[
    "Arabidopsis thaliana",
    "Bacillus subtilis",
    "Caenorhabditis elegans",
    "Chlamydomonas reinhardtii",
    "Danio rerio",
    "Drosophila melanogaster",
    "Homo sapiens",
    "Mus musculus",
    "Nicotiana tabacum",
    "Pseudomonas putida",
    "Saccharomyces cerevisiae",
    "Escherichia coli general",
]

# Styles
bt_style = {"align-items": "center", "background-color": "#F2F3F4", "border": "2px solid #000",
            "box-sizing": "border-box", "color": "#000", "cursor": "pointer", "display": "inline-flex",
             'padding':'0.3em 1.2em', 'margin':'0.5em 0.3em 0.3em 0',
            "font-size": "0.9em", 'font-weight':'500', 'border-radius':'2em', 'text-align':'center',
            'transition':'all 0.2s'}
dropdown_style={'background-color':'#F2F3F4',
                "color": "#000", "cursor": "pointer",'width':'260px',
                 'margin':'0.5em 0.3em 0.3em 0',
                "font-size": "1.2em", 'font-weight':'500', 'text-align':'center',
                'transition':'all 0.2s'}

table_viewer_div_style={'display': 'flex', 'justifyContent': 'flex-start'}
table_div_style={'marginRight': '30px'}
table_style={'overflowX': 'auto','width': '100%', 'margin': '0 auto'}
table_cell_style={'minWidth': '150px', 'width': '150px', 'maxWidth': '150px', 'overflow': 'hidden', 'textOverflow': 'ellipsis', 'padding': '5px',}
title_style = {"margin-left": "15px", "margin-top": "15px", "margin-bottom": "0em", "color": "Black", "font-family" : "Helvetica", "font-size":"2.5em"}
box_style3 = {"font-size":"0.9em",'padding':'0.3em 1.2em','width':"18%", "margin-left": "0%","margin-right": "1%", "color": "black","font-family" : "Helvetica", 'vertical-align': 'center', "margin-bottom":"2px"}
extraction_box_style={'padding': '20px','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'height':'100px'}
extraction_box_cl_style={'padding': '20px','max-height':'240px','overflow-y':'auto','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'height':'290px'}
cool_button_style={'font-size': '24px','padding': '30px 60px','background': 'linear-gradient(135deg, #4CAF50, #81C784)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex':'1','min-width': '300px', 'flex-basis': '350px','height': '290px'}
stop_button_style = {'font-size': '18px','padding': '30px 60px','background': 'linear-gradient(135deg, #FF0000, #FF6347)','color': 'white','border': 'none','border-radius': '10px','box-shadow': '0px 4px 6px rgba(0, 0, 0, 0.1)','cursor': 'pointer','transition': 'transform 0.2s, box-shadow 0.2s','flex': '1','min-width': '300px','flex-basis': '350px','height': '200px'}
extraction_cl_style={'margin-top':'50px','max-height':'750px','overflow-y':'auto','background-color': '#f9f9f9','border-radius': '10px','box-shadow': '0 4px 6px rgba(0, 0, 0, 0.1)','flex': '1','min-width': '250px', 'flex-basis': '300px', 'display':'flex', 'align-items':'flex-start'}
# Color definitions
color_points = 'cornflowerblue'

# Dash app initialization
# app = dash.Dash(__name__)

import dash_bootstrap_components as dbc

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
    assets_folder='BFmonitor/assets'
)
app.title = "BinderFlow Monitor"
server = app.server

# Define layout
def serve_layout():
    return html.Div(
        className="app-wrapper",
        children=[
            # Sidebar
            html.Div(
                className="sidebar",
                children=[
                    html.H1("BinderFlow Monitor", className="app-title"),
                    # Removed: html.Div(id='row-count'),
                    dcc.Dropdown(
                        options=[
                            {
                                'label': (
                                    f"{os.path.basename(d)}"
                                ),
                                'value': d
                            }
                            for d in directories_list
                        ],
                        value=directories_list[0],
                        placeholder='Select a folder',
                        id='directory-dropdown',
                        className='dropdown'
                    ),
                    html.Br(),
                    html.Button('STOP CAMPAIGN', id='stop-campaign', n_clicks=0, className="button button-danger"),
                    html.Br(),
                    html.Button('Filters & Axes', id='open-filters', n_clicks=0, className='button'),
                    dbc.Collapse(
                        html.Div([
                            html.H5('PAE Interaction'),
                            dcc.RangeSlider(id='pae_interaction_thres', min=0, max=30, value=[0,10], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('CUTRE'),
                            dcc.RangeSlider(id='CUTRE_thres', min=0, max=70, value=[0,10], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('pLDDT binder'),
                            dcc.RangeSlider(id='plddt_binder_thres', min=0, max=100, value=[80,100], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('dSASA'),
                            dcc.RangeSlider(id='dsasa_thres', min=0, max=10000, value=[1000,10000], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Shape Complementarity'),
                            dcc.RangeSlider(id='shape_complementarity_thres', min=0, max=1, value=[0.5,1], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Interface HBond'),
                            dcc.RangeSlider(id='interface_hbond_thres', min=0, max=15, value=[3,15], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Interface Unsatisfied HBond'),
                            dcc.RangeSlider(id='interface_unsat_hbond_thres', min=0, max=15, value=[0,4], tooltip={"placement":"bottom","always_visible":True}),
                            html.H5('Binder Surface Hydrophobicity'),
                            dcc.RangeSlider(id='binder_surf_hyd_thres', min=0, max=1, value=[0,0.35], tooltip={"placement":"bottom","always_visible":True}),
                            html.Hr(),
                            html.Div([
                                html.H4('Axis Values')]),
                            html.Div([
                                html.Span('X', className='axis-label'),
                                dcc.Dropdown([
                                                'plddt_binder',
                                                'pae_interaction',
                                                'CUTRE',
                                                'dG',
                                                'dSASA',
                                                'Shape_complementarity',
                                                'Packstat',
                                                'dG_SASA_ratio',
                                                'length',
                                                'SAP',
                                                'binder_int_hyd',
                                                'binder_surf_hyd',
                                                'interface_hbonds',
                                                'interface_unsat_hbonds',
                                                'ipSAE',
                                                'RMSD'
                                            ], 'pae_interaction', id='xaxis_value', className='dropdown', style = {'width': '100%'})
                            ], className='axis-control'),
                            html.Div([
                                html.Span('Y', className='axis-label'),
                                dcc.Dropdown( [
                                                'plddt_binder',
                                                'pae_interaction',
                                                'CUTRE',
                                                'dG',
                                                'dSASA',
                                                'Shape_complementarity',
                                                'Packstat',
                                                'dG_SASA_ratio',
                                                'length',
                                                'SAP',
                                                'binder_int_hyd',
                                                'binder_surf_hyd',
                                                'interface_hbonds',
                                                'interface_unsat_hbonds',
                                                'ipSAE',
                                                'RMSD'
                                            ], 'plddt_binder', id='yaxis_value', className='dropdown', style= {'width': '100%'})
                            ], className='axis-control'),
                        ], style={'padding':'10px'}),
                        id='filters-collapse',
                        is_open=False
                    ),
                    html.Br(),
                ]
            ),
            # Main content
            html.Div(
                className="main-content",
                children=[
                    html.Div([
                        dcc.Tabs([
                            dcc.Tab(label='Live Watcher', children=[
                                # Metrics cards row
                                html.Div(
                                    children=[
                                        dbc.Card(
                                            [
                                                html.Div(id='row-count', className='metric-value'),
                                                html.Div("Finished Designs", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                        dbc.Card(
                                            [
                                                html.Div(id='hit-count', className='metric-value'),
                                                html.Div("Hits", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                        dbc.Card(
                                            [
                                                html.Div(id='hit-efficiency', className='metric-value'),
                                                html.Div("Hit Efficiency", className='metric-label'),
                                            ],
                                            className="metric-card"
                                        ),
                                    ],
                                    className="metrics-row"
                                ),
                                html.Div([
                                    html.Div([
                                        html.Div([
                                            # Removed filter/axis toggles and panels; controls now in Offcanvas
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader("Scatter Plot"),
                                                    dbc.CardBody(
                                                        dcc.Graph(
                                                            id='scatter-plot',
                                                            className='graph-container',
                                                            config={'toImageButtonOptions': {'format': 'svg', 'filename': 'scatter_plot', 'height': 600, 'width': 800, 'scale': 1}, 'displaylogo': False}
                                                        )
                                                    ),
                                                ],
                                                className="graph-card"
                                            ),
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader("Radar Plot"),
                                                    dbc.CardBody(
                                                        dcc.Graph(
                                                            id='radar-plot',
                                                            className='graph-container',
                                                            config={'toImageButtonOptions': {'format': 'svg', 'filename': 'radar_plot', 'height': 600, 'width': 800, 'scale': 1}, 'displaylogo': False}
                                                        )
                                                    ),
                                                ],
                                                className="graph-card"
                                            ),
                                        ], className='plot-and-controls'),
                                    ], className='plot-and-controls'),
                                ]),
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                html.H4('Original Input Path (For Partial Diffusion comparison)'),
                                                dcc.Input(
                                                    id='input_pdb_path',
                                                    type='text',
                                                    placeholder='Input PDB file',
                                                    className="input-box"
                                                )
                                            ],
                                            style={'width':'100%', 'padding':'10px'}
                                        )
                                    ],
                                    style={'display': 'flex', 'flex-direction': 'row', 'align-items': 'center'}
                                ),
                            ]),
                            dcc.Tab(label='Pipeline Tracking', children=[
                                dbc.Row([
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader(id='job-status-counts', className='card-header-primary'),
                                            dbc.ListGroup(id='job-status-list', flush=True, className='mb-4', style={'width':'100%', 'overflowY': 'auto', 'maxHeight': '750px'}),
                                            dcc.Interval(id='interval-component', interval=60000, n_intervals=0),
                                            dcc.Store('filtered_df')], className='outfile_cards'
                                            ),
                                            width=4
                                    ),
                                    dbc.Col([
                                        dbc.Row([
                                            dbc.Col(
                                                dbc.Card([
                                                    dbc.CardHeader("Outfile Viewer", className='card-header-primary'),
                                                    dbc.CardBody([
                                                        html.Div(id='log-viewer', className='outfile_viewer', style={'maxHeight':'325px'}),
                                                        html.Div(id='log-viewer-content', style={'display': 'none'})
                                                ]),
                                                ], className='outfile_card',style={'height': 400})
                                            )                                      
                                        ]),
                                        dbc.Row([
                                            dbc.Col(
                                                dbc.Card([
                                                    dbc.CardHeader("Errorfile Viewer", className='card-header-primary'),
                                                    dbc.CardBody([
                                                        html.Div(id='error-viewer',className='outfile_viewer', style={'maxHeight':'325px'}),
                                                        html.Div(id='error-viewer-content', style={'display': 'none'})
                                                    ]),
                                                ], className='outfile_card',style={'height': 400})
                                            )
                                        ]),
                                    ], width=8)
                                ])
                            ]), 
                            dcc.Tab(label='Extraction', children=[
                                # Hit Viewer & Selection
                                dbc.Row([
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader("Hit PDB preview", className='card-header-primary'),
                                            dbc.CardBody([
                                                dcc.Dropdown(
                                                    options=[],
                                                    value=None,
                                                    placeholder='Select a hit',
                                                    id='extractions_molecule_dropdown',
                                                    className='dropdown'
                                                ),
                                                dashmolstar.MolstarViewer(
                                                                        id="extractions_molecule",
                                                                        style={'height': '600px', 'width': '100%'},
                                                                        )
                                            ])
                                        ], className='mb-4'),
                                        width=7
                                    ),
                                    dbc.Col(
                                        dbc.Card([
                                            dbc.CardHeader("Selection of hits for extraction", className='card-header-primary'),
                                            dbc.CardBody([
                                                dcc.Checklist(
                                                    options=[],
                                                    value=[],
                                                    id='extraction-selection',
                                                    style={'display':'flex','flexDirection':'column','gap':'5px', 'width':'100%','maxHeight':'600px', 'overflowY':'auto'},
                                                )
                                            ])
                                        ], className='mb-4'),
                                        width=5
                                    )
                                ]),
                                # Extraction Options Card
                                dbc.Card([
                                    dbc.CardHeader("Extraction Model", className='card-header-primary'),

                                    dbc.CardBody([
                                        dbc.Row([
                                            dbc.Col([
                                                dcc.Store(id='selected-model', data='PDB'),
                                                dbc.DropdownMenu(label='Select Model',
                                                                    children=[
                                                                        dbc.DropdownMenuItem("PDB", id={'type':'dna-model', 'name':'PDB'}, n_clicks=0),
                                                                        dbc.DropdownMenuItem("CodonTransformer", id={'type':'dna-model', 'name':'CT'}, n_clicks=0)
                                                                    ],
                                                                    id='dna-model-dropdown'),
                                                html.Div( id='model-specific-options'), # Division to put model-specific options
                                            ], width=6),
                                            dbc.Col([
                                                dbc.Card([
                                                    dbc.CardHeader("Extraction outfile"),
                                                    dbc.CardBody([
                                                        html.Div(id='extraction-viewer'),
                                                    ])
                                                ], className='outfile_viewer', style={'height':'100%'})
                                            ], width=6)
                                        ]),
                                        html.Div(
                                            dbc.Button('Extract Hits', id='execute-hits', n_clicks=0, className='button'),
                                            style={'textAlign':'center', 'margin-top':'10px'}
                                        )
                                    ])
                                ], className='mb-4')
                            ]),
                        ])
                    ])
                ]
            ),
            # Offcanvas removed
        ]
    )

app.layout = serve_layout

# Callback for molecule representation
@callback(
    Output("extractions_molecule", "data"),
    [
    Input("extractions_molecule_dropdown", "value"),
    Input("directory-dropdown", "value"),
    Input("input_pdb_path", "value")
    ]
)
def return_molecule(value, directory, input_path):
    if not value or not directory:
        raise PreventUpdate

    data_path, filename = get_design_file_path_and_name(value, directory, input_path)


    file_path = os.path.join(data_path, filename + ".pdb")


    print(f"Loading molecule from: {file_path}")

    chainA = molstar_helper.get_targets(chain="A")
    chainB = molstar_helper.get_targets(chain="B")

    chainA_representation = Representation(type='cartoon', color="uniform")
    chainA_representation.set_color_params({'value': 15577217})
    chainB_representation = Representation(type='gaussian-surface', color="uniform")
    chainB_representation.set_type_params({'alpha': 1})
    chainB_representation.set_color_params({'value': 8815233})

    component_A = molstar_helper.create_component("Chain A", chainA, chainA_representation)
    component_B = molstar_helper.create_component("Chain B", chainB, chainB_representation)

    data = molstar_helper.parse_molecule(file_path, "pdb", component=[component_A, component_B], preset={'kind': 'empty'})



    return data

# Callback to update graphs and table
@callback(
    Output('scatter-plot', 'figure'),
    Output('row-count', 'children'),
    Output('hit-count', 'children'),
    Output('hit-efficiency', 'children'),
    Output('job-status-counts', 'children'),
    Output('extractions_molecule_dropdown', 'options'),
    Output('filtered_df', 'data'),
    Output('extraction-selection', 'options'),
    Output('extraction-selection', 'value'),
    [
        Input('directory-dropdown', 'value'),
        Input('interval-component', 'n_intervals'),
        Input('stop-campaign', 'n_clicks'),
        Input('xaxis_value', 'value'),
        Input('yaxis_value','value'),
        Input('input_pdb_path', 'value'),
        Input('pae_interaction_thres', 'value'),
        Input('CUTRE_thres', 'value'),
        Input('plddt_binder_thres', 'value'),
        Input('dsasa_thres', 'value'),
        Input('shape_complementarity_thres', 'value'),
        Input('interface_hbond_thres', 'value'),
        Input('interface_unsat_hbond_thres', 'value'),
        Input('binder_surf_hyd_thres', 'value'),
    ]
)

def update_graph(working_dir, n, n_clicks_stop, xaxis_value, yaxis_value, input_pdb_path, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    directory = f'{working_dir}/output/'
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(working_dir, input_pdb_path)
    filtered_df = filtering_df(merged_df, pae_interaction_thres, CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres)
    status_df = pd.DataFrame()
    status_df = track_job_status(directory)
    status_df = status_df
    if not status_df.empty:
        status_df_records = status_df.to_dict('records')
    else:
        status_df_records = []

    # Make scatter plot
    scatter_plot, row_count_text = update_scatter_plot(working_dir, merged_df, filtered_df, xaxis_value, yaxis_value, input_pdb_path)

    ctx = dash.callback_context

    # STOP CAMPAIGN
    if ctx.triggered and ctx.triggered[0]['prop_id'] == 'stop-campaign.n_clicks':
        command = f'touch {working_dir}/campaign_done'
        subprocess.run(command, shell=True)

    # Job status counts
    status_counts = status_df['status'].value_counts().to_dict()
    status_counts_formatted = ", ".join([f"{status}: {count}" for status, count in status_counts.items()])
    job_status_counts_text = f"Job Status Counts: {status_counts_formatted}"
    # Dropdown HITS
    dropdown_options = get_hit_names(filtered_df, xaxis_value)

    # jasonified df for communication between callbacks
    jasonified_df = filtered_df.to_json(date_format='iso', orient='split')

    # Metrics for cards
    finished_models = len(merged_df)
    hit_count = len(dropdown_options)
    if dropdown_options == ["No hits found using current filters"]:
        hit_count = 0
        hit_efficiency = "0%"
    elif finished_models > 0:
        hit_efficiency = f"{(hit_count/finished_models*100):.1f}%"
    else:
        hit_efficiency = "N/A"

    return (
        scatter_plot,
        row_count_text,
        hit_count,
        hit_efficiency,
        job_status_counts_text,
        dropdown_options,
        jasonified_df,
        dropdown_options,
        dropdown_options
    )
# callback for job-status-list
@callback(
    Output('job-status-list', 'children'),
    Input('interval-component', 'n_intervals'),
    Input('directory-dropdown', 'value')
)
def update_job_list(n, working_dir):
    color = {
        'Waiting': '#898952',
        'WAITING': '#B6D369',
        'Finished':'#93C48B',
        'FAILED':'dangerous' 
    }

    df = track_job_status(f'{working_dir}/output/')
    # Extract run number and sort descending
    df_sorted = df.copy()
    df_sorted['run_number'] = df_sorted['job'].str.extract(r'(\d+)', expand=False).astype(int)
    df_sorted = df_sorted.sort_values(['run_number','gpu'], ascending=[False,True ])
    items = []
    for job in df_sorted['job'].unique():
        items.append(
            dbc.ListGroupItem(
                dbc.Row([
                    dbc.Col(
                        html.Span(job, className='me-2', style={'fontWeight': '600'}),
                    width=4),
                    dbc.Col(
                            [
                                dbc.Row(
                                    dbc.DropdownMenu(
                                        label=f'GPU {gpu}: \t {df_sorted["status"][(df_sorted["job"] == job) & (df_sorted["gpu"] == gpu)].values[0]}',
                                        children=[
                                            dbc.DropdownMenuItem('Backbone',  id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/rfd"}, n_clicks=0),
                                            dbc.DropdownMenuItem('Filtering', id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/aligning_filtering"}, n_clicks=0),
                                            dbc.DropdownMenuItem('Sequence',  id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/pmpnn"}, n_clicks=0),
                                            dbc.DropdownMenuItem('Scoring',   id={'type':'dropdown-item', 'value':f"{df_sorted['path'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]}/scoring"}, n_clicks=0),
                                        ], color=color[df_sorted['status'][(df_sorted['job'] == job) & (df_sorted['gpu'] == gpu)].values[0]],
                                    )
                                )
                                for gpu in df_sorted['gpu'][df_sorted['job'] == job].unique()
                            ])
                ])
            )
        )
    return items
#Callback to update the log and error viewer
@callback(
    Output('log-viewer', 'children'),
    Output('error-viewer', 'children'),
    Input({"type":'dropdown-item', "value":dash.ALL}, 'n_clicks'),
    prevent_initial_call=True
)

def update_log_viewer(n_clicks_list):
    ctx = dash.callback_context
    if not ctx.triggered:
        return " ", " "
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    value= eval(triggered_id)['value']

    outfile_path = f"{value}.out"
    errorfile_path = f"{value}.err"
    try:
        with open(outfile_path, 'r') as outfile:
            outfile_content = outfile.read()
    except FileNotFoundError:
        outfile_content = f'{outfile_path} cannot be accessed'
    try:
        with open(errorfile_path, 'r') as errorfile:
            errorfile_content = errorfile.read()
    except FileNotFoundError:
        errorfile_content = f'{errorfile_path} cannot be accessed'

    return outfile_content, errorfile_content

@callback(
    Output('model-specific-options', 'children'),
    Output('selected-model', 'data'),
    Output('dna-model-dropdown', 'label'),
    Input({'type':'dna-model', 'name': dash.ALL}, 'n_clicks'),
)

# Model specific options
# This is a little bit of a hell, probably we should think on moving it into other script and then exporting it
def update_model_specific_options(n_clicks_list):
    ctx = dash.callback_context
    if not ctx.triggered:
        selected_model='PDB'
        return (html.Div(), 'PDB', 'PDB')
    else:
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        selected_model = eval(triggered_id)['name']


    if selected_model == 'CT':
        return (html.Div([
                html.Br(),
                html.Div([
                    html.Span("CodonTransformer Options", className='model-options-title')
                ], style={'textAlign':'center'}),
                html.Br(),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Add Initial Methionine"),
                            dbc.CardBody(
                                dbc.RadioItems(
                                    [
                                        {'label': 'True', 'value': True},
                                        {'label': 'False', 'value': False}
                                    ],
                                    value=False,
                                    id={'type':'param', 'name':'add-met'},
                                    inputStyle={'margin-right':'5px'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    ),
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Organism"),
                            dbc.CardBody(
                                dcc.Dropdown(
                                    initial_organisms,
                                    'Escherichia coli general',
                                    id={'type':'param', 'name':'organism'},
                                    className='dropdown'
                                ), style={'width': '100%'}
                            )
                        ], className='mb-3'),
                        width=6
                    )
                ]),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("3' Overhang"),
                            dbc.CardBody(
                                dcc.Input(
                                    placeholder='Sequence',
                                    value='',
                                    id={'type':'param', 'name':'three_prime_overhang'},
                                    className='input-box',
                                    style={'width': '100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    ),
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("5' Overhang"),
                            dbc.CardBody(
                                dcc.Input(
                                    placeholder='Sequence',
                                    value='',
                                    id={'type':'param', 'name':'five_prime_overhang'},
                                    className='input-box',
                                    style = {'width': '100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    )
                ]),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Minimum Sequence Length"),
                            dbc.CardBody(
                                dcc.Input(
                                    type='number',
                                    value=300,
                                    id={'type':'param', 'name':'random_sequence'},
                                    className='input-box',
                                    style={'width':'100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    ),
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("GC % of Random Sequence"),
                            dbc.CardBody(
                                dcc.Input(
                                    type='number',
                                    value=50,
                                    id={'type':'param', 'name':'GC_content'},
                                    className='input-box',
                                    style={'width':'100%'}
                                )
                            )
                        ], className='mb-3'),
                        width=6
                    )
                ]),
                dbc.Row([
                    dbc.Col(
                        dbc.Card([
                            dbc.CardHeader("Check Restriction Sites"),
                            dbc.CardBody(
                                dcc.Checklist(
                                    options=[{'label': enz, 'value': enz} for enz in [
                                        'EcoRI','BamHI','HindIII','NotI','XhoI','PstI','SacI','KpnI',
                                        'SmaI','XbaI','SpeI','NcoI','SalI','ApaI','HaeIII','AluI',
                                        'TaqI','BglII','ClaI','MluI','BsaI'
                                    ]],
                                    value=[],
                                    id={'type':'param', 'name':'enzyme'},
                                    inputStyle={'margin-right':'5px'},
                                    className='enzyme-grid'
                                )
                            )
                        ], className='mb-3', style={'overflowY':'auto', 'maxHeight':'300px'}),
                        width=12
                    )
                ]),
            ]), 'CT', 'CodonTransformer')
    elif selected_model == 'PDB':
        return (html.Div(), 'PDB', 'PDB')


#callback to update the radar plot
@callback(
    Output('radar-plot', 'figure'),
    [
    Input('scatter-plot', 'clickData'),
    Input('interval-component', 'n_intervals'),
    Input('directory-dropdown', 'value'),
    Input('input_pdb_path', 'value'),
    Input('pae_interaction_thres', 'value'),
    Input('CUTRE_thres', 'value'),
    Input('plddt_binder_thres', 'value'),
    Input('dsasa_thres', 'value'),
    Input('shape_complementarity_thres', 'value'),
    Input('interface_hbond_thres', 'value'),
    Input('interface_unsat_hbond_thres', 'value'),
    Input('binder_surf_hyd_thres', 'value'),
    ])

#This prints the radar plot
def update_radar_plot(design_to_plot, n, directory, input_pdb_path, pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres):
    #This is meant to store the original designs for the following plotting; is a little bit cutre 
    update_designs_list(designs_list, design_to_plot)
    merged_df = pd.DataFrame()
    merged_df = merge_csv_files(directory, input_pdb_path)
    #Most of the heavy work is carry in the utils file, go there for further info
    radar_figure=radar_plot(designs_list, merged_df,pae_interaction_thres,CUTRE_thres, plddt_binder_thres, dsasa_thres, shape_complementarity_thres, interface_hbond_thres, interface_unsat_hbond_thres, binder_surf_hyd_thres )
    # Override radar area colors with transparency
    fill_colors = ['rgba(134,48,113,0.5)', 'rgba(237,176,129,0.5)']
    line_colors = ['#863071', '#edb081']
    for i, trace in enumerate(radar_figure.data):
        if hasattr(trace, 'fillcolor'):
            trace.fillcolor = fill_colors[i % len(fill_colors)]
            trace.line.color = line_colors[i % len(line_colors)]
        elif hasattr(trace, 'marker'):
            trace.marker.color = line_colors[i % len(line_colors)]
    return radar_figure

# Callback to extract hits
@callback(
    Output('extraction-viewer', 'children'),
    [
        Input('directory-dropdown', 'value'),
        Input('interval-component', 'n_intervals'),
        Input('execute-hits', 'n_clicks'),
        Input({'type':'param', 'name':dash.ALL}, 'value'),
        Input('extraction-selection', 'value')
    ],
    State('selected-model', 'data'),
)

def extract_hits(working_dir, n, clicks, param_values, extraction_list, model_clicks):

    #Create the hits folder
    hits_folder = os.path.join(working_dir, 'hits')
    fastas_folder = os.path.join(hits_folder, 'fastas')
    dnaseq_folder=os.path.join(hits_folder, 'dna_seqs')
    pdbs_folder = os.path.join(hits_folder, 'pdbs')
    if not os.path.exists(hits_folder):
        os.makedirs(hits_folder)
    if not os.path.exists(fastas_folder) or not os.path.exists(dnaseq_folder) or not os.path.exists(pdbs_folder):
        os.makedirs(fastas_folder, exist_ok=True)
        os.makedirs(dnaseq_folder, exist_ok=True)
        os.makedirs(pdbs_folder, exist_ok=True)

    dnaseq_folder = dnaseq_folder + '/'


    #Name of the multifastas output
    multifastas_output = 'all_sequences.fasta'
    multifastas_path = os.path.join(fastas_folder, multifastas_output)

    ctx = dash.callback_context

    if ctx.triggered_id != "execute-hits":
        raise PreventUpdate

    else:
        #Getting the list of clicks
        model_states = ctx.states_list[0]

        chosen_model = model_states['value']

        print('CHOSEN MODEL:', chosen_model )

        # Getting all the parameters in a dict object, with name as key and value as value

        inputs_with_ids = ctx.inputs_list

        # Generator for dictionary
        params = (
            (str(param['id']['name']), param['value'])
            for item in inputs_with_ids
            if isinstance(item, list)
            for param in item
            if 'id' in param and 'name' in param['id'] and 'type' in param['id'] and 'param' in param['id']['type'] and 'value' in param
        )

        # Convert generator to dict when needed
        params_dict = dict(params)
        # Whatever has been selected, the pdbs are extracted
        
        if extraction_list:
            for description in extraction_list:
        
                # Get the PDBs

                extract_pdbs(description,working_dir,pdbs_folder)
        
                #Get the fastas of the hits
        
                extract_fasta_seq(os.path.join(pdbs_folder, description), fastas_folder)
        
        # Make the multifasta file

        multifastas(fastas_folder, multifastas_output)

        # Extract the DNA sequence based on the model selected
        
        if chosen_model == 'CT':
            stdout = extract_dna_seq_CT(multifastas_path, dnaseq_folder, params_dict)
        
        else:
            stdout = 'PDB OF THE HITS EXTRACTED. LIST OF PDBS EXTRACTED:\n' \
                     f'{os.listdir(pdbs_folder)}\n' \

        # Make a csv file with the hits metrics
        
        metrics_df = pd.read_csv(os.path.join(working_dir, 'Scoring_Stats.csv'))
        filtered_df = metrics_df[metrics_df['description'].isin(extraction_list)]
        filtered_df.to_csv(f'{working_dir}/hits/hits_metrics.csv', index=False)
    
    return stdout

# Callback to toggle the Collapse for filters and axes in the sidebar
@callback(
    Output('filters-collapse', 'is_open'),
    Input('open-filters', 'n_clicks'),
    State('filters-collapse', 'is_open')
)

#

def toggle_filters_collapse(n, is_open):
    if n:
        return not is_open
    return is_open
# Run the app
if __name__ == '__main__':
    app.run(debug=False, dev_tools_hot_reload = False, use_reloader=True,
                   host=hostname, port=port_number)


