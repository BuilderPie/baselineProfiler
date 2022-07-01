import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.io as pio
import plotly.express as px

# define plotly_heatmap function
def plotly_heatmap(df, colorbar_title = "log2(TPM)"):
    # get data
    # df = pd.read_csv('../RNA_seq_V3/plots/heatmap/ccle/heatmap_CCLE_Blood.csv', index_col = 0)
    data_array = df.to_numpy()
    labels_y = df.index.to_list()
    labels_x = df.columns.to_list()

    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(data_array.T, orientation='bottom', labels=labels_x)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'


    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array, orientation='right', labels=labels_y)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'


    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)


    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves_ind = list(map(labels_y.index, dendro_leaves))

    fig_leaves = fig['layout']['xaxis']['ticktext']
    fig_leaves_ind = list(map(labels_x.index, fig_leaves))

    heat_data = data_array.copy()
    heat_data = heat_data[dendro_leaves_ind, :]
    heat_data = heat_data[:, fig_leaves_ind]

    heatmap = [
        go.Heatmap(
            x = fig_leaves_ind,
            y = dendro_leaves_ind,
            z = heat_data,
            colorscale = 'RdYlBu',
            reversescale=True,
            opacity = 0.9,
            colorbar = dict(len = 0.4, lenmode = 'fraction', title = colorbar_title)
        )
    ]


    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)

    # calculate fig height
    if heat_data.shape[0] <= 30:
        fig_height = 400+heat_data.shape[0]*7
    else:
        fig_height = 200+heat_data.shape[0]*7

    # Edit Layout
    fig.update_layout({
                        # 'width':400+heat_data.shape[1]*25, 
                        'height':fig_height,
                        'autosize':True,
                        'showlegend':False, 'hovermode': 'closest',
                        'plot_bgcolor': 'rgba(0,0,0,0)',
                        'paper_bgcolor': 'rgba(0,0,0,0)',
                        'margin_t': 20
                        })
    # Edit xaxis
    fig.update_layout(xaxis={'domain': [.15, 1],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'ticks':""})
    # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .145],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'showticklabels': False,
                                    'ticks':""})
    # add y-axis labels to heatmap
    fig['layout']['yaxis']['ticktext'] = dendro_side['layout']['yaxis']['ticktext']
    fig['layout']['yaxis']['tickvals'] = dendro_side['layout']['yaxis']['tickvals']

    # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .85],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'showticklabels': False,
                                    'ticks': ""
                            })
    # Edit yaxis2
    fig.update_layout(yaxis2={'domain':[.85, .975],
                                    'mirror': False,
                                    'showgrid': False,
                                    'showline': False,
                                    'zeroline': False,
                                    'showticklabels': False,
                                    'ticks':""})

    # reverse y-axis
    fig['layout']['yaxis']['autorange'] = "reversed"

    # Plot!
    # fig.show()

    # return the fig object
    return fig


# define plotly_heatmap function without dendrogram
def plotly_heatmap_wo_dendrogram(df, colorbar_title = "log2(TPM)"):
    # get data
    # df = pd.read_csv('../RNA_seq_V3/plots/heatmap/ccle/heatmap_CCLE_Blood.csv', index_col = 0)
    data_array = df.to_numpy()
    labels_y = df.index.to_list()
    labels_x = df.columns.to_list()

    fig = go.Figure(data=go.Heatmap(
        z=data_array,
        x=labels_x,
        y=labels_y,
        colorscale = 'RdYlBu',
        reversescale=True,
        opacity = 0.9,
        colorbar = dict(len = 0.4, lenmode = 'fraction', title = colorbar_title)
        )
    )

    # fig.update_layout(
    # title='GitHub commits per day',
    # xaxis_nticks=36)

    # # Add Heatmap Data to Figure
    # for data in heatmap:
    #     fig.add_trace(data)

    # calculate fig height
    if data_array.shape[0] <= 30:
        fig_height = 400+data_array.shape[0]*7
    else:
        fig_height = 200+data_array.shape[0]*7

    # Edit Layout
    fig.update_layout({
                        # 'width':400+heat_data.shape[1]*25, 
                        'height':fig_height,
                        'autosize':True,
                        'showlegend':False, 'hovermode': 'closest',
                        'plot_bgcolor': 'rgba(0,0,0,0)',
                        'paper_bgcolor': 'rgba(0,0,0,0)',
                        'margin_t': 20
                        })
    # # Edit xaxis
    # fig.update_layout(xaxis={'domain': [.15, 1],
    #                                 'mirror': False,
    #                                 'showgrid': False,
    #                                 'showline': False,
    #                                 'zeroline': False,
    #                                 'ticks':""})
    # # Edit xaxis2
    # fig.update_layout(xaxis2={'domain': [0, .145],
    #                                 'mirror': False,
    #                                 'showgrid': False,
    #                                 'showline': False,
    #                                 'zeroline': False,
    #                                 'showticklabels': False,
    #                                 'ticks':""})
    # # add y-axis labels to heatmap
    # fig['layout']['yaxis']['ticktext'] = dendro_side['layout']['yaxis']['ticktext']
    # fig['layout']['yaxis']['tickvals'] = dendro_side['layout']['yaxis']['tickvals']

    # # Edit yaxis
    # fig.update_layout(yaxis={'domain': [0, .85],
    #                                 'mirror': False,
    #                                 'showgrid': False,
    #                                 'showline': False,
    #                                 'zeroline': False,
    #                                 'showticklabels': False,
    #                                 'ticks': ""
    #                         })
    # # Edit yaxis2
    # fig.update_layout(yaxis2={'domain':[.85, .975],
    #                                 'mirror': False,
    #                                 'showgrid': False,
    #                                 'showline': False,
    #                                 'zeroline': False,
    #                                 'showticklabels': False,
    #                                 'ticks':""})

    # reverse y-axis
    # fig['layout']['yaxis']['autorange'] = "reversed"

    # Plot!
    # fig.show()

    # return the fig object
    return fig

# plotly_heatmap(df = pd.read_csv('../RNA_seq_V3/plots/heatmap/heatmap_CCLE_Blood.csv', index_col = 0))