import plotly.graph_objs as go
import plotly.io as pio

from pamogk import config


def plot(nx_g, title, auto_open=False):
    # import pdb; pdb.set_trace()
    nodes = [nx_g.nodes[n] for n in nx_g.nodes()]
    marker = None
    types = [n['type'] for n in nodes if 'type' in n]
    # if all nodes have type show color depending on type
    if len(types) == len(nodes):
        types = list(set(types))
        print('Type set for plot:', types)
        marker = go.scatter.Marker(
            color=[types.index(n['type']) for n in nodes],
            colorscale='Viridis',
            showscale=True)
    fig = go.Figure(
        data=[
            go.Scatter(
                x=[n['x'] for n in nodes],
                y=[n['y'] for n in nodes],
                mode='markers+text',
                text=[n['n'] + '-' + n['type'] for n in nodes],
                textposition='bottom center',
                marker=marker)
        ],
        layout=go.Layout(
            title=title,
            hovermode='closest',
            xaxis=dict(
                ticklen=5,
                zeroline=False,
                gridwidth=2,
            ),
            yaxis=dict(
                ticklen=5,
                gridwidth=2,
            ),
            showlegend=False,
            width=1200,
            height=900,
        ))
    pio.write_html(fig,
                   config.DATA_DIR / f'{title}.html',
                   include_plotlyjs='directory',  # adds plotly.js file separately
                   auto_open=auto_open)
