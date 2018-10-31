import plotly.offline as pyoff
import plotly.plotly as py
import plotly.graph_objs as go

def plot(nx_G):
    colors = list(set(nx_G[n]['type'] for n in nx_G.nodes()))
    fig = go.Figure(
        data = [
            go.Scatter(
            x=[nx_G[n]['x'] for n in nx_G.nodes()],
            y=[nx_G[n]['y'] for n in nx_G.nodes()],
            mode='markers+text',
            text=[nx_G[n]['n'] for n in nx_G.values()],
            textposition='bottom center',
            marker=go.scatter.Marker(color=[colors.index(n['type']) for n in nx_G.nodes()], colorscale='Viridis', showscale=True))
        ],
        layout = go.Layout(
            title= 'TSNE Representation of {} p={:0.2f} q={:0.2f} run={}'.format(pathway_id, args.p, args.q, args.rid),
            hovermode= 'closest',
            xaxis= dict(
                ticklen= 5,
                zeroline= False,
                gridwidth= 2,
            ),
            yaxis = dict(
                ticklen=5,
                gridwidth=2,
            ),
            showlegend=False,
            width=1200,
            height=900,
        ))
    pyoff.plot(fig)
