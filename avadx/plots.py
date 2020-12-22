from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.offline import plot
from .logger import Logger


LOG_LEVEL = 'INFO'
LOG_FILE = None


class Figure:

    def __init__(self, vis, dataset, genescore_normalize, wd_folder, out_folder):

        self.vis = vis
        self.dataset = dataset
        self.genescore_normalize = genescore_normalize
        self.in_file = None
        self.out_file = None
        self.data = None
        self.log = self.get_logger()

        if vis == 'heatmap_clust_dend':
            if dataset == 'genescores':
                self.in_file = out_folder / Path(f'genescores/GeneScoreTable_{"normalized" if self.genescore_normalize else "raw"}_variation_filtered.csv')
                selector_file = out_folder / Path('results/model_features.txt')
                cv_scheme_path = wd_folder / 'tmp' / 'cv-scheme.csv'
                self.out_file = out_folder / Path('results/genescores_heatmap.html')
                if (self.in_file and self.in_file.exists()) and (selector_file and selector_file.exists()):
                    data_raw = pd.read_csv(self.in_file)
                    cv_scheme = pd.read_csv(cv_scheme_path, header=None, names=['sampleid', 'fold', 'status'])
                    data_raw = data_raw.merge(cv_scheme[['sampleid', 'status']], how='left', on=['sampleid'])
                    data_raw = data_raw.T
                    selector = pd.read_csv(selector_file, header=None, names=['gene', 'count'])
                    data_subset = data_raw.loc[data_raw.index.isin(selector.gene.values)]
                    self.data = data_subset.append(data_raw.loc['status']).T
                else:
                    print(f'Data file not found: {self.in_file}')
        elif vis == 'piechart':
            if dataset == 'ethseq':
                reports_path = out_folder / 'reports'
                self.in_file = reports_path / '2_2-EthSEQ_summary_Report.csv'
                self.out_file = reports_path / '2_2-EthSEQ_summary_Report.html'
                report_splits = []
                for report_split in reports_path.glob('2_2-EthSEQ_split_*_Report.txt'):
                    report_splits += [pd.read_csv(report_split, sep='\t')]
                df = pd.concat(report_splits)
                df.to_csv(self.in_file, index=False)
                self.data = df

    def get_logger(self):
        logger = Logger(self.__class__.__name__, level=LOG_LEVEL)
        logger.addConsoleHandler()
        log = logger.getLogger()
        return log

    def run(self, func_name):
        func = getattr(self, func_name, None)
        if func:
            return func()
        else:
            self.log.warning(f'No visualization: {func_name}')

    def create(self):
        try:
            import scipy  # noqa: F401
            self.run(self.vis)
        except ImportError:
            self.log.error('|6.00| Visualizations require the "scipy" package: pip install scipy')

    def save(self, fig):
        with self.out_file.open('w') as fout:
            if self.out_file.suffix == '.html':
                fout.write(plot(fig, output_type='div'))

    def set_tight_margin(self, fig):
        fig.update_layout(margin=dict(t=0, l=0, r=0, b=0))

    def heatmap_simple(self):
        fig = go.Figure(
            data=go.Heatmap(
                z=self.data,
                x=self.data.columns,
                y=self.data.index,
                colorscale='Viridis'
            )
        )
        fig.update_xaxes(side="top")
        self.save(fig)

    def heatmap_clust_dend(self):
        df = self.data.copy()
        class_labels = df.status.values
        df.drop('status', axis=1, inplace=True)
        df_array = df.to_numpy()
        genes = df.columns.values
        samples = df.index.values

        # Initialize figure by creating upper dendrogram
        fig = ff.create_dendrogram(df_array.T, orientation='bottom', labels=genes)
        dendro_leaves_genes = list(fig['layout']['xaxis']['ticktext'])
        for i in range(len(fig['data'])):
            fig['data'][i]['yaxis'] = 'y2'

        # Create Side Dendrogram
        dendro_side = ff.create_dendrogram(df_array, orientation='right')
        dendro_leaves_samples = dendro_side['layout']['yaxis']['ticktext']
        c = [int(_) for _ in dendro_leaves_samples]
        dendro_side['layout']['yaxis']['ticktext'] = samples[c]
        tickvals_samples = np.asarray(dendro_side['layout']['yaxis']['tickvals'])
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'

        # Add Side Dendrogram Data to Figure
        for data in dendro_side['data']:
            fig.add_trace(data)

        # Create Heatmap
        dendro_leaves_samples = list(map(int, dendro_leaves_samples))
        df = df.iloc[dendro_leaves_samples, :]
        df = df.loc[:, dendro_leaves_genes]
        heatmap = [
            go.Heatmap(
                x=dendro_leaves_genes,
                y=dendro_leaves_samples,
                z=df,
                colorscale='Viridis'
            )
        ]

        heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
        heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

        # Add Heatmap Data to Figure
        for data in heatmap:
            fig.add_trace(data)

        # Edit Layout
        fig.update_layout({
            'width': 1400,
            'height': 1000,
            'showlegend': False,
            'hovermode': 'closest',
        })

        # Edit xaxis
        fig.update_layout(
            xaxis={
                'domain': [.15, 1],
                'mirror': False,
                'showgrid': False,
                'showline': False,
                'zeroline': False,
                'ticks': '',
            }
        )

        # Edit xaxis2
        fig.update_layout(
            xaxis2={
                'domain': [0, .15],
                'mirror': False,
                'showgrid': False,
                'showline': False,
                'zeroline': False,
                'showticklabels': False,
                'ticks': '',
            }
        )

        # Edit yaxis
        fig.update_layout(
            yaxis={
                'domain': [0, .85],
                'mirror': False,
                'showgrid': False,
                'showline': False,
                'zeroline': False,
                'showticklabels': True,
                'ticks': "outside",
                'ticktext': [f'{s} [{int(c)}]' for s, c in zip(samples[c], np.array(class_labels)[c])],
                'tickvals': tickvals_samples,
                'ticklen': 200,
                'tickcolor': 'rgba(255,255,255,0)',
            }
        ),

        # Edit yaxis2
        fig.update_layout(
            yaxis2={
                'domain': [.825, .975],
                'mirror': False,
                'showgrid': False,
                'showline': False,
                'zeroline': False,
                'showticklabels': False,
                'ticks': ''
            }
        )

        # Save
        self.save(fig)

    def piechart(self):
        self.data['share'] = self.data.groupby('type')['type'].transform('count')
        fig = px.sunburst(
            self.data, path=['pop', 'type'], color='pop',
        )
        fig.update_layout(
            title_text='EthSEQ ethnicity report'
        )
        self.save(fig)
