import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
from pathlib import Path
from multiprocessing import Pool
from joblib import dump, load
from sklearn import metrics
from sklearn.feature_selection import VarianceThreshold
from genescore import ScoringFunction


def flatten(list):
    return [item for sublist in list for item in sublist]


def get_fselection(featureselection, kwargs_dict, cvscheme, maxgenes, progress=False):
    from feature_selections.fselection_avadx import FSelection
    return FSelection(featureselection, kwargs_dict, cvscheme, max_genes=maxgenes, progress=not progress)


def get_model(model, kwargs_dict, fselection):
    from models.models_avadx import Model
    return Model(model, kwargs_dict, fselection)


def load_model(filename):
    return load(filename)


def save_model(model, filename):
    dump(model, filename)


def load_features(filename):
    return pd.read_csv(filename, header=None)


def save_features(features, filename):
    features.to_csv(filename, columns=[], header=False)


def get_NA_score(variantfn, genefn):
    scoring_functions = ScoringFunction(variantfn, genefn)
    return scoring_functions.gene.NA_SCORE


def plot_curve(x, y, save_as, color='darkorange', label='Label', x_lab='x', y_lab='y', title='Plot', x_lim=[0.0, 1.0], y_lim=[0.0, 1.05], lw=2, legend_loc='lower right', show_diag_col='navy', show_diag_ln='--', diag_x=[0, 1], diag_y=[0, 1]):
    plt.figure()
    plt.plot(x, y, color=color, lw=lw, label=label)
    if show_diag_col not in [False, None]:
        plt.plot(diag_x, diag_y, color=show_diag_col, lw=lw, linestyle=show_diag_ln)
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title(title)
    plt.legend(loc=legend_loc)
    plt.savefig(save_as)
    plt.close()


def plot_jitter(data_raw, x='class', y='score_1', colors=['#356288', '#fe1100'], jitter=0.2, size=6, alpha=.75, label='Label', x_lab='x', y_lab='y', title='Plot', y_lim=[-0.1, 1.05], legend=False, data_overlay=None, save_as=None):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(x_lab)
    ax1.set_ylabel(y_lab)
    ax2 = ax1.twinx()
    classes = data_raw['class'].unique()
    medians = data_raw.groupby('class').median()[['score_0', 'score_1']].score_1.to_list()
    if len(classes) == 2:
        median_0, median_1 = medians
        mean_of_medians = np.mean([median_0, median_1])
    elif 0 not in classes:
        median_0, median_1 = None, medians[0]
        mean_of_medians = None
    elif 1 not in classes:
        median_0, median_1 = medians[0], None
        mean_of_medians = None
    plt.ylim(y_lim)
    plt.title(title)
    if data_overlay is None:
        data = data_raw.copy()
    else:
        data = data_raw.copy()
        data = data.append(pd.DataFrame(zip(data_overlay.score_1.to_list(), [0.5] * data_overlay.shape[0]), columns=['score_1', 'class']), ignore_index=True)
        colors = [colors[0], '#000000', colors[-1]]
    ax1 = sns.boxplot(x=x, y=y, palette=colors, data=data)
    for patch in ax1.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .3))
    if data_overlay is None:
        ax1 = sns.stripplot(x=x, y=y, data=data, palette=colors, jitter=jitter, size=size, alpha=alpha)
    else:
        ax1 = sns.stripplot(x=x, y=y, data=data, palette=colors, jitter=jitter, size=size, alpha=alpha)
        training_0 = mlines.Line2D([], [], color=colors[0], marker='o', linestyle='None', markersize=size, label='training (0)')
        prediction = mlines.Line2D([], [], color=colors[1], marker='o', linestyle='None', markersize=size, label='prediction')
        training_1 = mlines.Line2D([], [], color=colors[2], marker='o', linestyle='None', markersize=size, label='training (1)')
        if legend:
            box = ax1.get_position()
            ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax1.legend(handles=[training_0, prediction, training_1], loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xticks(ticks=[0, 1, 2], labels=['0', 'prediction', '1'])
    if mean_of_medians is not None:
        ax1.axhline(mean_of_medians, ls='--', color='grey')
    y_ticks = [0]
    y_ticks += [round(median_0, 2)] if median_0 is not None else [round(median_1, 2)]
    y_ticks += [round(mean_of_medians, 2)] if mean_of_medians is not None else [y_ticks[1]]
    y_ticks += [round(median_1, 2)] if median_1 is not None else [round(median_0, 2)]
    y_ticks += [1]
    ax2.set_yticks(y_ticks)
    ax2.set_ylabel('')
    mean_0_ytick, cutoff_ytick, mean_1_ytick = ax2.get_yticklabels()[1:4]
    mean_0_ytick.set_color('grey')
    mean_1_ytick.set_color('grey')
    cutoff_ytick.set_fontsize(10)
    cutoff_ytick.set_color('red')
    cutoff_ytick.set_weight('bold')
    if save_as:
        plt.savefig(save_as)
        plt.close()
    else:
        plt.show()


def main(args, kwargs):
    if args.makepredictions:
        pred_id, model, genescores, features, variantfn, genefn = kwargs
        predict(pred_id, model, genescores, features, variantfn, genefn, args.out)
    else:
        run(args.genescores, args.featureselection, args.featurelist, args.model, args.cvscheme, args.protlength, args.kfold, args.variance, args.variation, args.pvalue, args.maxgenes, args.topgenes, args.stepsize, args.out, args.wd, args.cores, args.quiet, kwargs)


def run(genescores_path, featureselection, featurelist, model, cvscheme_path, protlength_path, kfold, variance_cutoff, variation_cutoff, pval_cutoff, maxgenes, topgenes_ratio, stepsize, out_path, wd_path, cores, quiet, kwargs):
    kwargs_dict = {_.split('=')[0]: _.split('=')[1] for _ in kwargs}
    kwargs_dict['pval_cutoff'] = pval_cutoff
    bar_prefix = '[     INFO ] --- |5.10| '
    logs = []

    # create output directories
    out_path.mkdir(exist_ok=True)
    wd_path.mkdir(exist_ok=True)
    model_path = out_path / 'model.joblib'
    features_path = out_path / 'model_features.txt'
    crossval_bestauc_out = out_path / '1-crossval_genesets_models'
    crossval_bestauc_out.mkdir(exist_ok=True)
    crossval_final_out = out_path / '2-crossval_final_model'
    crossval_final_out.mkdir(exist_ok=True)

    # load datasets
    cvscheme = pd.read_csv(cvscheme_path, header=None, names=['sampleid', 'fold', 'class'], dtype={'sampleid': str, 'fold': int, 'class': int})
    cvscheme = cvscheme.set_index('sampleid')
    kfold_is_set = False if kfold is None else True
    kfold = kfold if kfold_is_set else len(cvscheme.fold.unique())
    kfold_steps = range(1, (cvscheme.shape[0] // kfold)) if kfold_is_set else sorted(cvscheme.fold.unique())
    # NOTE folds auto generation here is deprecated and was moved to the pipeline preprocessing step
    if kfold_is_set:
        folds_new = flatten([[_] * (cvscheme.shape[0] // kfold) for _ in range(1, (cvscheme.shape[0] // kfold))])
        remaining = cvscheme.shape[0] % kfold
        if remaining:
            # NOTE distribution of remaining instances among all folds is not optimally solved here (added to the last fold)
            folds_new += [folds_new[-1]] * remaining
        cvscheme.fold = folds_new

    genescores = pd.read_csv(genescores_path)
    logs += [f'Number of features from GeneScoresTable: {genescores.shape[0]}']
    protlength = pd.read_csv(protlength_path)

    # in case of multiple transcripts for the same gene only keep the longest one
    if sum(genescores.Gene.duplicated()):
        genescores = genescores.merge(protlength, how='left', on=['Gene', 'Transcript'])
        genescores = genescores.loc[genescores.Prot_length.eq(genescores.groupby('Gene').Prot_length.transform(max), fill_value=0)]
        genescores.drop('Prot_length', axis=1, inplace=True)
        genescores.reset_index(drop=True, inplace=True)
    logs += [f'Number of features after removal of duplicate transcripts: {genescores.shape[0]}']

    # select relevant samples and transpose dataframe to create the final dataset
    genescores = genescores[['Gene', 'Transcript'] + cvscheme.index.tolist()]
    dataset = genescores.drop('Transcript', axis=1).set_index('Gene').T.rename_axis('sampleid', axis='rows')

    # variance filtering
    var_selector = VarianceThreshold(variance_cutoff)
    try:
        var_selector.fit(dataset)
        dataset = dataset[dataset.columns[var_selector.get_support(indices=True)]]
    except ValueError as err:
        print(f'Error at variance filter step: {err}')
    logs += [f'Number of features after variance filtering: {dataset.shape[1]}']

    # variation filtering
    gene_ratio_most_common_score = dataset.apply(pd.Series.value_counts, axis=0, normalize=True).apply(pd.Series.max, axis=0)
    gene_variation_subset = gene_ratio_most_common_score[gene_ratio_most_common_score <= (variation_cutoff / 100)].index
    if len(gene_variation_subset):
        dataset = dataset[gene_variation_subset]
    else:
        print(f'Error at variation filter step: no feature met the required threshold of <= {variation_cutoff}%')
    logs += [f'Number of features after variation filtering: {dataset.shape[1]}']

    # save filtered dataset and report
    dataset.to_csv(genescores_path.parent / f'{genescores_path.stem}_variation_filtered.csv')
    with (wd_path / 'genescore_filtering_report.txt').open('w') as fout:
        fout.writelines([f'{msg}\n' for msg in logs])

    # change maxgenes after filtering
    maxgenes = min(maxgenes, dataset.shape[1])

    # add class labels
    dataset = dataset.merge(cvscheme['class'], how='left', on='sampleid')

    # run feature selection or load list of features (genes) to use
    if featurelist and featurelist.exists():
        fselection_df = pd.read_csv(featurelist, header=None, index_col=0)
        index_overlap = fselection_df.index.intersection(dataset.columns)
        fselection_df = fselection_df.loc[index_overlap]
        genes_manual = fselection_df.shape[0]
        fselection = get_fselection('manual', kwargs_dict, cvscheme, genes_manual, quiet)
        fselection.selected = {k: fselection_df for k in kfold_steps}
        genes_considered = [genes_manual]
    else:
        fselection = get_fselection(featureselection, kwargs_dict, cvscheme, maxgenes, quiet)
        with Pool(processes=cores) as fselection_pool:
            if not quiet:
                print(f'progress:start:{len(kfold_steps)}:{bar_prefix}{fselection.name}')
            fselection_pooled = [fselection_pool.apply_async(fselection.fn, args=(dataset, k)) for k in kfold_steps]
            fselection_res = [p.get() for p in fselection_pooled]
            if not quiet:
                print(f'progress:end:{fselection.name}')
            fselection.selected = {genes.name: genes for res, genes in fselection_res}
            fselection_df = pd.DataFrame([res for res, genes in fselection_res]).rename_axis('folds', axis='rows')
            fselection_df.to_csv(wd_path / f'{kfold}F-CV-{fselection.name}-selectedGenes.csv')
            fselection_df.to_csv(crossval_bestauc_out / 'fselection_all_genes_performance_per_fold.csv')
        genes_considered = list(range(stepsize, maxgenes + stepsize, stepsize))
        remaining = maxgenes - genes_considered[-1]
        if remaining:
            genes_considered += [maxgenes + remaining]

    # save order list (descending) of best scoring features (genes) as determimed by Feature Selection step
    fselection_performance = pd.DataFrame(dict([(k, v.index.to_series(index=range(1, v.shape[0] + 1))) for k, v in fselection.selected.items()])).rename_axis('pos', axis='rows')
    fselection_performance.to_csv(crossval_bestauc_out / 'fselection_best_genes_performance_per_fold.csv')
    fselection_performance.to_csv(wd_path / 'fselection_performance.csv')

    # run model training and performance evaluation
    model_eval = get_model(model, kwargs_dict, fselection)
    with Pool(processes=cores) as model_pool:
        predictions_all, performances_roc_data, performances_roc_auc, performances_prc_data, performances_prc_avg = {}, {}, {}, {}, {}
        for max_genes in genes_considered:
            if not quiet:
                print(f'progress:start:{len(kfold_steps)}:{bar_prefix}{model_eval.name}/{model_eval.fselection.name} {max_genes} genes')
            model_pooled = [model_pool.apply_async(model_eval.train, args=(dataset, max_genes, k)) for k in kfold_steps]
            model_res = [p.get() for p in model_pooled]
            model_predictions = pd.DataFrame.from_dict({k: v for d in model_res for k, v in d.items()}, orient='index', columns=['0', '1']).loc[dataset.index]
            y_true, y_scores = dataset['class'], model_predictions['1']
            roc_data = metrics.roc_curve(y_true, y_scores)
            roc_auc = metrics.roc_auc_score(y_true, y_scores)
            prc_data = metrics.precision_recall_curve(y_true, y_scores)
            prc_auc = metrics.average_precision_score(y_true, y_scores)
            predictions_all[max_genes] = model_predictions
            performances_roc_data[max_genes] = roc_data
            performances_roc_auc[max_genes] = roc_auc
            performances_prc_data[max_genes] = prc_data
            performances_prc_avg[max_genes] = prc_auc
            if not quiet:
                print(f'progress:end:{model_eval.name}')
        max_auc_genes = max(performances_roc_auc, key=lambda key: performances_roc_auc[key])
        print(f'|5.10| {model_eval.name}/{model_eval.fselection.name}: best-geneset model ({max_auc_genes} genes) cross-validation performance = {performances_roc_auc[max_auc_genes]:.2f} ROCauc')
        predictions_all_best = pd.DataFrame(predictions_all[max_auc_genes]).rename_axis('sampleid', axis='rows')
        predictions_all_best_samples = predictions_all_best.merge(cvscheme[['class']], how='left', left_index=True, right_index=True).rename(columns={'0': 'score_0', '1': 'score_1'})
        predictions_all_best_samples.to_csv(crossval_bestauc_out / 'best-geneset_model_predictions.csv')
        predictions_all_best_stats = predictions_all_best_samples[['class', 'score_1']].groupby('class').describe()
        predictions_all_best_stats.columns = predictions_all_best_stats.columns.droplevel(0)
        predictions_all_best_stats.to_csv(crossval_bestauc_out / 'best-geneset_model_statistics.csv')
        plot_jitter(predictions_all_best_samples, title='best geneset model predictions (crossval)', x_lab='class', y_lab='prediction score', save_as=out_path / '1-best-geneset_model_predictions.png')
        roc_auc_df = pd.DataFrame.from_dict(performances_roc_auc, orient='index', columns=['AUC']).rename_axis('selected_genes', axis='rows').sort_values(by='AUC', ascending=False)
        roc_auc_df.to_csv(wd_path / f'{kfold}F-CV-{model_eval.name}-performance.csv')
        roc_auc_df.to_csv(crossval_bestauc_out / 'all_performances_ROC-AUC.csv')
        prc_avg_df = pd.DataFrame.from_dict(performances_prc_avg, orient='index', columns=['AVGpr']).rename_axis('selected_genes', axis='rows').sort_values(by='AVGpr', ascending=False)
        prc_avg_df.to_csv(crossval_bestauc_out / 'all_performances_PRC-AVGpr.csv')
        roc_df = pd.DataFrame(performances_roc_data[max_auc_genes], index=['fpr', 'tpr', 'thresholds']).T
        roc_df.to_csv(crossval_bestauc_out / 'best-geneset_model_ROC.csv', index=False)
        plot_curve(roc_df.dropna().fpr, roc_df.dropna().tpr, out_path / '1-best-geneset_model_ROC.png', x_lab='fpr', y_lab='tpr', label=f'Area Under ROC curve (AUC) = {performances_roc_auc[max_auc_genes]:.2f}', title=f'Receiver Operating Characteristic (ROC) curve for top {max_auc_genes} genes [{model_eval.name}/{model_eval.fselection.name}]')
        prc_df = pd.DataFrame(performances_prc_data[max_auc_genes], index=['precision', 'recall', 'thresholds']).T
        prc_df.to_csv(crossval_bestauc_out / 'best-geneset_model_PRC.csv', index=False)
        plot_curve(prc_df.dropna().recall, prc_df.dropna().precision, out_path / '1-best-geneset_model_PRC.png', x_lab='recall', y_lab='precision', label=f'Average precision (AP) = {performances_prc_avg[max_auc_genes]:.2f}', title=f'Precision-Recall curve (PRC) for top {max_auc_genes} genes [{model_eval.name}/{model_eval.fselection.name}]', diag_x=[0, 1], diag_y=[1, 0])

    # get list of selected genes for best AUC over all folds (merge)
    genes_best_merged = {}
    rank1_df = None
    for rank, max_genes in enumerate(roc_auc_df.index, 1):
        genes_best_folds = {}
        for k in kfold_steps:
            genes_selected = fselection.selected[k]
            genes_count = min(max_genes, genes_selected.shape[0])
            for gene in genes_selected.index[0:genes_count]:
                genes_best_folds[gene] = genes_best_folds.get(gene, 0) + 1
        genes_best_merged[max_genes] = genes_best_folds
        genes_best_df = pd.DataFrame.from_dict(genes_best_folds, orient='index', columns=['frequency']).sort_values(by='frequency', ascending=False)
        genes_best_df.to_csv(wd_path / f'AUC_rank.{rank}_top.{max_genes}_{kfold}F-CV-{fselection.name}-selectedGenes.csv')
        if rank == 1:
            rank1_df = genes_best_df
            rank1_df.to_csv(wd_path / 'crossval_genes.csv', columns=[], header=False)
            rank1_df.to_csv(crossval_bestauc_out / 'best-geneset_model_genes.csv', header=False)

    # evaluate final model
    if featurelist and featurelist.exists():
        maxgenes_final = fselection_df.shape[0]
        fselection_final_eval = fselection
        fselection_final_eval.selected = {k: fselection_df for k in kfold_steps}
    else:
        rank1_df_final = rank1_df[((rank1_df / len(kfold_steps)) > topgenes_ratio).frequency]
        maxgenes_final = rank1_df_final.shape[0]
        fselection_final_eval = get_fselection(featureselection, kwargs_dict, cvscheme, maxgenes_final, quiet)
        fselection_final_eval.selected = {k: rank1_df_final for k in kfold_steps}
    model_final_eval = get_model(model, kwargs_dict, fselection_final_eval)
    with Pool(processes=cores) as model_pool:
        if not quiet:
            print(f'progress:start:{len(kfold_steps)}:{bar_prefix}{model_final_eval.name}/{model_final_eval.fselection.name} {maxgenes_final} genes')
        model_final_eval_pooled = [model_pool.apply_async(model_final_eval.train, args=(dataset, maxgenes_final, k)) for k in kfold_steps]
        model_final_eval_res = [p.get() for p in model_final_eval_pooled]
        model_final_eval_predictions = pd.DataFrame.from_dict({k: v for d in model_final_eval_res for k, v in d.items()}, orient='index', columns=['0', '1']).loc[dataset.index]
        y_final_eval_true, y_final_eval_scores = dataset['class'], model_final_eval_predictions['1']
        roc_final_eval_data = metrics.roc_curve(y_final_eval_true, y_final_eval_scores)
        roc_final_eval_auc = metrics.roc_auc_score(y_final_eval_true, y_final_eval_scores)
        prc_final_eval_data = metrics.precision_recall_curve(y_final_eval_true, y_final_eval_scores)
        prc_final_eval_auc = metrics.average_precision_score(y_final_eval_true, y_final_eval_scores)
        if not quiet:
            print(f'progress:end:{model_final_eval.name}')
        print(f'|5.10| {model_final_eval.name}/{model_final_eval.fselection.name}: final model ({maxgenes_final} genes) cross-validation performance = {roc_final_eval_auc:.2f} ROCauc')
        predictions_final_eval = pd.DataFrame(model_final_eval_predictions).rename_axis('sampleid', axis='rows')
        predictions_final_eval_samples = predictions_final_eval.merge(cvscheme[['class']], how='left', left_index=True, right_index=True).rename(columns={'0': 'score_0', '1': 'score_1'})
        predictions_final_eval_samples.to_csv(crossval_final_out / 'final_model_predictions.csv')
        predictions_final_eval_stats = predictions_final_eval_samples[['class', 'score_1']].groupby('class').describe()
        predictions_final_eval_stats.columns = predictions_final_eval_stats.columns.droplevel(0)
        predictions_final_eval_stats.to_csv(crossval_final_out / 'final_model_statistics.csv')
        plot_jitter(predictions_final_eval_samples, title='final model predictions (crossval)', x_lab='class', y_lab='prediction score', save_as=out_path / '2-final_model_predictions.png')
        roc_df = pd.DataFrame(roc_final_eval_data, index=['fpr', 'tpr', 'thresholds']).T
        roc_df.to_csv(crossval_final_out / 'final_model_ROC.csv', index=False)
        plot_curve(roc_df.dropna().fpr, roc_df.dropna().tpr, out_path / '2-final_model_ROC.png', x_lab='fpr', y_lab='tpr', label=f'Area Under ROC curve (AUC) = {roc_final_eval_auc:.2f}', title=f'Receiver Operating Characteristic (ROC) curve for top {max_auc_genes} genes [{model_eval.name}/{model_eval.fselection.name}]')
        prc_df = pd.DataFrame(prc_final_eval_data, index=['precision', 'recall', 'thresholds']).T
        prc_df.to_csv(crossval_final_out / 'final_model_PRC.csv', index=False)
        plot_curve(prc_df.dropna().recall, prc_df.dropna().precision, out_path / '2-final_model_PRC.png', x_lab='recall', y_lab='precision', label=f'Average precision (AP) = {prc_final_eval_auc:.2f}', title=f'Precision-Recall curve (PRC) for top {max_auc_genes} genes [{model_eval.name}/{model_eval.fselection.name}]', diag_x=[0, 1], diag_y=[1, 0])

    # build and save final model
    if featurelist and featurelist.exists():
        fselection_final = fselection_final_eval
        fselection_final.final = True
        fselection_final.selected = {'all': fselection_df}
    else:
        fselection_final = fselection_final_eval
        fselection_final.final = True
        fselection_final.selected = {'all': rank1_df_final}
    model_final = get_model(model, kwargs_dict, fselection_final)
    model_final.final = True
    model_final_res = model_final.train(dataset, maxgenes_final, 'all')
    model_final.crosscal_pred = predictions_final_eval_samples
    model_final_predictions = pd.DataFrame.from_dict({k: v for k, v in model_final_res.items()}, orient='index', columns=['0', '1']).loc[dataset.index]
    pd.DataFrame(model_final_predictions).rename_axis('sampleid', axis='rows').to_csv(wd_path / 'complete_model_reprediction_predictions.csv')
    y_final_true, y_final_scores = dataset['class'], model_final_predictions['1']
    roc_final_data = metrics.roc_curve(y_final_true, y_final_scores)
    roc_final_auc = metrics.roc_auc_score(y_final_true, y_final_scores)
    prc_final_data = metrics.precision_recall_curve(y_final_true, y_final_scores)
    prc_final_auc = metrics.average_precision_score(y_final_true, y_final_scores)
    save_model(model_final, model_path)
    save_features(model_final.get_selected_genes(maxgenes_final, 'all'), features_path)
    roc_final_df = pd.DataFrame(roc_final_data, index=['fpr', 'tpr', 'thresholds']).T
    roc_final_df.to_csv(wd_path / 'complete_model_reprediction_ROC.csv', index=False)
    plot_curve(roc_final_df.dropna().fpr, roc_final_df.dropna().tpr, wd_path / 'complete_model_reprediction_ROC.png', x_lab='fpr', y_lab='tpr', label=f'Area Under ROC curve (AUC) = {roc_final_auc:.2f}', title=f'Receiver Operating Characteristic (ROC) curve for final model using {maxgenes_final} genes [{model_eval.name}/{model_eval.fselection.name}]')
    prc_final_df = pd.DataFrame(prc_final_data, index=['precision', 'recall', 'thresholds']).T
    prc_final_df.to_csv(wd_path / 'complete_model_reprediction_PRC.csv', index=False)
    plot_curve(prc_final_df.dropna().recall, prc_final_df.dropna().precision, wd_path / 'complete_model_reprediction_PRC.png', x_lab='recall', y_lab='precision', label=f'Average precision (AP) = {prc_final_auc:.2f}', title=f'Precision-Recall curve (PRC) for final model using {maxgenes_final} genes [{model_eval.name}/{model_eval.fselection.name}]', diag_x=[0, 1], diag_y=[1, 0])


def predict(pred_id, model_file, genescores, features, variantfn, genefn, outfolder):
    model_file = Path(model_file)
    model = load_model(model_file)
    features = Path(features)
    genescores = Path(genescores)
    genescores_df = pd.read_csv(genescores)
    missing_features_list = []
    if not features.exists():
        print(
            'No features file found, omitting features validation and using all features supplied in dataset. '
            'Make sure all features used in training are provided and sorted accordingly.'
        )
    else:
        features_s = pd.read_csv(features, header=None)[0]
        features_list = features_s.to_list()
        features_availability = features_s.isin(genescores_df.Gene)
        if sum(features_availability) != features_s.shape[0]:
            missing = features_s[~features_availability].values
            print(f'|8.00| {len(missing)} missing feature(s) compared with features file: {missing}. Using default NA score.')
            missing_features_list = missing
        genescores_df = genescores_df[genescores_df.Gene.isin(features_list)]

    features_s = pd.Series(model.features)
    features_availability = features_s.isin(genescores_df.Gene)
    if sum(features_availability) != features_s.shape[0]:
        missing = features_s[~features_availability].values
        missing_non_overlap = list(set(missing) - set(missing_features_list))
        if missing_non_overlap:
            missing_features_list = missing
            print(f'|8.00| {len(missing_non_overlap)} additional missing feature(s) compared with model features: {missing_non_overlap}. Using default NA score.')

    dataset = genescores_df.drop('Transcript', axis=1).set_index('Gene').T.rename_axis('sampleid', axis='rows')

    # add missing features (genes) using the default NA value from the gene_score class used to generate genescores
    nan = get_NA_score(variantfn, genefn)
    dataset = dataset.assign(**{col: nan for col in missing_features_list})

    # select features subset in relevant order
    dataset = dataset[features_s]

    # save genescores used for prediction
    dataset.to_csv(outfolder / f'{pred_id}_genescores_prediction.csv')

    classes = list(model.model.classes_)
    try:
        y_pred = model.predict(dataset)
    except ValueError as err:
        print(f'|8.00| ERROR: {err}')
        exit(1)
    y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
    y_pred_ordered_dict = dict(zip(dataset.index.tolist(), list(y_pred_ordered)))
    predictions_df = pd.DataFrame.from_dict(y_pred_ordered_dict, orient='index', columns=['0', '1']).rename_axis('sampleid', axis='rows').rename(columns={'0': 'score_0', '1': 'score_1'})
    predictions_df.to_csv(outfolder / f'{pred_id}_predictions.csv')

    # plot predictions
    plot_jitter(model.crosscal_pred, data_overlay=predictions_df, title='Predictions', x_lab='class', y_lab='prediction score', save_as=outfolder / f'{pred_id}_predictions.png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-M', '--makepredictions', action='store_true',
        help='predict samples using specified model'
    )
    parser.add_argument(
        '-g', '--genescores', type=Path,
        help='path to the genescores file'
    )
    parser.add_argument(
        '-f', '--featureselection', type=str, default='ks',
        help='feature selection method'
    )
    parser.add_argument(
        '-m', '--model', type=str, default='rf',
        help='machine learning model'
    )
    parser.add_argument(
        '-c', '--cvscheme', type=Path,
        help='path to cross-validation scheme file'
    )
    parser.add_argument(
        '-p', '--protlength', type=Path,
        help='path to protein length file'
    )
    parser.add_argument(
        '-K', '--kfold', type=int,
        help='use specified k-fold cross-validation overwriting any provided scheme'
    )
    parser.add_argument(
        '-v', '--variation', type=float, default=80,
        help='cutoff for variation pre-filter (%%)'
    )
    parser.add_argument(
        '-V', '--variance', type=float, default=0,
        help='cutoff for variance pre-filter'
    )
    parser.add_argument(
        '-P', '--pvalue', type=float, default=0,
        help='set p-value cutoff for methods like ks'
    )
    parser.add_argument(
        '-G', '--maxgenes', type=int, default=200,
        help='number of top-ranked genes to use for model evaluation'
    )
    parser.add_argument(
        '-T', '--topgenes', type=float, default=0.5,
        help='required ratio for top-ranked genes to include into final model training'
    )
    parser.add_argument(
        '-S', '--stepsize', type=int, default=5,
        help='stepsize for maxgenes'
    )
    parser.add_argument(
        '-F', '--featurelist', type=Path,
        help='list of features (genes) to use instead of feature selection'
    )
    parser.add_argument(
        '-o', '--out', type=Path,
        help='path to the output folder'
    )
    parser.add_argument(
        '-w', '--wd', type=Path,
        help='path to the wd folder'
    )
    parser.add_argument(
        '-C', '--cores', type=int, default=1,
        help='number of cores to use for computation'
    )
    parser.add_argument(
        '-q', '--quiet', action='store_true',
        help='do not show realtime progress'
    )
    main(*parser.parse_known_args())
