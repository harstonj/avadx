import argparse
import pandas as pd
# import matplotlib.pyplot as plt
# from timeit import default_timer as timer
from pathlib import Path
from multiprocessing import Pool
from joblib import dump, load
from sklearn import metrics
from sklearn.feature_selection import VarianceThreshold


def flatten(list):
    return [item for sublist in list for item in sublist]


def get_fselection(featureselection, kwargs_dict, cvscheme, maxgenes):
    from feature_selections.fselection_avadx import FSelection
    return FSelection(featureselection, kwargs_dict, cvscheme, maxgenes)


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


def main(args, kwargs):
    run(args.genescores, args.featureselection, args.featurelist, args.model, args.cvscheme, args.protlength, args.kfold, args.variance, args.variation, args.pvalue, args.maxgenes, args.stepsize, args.out, args.wd, args.cores, kwargs)


def run(genescores_path, featureselection, featurelist, model, cvscheme_path, protlength_path, kfold, variance_cutoff, variation_cutoff, pval_cutoff, maxgenes, stepsize, out_path, wd_path, cores, kwargs):
    kwargs_dict = {_.split('=')[0]: _.split('=')[1] for _ in kwargs}
    kwargs_dict['pval_cutoff'] = pval_cutoff

    # create output directories
    out_path.mkdir(exist_ok=True)
    wd_path.mkdir(exist_ok=True)
    model_path = out_path / 'model.joblib'
    features_path = out_path / 'model_features.txt'

    # load datasets
    cvscheme = pd.read_csv(cvscheme_path, header=None, names=['sampleid', 'fold', 'class'], dtype={'sampleid': str, 'fold': int, 'class': int})
    cvscheme = cvscheme.set_index('sampleid')
    kfold_is_set = False if kfold is None else True
    kfold = kfold if kfold_is_set else len(cvscheme.fold.unique())
    kfold_steps = range(1, (cvscheme.shape[0] // kfold)) if kfold_is_set else sorted(cvscheme.fold.unique())
    if kfold_is_set:
        folds_new = flatten([[_] * (cvscheme.shape[0] // kfold) for _ in range(1, (cvscheme.shape[0] // kfold))])
        remaining = cvscheme.shape[0] % kfold
        if remaining:
            # TODO improve distribution of remaining among all folds instead of adding them to the last one
            folds_new += [folds_new[-1]] * remaining
        cvscheme.fold = folds_new

    genescores = pd.read_csv(genescores_path)
    protlength = pd.read_csv(protlength_path)

    # in case of multiple transcripts for the same gene only keep the longest one
    if sum(genescores.Gene.duplicated()):
        genescores = genescores.merge(protlength, how='left', on=['Gene', 'Transcript'])
        genescores = genescores.loc[genescores.Prot_length.eq(genescores.groupby('Gene').Prot_length.transform(max).fillna(0), fill_value=0)]
        genescores.drop('Prot_length', axis=1, inplace=True)
        genescores.reset_index(drop=True, inplace=True)

    # set all NaN to 0 and transpose to create the final dataset
    genescores = genescores[['Gene', 'Transcript'] + cvscheme.index.tolist()].fillna(0)
    dataset = genescores.drop('Transcript', axis=1).set_index('Gene').T.rename_axis('sampleid', axis='rows')

    # variance filtering
    var_selector = VarianceThreshold(variance_cutoff)
    try:
        var_selector.fit(dataset)
        dataset = dataset[dataset.columns[var_selector.get_support(indices=True)]]
    except ValueError as err:
        print(f'Error at variance filter step: {err}')

    # variation filtering
    gene_variation_across_samples = (dataset.nunique() / dataset.shape[0]) * 100
    gene_variation_subset = gene_variation_across_samples[gene_variation_across_samples <= variation_cutoff].index
    if len(gene_variation_subset):
        dataset = dataset[gene_variation_subset]
    else:
        print(f'Error at variation filter step: no feature met the required threshold of <= {variation_cutoff}%')

    # save filtered dataset
    dataset.to_csv(genescores_path.parent / f'{genescores_path.stem }_variation_filtered.csv')

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
        fselection = get_fselection('manual', kwargs_dict, cvscheme, genes_manual)
        fselection.selected = {k: fselection_df for k in kfold_steps}
        genes_considered = [genes_manual]
    else:
        fselection = get_fselection(featureselection, kwargs_dict, cvscheme, maxgenes)
        with Pool(processes=cores) as fselection_pool:
            fselection_pooled = [fselection_pool.apply_async(fselection.fn, args=(dataset, k)) for k in kfold_steps]
            fselection_res = [p.get() for p in fselection_pooled]
            fselection.selected = {genes.name: genes for res, genes in fselection_res}
            fselection_df = pd.DataFrame([res for res, genes in fselection_res]).rename_axis('folds', axis='rows')
            fselection_df.to_csv(wd_path / f'{kfold}F-CV-{fselection.name}-selectedGenes.csv')
            fselection_df.to_csv(out_path / 'performance_folds_genes.csv')
        genes_considered = list(range(stepsize, maxgenes + stepsize, stepsize))
        remaining = maxgenes - genes_considered[-1]
        if remaining:
            genes_considered += [maxgenes + remaining]

    # save order list (descending) of best scoring features (genes) as determimed by Feature Selection step
    pd.DataFrame(dict([(k, v.index.to_series(index=range(1, v.shape[0] + 1))) for k, v in fselection.selected.items()])).rename_axis('pos', axis='rows').to_csv(out_path / 'featureselection_genes_per_fold_ordered.csv')

    # run model training and performance evaluation
    model_eval = get_model(model, kwargs_dict, fselection)
    with Pool(processes=cores) as model_pool:
        performances_roc = {}
        performances_auc = {}
        for max_genes in genes_considered:
            print(f'{model_eval.name}/{model_eval.fselection.name}: {max_genes} genes...')
            model_pooled = [model_pool.apply_async(model_eval.fn, args=(dataset, max_genes, k)) for k in kfold_steps]
            model_res = [p.get() for p in model_pooled]
            model_predictions = pd.DataFrame.from_dict({k: v for d in model_res for k, v in d.items()}, orient='index', columns=['0', '1']).sort_index()
            roc_all = metrics.roc_curve(dataset['class'], model_predictions['1'])
            roc_auc = metrics.roc_auc_score(dataset['class'], model_predictions['1'])
            performances_roc[max_genes] = roc_all
            performances_auc[max_genes] = roc_auc
            print(f'{model_eval.name}/{model_eval.fselection.name}: {max_genes} genes AUC: {roc_auc}')
        max_auc_genes = max(performances_auc, key=lambda key: performances_auc[key])
        print(f'{model_eval.name}/{model_eval.fselection.name}: max AUC for {max_auc_genes} genes = {performances_auc[max_auc_genes]}')
        model_df = pd.DataFrame.from_dict(performances_auc, orient='index', columns=['AUC']).rename_axis('selected_genes', axis='rows').sort_values(by='AUC', ascending=False)
        model_df.to_csv(wd_path / f'{kfold}F-CV-{model_eval.name}-performance.csv')
        model_df.to_csv(out_path / 'performance.csv')

    # get list of selected genes for best AUC over all folds (merge)
    genes_best_merged = {}
    rank1_df = None
    for rank, max_genes in enumerate(model_df.index, 1):
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
            rank1_df.to_csv(wd_path / 'AUC_rank.1-genes-list.csv', columns=[], header=False)
            rank1_df.to_csv(out_path / 'AUC_rank.1-genes.csv', header=False)

    # build and save final model
    maxgenes_final = rank1_df.shape[0]
    fselection_final = get_fselection(featureselection, kwargs_dict, cvscheme, maxgenes_final)
    fselection_final.final = True
    fselection_final_res = fselection_final.fn(dataset, 'all')
    fselection_final.selected = {fselection_final_res[1].name: fselection_final_res[1]}
    model_final = get_model(model, kwargs_dict, fselection_final)
    model_final.final = True
    model_final_res = model_final.fn(dataset, maxgenes_final, 'all')
    model_final_predictions = pd.DataFrame.from_dict({k: v for k, v in model_final_res.items()}, orient='index', columns=['0', '1']).sort_index()
    # roc_final_all = metrics.roc_curve(dataset['class'], model_final_predictions['1'])
    roc_final_auc = metrics.roc_auc_score(dataset['class'], model_final_predictions['1'])
    save_model(model_final.model, model_path)
    save_features(model_final.get_selected_genes(maxgenes_final, 'all'), features_path)
    print(f'FINAL - {model_final.name}/{model_eval.fselection.name}: {maxgenes_final} genes AUC: {roc_final_auc}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        '-v', '--variation', type=float, default=99,
        help='cutoff for variation pre-filter (%)'
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
        help='number of top-ranked genes to use for model building'
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
        help='path to the output folder'
    )
    parser.add_argument(
        '-C', '--cores', type=int, default=1,
        help='number of cores to use for computation'
    )
    main(*parser.parse_known_args())
