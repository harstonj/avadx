import pandas as pd
from scipy.stats import ks_2samp
from sklearn.feature_selection import RFE, SelectKBest, f_classif
from skrebate import ReliefF
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import VarianceThreshold


class FSelection:
    def __init__(self, name, kwargs, cvscheme, max_genes=200, fold_variance_filter=0, final=False, seed=42):
        self.name = name.upper()
        self.kwargs = kwargs
        self.cvscheme = cvscheme
        self.max_genes = max_genes
        self.fold_variance_filter = fold_variance_filter
        self.final = final
        self.seed = seed
        self.fn = getattr(self, self.name, None)
        self.selected = None

    def run(self, kwargs):
        if self.fn:
            self.fn(**kwargs)
        else:
            print(f'ERROR: {self.name}() not available.')

    def get_arg(self, arg):
        return self.kwargs.get(arg, None)

    def split_train_test(self, dataset, k):
        if self.final:
            dataset_ = dataset.copy()
            if self.fold_variance_filter is not None:
                var_selector = VarianceThreshold(self.fold_variance_filter)
                try:
                    var_selector.fit(dataset_)
                    dataset_ = dataset_[dataset_.columns[var_selector.get_support(indices=True)]]
                except ValueError as err:
                    print(f'Error at fold variance filter step: {err}')
            return dataset_, None
        else:
            train = dataset.drop(self.cvscheme[self.cvscheme.fold == k].index)
            test = dataset.loc[self.cvscheme[self.cvscheme.fold == k].index]
            if self.fold_variance_filter is not None:
                var_selector = VarianceThreshold(self.fold_variance_filter)
                try:
                    var_selector.fit(train)
                    train = train[train.columns[var_selector.get_support(indices=True)]]
                    test = test[test.columns[var_selector.get_support(indices=True)]]
                except ValueError as err:
                    print(f'Error at fold variance filter step: {err}')
            return train, test

    def KS(self, dataset, k):
        class_0_idx = dataset[dataset['class'] == 0].index
        class_1_idx = dataset[dataset['class'] == 1].index
        class_0_idx_k = class_0_idx.drop(self.cvscheme[self.cvscheme.fold == k].index, errors='ignore')
        class_1_idx_k = class_1_idx.drop(self.cvscheme[self.cvscheme.fold == k].index, errors='ignore')
        res = dataset.drop(self.cvscheme[self.cvscheme.fold == k].index).drop('class', axis=1)\
            .apply(lambda x: ks_2samp(x[class_1_idx_k], x[class_0_idx_k]), axis=0)\
            .rename(index={0: 'statistic', 1: 'pvalue'})\
            .rename_axis(k, axis='rows')\
            .loc['pvalue'].rename(k)
        selected = res[res < self.get_arg('pval_cutoff')].sort_values() if self.get_arg('pval_cutoff') is not None else res.sort_values()
        print(f'progress:update:{self.name} - fold {k}...')
        return (res, selected)

    def KBEST(self, dataset, k):
        train, _ = self.split_train_test(dataset, k)
        y_train = train.pop('class')
        X_train = train
        kbest = SelectKBest(score_func=f_classif, k=self.max_genes)
        kbest.fit(X_train, y_train)
        res = pd.DataFrame(kbest.scores_, X_train.columns.tolist(), columns=[k])[k]
        selected = res[kbest.get_support(indices=True)].sort_values(ascending=False)
        print(f'progress:update:{self.name} - fold {k}...')
        return (res, selected)

    def RFE(self, dataset, k):
        train, _ = self.split_train_test(dataset, k)
        y_train = train.pop('class')
        X_train = train
        model = LogisticRegression(solver='lbfgs')
        rfe = RFE(model, self.max_genes)
        rfe.fit(X_train.to_numpy(), y_train.to_numpy())
        res = pd.DataFrame(rfe.ranking_, X_train.columns, columns=[k])[k]
        selected = res[rfe.get_support(indices=True)].sort_values()
        print(f'progress:update:{self.name} - fold {k}...')
        return (res, selected)

    def RELIEFF(self, dataset, k):
        train, _ = self.split_train_test(dataset, k)
        y_train = train.pop('class')
        X_train = train
        fs_relieff = ReliefF(n_features_to_select=2, n_neighbors=100, verbose=True)
        fs_relieff.fit(X_train.to_numpy(), y_train.to_numpy())
        res = pd.DataFrame(fs_relieff.feature_importances_, X_train.columns, columns=[k])[k]
        selected = res.sort_values()[0:self.max_genes]
        print(f'progress:update:{self.name} - fold {k}...')
        return (res, selected)
