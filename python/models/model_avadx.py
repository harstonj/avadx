import numpy as np
import pandas as pd
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D


class Model:
    def __init__(self, name, kwargs, fselection, final=False, seed=42):
        self.name = name.upper()
        self.kwargs = kwargs
        self.fselection = fselection
        self.final = final
        self.seed = seed
        self.train = getattr(self, f'{self.name}_train', None)
        self.predict = getattr(self, f'{self.name}_predict', None)
        self.eval = getattr(self, f'{self.name}_eval', None)
        self.model = None
        self.features = None

    def run(self, kwargs):
        if self.train:
            self.train(**kwargs)
        else:
            print(f'ERROR: {self.name}() not available.')

    def update_progress(self, k):
        if self.fselection.progress:
            print(f'progress:update:{self.name}')

    def get_arg(self, arg):
        return self.kwargs.get(arg, None)

    def split_train_test(self, dataset, k):
        if self.final:
            return dataset.copy(), dataset.copy()
        else:
            train = dataset.drop(self.fselection.cvscheme[self.fselection.cvscheme.fold == k].index)
            test = dataset.loc[self.fselection.cvscheme[self.fselection.cvscheme.fold == k].index]
            return train, test

    def get_class_weights(self, train):
        train_class_weights = 1 / train.groupby('class')['class'].count()
        return train_class_weights.div(train_class_weights.sum())

    def get_selected_genes(self, max_genes, k):
        genes_selected = self.fselection.selected[k]
        genes_count = min(max_genes, genes_selected.shape[0])
        return genes_selected[0:genes_count]

    def resample(self, df):
        from sklearn.utils import resample
        class_0 = df[df['class'] == 0]
        class_1 = df[df['class'] == 1]
        df_minority, df_majority = (class_0, class_1) if class_0.shape[0] < class_1.shape[0] else (class_1, class_0)
        upsample_cnt = df_majority.shape[0] * 2
        df_minority_upsampled = resample(df_minority, replace=True, n_samples=upsample_cnt, random_state=self.seed)
        df_majority_upsampled = resample(df_majority, replace=True, n_samples=upsample_cnt, random_state=self.seed)
        df_upsampled = pd.concat([df_minority_upsampled, df_majority_upsampled])
        return df_upsampled

    def eval_rf_max_depth(self, rf_classifier, max_genes, k, y_train, y_test, X_train, X_test, depth_from=1, depth_to=None, min_depth=5, plot=False):
        eval_out = self.get_arg('out_path') / self.get_arg('eval')
        eval_out.mkdir(exist_ok=True)
        depth_to = max_genes if depth_to is None else depth_to
        max_depths = np.linspace(depth_from, depth_to, (depth_to - depth_from) + 1, endpoint=True)
        train_results, test_results, max_test_idx = [], [], 0
        max_depth_original = rf_classifier.get_params()['max_depth']
        for max_depth in max_depths:
            rf_classifier.set_params(**{'max_depth': max_depth})
            rf_classifier.fit(X_train, y_train)
            train_pred = rf_classifier.predict(X_train)
            false_positive_rate_train, true_positive_rate_train, thresholds_train = roc_curve(y_train, train_pred)
            roc_auc_train = auc(false_positive_rate_train, true_positive_rate_train)
            train_results.append(roc_auc_train)
            y_pred = rf_classifier.predict(X_test)
            false_positive_rate_test, true_positive_rate_test, thresholds_test = roc_curve(y_test, y_pred)
            roc_auc_test = auc(false_positive_rate_test, true_positive_rate_test)
            test_results.append(roc_auc_test)
            max_test_idx = int(max_depth) - 1 if roc_auc_test > test_results[max_test_idx] else max_test_idx
            if max_depth > (max_depth + int(np.round((max_test_idx + 1) / 2))):
                if max_depth >= min_depth and roc_auc_test <= test_results[max_test_idx]:
                    break
            max_test_auc = test_results[max_test_idx]
            while max_test_idx > 0 and max_test_auc == test_results[max_test_idx - 1]:
                max_test_idx -= 1
        rf_classifier.set_params(**{'max_depth': max_depth_original})
        if plot:
            line1, = plt.plot(max_depths, train_results, 'b', label='Train AUC')
            line2, = plt.plot(max_depths, test_results, 'r', label='Test AUC')
            plt.legend(handler_map={line1: HandlerLine2D(numpoints=2)})
            plt.ylabel('AUC score')
            plt.xlabel('Tree depth')
            plt.ylim((0, 1))
            plt.xlim((1, len(train_results)))
            plt.xticks([_ for _ in range(1, len(train_results) + 1)])
            plt.savefig(eval_out / f'eval_max_depth_{max_genes}_{k}.png')
            plt.close()
        else:
            with (eval_out / f'eval_max_depth_{max_genes}_{k}.tmp').open('w') as fout:
                fout.write(f'{max_test_idx + 1},{test_results[max_test_idx]}\n')
        return max_test_idx + 1

    def RF_train(self, dataset, max_genes, k, predict=False):
        train, test = self.split_train_test(dataset, k)
        train_class_weights = self.get_class_weights(train)
        genes_selected = self.get_selected_genes(max_genes, k)
        y_train, y_test = train.pop('class'), test.pop('class')
        X_train, X_test = train[genes_selected.index], test[genes_selected.index]
        rf_classifier = RandomForestClassifier(
            random_state=self.seed,
            class_weight=train_class_weights.to_dict(),
            n_estimators=min(round(X_train.shape[0] * 1.5), 500),
            bootstrap=False
        )
        rf_classifier.fit(X_train, y_train)
        if self.final:
            self.model = rf_classifier
            self.features = list(X_train.columns)
        classes = list(rf_classifier.classes_)
        y_pred = rf_classifier.predict_proba(X_test)
        y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
        self.update_progress(k)
        return dict(zip(y_test.index.tolist(), list(y_pred_ordered)))

    def RF_predict(self, dataset):
        y_pred = self.model.predict_proba(dataset)
        return y_pred

    def RF_RESAMPLE_train(self, dataset, max_genes, k, predict=False):
        train, test = self.split_train_test(dataset, k)
        train_class_weights = self.get_class_weights(train)
        genes_selected = self.get_selected_genes(max_genes, k)
        train = self.resample(train)
        y_train, y_test = train.pop('class'), test.pop('class')
        X_train, X_test = train[genes_selected.index], test[genes_selected.index]
        rf_classifier = RandomForestClassifier(
            random_state=self.seed,
            class_weight=train_class_weights.to_dict(),
            n_estimators=min(round(X_train.shape[0] * 1.5), 500),
            bootstrap=False
        )
        rf_classifier.fit(X_train, y_train)
        if self.final:
            self.model = rf_classifier
            self.features = list(X_train.columns)
        classes = list(rf_classifier.classes_)
        y_pred = rf_classifier.predict_proba(X_test)
        y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
        self.update_progress(k)
        return dict(zip(y_test.index.tolist(), list(y_pred_ordered)))

    def RF_RESAMPLE_predict(self, dataset):
        y_pred = self.model.predict_proba(dataset)
        return y_pred

    def RFE_train(self, dataset, max_genes, k, predict=False):
        train, test = self.split_train_test(dataset, k)
        train_class_weights = self.get_class_weights(train)
        genes_selected = self.get_selected_genes(max_genes, k)
        y_train, y_test = train.pop('class'), test.pop('class')
        X_train, X_test = train[genes_selected.index], test[genes_selected.index]
        rf_classifier = RandomForestClassifier(
            random_state=self.seed,
            class_weight=train_class_weights.to_dict(),
            n_estimators=min(round(X_train.shape[0] * 1.5), 500),
            bootstrap=False
        )
        val_fraction = 0.10
        y_val = y_train.sample(frac=val_fraction, weights=y_train.groupby(y_train).transform('count'))
        X_val = X_train.loc[y_val.index]
        y_trainsub = y_train.drop(y_val.index)
        X_trainsub = X_train.drop(y_val.index)
        max_depth = self.eval_rf_max_depth(rf_classifier, max_genes, k, y_trainsub, y_val, X_trainsub, X_val, depth_from=1)
        rf_classifier.set_params(**{'max_depth': max_depth})
        rf_classifier.fit(X_train, y_train)
        if self.final:
            self.model = rf_classifier
            self.features = list(X_train.columns)
        classes = list(rf_classifier.classes_)
        y_pred = rf_classifier.predict_proba(X_test)
        y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
        self.update_progress(k)
        return dict(zip(y_test.index.tolist(), list(y_pred_ordered)))

    def RFE_eval(self, dataset, max_genes):
        eval_out = self.get_arg('out_path') / self.get_arg('eval')
        data = list(eval_out.glob(f'eval_max_depth_{max_genes}_*.tmp'))
        data_depth, data_auc = [], []
        for d in data:
            with d.open() as fin:
                depth, auc = fin.readline().strip().split(',')
                data_depth += [int(depth)]
                data_auc += [float(auc)]
            d.unlink()
        # selected_depth = max(data_depth, key=data_depth.count)
        selected_depth = int(round(np.median(data_depth)))
        y, x, _ = plt.hist(data_depth, density=False, bins=len(set(data_depth)))
        plt.ylabel('counts')
        plt.xlabel('RF max depth')
        plt.savefig(eval_out / f'eval_max_depth_{max_genes}_hist_{selected_depth}.png')
        plt.close()

    def SVM_train(self, dataset, max_genes, k):
        train, test = self.split_train_test(dataset, k)
        train_class_weights = self.get_class_weights(train)
        genes_selected = self.get_selected_genes(max_genes, k)
        y_train, y_test = train.pop('class'), test.pop('class')
        X_train, X_test = train[genes_selected.index], test[genes_selected.index]
        svm_classifier = SVC(
            random_state=self.seed,
            probability=True,
            class_weight=train_class_weights.to_dict()
        )
        svm_classifier.fit(X_train, y_train)
        if self.final:
            self.model = svm_classifier
            self.features = list(X_train.columns)
        classes = list(svm_classifier.classes_)
        y_pred = svm_classifier.predict_proba(X_test)
        y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
        self.update_progress(k)
        return dict(zip(y_test.index.tolist(), list(y_pred_ordered)))

    def SVM_predict(self, dataset):
        y_pred = self.model.predict_proba(dataset)
        return y_pred
