import pandas as pd
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier


class Model:
    def __init__(self, name, kwargs, fselection, final=False, seed=42):
        self.name = name.upper()
        self.kwargs = kwargs
        self.fselection = fselection
        self.final = final
        self.seed = seed
        self.train = getattr(self, f'{self.name}_train', None)
        self.predict = getattr(self, f'{self.name}_predict', None)
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
