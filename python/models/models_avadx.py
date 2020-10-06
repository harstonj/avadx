from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier


class Model:
    def __init__(self, name, kwargs, fselection, seed=42):
        self.name = name.upper()
        self.kwargs = kwargs
        self.fselection = fselection
        self.seed = seed
        self.fn = getattr(self, self.name, None)

    def run(self, kwargs):
        if self.fn:
            self.fn(**kwargs)
        else:
            print(f'ERROR: {self.name}() not available.')

    def get_arg(self, arg):
        return self.kwargs.get(arg, None)

    def split_train_test(self, dataset, k):
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

    def RF(self, dataset, max_genes, k):
        train, test = self.split_train_test(dataset, k)
        train_class_weights = self.get_class_weights(train)
        genes_selected = self.get_selected_genes(max_genes, k)
        y_train, y_test = train.pop('class'), test.pop('class')
        X_train, X_test = train[genes_selected.index], test[genes_selected.index]
        rf_classifier = RandomForestClassifier(
            random_state=self.seed,
            class_weight=train_class_weights.to_dict()
        )
        rf_classifier.fit(X_train, y_train)
        classes = list(rf_classifier.classes_)
        y_pred = rf_classifier.predict_proba(X_test)
        y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
        return dict(zip(y_test.index.tolist(), list(y_pred_ordered)))

    def SVM(self, dataset, max_genes, k):
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
        classes = list(svm_classifier.classes_)
        y_pred = svm_classifier.predict_proba(X_test)
        y_pred_ordered = [[y_pred_instance[classes.index(0)], y_pred_instance[classes.index(1)]] for y_pred_instance in y_pred]
        return dict(zip(y_test.index.tolist(), list(y_pred_ordered)))
