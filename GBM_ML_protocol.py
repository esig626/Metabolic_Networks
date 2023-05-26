import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# Load the microarray data
data = np.loadtxt("microarray_data.txt", delimiter="\t")

# Separate the features (gene expression values) and labels (hypoxic phenotypes)
X = data[:, :-1]
y = data[:, -1]

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize and fit the gradient boosting classifier
gbm = GradientBoostingClassifier()
gbm.fit(X_train, y_train)

# Predict the hypoxic phenotypes on the test set
y_pred = gbm.predict(X_test)

# Calculate the accuracy of the predictions
accuracy = accuracy_score(y_test, y_pred)

# Get the feature importances from the GBM model
feature_importances = gbm.feature_importances_

# Sort the feature importances in descending order
sorted_indices = np.argsort(feature_importances)[::-1]

# Get the top N genes based on feature importances
top_genes = sorted_indices[:N]

# Print the top N genes
print("Top {} genes:".format(N))
for gene_index in top_genes:
    print("Gene {}: Importance = {}".format(gene_index, feature_importances[gene_index]))
