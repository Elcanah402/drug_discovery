import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data
file_path = "swissadme.csv"
data = pd.read_csv(file_path)

# Set the 'Ligands' column as index for a proper heatmap layout
data.set_index('Ligands', inplace=True)

# Transpose the DataFrame so properties are on the y-axis and ligands on the x-axis
data_transposed = data.T

# Create the heatmap
plt.figure(figsize=(10, 10))
sns.heatmap(data_transposed, annot=True, cmap='Reds', cbar=True)

# Customize the heatmap
plt.title("Heatmap of ADME Properties", fontsize=16)
plt.xlabel("Ligands")
plt.ylabel("ADME Properties")

# Display the heatmap
plt.tight_layout()
plt.show()
