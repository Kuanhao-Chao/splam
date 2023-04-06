import numpy as np
import matplotlib.pyplot as plt

# Generate sample data
data = np.random.randn(800)

# Reshape the data into a 2D array with 20 rows and 40 columns
data = data.reshape((1, 800))

# Set the figure size
fig = plt.figure(figsize=(12, 8))

# Plot the heatmap with linear axis and increased bar height
plt.imshow(data, cmap='coolwarm', aspect='auto', interpolation='nearest', vmin=-3, vmax=3)
plt.colorbar()

# Show the plot
plt.show()
