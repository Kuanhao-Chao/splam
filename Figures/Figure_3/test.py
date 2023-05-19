import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

# Generate random data
np.random.seed(1)
x = np.random.randn(1000)
y = np.random.randn(1000)

# Calculate density using Gaussian kernel density estimation
density = gaussian_kde(np.vstack([x, y]))

# Calculate density values for each data point
density_values = density(np.vstack([x, y]))

# Define colormap and normalize density values
cmap = 'Purples'
norm = plt.Normalize(vmin=density_values.min(), vmax=density_values.max())

# Create scatter plot with varying dot darkness based on density
plt.scatter(x, y, c=density_values, cmap=cmap, alpha=0.5, norm=norm)
plt.colorbar()

# Set plot title and labels
plt.title('Scatter Plot with Varying Dot Darkness')
plt.xlabel('X')
plt.ylabel('Y')

# Show the plot
plt.show()
