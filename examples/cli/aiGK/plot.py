import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('velocities.csv', index_col=0)

df.plot()
plt.savefig('plot.pdf')

