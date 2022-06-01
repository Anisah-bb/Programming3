import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('output/timings.txt', header = None)
time = df[0]
x = np.linspace(1, len(time), len(time))


plt.plot(x, time)
plt.xlabel('num_of_threads')
plt.ylabel('time_taken')
plt.title('Time Graph')
plt.show()
plt.savefig('output/timings.png')