import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

File = sys.argv[1]
Hour = sys.argv[2]
Name = sys.argv[3]
Expected = sys.argv[4]

df = pd.read_csv(File,sep="\t",header=None,index_col=None)
df.rename(columns={0: 'Base Pairs'}, inplace=True)

avg = df["Base Pairs"].mean()
stdev = df["Base Pairs"].std()
score = (float(Expected) - avg)/stdev
#print(File)
#print("STDEV: " + str(stdev))
sns.set(style="darkgrid")
plt.title("Probability Distribution for " + str(Hour) + " Intersected with " + str(Name),size=10)
sns.histplot(data=df, x="Base Pairs", kde=True)

plt.axvline(avg,linestyle='--',color="black")
if(score<=4 and score>=-4):
    plt.axvline(int(Expected),linestyle='--',color="black")
    plt.text(int(Expected)+2,plt.ylim()[1]-1,'Observed BP',size=10)
plt.text(int(avg)+2,plt.ylim()[1]-1,'Expected BP',size=10)
plt.savefig(Hour+"."+Name+".png")

print(str(stdev))
print(str(avg) +"\t"+ str(score))
