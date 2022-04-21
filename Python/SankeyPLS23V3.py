import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import plotly
import plotly.express as px
import plotly.graph_objects as go
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import LeaveOneOut
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import MaxAbsScaler

#inputs=['HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11','SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7']
inputs=['HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HDAC8','HDAC9','HDAC10','HDAC11','SIRT1','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7', 'KAT2A', 'KAT2B', 'HAT1', 'ATF2', 'KAT5', 'KAT6A', 'KAT6B', 'KAT7', 'EP300', 'CREBBP', 'NCOA1', 'NCOA3', 'TAF1','GTF3C1', 'CLOCK']
midputs=['FOXO1', 'FOXO3', 'GATA4', 'GATA6', 'HIF1A', 'KLF4', 'KLF5', 'MEF2A', 'NFAT5', 'NFKB1', 'NKX25', 'NOTCH1', 'RUNX1', 'SHMT2', 'SOD1', 'TBX5', 'TGFB1','TRIM28', 'YY1']
outputs=['SCN5A', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'ATP1A1', 'SLC8A1', 'ATP2A2', 'RYR2', 'GJA1']

MaleAndFemale=pd.read_csv('Male+Female.csv')
MaleAndFemale.dropna()

ID = MaleAndFemale['ID']
Sex = MaleAndFemale['Sex']
MaleAndFemale = MaleAndFemale.drop('ID', axis = 1)
MaleAndFemale = MaleAndFemale.drop('Sex', axis = 1)

Male_Data = MaleAndFemale[Sex==0]
Female_Data = MaleAndFemale[Sex==1]

s = 'female'
scaled = True
threshold = 0.00
cutoff = 15 #If you have to many input varaibles, the amount of variables to not display so that it fits in the figure

X_data_M = Male_Data[inputs]
X_data_F = Female_Data[inputs]

y_data_M = Male_Data[midputs]
y_data_F = Female_Data[midputs]

z_data_M = Male_Data[outputs]
z_data_F = Female_Data[outputs]

if s =='female':
    print("Female Selected")
    X_data = X_data_F
    y_data = y_data_F
    z_data = z_data_F
elif s == 'male':
    print("Male Selected")
    X_data = X_data_M
    y_data = y_data_M
    z_data = z_data_M

    
else:
    print("Both Selected")
    X_data = MaleAndFemale[inputs]
    y_data = MaleAndFemale[midputs]
    z_data = MaleAndFemale[outputs]

pls_reg1 = PLSRegression(n_components=4, scale=scaled)
pls_reg1.fit(X_data,z_data)
pls_reg2 = PLSRegression(n_components=4, scale=scaled)
pls_reg2.fit(X_data,y_data)
pls_reg3 = PLSRegression(n_components=4, scale=scaled)
pls_reg3.fit(y_data,z_data)
imp_1 = pls_reg1.coef_
imp_2 = pls_reg2.coef_
imp_3 = pls_reg3.coef_

imp_1s = MaxAbsScaler().fit_transform(imp_1)
imp_2s = MaxAbsScaler().fit_transform(imp_2)
imp_3s = MaxAbsScaler().fit_transform(imp_3)
pd.DataFrame(imp_1s, columns = outputs, index = inputs).to_csv('imp_1sf.csv')
pd.DataFrame(imp_2s, columns = midputs, index = inputs).to_csv('imp_2sf.csv')
pd.DataFrame(imp_3s, columns = outputs, index = midputs).to_csv('imp_3sf.csv')

pls1 = pd.DataFrame(imp_1, columns = outputs, index = inputs).to_csv('Beta_1f.csv')
pls2 = pd.DataFrame(imp_2, columns = midputs, index = inputs).to_csv('Beta_2f.csv')
pls3 = pd.DataFrame(imp_3, columns = outputs, index = midputs).to_csv('Beta_3f.csv')



labels = []
for i in range(len(inputs)-cutoff):
    labels.append(inputs[i])
    
for i in range(len(midputs)):
    labels.append(midputs[i])

for i in range(len(outputs)):
    labels.append(outputs[i])
    
source = []
target = []
values = []
for i in range(len(midputs)):
    for j in range(len(inputs)-cutoff): # We do not want to view HATs so we subtract total by 15
        if abs(imp_2[j][i]) > threshold:
            source.append(j)
            target.append(i+len(inputs)-cutoff)
            values.append(imp_2s[j][i])

for i in range(len(outputs)):
    for j in range(len(midputs)):
        if abs(imp_3[j][i]) > threshold:
            source.append(j+len(midputs)-1)
            target.append(i+len(inputs)-cutoff + len(midputs))
            values.append(imp_3s[j][i])
        
v = np.array(values)
values_n = np.interp(v, (-1, 0, 1), (0,+.5, +1))
c_colorscale = [[0.0,'rgb(248, 105, 107)'],
               [0.5, 'rgb(255, 255, 255)'],
               [1.0, 'rgb(99, 190, 123)']]

colors =  px.colors.sample_colorscale(c_colorscale,values_n, low = 0, high = 1)

y1 = np.linspace(0.1, .9, len(inputs)-cutoff)
y2 = np.linspace(0.1, .9, len(midputs))
y3 = np.linspace(0.3, .7, len(outputs))
y = np.append(y1,y2)
y = np.append(y,y3)

x1 = np.ones((1, len(inputs)-cutoff))*.2
x2 = np.ones((1, len(midputs)))*.50
x3 = np.ones((1, len(outputs)))*.8
x = np.append(x1,x2)
x = np.append(x,x3)

fig = go.Figure(
    data=[
        go.Sankey(
            node=dict(
                thickness=10,
                line=dict(color="black", width=0.5),
                label=labels,
                x = x,
                y = y,
                
            ),
            link=dict(
                # flow source node index
                source=source,
                # flow target node index
                target=target,
                # flow quantity for each source/target pair
                value=np.abs(v),
                color=colors,
            ),
            arrangement = 'snap',
        )
    ] 
)

fig.update_layout(width=2400, height=900, font_size = 18)
fig.show()