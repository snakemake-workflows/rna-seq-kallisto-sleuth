import altair as alt
from altair_saver import save
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('aligned_POS.txt', sep="\t")
#TransID_sorted = data.groupby(['Transcript_ID'])['Start'].transform(max)
data['max_rep_pos']=data.groupby('Transcript_ID')['Start'].transform(max)
#data['read_fract_num'] = data['Start']/data['max_rep_pos'].replace({ 0 : np.inf })
data['read_fract_num'] = data['Start']/data['max_rep_pos']
#data.to_csv('output_File2.csv')
#alt.data_transformers.disable_max_rows()
#alt.data_transformers.enable('data_server')
#alt.renderers.enable('altair_saver', fmts=['vega-lite', 'png'])   
#chart = alt.Chart(data).mark_bar().encode(x = 'read_fract_num', y = 'count()')
plt.hist(data['read_fract_num'], bins=150)
plt.xlabel('reads fractional number')
plt.ylabel("Frequency")
plt.tight_layout()
plt.savefig('hist_full.png')
"""
chart = alt.Chart(data).mark_bar().encode(x = alt.X('read_fract_num',
                                                       bin = alt.BinParams(maxbins = 60)),
                                             y = 'count()')
"""
#chart.save('chart4.html')