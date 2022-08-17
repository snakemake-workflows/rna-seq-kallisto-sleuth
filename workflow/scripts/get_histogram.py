from turtle import title
import altair as alt
from altair_saver import save
import pandas as pd
import numpy as np
import pysam
from altair_transform import transform_chart
from scipy.stats import gaussian_kde
from scipy import stats
import sys
import json
sys.stderr = open(snakemake.log[0], "w")


samples = snakemake.params["samples"]
f = open(snakemake.params["read_length"])
read_length = json.load(f)
f.close()


allsamp_disdata = pd.DataFrame([])
for each_sample in samples:
    sample_name = each_sample.split("/")[2]
    print(each_sample)
    print(sample_name)
    #Bam sorting and Indexing
    pysam.sort(each_sample + "/pseudoalignments.bam", "-o", each_sample + "/pseudoalignments_sort.bam")
    pysam.index("%s/pseudoalignments_sort.bam" % each_sample)
    #Bam file reading
    bam_file = pysam.AlignmentFile(each_sample + '/pseudoalignments_sort.bam',"rb")
    bam_header = bam_file.header.to_dict()
    trans_length_data = pd.DataFrame(bam_header.get('SQ'))
    trans_length_data.rename(columns={'SN': 'Transcript_ID'}, inplace=True)
    #Aligned text file reading
    align_bam_txt = pd.read_csv('results/QC/' + sample_name + '.aligned.txt', 
        sep="\t",names=["read_Name", "Transcript_ID", "Start"])
    #Both transcript len and start postion are merged based on same transcript ID
    merge_data = align_bam_txt.merge(trans_length_data, on='Transcript_ID')
    #Each read postion is calcuated
    merge_data[sample_name + '_read_dist'] = merge_data['LN'] - merge_data['Start']
    #Gaussian calcuation
    kde_dist = gaussian_kde(merge_data[sample_name + '_read_dist'])
    data_space = np.linspace(merge_data[sample_name + '_read_dist'].min(), 10000, num=100)
    evaluated_dist = kde_dist.evaluate(data_space)
    #import pdb; pdb.set_trace()
    #Python will evaluate each step

    eachsample_dist = pd.DataFrame({'sample_name': sample_name, 'data_points': data_space,
        'kde_dist': evaluated_dist })
    allsamp_disdata = pd.concat([allsamp_disdata, eachsample_dist])

hist = alt.Chart().mark_line().encode(x = alt.X('data_points', 
    title="difference between transcript length and read start"), 
        y =alt.Y('kde_dist', title = "density"), color = 'sample_name')

allsamp_disdata['read_length'] = read_length
chart_read_length = alt.Chart().mark_rule().encode(
    x='read_length')

final_chart = (hist+chart_read_length).facet(row='site', data=allsamp_disdata)

final_chart.save(snakemake.output["histogram"])

