import matplotlib.pyplot as plt
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import plotly.express as px
import csv
import sys
import pandas as pd
import numpy as np
import os
import subprocess as sp

BED = sys.argv[1]
META = sys.argv[2]
SIZE = sys.argv[3]
GENE = sys.argv[4]

meta_df = pd.read_csv(META,header=None, sep ="\t", index_col=3) #Read meta features file with 4 column as index

#######################################################################################################
#Determine if working with genes or other type of feature
if(GENE=='True'):
    if(len(meta_df.axes[1]) == 4): #Determine if strand column is present
        strand=True
        strand_list = meta_df[4].tolist()
    else:
        strand=False
    
    print(strand)
    if(strand==True): #Flank TSS taking strand into account
        args = ["awk '{OFS=\"\t\" ; if($5==\"+\") {print $1,$2-" + str(SIZE) + ",$2+" + str(SIZE) + ",$4} else {print $1,$3-" + str(SIZE) + ",$3+" + str(SIZE) + ",$4}}' " + META + " > flanks.bed"]
        process = sp.Popen(args, shell=True)
        process.wait()

    else: #No strand info so assume everything is + strand
        args = ["awk '{OFS=\"\t\" ; print $1,$2-" + str(SIZE) + ",$2+" + str(SIZE) + ",$4}' " + META + " > flanks.bed"]
        process = sp.Popen(args, shell=True)
        process.wait()

else: #With any other feature other than gene, strand does not matter since the flanks are taken from the midpoint

    strand=False
    args = ["awk '{OFS=\"\t\" ; if(($2+$3) % 2==0) {a=int(($2+$3)/2); print $1,a-" + str(SIZE) + ",a+" + str(SIZE) + ",$4} else {a=int(($2+$3)/2+rand()); print $1,a-" + str(SIZE) + ",a+" + str(SIZE) + ",$4}}' " + META + " > flanks.bed" ]   
    process = sp.Popen(args, shell=True)
    process.wait()


#Retain Coordinates from META source file after removing features at beginning of chr
flanks_df = pd.read_csv("flanks.bed",header=None, sep ="\t", index_col=None) #Read in flanks file 
temp = flanks_df.merge(meta_df,on=flanks_df.index) #Merge source META file and flanks. They are in same order
temp = temp[temp['1_x'] > 0] #Remove rows that have a negative flank coordinate (beginning of chr)
temp2 = temp[['0_x','1_x','2_x',3]]
temp2.to_csv("flanks.bed", header=False, index=False, sep="\t")
print(len(temp2))
feature = temp[3].tolist() #Retain feature order in a list
feature_arr = np.array(feature) #Convert feature list to array
total_features = len(temp.axes[0]) #Retain total number of individual features
coord_df = temp[['0_y','1_y','2_y']] #Retain source file coordinates
feature_exclusive = list(set(feature)) #Determine number of unique feature types by removing duplicates from total occurence feature list
feature_exclusive = sorted(feature_exclusive) #Sort alphabetically
del temp
del temp2
del feature
del flanks_df

if(strand==True):
    feature_exclusive = []
    feature_exclusive.append("Genocode_Annotations")

args = ["bedtools coverage -a flanks.bed -b " + BED + " -d > flanks_1bp_cov.bed"] #Map single base coverage
process = sp.Popen(args, shell=True)
process.wait()

df_scale = pd.read_csv(BED,header=None, sep ="\t", index_col=None) 
print(len(df_scale))
scale_factor = 1000000/float(len(df_scale))
print(scale_factor)
del df_scale
###########################################################################################################
#Read coverage data of flanks windows of features and tranpose

df = pd.read_csv("flanks_1bp_cov.bed",header=None, sep ="\t", index_col=False)
df.drop([0,1,2,3,4], axis=1,inplace=True) #Remove all other columns but coverage counts

arr = np.array(df[5]) #Make array of raw counts
arr2 = np.reshape(arr, (int(total_features), int(SIZE)*2)) # Transpose raw counts where (rows=total feature count,columns=flanking window size i.e. total positions)
df = pd.DataFrame(arr2)


###########################################################################################################
#bedtools calculates coverage upstream to into the gene body for +
#bedtools calculates coverage gene body to upstream for -
#Reverse rows for negative strand genes so positions are uniform
if(strand==True):
    df["strand"] = strand_list
    df["gene"] = feature_arr
    df = df.sort_values('strand')

    plus = df[df["strand"] == "+"] #Separate positive strand and negative strand genes

    minus = df[df["strand"] == "-"]
    minus_feature = minus["gene"].tolist() #Retain gene name order in list
    minus = minus[minus.columns[::-1]] #Reverse all rows in array so position is uniform with positive strand
    minus.drop(["strand","gene"], axis=1,inplace=True) #Remove other columns
    minus.columns = range(minus.shape[1]) #Reset header of reversed array
    minus["gene"] = minus_feature #Add gene order back to negative strand array
    frames = [plus,minus] 
    df = pd.concat(frames) #Rejoin positve and negative strand genes
    feature = df["gene"].tolist() #Retain total gene order in list
    feature_arr = np.array(feature) #Convert gene list to array
    df.drop(["strand","gene"], axis=1,inplace=True) #Remove strand and gene columns. Not needed anymore


##########################################################################################################
#Feature Counts. Genes for TSS and Cooridnates for other features
df['Sum'] = df.sum(axis=1) #Caluclate the sum of raw reads for each row. (total count for a feature)
df.set_index(feature_arr,inplace=True) #Insert the features array as index so rows are labeled

counts_df = df[['Sum']] #Create new DF with only sums and features
counts_df.reset_index(inplace=True) #Keep index column (features) but remove it as index for further manipulation
if(strand!=True): #Add coordinates to corresponding features if not genes and then sort
    feature_counts_df = coord_df.join(counts_df) #add genomic coordinates for features if not genes
    feature_counts_df.sort_values('Sum',ascending=False,inplace=True)
else:
    feature_counts_df = counts_df.sort_values('Sum',ascending=False,inplace=False) #For genes keep everyhting the same and sort
feature_counts_df.to_csv("feature_counts.csv", header=False, index=False, sep="\t") #Output the total raw counts for each feature


###########################################################################################################
#Position counts genome wide for all features, and individual features if features file is not genes
position = []
num = -int(SIZE)
for i in range(int(SIZE)*2): #Create position list (-size to +size with 0 being either tss or midpoint of feature)
    position.append(num)
    num = num + 1 
total = []
total = df.mean() #Sum each column (position) and put in list
combined = list(zip(position, total)) #Create 2D list joining the position and its corresponding counts
feature_pos_overlay = pd.DataFrame(combined)

with open("position_counts.csv", "w") as f: #Output position and count mapping for all features
    writer = csv.writer(f)
    writer.writerows(combined)

top = feature_pos_overlay[1].head(100)
bottom = feature_pos_overlay[1].tail(100)
edges = top.append(bottom)
avg = float(round(edges.mean(),2))

if(avg>0):
    score = float(round(feature_pos_overlay[1].max() / float(avg),2))
    score_min = float(round(feature_pos_overlay[1].min() / float(avg),2))
else:
    score = 0
    score_min = 0

if(int(len(feature_exclusive)) == 1):
    if(strand!=True):
        with open("feature_enrichment_scores.csv", "w") as f:
            f.write("Feature,100_Flank_Avg,Count_Max,Position_Max,Score_Max,Count_Min,Position_Min,Score_Min\n")
            f.write(str(feature_exclusive[0]) + ",")
            f.write(str(avg) + ",")
            f.write(str(feature_pos_overlay[1].max()) + ",")
            f.write(str(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]) + ",")
            f.write(str(score) + ",")
            f.write(str(feature_pos_overlay[1].min()) + ",")
            f.write(str(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]) + ",")
            f.write(str(score_min) + "\n")

    else:
        with open("feature_enrichment_scores.csv", "w") as f:
            f.write("Feature,100_Flank_Avg,Count_Max,Position_Max,Score_Max,Count_Min,Position_Min,Score_Min\n")
            f.write("Gencode_Annotations,")
            f.write(str(avg) + ",")
            f.write(str(feature_pos_overlay[1].max()) + ",")
            f.write(str(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]) + ",")
            f.write(str(score) + ",")
            f.write(str(feature_pos_overlay[1].min()) + ",")
            f.write(str(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]) + ",")
            f.write(str(score_min) + "\n")

if(strand!=True and int(len(feature_exclusive)) > 1):
    feature_pos_overlay[2] = "All"
    with open("feature_enrichment_scores.csv", "w") as f:
        f.write("Feature,100_Flank_Avg,Count_Max,Position_Max,Score_Max,Count_Min,Position_Min,Score_Min\n")
        f.write("All Features - ")
        f.write(str(', '.join(feature_exclusive)) + ",")
        f.write(str(avg) + ",")
        f.write(str(feature_pos_overlay[1].max()) + ",")
        f.write(str(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]) + ",")
        f.write(str(score) + ",")
        f.write(str(feature_pos_overlay[1].min()) + ",")
        f.write(str(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]) + ",")
        f.write(str(score_min) + "\n")

    for i in range(len(feature_exclusive)): #Output position and count mapping for individual features
        feature_temp = df[df.index == str(feature_exclusive[i])] #isolate feature from feature list with no duplicates
        total = []
        total = feature_temp.sum() #Sum counts for each position i.e. sum counts along each column
        combined = list(zip(position, total)) #Create 2D list joining the position and its corresponding counts
        temporary = pd.DataFrame(combined) #convert to DF

        #Compute enrichment score for each feature
        top = temporary[1].head(100) #get first and last 100 base pairs counts
        bottom = temporary[1].tail(100)
        edges = top.append(bottom)
        avg = float(round(edges.mean(),2)) #get average of those 200 basepairs
        if(avg>0):
            score = float(round(temporary[1].max() / float(avg),2)) #divide the peak count by the 200 base pair average
            score_min = float(round(temporary[1].min() / float(avg),2))
        else:
            score = 0
            score_min = 0
        with open("feature_enrichment_scores.csv", "a") as f: #Output enrichment score
            f.write(str(feature_exclusive[i]) + ",")
            f.write(str(avg) + ",")
            f.write(str(temporary[1].max()) + ",")
            f.write(str(temporary[0][temporary[1].idxmax()]) + ",")
            f.write(str(score) + ",")
            f.write(str(temporary[1].min()) + ",")
            f.write(str(temporary[0][temporary[1].idxmin()]) + ",")
            f.write(str(score_min) + "\n")

        temporary[2] = str(feature_exclusive[i]) #Add column identifier for feature for plotly
        frames = [feature_pos_overlay,temporary] #append individual feature counts to overall DF
        feature_pos_overlay = pd.concat(frames) 
        with open(str(feature_exclusive[i]) + "_position_counts.csv", "w") as f: #Output file with feature as identifier
            writer = csv.writer(f)
            writer.writerows(combined)


#############################################################################################################
#Plots
if(strand!=True and int(len(feature_exclusive)) > 1):
    fig = px.scatter(feature_pos_overlay, y=1, x=0,color=2,title="Multiple Feature Comparison")
    fig.update_yaxes(title_text = "Raw Reads")
    fig.update_xaxes(title_text = "Position")
    fig.write_image("Feature.png")

    feature_pos_overlay[1] = feature_pos_overlay[1] * scale_factor

    fig = px.scatter(feature_pos_overlay, y=1, x=0,color=2,title="Multiple Feature Comparison")
    fig.update_yaxes(title_text = "CPM Reads")
    fig.update_xaxes(title_text = "Position")
    fig.write_image("Feature_cpm.png")

else:
    temp = pd.read_csv("position_counts.csv",header=None,index_col=None)
    fig = px.scatter(temp, y=1, x=0, title="Single Feature")
    fig.update_yaxes(title_text = "Raw Reads")
    fig.update_xaxes(title_text = "Position")
    fig.add_vline(x=int(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]), line_dash="dash" , line_color="black",annotation_text="MAX " + str(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]), annotation_position="top right")
    fig.add_vline(x=int(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]), line_dash="dash" , line_color="grey",annotation_text="MIN " + str(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]), annotation_position="top right")
    fig.write_image("Feature.png")

    
    print(temp)
    temp[1] = temp[1] * scale_factor
    print(temp)
    fig = px.scatter(temp, y=1, x=0, title="Single Feature")
    fig.update_yaxes(title_text = "CPM Reads")
    fig.update_xaxes(title_text = "Position")
    fig.add_vline(x=int(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]), line_dash="dash" , line_color="black",annotation_text="MAX " + str(feature_pos_overlay[0][feature_pos_overlay[1].idxmax()]), annotation_position="top right")
    fig.add_vline(x=int(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]), line_dash="dash" , line_color="grey",annotation_text="MIN " + str(feature_pos_overlay[0][feature_pos_overlay[1].idxmin()]), annotation_position="top right")
    fig.write_image("Feature_cpm.png")
