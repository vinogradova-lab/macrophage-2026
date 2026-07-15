import pandas as pd 
import seaborn as sb
import numpy as np
import matplotlib.pylab as plt
from scipy.stats.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')
from sklearn.decomposition import PCA
import plotly.express as px
from pathlib import Path
import plotly.graph_objects as go
from plotly.graph_objects import Layout
from matplotlib.patches import Rectangle
from scipy.stats import linregress
import scipy.stats as stat
from matplotlib.ticker import MaxNLocator
# goatools (GODag / GoSearch / Gene2GoReader) was removed when GO curation moved to MSigDB
# 2026.1.Hs; see query_gene_ontology_msigdb below and whole_proteome/export_msigdb_go.R.
import re
from functools import reduce
import math
import seaborn as sns
import sys
sys.path.append("..")
from src.pca_utils import run_pca, percent_explained_df

FC_CUTOFF = np.log2(1.5)

def get_pca_plot(df, title_string, with_loadings = False): #df needs to be without log

    pepnum_list = df.filter(regex="pepNum").columns.tolist()
    df = df.drop(["uniprot", "description"] + pepnum_list, axis = 1)

    #drop na
    df = df.dropna()
    #number_of_proteins_in_common = df.shape[0]
    df = df.set_index("protein")

    #drop TLR9 channels before PCA (channels are named "<condition>_...")
    df = df.drop(columns=df.filter(regex=r"^TLR9").columns)

    #run centralized PCA core (log2 -> drop non-finite rows -> transpose -> PCA(2))
    pcs_df, loadings_df, explained_variance_ratio = run_pca(df, n_components=2)
    number_of_proteins_in_common = loadings_df.shape[0]

    #rebuild finalDf with the historical CSV schema (read by the R pca_plot())
    finalDf = pcs_df.rename(columns={
        "PC1": "principal component 1",
        "PC2": "principal component 2",
        "sample_name": "channel_name",
    })

    #get PCA %
    pca_1_percent = round((explained_variance_ratio[0] * 100),2)
    pca_2_percent = round((explained_variance_ratio[1] * 100),2)

    #add condition
    finalDf["condition"] = finalDf["channel_name"]
    #finalDf["condition"] = finalDf["condition"].str[:-7]
    #finalDf["stimulation"] = finalDf["condition"].str.split("_").str[1]
    
    
    finalDf["condition"] = finalDf["condition"].str.split("_").str[0]
    color_discrete_sequence_list = ["salmon", "#A3A533", "#56BC82", "#4EADF0", "#D772EC"]
    title_info = title_string+" ("+str(number_of_proteins_in_common)+" Proteins)"
    arrow_increase = 120
        
    finalDf = finalDf.sort_values(by=['condition'])
    
    #get range
    range_list = []

    range_list.append(finalDf["principal component 1"].max())
    range_list.append(finalDf["principal component 1"].min())

    range_list.append(finalDf["principal component 2"].max())
    range_list.append(finalDf["principal component 2"].min())

    max_number = max(range_list)+2
    min_number = min(range_list)-2
    
    #print(finalDf)
    #symbols = ["circle", "triangle-up"]
    
    fig = px.scatter(finalDf, 
                     x='principal component 1', 
                     y='principal component 2', 
                     hover_data=['channel_name'], 
                     color=finalDf['condition'],
                     #symbol_sequence = symbols,
                     #symbol=finalDf['stimulation'],
                     #color_discrete_sequence=color_discrete_sequence_list,
                     labels={
                         "principal component 1": "PC1 ({}%)".format(pca_1_percent),
                         "principal component 2": "PC2 ({}%)".format(pca_2_percent),
                         "condition": "Condition"},
                     title=title_info,
                     #color_discrete_map=color_discrete_map,
                     template="plotly_white",
                     )
    

    fig.update_xaxes(dtick=5, range=[min_number, max_number]) #, range=[min_number, max_number]) dtick=4
    fig.update_yaxes(dtick=5, range=[min_number, max_number]) #, range=[min_number, max_number]) dtick=4
    
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    fig.update_xaxes(showgrid=True, gridwidth=2, gridcolor='#E8E8E8')
    fig.update_yaxes(showgrid=True, gridwidth=2, gridcolor='#E8E8E8')
    
    fig.update_xaxes(zerolinewidth=2, zerolinecolor='#E8E8E8')
    fig.update_yaxes(zerolinewidth=2, zerolinecolor='#E8E8E8')
    
    fig.update_traces(marker=dict(size=12),
                      selector=dict(mode='markers')) #line=dict(width=2,color='DarkSlateGrey')),
    
    
    fig.update_layout(height=500, width=600, showlegend=True)
    
    fig.update_layout(legend_title_text='')
    
    #loadings_df (PC1, PC2, variable) is returned by run_pca above
    #get one laodings worth, rememeber sum of squares per PC loading = 1
    one_loading_value = math.sqrt(1/len(loadings_df))
    #then filter out ones that contribute more than one variables worth 
    filtered_loadings = loadings_df[(loadings_df["PC1"] > one_loading_value) | (loadings_df["PC1"] < -one_loading_value) | (loadings_df["PC2"] > one_loading_value) | (loadings_df["PC2"] < -one_loading_value)]
    
    #filtered_loadings = filtered_loadings.head(10)
    
    top_pos_PC1_df = filtered_loadings.sort_values('PC1', ascending=False).head(20)
    top_neg_PC1_df = filtered_loadings.sort_values('PC1').head(20)
    
    top_pos_PC2_df = filtered_loadings.sort_values('PC2', ascending=False).head(20)
    top_neg_PC2_df = filtered_loadings.sort_values('PC2').head(20)
    
    df_list = [top_pos_PC1_df, top_neg_PC1_df, top_pos_PC2_df, top_neg_PC2_df]
    
    merged_loadings_df = reduce(lambda x, y: pd.merge(x, y, on = ["variable", "PC1", "PC2"],  how = "outer"), df_list)
    
    if with_loadings == True:
        n = merged_loadings_df.shape[0]
        for i in range(n):
            fig.add_annotation(x= 0,
                               y= 0,
                               ax=merged_loadings_df.iloc[i,0] * arrow_increase,
                               ay=merged_loadings_df.iloc[i,1] * arrow_increase,
                               xref='x',
                               yref='y',
                               axref='x',
                               ayref='y',
                               text=merged_loadings_df.iloc[i,2],
                               showarrow=True,
                               arrowhead=3,
                               arrowsize=1,
                               arrowwidth=1,
                               arrowcolor='red',
                               opacity=0.6,
                               arrowside='start')
    
    #plt.figure(figsize=(5,20))
    #sb.heatmap(loadings_df.set_index("variable"),
    #           cmap='YlGnBu',
    #           #linewidths=0.7,
    #           linecolor="black").set(title='Loadings', ylabel=None)
    #plt.show()
    
    loadings_df = loadings_df.set_index("variable")

    percent_df = percent_explained_df(explained_variance_ratio)

    return(fig, loadings_df, finalDf, percent_df)

color_discrete_map = {'Significant Up': '#ff8080', 
                      'Not Significant Up': '#ffcccc', 
                      'Significant but <1.5 FC': 'darkgrey', 
                      'Not Significant': 'lightgrey',
                      'Significant Down': '#71A0C6',
                      'Not Significant Down': '#CDDEFA' 
                    }

color_heatmap_up = 240
color_heatmap_down = 15 

def get_p_value(row, cond_1, cond_2):
    ttest_result = stat.ttest_ind(row[cond_1],  row[cond_2], nan_policy='omit')
    return ttest_result[1]

def get_expr(row):
    p_value_column = '-log10_pval'
    p_value_cutoff = -1*math.log10(0.05)
    FC_CUTOFF = math.log2(1.5)
    if (row['log2_FC'] > FC_CUTOFF) & (row[p_value_column] > p_value_cutoff):
        return "Significant Up"
    if (row['log2_FC'] > FC_CUTOFF) & (row[p_value_column] < p_value_cutoff):
        return "Not Significant Up"
    if (row['log2_FC'] < -FC_CUTOFF) & (row[p_value_column] > p_value_cutoff):
        return "Significant Down"
    if (row['log2_FC'] < -FC_CUTOFF) & (row[p_value_column] < p_value_cutoff):
        return "Not Significant Down"
    if (row['log2_FC'] > -FC_CUTOFF) & (row['log2_FC'] < FC_CUTOFF) & (row[p_value_column] > p_value_cutoff):
        return "Significant but <1.5 FC"
    else:
        return "Not Significant"
def annotate_functional_group(row, column_name, protein_set):
    row[column_name] = row.name in protein_set
    return row

def annotate_functional_categories(volcano_df):
    """ Annotate the volcano_data plot with functional categories. See 'protein_lists/README.md' for 
    further documentation
    """
    fun_groups = pd.read_csv("protein_lists/group_annotations.csv")
    for i,row in fun_groups.iterrows():
        column_name = row["column_name"]
        file_name = row["file_name"]
        protein_table = pd.read_csv(file_name)
        functional_proteins = set(protein_table["uniprot"])
        volcano_df = volcano_df.apply(
            annotate_functional_group, 
            column_name=column_name, 
            protein_set = functional_proteins, axis=1)
    return volcano_df

def add_peroxisome_annotation(volcano_df):

    obodag = GODag("go-basic.obo") # go-basic.obo
    obodag

    gene2go = download_ncbi_associations()

    objanno = Gene2GoReader("gene2go", taxids=[9606])

    go2geneids_human = objanno.get_id2gos(namespace='CC', go2geneids=True)

    print("{N:} GO terms associated with human NCBI Entrez GeneIDs".format(N=len(go2geneids_human)))

    srchhelp = GoSearch("go-basic.obo", go2items=go2geneids_human)

    # Compile search pattern for 'peroxisome'
    peroxisome_all = re.compile(r'peroxisom', flags=re.IGNORECASE)
    peroxisome_not = re.compile(r'peroxisome.independent', flags=re.IGNORECASE)

    # Find ALL GOs and GeneIDs associated with 'peroxisome'.

    # Details of search are written to a log file
    fout_allgos = "peroxisome_gos_human.log" 
    with open(fout_allgos, "w") as log:
        # Search for 'cell cycle' in GO terms
        gos_cc_all = srchhelp.get_matching_gos(peroxisome_all, prt=log)
        # Find any GOs matching 'peroxisome-independent'
        gos_no_cc = srchhelp.get_matching_gos(peroxisome_not, gos=gos_cc_all, prt=log)
        # Remove GO terms that are not "cell cycle" GOs
        gos = gos_cc_all.difference(gos_no_cc)
        # Add children GOs of cell cycle GOs
        gos_all = srchhelp.add_children_gos(gos)
        # Get Entrez GeneIDs for cell cycle GOs
        geneids = srchhelp.get_items(gos_all)
    print("{N} human NCBI Entrez GeneIDs related to 'peroxisome' found.".format(N=len(geneids)))

    from genes_ncbi_9606_proteincoding import GENEID2NT
    peroxisome_genes = []
    for geneid in geneids: # geneids associated with cell-cycle
        nt = GENEID2NT.get(geneid, None)
        if nt is not None:
            peroxisome_genes.append({"gene":nt.Symbol,
                                    "description":nt.description})

    peroxisome_genes = pd.DataFrame(peroxisome_genes)

    protein_col = [col for col in volcano_df.columns if 'protein' in col][0]

    volcano_df["peroxisome"] = volcano_df[protein_col].isin(peroxisome_genes["gene"])
    peroxisome_genes.to_csv("protein_lists/peroxisome_genes.csv")
    
    return volcano_df

def get_volcano_plot_treatment_vs_control(conditions_list, control_labelling, df, file_name, folder_path):
    volcano_plot_list = []
    volcano_df_list = []
    long_format_list = []
    index_cols = ["uniprot", "protein", "description"]
    df = df.set_index(index_cols)

    for condition in conditions_list:
        copy_df = df.copy()
        list_cond_1 = copy_df.filter(like=condition).columns.tolist()
        list_cond_2 = copy_df.filter(like=control_labelling).columns.tolist()
        if len(list_cond_1) <= 1 or len(list_cond_2) <= 1:
            print(f"Condition {condition} does not have enough replicates to be shown in volcano plot!")
            print(list_cond_1)
            print(list_cond_2)
            continue

        for labelling in [control_labelling, condition]:
            list_columns_labelling = copy_df.filter(like=labelling).columns.tolist()
            labelling_df = copy_df[list_columns_labelling]
            copy_df["mean_"+labelling] = labelling_df.mean(axis=1)
            
        copy_df["FC"] = copy_df["mean_"+condition] / copy_df["mean_"+control_labelling]
        x_axis_name = "log2(" + condition + "/" + control_labelling + ")"

        idx_cond_1 = copy_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = copy_df.columns.get_indexer(list_cond_2)
        copy_df["p_value"] = copy_df.apply(get_p_value, axis=1, args=(idx_cond_1, idx_cond_2), )
        
        copy_df["log2_FC"] = np.log2(copy_df["FC"])
        volcano_df = copy_df[["p_value", "log2_FC"]]
            
        volcano_df = volcano_df.dropna() 
        volcano_df["-log10_pval"] = -1*np.log10(volcano_df["p_value"])
        # Using BH correction here (previous versions used Boneferroni)
        volcano_df["-log10_pval_adj"] = -1*np.log10(stat.false_discovery_control(volcano_df["p_value"]))
        volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)
            
        volcano_df["Regulation"] = volcano_df["Regulation"].astype('category')
        
        volcano_df = volcano_df.reset_index()
        #volcano_df = volcano_df.set_index(["annotation"])
        #volcano_df["name"] = volcano_df['description'].str.split(" ").str[0]
        volcano_df["name"] = volcano_df["protein"]

        title_name = file_name + " - " + condition + " vs. " + control_labelling + " (" + str(len(volcano_df)) + " Proteins)"

        fig = px.scatter(volcano_df, 
                         x='log2_FC', 
                         y='-log10_pval',
                         color='Regulation',
                         color_discrete_map=color_discrete_map,
                         #color_discrete_sequence=colors_volcano,
                         hover_data=['name', 'uniprot'],
                         title = title_name.replace("processed_census-out_", ""),
                         labels = {"log2_FC": x_axis_name},
                         template = "simple_white",
                         category_orders={'Regulation': np.sort(volcano_df['Regulation'].unique())})
        fig.add_vline(x=FC_CUTOFF, line_width=2, line_dash="dash", line_color="grey")
        fig.add_vline(x=-FC_CUTOFF, line_width=2, line_dash="dash", line_color="grey")
        fig.add_hline(y=1.3, line_width=2, line_dash="dash", line_color="grey")
        fig.update_layout(legend=dict(title=""), title_x=0.5, font_family="Arial")
        
        fig.write_image(folder_path / "volcano_plots" / ("volcano_plot_" + title_name.split(" - ")[1].split(" (")[0] + ".svg"), engine="kaleido")

        sign_up_df = volcano_df.loc[volcano_df["Regulation"] == "Significant Up"].sort_values(by="-log10_pval", ascending=False).head(30)
        sign_up_df_fc = volcano_df.loc[volcano_df["Regulation"] == "Significant Up"].sort_values(by="log2_FC", ascending=False).head(30)
        sign_down_df = volcano_df.loc[volcano_df["Regulation"] == "Significant Down"].sort_values(by="-log10_pval", ascending=False).head(30)
        sign_down_df_fc = volcano_df.loc[volcano_df["Regulation"] == "Significant Down"].sort_values(by="log2_FC", ascending=True).head(30)
        labels_df = pd.concat([sign_up_df, sign_down_df, sign_up_df_fc, sign_down_df_fc], axis=0)
        labels_df = labels_df.drop_duplicates()
        labels_df = labels_df.loc[~np.isinf(labels_df["log2_FC"])]
        for i,r in labels_df.iterrows():
            if r['Regulation'] == 'Significant Down':
                color = color_discrete_map['Significant Down'] 
            elif r['Regulation'] == 'Significant Up':
                color = color_discrete_map['Significant Up'] 

        
            fig.add_annotation(x=r['log2_FC'],
                               y=r["-log10_pval"],
                               text= r['name'], 
                               showarrow=False,
                               xanchor='center',
                               yanchor='bottom',
                               font=dict(size=10, color=color))

        volcano_plot_list.append(fig)
        volcano_df = volcano_df.reset_index().set_index(index_cols)
        #volcano_df = volcano_df.drop(["description", "name"], axis=1)
        long_df = volcano_df.copy()
        long_df["condition"] = condition
        long_format_list.append(long_df)
        volcano_df = volcano_df.add_suffix("_"+title_name)
        volcano_df_list.append(volcano_df)
        
    with open(folder_path / "volcano_plots" / 'volcano_plots.html' , 'w') as f:
        for fig in volcano_plot_list:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

    volcano_df = pd.concat(volcano_df_list, axis=1)
    long_df = pd.concat(long_format_list, axis=0)
            
    return volcano_df, copy_df, long_df

def get_heatmap(df, volcano_df, file_name, folder_path):
    df = df.set_index(["uniprot", "protein", "description"])
    cols_list = df.columns.tolist()
    df = df.reset_index() #.drop(["pep_num", "annotation"], axis=1) #.set_index("description")
    df.description = df.description.str.split(" ").str[0]
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)

    volcano_cols = []

    # exclude some metadata columns from this analysis
    for col in volcano_df.columns.tolist():
        if " - " in col:
            volcano_cols.append(col)
    cols_list = [col.split(" - ")[1].split(" (")[0] for col in volcano_cols]
    unique_col_list = list(set(cols_list))

    uniprot_list = []

    volcano_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    volcano_df.dropna(inplace=True)

    for comparison in unique_col_list: 
        #print(comparison)
        test_df = volcano_df.filter(like=comparison)
        reg_col = [col for col in test_df if col.startswith('Regulation_')]
        fc_col = [col for col in test_df if col.startswith('log2_FC_')]
        pval_col = [col for col in test_df if col.startswith('-log10_pval_adj_')]

        sign_up_df = test_df.loc[test_df[reg_col[0]] == "Significant Up"].sort_values(by=pval_col[0], ascending=False).head(30)
        sign_up_df_fc = test_df.loc[test_df[reg_col[0]] == "Significant Up"].sort_values(by=fc_col[0], ascending=False).head(30)
        sign_down_df = test_df.loc[test_df[reg_col[0]] == "Significant Down"].sort_values(by=pval_col[0], ascending=False).head(30)
        sign_down_df_fc = test_df.loc[test_df[reg_col[0]] == "Significant Down"].sort_values(by=fc_col[0], ascending=True).head(30)
        labels_df = pd.concat([sign_up_df, sign_down_df, sign_up_df_fc, sign_down_df_fc], axis=0)
        labels_df = labels_df.drop_duplicates()
        #print(len(labels_df))
        for uniprot_id in labels_df.index: 
            if uniprot_id not in uniprot_list: 
                uniprot_list.append(uniprot_id)

    df = df.loc[df.uniprot.isin(uniprot_list)]
    df = df.set_index("uniprot")
    subset_for_heatmap_merged = df.join(volcano_df, how="left")

    df = df.reset_index("uniprot")
    df = df.drop(["uniprot", "description"], axis=1).set_index("protein")
    df.columns = df.columns.str.split("_processed").str[0]

    g = sns.clustermap(df, z_score=0, cmap=sns.diverging_palette(color_heatmap_up, color_heatmap_down, s=60, as_cmap=True), center=0, figsize=(18,int(len(df) / 5 )), yticklabels=True)
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
    g.ax_heatmap.set_title(file_name) 
    g.ax_cbar.set_ylabel("z-score",size=15)
    g.ax_cbar.set_position((0.1, .2, .03, .4))

    heatmap_path = folder_path / "heatmap" 
    Path(heatmap_path).mkdir(exist_ok=True)
    g.savefig(heatmap_path/ "heatmap_signupdown_cond.pdf", dpi=400)
    subset_for_heatmap_merged.to_csv(folder_path / "heatmap" / "heatmap_data.csv") 
    return "done"

def corr_plot_wp(df, condition_1, condition_2, control, corr_folder):
    df = df.reset_index()
    #get fold change 
    df["FC_"+condition_1] = df[condition_1] / df[control]
    df["FC_"+condition_2] = df[condition_2] / df[control]
    
    #get log fold change
    df["logFC"+condition_1] = np.log2(df["FC_"+condition_1])
    df["logFC"+condition_2] = np.log2(df["FC_"+condition_2])

    #this will be the x and y labeling
    condition_1_string = "LFC (" + condition_1 + " vs. " + control + ")"
    condition_2_string = "LFC (" + condition_2 + " vs. " + control + ")"
    
    #add categorial column to identify values that are >=2 fold up or down regulated
    #we call the new column "fold_change"
    conditions = [
                (df["logFC"+condition_1] <= -FC_CUTOFF) & (df["logFC"+condition_2] <= -FC_CUTOFF),
                (df["logFC"+condition_1] >=FC_CUTOFF) & (df["logFC"+condition_2] >=FC_CUTOFF),
                (df["logFC"+condition_1] <=FC_CUTOFF) & (df["logFC"+condition_2] >=FC_CUTOFF),
                (df["logFC"+condition_1] >= -FC_CUTOFF) & (df["logFC"+condition_2] <= -FC_CUTOFF),
                (df["logFC"+condition_1] >=FC_CUTOFF) & (df["logFC"+condition_2] <=FC_CUTOFF),
                (df["logFC"+condition_1] <= -FC_CUTOFF) & (df["logFC"+condition_2] >= -FC_CUTOFF)]
                

    choices = [">= 1.5 fold down reg",
               ">= 1.5 fold up reg", 
               "interesting up reg "+condition_2, 
               "interesting down reg "+condition_2,
               "interesting up reg "+condition_1, 
               "interesting down reg "+condition_1,
              ]
    
    color_discrete_map = {
        ">= 1.5 fold down reg":'#fc8961',
        ">= 1.5 fold up reg":'#1d1147', 
        "interesting up reg "+condition_2:'#bade28', 
        "interesting down reg "+condition_2:'#832681',
        "interesting up reg "+condition_1:'#e75263', 
        "interesting down reg "+condition_1:'#fec488',
        "nc":'darkgray'

    }

    df['fold_change'] = np.select(conditions, choices, default="nc")
    
    #save df 
    df.to_csv(corr_folder / ("corr_" + condition_1+" vs. "+condition_2+".csv"))
    
    #remove inf values 
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=["logFC"+condition_1, "logFC"+condition_2])

    r_pval = pearsonr(df["logFC"+condition_1].tolist(), df["logFC"+condition_2].tolist())

    if round(r_pval[1],3) == 0.0:
        p_val = 0.001
    else:
        p_val = round(r_pval[1], 3)
    

    x_max = np.max(df["logFC"+condition_1]) + 1
    x_min = np.min(df["logFC"+condition_1]) - 1
    y_max = np.max(df["logFC"+condition_2]) + 1
    y_min = np.min(df["logFC"+condition_2]) - 1
    
    fig = px.scatter(data_frame=df,
        x="logFC"+condition_1, 
        y="logFC"+condition_2, 
        color = 'fold_change',
        hover_data=["uniprot","protein","description"],
        template = "simple_white",
        color_discrete_map=color_discrete_map ,
        labels={"logFC"+condition_1:condition_1_string,
                "logFC"+condition_2:condition_2_string})



    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    fig.add_hline(FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")
    fig.add_hline(-FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")
    fig.add_vline(FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")
    fig.add_vline(-FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")

    
    fig.add_trace(go.Scatter(
        x = [x_min * 0.9, 
             x_min * 0.9
             ],
        y = [y_max * 0.9, 
             y_max * 0.8
             ],
        mode = "text",
        text=[r"PCC = {}".format(round(r_pval[0],3)), 
              r"p-value < {}".format(p_val)
              ],
        showlegend=False    
    ))

    fig.show()
    fig.write_html(corr_folder / ("corr_" + condition_1+" vs. "+condition_2+".html"))
    
    return

def corr_plot_wp_rnaseq(file_path, condition, control, folder_path):
    #read in dataframe
    df = pd.read_csv(file_path).drop("Unnamed: 0", axis=1)
    df = (df.drop(list(df.filter(regex = 'rep')), axis = 1)
          .dropna(subset = ['log2FoldChange', condition, control]))
    
    plot_title = condition + " vs " + control + " (" + str(df.shape[0]) + " Proteins)"
    
    #get fold change
    df["FC"] = df[condition] / df[control]
    #get log2 fold change
    df["logFC_wp"] = np.log2(df["FC"])

    
    #this will be the x and y axis labelling 
    wp_string = "Whole proteome LFC (" + condition + " vs. " + control + ")"
    rnaseq_string = "RNASeq LFC (" + condition + " vs. " + control + ")"
    
    #add categorial column to identify values that are >=cutoff fold up or down regulated
    #we call the new column "fold_change"
    conditions = [
                (df['log2FoldChange'] <= -FC_CUTOFF) & (df['logFC_wp'] <= -FC_CUTOFF),
                (df['log2FoldChange'] >= FC_CUTOFF) & (df['logFC_wp'] >= FC_CUTOFF)
                ]

    choices = [">= 1.5 fold down reg",">= 1.5 fold up reg"]
    df['fold_change'] = np.select(conditions, choices, default="nc")
    
    #add another categorial column to identify proteins which dont change in gene expression but in protein expression
    conditions = [
                (abs(df['log2FoldChange']) <= FC_CUTOFF) & (abs(df['logFC_wp']) >= FC_CUTOFF),
                (abs(df['logFC_wp']) <= FC_CUTOFF) & (abs(df['log2FoldChange']) >= FC_CUTOFF),
                ((df['logFC_wp'] >= FC_CUTOFF) & (df['log2FoldChange'] <= -FC_CUTOFF)) | ((df['logFC_wp'] <= -FC_CUTOFF) & (df['log2FoldChange'] >= -FC_CUTOFF)),
                ((df['logFC_wp'] >= FC_CUTOFF) & (df['log2FoldChange'] >= FC_CUTOFF)) | ((df['logFC_wp'] <= -FC_CUTOFF) & (df['log2FoldChange'] <= -FC_CUTOFF))
                ]

    choices = ["yellow", "violet", "teal", "dark_grey"]

    df['expression'] = np.select(conditions, choices, default="light_grey")
    
    
    #save df 
    df.to_csv(folder_path / (condition+" vs. "+control+".csv"))   

    r_pval = pearsonr(df["logFC_wp"].tolist(), df["log2FoldChange"].tolist())
    if round(r_pval[1],3) == 0.0:
        p_val = 0.001
    else:
        p_val = round(r_pval[1],3)
    
    #remove inf values 
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=['logFC_wp', 'log2FoldChange'])

    x_max = np.max(df["logFC_wp"]) + 1
    x_min = np.min(df["logFC_wp"]) - 1
    y_max = np.max(df["log2FoldChange"]) + 1
    y_min = np.min(df["log2FoldChange"]) - 1
    
    #print rows which are >=cutoff up or down regulated
    subdf = df[df["fold_change"] != "nc"]
    subdf = subdf[["protein", "logFC_wp", "log2FoldChange", "fold_change"]]
    subdf.columns = ["protein", "logFC_wp", "logFC_rnaseq", "fold_change"]
    
    color_discrete_map = {'yellow': '#bade28', 
                      'violet': '#440154', 
                      "teal" : "#21918c",
                      'light_grey': 'lightgrey',
                      'dark_grey': 'dimgrey'
                    }
    
    
    fig = (px.scatter(data_frame=df, 
                    x="logFC_wp", 
                    y="log2FoldChange", 
                    color='expression', 
                    template = "simple_white",
                    color_discrete_map=color_discrete_map,
                    title = plot_title, 
                    hover_data=["protein","uniprot","description"],
                    labels={"logFC_wp":wp_string,
                            "log2FoldChange":rnaseq_string})
                            .update_layout(showlegend=False))

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)



    fig.add_hline(FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")
    fig.add_hline(-FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")
    fig.add_vline(FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")
    fig.add_vline(-FC_CUTOFF, line_width=1, line_color="black",line_dash="dash")

    #labelling of data in scatterplot
    subset_df = df[df["expression"] != "nc"]
    subset_df = subset_df.reset_index()
    
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

    
    fig.add_trace(go.Scattergl(
        x = [x_min * 0.9, 
             x_min * 0.9
             ],
        y = [y_max * 0.9, 
             y_max * 0.8
             ],
        mode = "text",
        text=[r"PCC = {}".format(round(r_pval[0],3)), 
              r"p-value < {}".format(p_val)
              ],
        textposition="top right",
        showlegend=False    
    ))

    x_lbl_sprd = min(abs(x_min), abs(x_max)) * 0.98
    y_lbl_sprd = min(abs(y_min), abs(y_max)) * 0.98 

    fig.add_trace(go.Scattergl(
        x = [-x_lbl_sprd,x_lbl_sprd,-x_lbl_sprd,x_lbl_sprd,-x_lbl_sprd,x_lbl_sprd,0,0,0],
        y = [0,0,y_lbl_sprd,-y_lbl_sprd,-y_lbl_sprd,y_lbl_sprd,-y_lbl_sprd,y_lbl_sprd,0],
        mode = "text",
        text = [
            len(df[(df["expression"] == "yellow") & (df["logFC_wp"] < 0)]),
            len(df[(df["expression"] == "yellow") & (df["logFC_wp"] > 0)]),
            len(df[(df["expression"] == "teal") & (df["logFC_wp"] < 0)]),
            len(df[(df["expression"] == "teal") & (df["logFC_wp"] > 0)]),
            len(df[(df["expression"] == "dark_grey") & (df["logFC_wp"] < 0)]),
            len(df[(df["expression"] == "dark_grey") & (df["logFC_wp"] > 0)]),
            len(df[(df["expression"] == "violet") & (df["log2FoldChange"] < 0)]),
            len(df[(df["expression"] == "violet") & (df["log2FoldChange"] > 0)]),
            len(df[df["expression"] == "light_grey"])
        ],
        textfont = dict(color = ["#bade28","#bade28","#21918c","#21918c","dimgrey","dimgrey","#440154","#440154", "black"])
    ))

    
    fig.write_html(folder_path / (condition+" vs. "+control+ "_wp_vs_rna.html"))
    
    
    return

def sub_mito_local(protein, localization_dict):
    if protein in localization_dict:
        return localization_dict[protein]
    return "NA"

def add_mitocarta_localization(df):
    mitocarta = pd.read_csv("protein_lists/Human.MitoCarta3.0.csv")
    prot_col = [col for col in df if col.startswith('protein')][0]
    mitocarta_locals = {}
    for i,r in mitocarta.iterrows():
        mitocarta_locals[r["Symbol"]] = r["MitoCarta3.0_SubMitoLocalization"]
    df["MitoCarta3.0_SubMitoLocalization"] = df[prot_col].apply(sub_mito_local, localization_dict=mitocarta_locals)
    return df


import functools


@functools.lru_cache(maxsize=4)
def _load_msigdb_go(csv_path):
    """Load and cache the flat MSigDB C5 GO table exported from R
    (reference_dbs/msigdb_go_2026.csv, columns: ontology, GO_id, term, gene).
    This is the SAME 2026.1.Hs membership .go_msigdb() / go_enrich() use in R."""
    return pd.read_csv(csv_path)


def query_gene_ontology_msigdb(
    substrings,
    domains,
    output_path,
    msigdb_csv=None,
    map_genes_to_uniprot=None,
):
    """MSigDB 2026.1.Hs drop-in for query_gene_ontology.

    Each token in ``substrings`` resolves to GO ids: a literal GO id (``GO:\\d+``) is used
    directly; any other token is a case-insensitive substring matched against the term name
    OR GO definition. ``domains`` filters ontology in {BP, CC, MF}. MSigDB C5 sets are
    already DAG-propagated, so a GO id yields its full descendant membership directly.

    Writes ``01_go_terms_matching_query.csv`` / ``02_genes_matching_query.csv`` and returns a
    one-row-per-gene frame; ``map_genes_to_uniprot``, if given, adds the ``uniprot`` column.
    """
    output_path.mkdir(exist_ok=True, parents=True)
    if msigdb_csv is None:
        msigdb_csv = str(Path(__file__).resolve().parent.parent / "reference_dbs" / "msigdb_go_2026.csv")
    df = _load_msigdb_go(str(msigdb_csv))
    df = df[df["ontology"].isin(list(domains))]

    go_re = re.compile(r"^GO:\d+$")
    go_cols = ["search_substring", "search_domain", "id", "name", "genes", "n_genes"]
    rows = []  # one row per (search token, matched GO term)
    for substr in substrings:
        if go_re.match(str(substr)):
            sel = df[df["GO_id"] == substr]
        else:
            # name OR definition; synonyms/xrefs are the only goatools fields MSigDB lacks
            hit = df["term"].str.contains(str(substr), case=False, regex=False, na=False)
            if "description" in df.columns:
                hit = hit | df["description"].str.contains(
                    str(substr), case=False, regex=False, na=False
                )
            sel = df[hit]
        for (go_id, term, ont), grp in sel.groupby(["GO_id", "term", "ontology"]):
            genes = sorted(grp["gene"].dropna().unique().tolist())
            rows.append(
                {
                    "search_substring": substr,
                    "search_domain": ont,
                    "id": go_id,
                    "name": term,
                    "genes": genes,
                    "n_genes": len(genes),
                }
            )

    go_df = pd.DataFrame(rows, columns=go_cols)

    # one row per gene (union across matched terms)
    gene_df = (
        go_df.explode("genes")
        .dropna(subset=["genes"])
        .groupby("genes")
        .agg(list)
        .reset_index()
    )
    if map_genes_to_uniprot is not None and not gene_df.empty:
        gene_df = map_genes_to_uniprot(gene_df)

    gene_df.to_csv(output_path / "02_genes_matching_query.csv")
    go_df.to_csv(output_path / "01_go_terms_matching_query.csv")
    return gene_df