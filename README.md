## Step 0: Install all the needed packages listed below
```r
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("caret")
install.packages("catch")
install.packages("openxlsx")
install.packages("MASS")
install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("UpSetR")
install.packages("grid")
install.packages("gridSVG")
install.packages("ggplotify")
install.packages("factoextra")
install.packages("Rtsne")
install.packages("PupillometryR")
install.packages("cowplot")
install.packages("ggdist")
install.packages("remotes")
install.packages("patchwork")
remotes::install_github("davidsjoberg/ggsankey")
```

## Step 1: Load the libraries into your R
```r
library("pheatmap")
library("ggplot2")
library("caret")
library("catch")
library("openxlsx")
library("MASS")
library("ggsankey")
library("dplyr")
library("RColorBrewer")
library("UpSetR")
library("grid")
library("gridSVG")
library("ggplotify")
library("factoextra")
library("Rtsne")
library("PupillometryR")
library("cowplot")
library("ggdist")
library("patchwork")
```

## Step 2: Source the R scripts
```r
source("model/catch_demo.R")
```

## Step 3: Load the datasets  
```r
X_GI=read.csv(file = "data/Tensor_GI.csv",check.names = FALSE)
U_GI=read.csv(file = "data/Covariant_GI.csv")
Y_GI=read.csv(file = "data/Response_GI.csv", header = T)
X_MSI=read.csv(file = "data/Tensor_MSI.csv",check.names = FALSE)
U_MSI=read.csv(file = "data/Covariant_MSI.csv")
Y_MSI=read.csv(file = "data/Response_MSI.csv", header = T)
X_Lung=read.csv(file = "data/Tensor_Lung.csv",check.names = FALSE)
U_Lung=read.csv(file = "data/Covariant_Lung.csv")
Y_Lung=read.csv(file = "data/Response_Lung.csv", header = T)
pathway_df=read.xlsx("Pathway/Pathway.xlsx",sheet = 1)
Metabolite_pathway=read.xlsx("Pathway/Pathway.xlsx",sheet = 2)
pathway_GI=data.frame(pathway=pathway_df[pathway_df$gene %in% colnames(U_GI),]$pathway,row.names = pathway_df[pathway_df$gene %in% colnames(U_GI),]$gene)
pathway_LUNG=data.frame(pathway=pathway_df[pathway_df$gene %in% colnames(U_Lung),]$pathway,row.names = pathway_df[pathway_df$gene %in% colnames(U_Lung),]$gene)
pathway_MSI=data.frame(pathway=pathway_df[pathway_df$gene %in% colnames(U_MSI),]$pathway,row.names = pathway_df[pathway_df$gene %in% colnames(U_MSI),]$gene)
pathway_GI=pathway_GI[colnames(U_GI)[-1],,drop=F]
pathway_LUNG=pathway_LUNG[colnames(U_Lung)[-1],,drop=F]
pathway_MSI=pathway_MSI[colnames(U_MSI)[-1],,drop=F]
```

Perform data preprocessing and CATCH model  
```r
MSI_brain = fit_and_catch_model(X_MSI[,-1], Y_MSI$brain, U_MSI[,-1], "Metastasis", pathway_MSI)
MSI_lung = fit_and_catch_model(X_MSI[,-1], Y_MSI$lung, U_MSI[,-1], "Metastasis", pathway_MSI)
MSI_liver = fit_and_catch_model(X_MSI[,-1], Y_MSI$liver, U_MSI[,-1], "Metastasis", pathway_MSI)
MSI_bone = fit_and_catch_model(X_MSI[,-1], Y_MSI$bone, U_MSI[,-1], "Metastasis", pathway_MSI)
MSI_kidney = fit_and_catch_model(X_MSI[,-1], Y_MSI$kidney, U_MSI[,-1], "Metastasis", pathway_MSI)
Lung_brain = fit_and_catch_model(X_Lung[,-1], Y_Lung$brain, U_Lung[,-1], "Metastasis", pathway_LUNG)
Lung_lung = fit_and_catch_model(X_Lung[,-1], Y_Lung$lung, U_Lung[,-1], "Metastasis", pathway_LUNG)
Lung_liver = fit_and_catch_model(X_Lung[,-1], Y_Lung$liver, U_Lung[,-1], "Metastasis", pathway_LUNG)
Lung_bone = fit_and_catch_model(X_Lung[,-1], Y_Lung$bone, U_Lung[,-1], "Metastasis", pathway_LUNG)
Lung_kidney = fit_and_catch_model(X_Lung[,-1], Y_Lung$kidney, U_Lung[,-1], "Metastasis", pathway_LUNG)
GI_brain = fit_and_catch_model(X_GI[,-1], Y_GI$brain, U_GI[,-1], "Metastasis", pathway_GI)
GI_lung = fit_and_catch_model(X_GI[,-1], Y_GI$lung, U_GI[,-1], "Metastasis", pathway_GI)
GI_liver = fit_and_catch_model(X_GI[,-1], Y_GI$liver, U_GI[,-1], "Metastasis", pathway_GI)
GI_bone = fit_and_catch_model(X_GI[,-1], Y_GI$bone, U_GI[,-1], "Metastasis", pathway_GI)
GI_kidney = fit_and_catch_model(X_GI[,-1], Y_GI$kidney, U_GI[,-1], "Metastasis", pathway_GI)
```

Generate supplementary materials: Each R function is used to generate a supplementary file in CSV or XLSX format.
```r
make_supplementary1(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney)
make_supplementary2(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney)
make_supplementary4(pathway_MSI=pathway_MSI,pathway_GI=pathway_GI,pathway_LUNG=pathway_LUNG)
make_supplementary7(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,pathway_df=pathway_df)
make_supplementary5(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI)
make_table_3(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney)
before_adj_matrix=Before_adjusted_matrix(X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,metabolite_name_trans=Metabolite_pathway$Metabolite_ori)

```

Generate the figures: Each R function is used to generate a figure.
```r
Sankey(df_table2_metastasis_4 = read.csv("supplementary/Figure_1_Sankey_Supplementary_Data_4.csv"))
make_heatmap(df_table2 = read.csv("supplementary/Figure_2and3_Upset_Plot_Supplementary_Data_4.csv"))
UpSet_plot(df_table2 = read.csv("supplementary/Figure_2and3_Upset_Plot_Supplementary_Data_4.csv"))
make_heatmap2(metabolite_225=Metabolite_pathway$Metabolite_trans,meta_225_col=Metabolite_pathway,metabolite_name_trans=Metabolite_pathway$Metabolite_ori,supplementary_7="supplementary/Figure_5_Supplementary_Data_7_Metabolic_Pathway_(Metabolite)_versus_Metabolic_Metapathway (Gene).xlsx")
dumbbell_overall(before_adj_matrix=before_adj_matrix,S5_df = read.csv("supplementary/Figure_4_Vertical_Dumbbell_Supplementary_Data_5.csv"))
dumbbell_with_sankey_significant(before_adj_matrix=before_adj_matrix,S5_df = read.csv("supplementary/Figure_4_Vertical_Dumbbell_Supplementary_Data_5.csv"),df_table2_metastasis_4 = read.csv("supplementary/Figure_1_Sankey_Supplementary_Data_4.csv"))
PC_plot(before_adj_matrix=before_adj_matrix,S5_df = read.csv("supplementary/Figure_4_Vertical_Dumbbell_Supplementary_Data_5.csv"))
tSNE_plot(before_adj_matrix=before_adj_matrix,S5_df = read.csv("supplementary/Figure_4_Vertical_Dumbbell_Supplementary_Data_5.csv"))
dumbbell_each(df_table2 = read.csv("supplementary/Figure_2and3_Upset_Plot_Supplementary_Data_4.csv"),S5_df = read.csv("supplementary/Figure_4_Vertical_Dumbbell_Supplementary_Data_5.csv"),before_adj_matrix=before_adj_matrix,GI_brain=GI_brain)
Generate_importance(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney)
make_meta_vs_non_meta(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,U_GI=U_GI,U_Lung=U_Lung,U_MSI=U_MSI,Y_GI=Y_GI,Y_LUNG=Y_Lung,Y_MSI=Y_MSI)
raincloud_plot_generate(GI_brain=GI_brain,GI_lung=GI_lung,GI_liver=GI_liver,GI_bone=GI_bone,GI_kidney=GI_kidney,Lung_brain=Lung_brain,Lung_lung=Lung_lung,Lung_liver=Lung_liver,Lung_bone=Lung_bone,Lung_kidney=Lung_kidney,MSI_brain=MSI_brain,MSI_lung=MSI_lung,MSI_liver=MSI_liver,MSI_bone=MSI_bone,MSI_kidney=MSI_kidney,X_GI=X_GI,X_Lung=X_Lung,X_MSI=X_MSI,U_GI=U_GI,U_Lung=U_Lung,U_MSI=U_MSI,Y_GI=Y_GI,Y_LUNG=Y_Lung,Y_MSI=Y_MSI,metabolite_name_trans=Metabolite_pathway$Metabolite_ori)
```
