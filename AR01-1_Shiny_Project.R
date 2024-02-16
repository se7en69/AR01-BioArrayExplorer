
# Settings ----------------------------------------------------------------

library(shiny)
library(tidyverse)
library(genefilter)
library(fgsea)
library(colourpicker)
library(RColorBrewer)
library(data.table)
library(glue)
library(pheatmap)
library(bslib)
library(ggbeeswarm)
library(DT)

# this is to allow the upload of the counts csv file as it is too large for default parameters
options(shiny.maxRequestSize = 100*1024^2)


# Definition of columns for selection within the app ----------------------

meta_columns <- c('Age of Death',
                  'PMI',
                  'RIN',
                  'mRNA-seq Reads',
                  'Age of Onset',
                  'Vonsattel Grade')

de_columns <- c(#'symbol',
                'baseMean',
                'HD.mean',
                'Control.mean',
                'log2FoldChange',
                'lfcSE',
                'stat',
                'pvalue',
                'padj')
gmt <- './data/c2.cp.v7.5.1.symbols.gmt'

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  # using `bslib` package to clean things up for me
  theme = bs_theme(version = 5, bootswatch = 'solar'),
  # title
  titlePanel('AR01-1—; BioArrayExplorer Rehman et al (2023)'),
  
  # disclaimer about app and data needed to run
  tags$h4('Disclaimer'),
  HTML('<p> Welcome to the Microarray Dataset Analysis App! Here what you can expect
<br>
1. Metadata Exploration:
   (I) Summary table providing an overview of your dataset metadata.
   (II) Two visualization options to analyze variables of interest.
<br>
2. Counts Data Analysis:
   (I) Summary tables and diagnostic scatter plots for insights into the counts data.
   (II) Heatmap visualization for gene expression patterns.
   - PCA plot for dimensionality reduction and clustering analysis.
<br>
3. Differential Expression Analysis:
   (I) Filtered tables of DESeq2 results and volcano plot for up- and down-regulated genes.
<br>
4. Gene Set Enrichment Analysis:
  (I)  `fgsea`-based enrichment analysis.
  (II) Bar plot of top enriched pathways.
  (III) Scatterplot of NES vs. -log(padj).
  (IV) Option to download filtered results as CSV.
<br>
(I) User-friendly interface and downloadable results for convenient analysis and sharing.
(II) Simplify your microarray dataset analysis with our app and uncover valuable insights effortlessly. 
(III) Enjoy exploring your data and feel free to contact us with any questions or feedback. 
<br>
Happy analyzing!</p>'),
  HTML('<p>To use this application, please download the Sample data from <a href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810">here</a>.</p>'),
  
  # let's begin building the UI in earnest here
  mainPanel(
    tabsetPanel(
      # Samples tab -------------------------------------------------------------
      
      tabPanel('Samples',
               tabsetPanel(
                 fileInput('file', label = 'Please upload your data', accept = '.csv', placeholder = 'Metadata.csv'),
                 tabPanel('Summary', hr(), tableOutput('summary_table')),
                 tabPanel('Data Table', p('Note: You can search by column at the bottom of the page.'),
                          p("In the search box you can filter for what ever you'd like.", br(),
                            "You could restrict data only to those in the Diagnosis column who are Neurologically normal, or who have Huntington's Disease."),
                          hr(), DTOutput('filtered_data_table'), width = 6),
                 tabPanel('Plots',
                          tabsetPanel(
                            
                            # histogram-----------------------
                            tabPanel('Density Plot', br(),
                                     radioButtons(inputId = 'x_axis',
                                                  label = 'Please select a column from the data that you would like to visualize as a histogram.',
                                                  choices = meta_columns,
                                                  selected = 'Age of Death'), HTML("<p><b>Choose colors you would like below</b></p>"),
                                     colourInput(
                                       inputId = 'fill_color1',
                                       label = 'Fill Color 1', 
                                       value = '#fa0000', 
                                       closeOnClick = T),
                                     colourInput(
                                       inputId = 'fill_color2',
                                       label = 'Fill Color 2', 
                                       value = '#6df900', 
                                       closeOnClick = T),
                                     submitButton(text = 'Change Density Plot', icon('heart')), plotOutput('histogram_plot')),
                            
                            # violin------------------------------
                            tabPanel('Violin Plot', br(),
                                     radioButtons(inputId = 'y_axis',
                                                  label = 'Please select a column from the data that you would like to visualize as a violin plot.',
                                                  choices = meta_columns,
                                                  selected = 'Age of Death'),
                                     colourInput(
                                       inputId = 'fill_color3',
                                       label = 'Fill Color 1', 
                                       value = '#fa0000', 
                                       closeOnClick = T),
                                     colourInput(
                                       inputId = 'fill_color4',
                                       label = 'Fill Color 2', 
                                       value = '#6df900', 
                                       closeOnClick = T), 
                                     submitButton(text = 'Change Violin Plot', icon('heart')), plotOutput('violin_plot')))))),
      
      
      
      
      
      
      
      
      # Counts tab --------------------------------------------------------------
      
      
      tabPanel('Counts',
               ## filter sliders-----------------------
               tabsetPanel(
                 fileInput('file1', label = 'Please upload your counts data', accept = '.csv', placeholder = 'Counts_mat.csv'), 
                 
                 sliderInput('slider_1', 
                             "Slider to include genes with at least X percent of variance.",  
                             min = 0, 
                             max = 100, 
                             value = 10, 
                             step = 1),
                 
                 sliderInput('slider_2', 
                             "Slider to filter by number of samples that are non-zero.", 
                             min = 0,
                             max = 100,
                             value = 50, 
                             step = 1,
                             ticks = FALSE),
                 tags$head(tags$style(HTML('.irs-from, .irs-to, .irs-min, .irs-max, .irs-single {
                                           visibility: hidden !important;}'))),
                 
                 submitButton(text = 'Submit', icon = icon('heart')),
                 
                 
                 ## table outputs----------------------------
                 tabPanel('Filter',
                          hr(HTML('<p><b>Summary table for X percent of variance: </b><p>')), br(),
                          tableOutput('variance_table'), 
                          hr(HTML('<p><b>Summary table for X percent of genes that are non-zero: </b><p>')), br(),
                          tableOutput('filter_zero')),
                 
                 ## diagnostic plots---------------------------
                 tabPanel('Diagnostic Plots', br(),
                          
                          colourInput(
                            inputId = 'c1',
                            label = 'Passing threshold', 
                            value = '#333233', 
                            closeOnClick = T),
                          colourInput(
                            inputId = 'c2',
                            label = 'Not passing threshold', 
                            value = '#960000', 
                            closeOnClick = T),
                          
                          hr(), br(),
                          plotOutput('scat1'),
                          hr(), br(),
                          plotOutput('scat2')),
                  
                 
                 ## clustered heatmap--------------------------
                 tabPanel('Clustered Heatmap', 
                          hr('Heatmap based on percent variance.
                             The legend is the log2(counts).', br(), 
                             '(The higher the percent variance, the longer the heatmap will take to render!)'), 
                          plotOutput('heatmap_vis')), 
                 
                 ## PCA plot-------------------------------
                 tabPanel('PCA Plots (biplot & projection plot)', br(),
                          numericInput(inputId = 'PCx',
                                       label = 'Select which PC to plot on the x-axis',
                                       value = 1,
                                       min = 1,
                                       max = 69),
                          numericInput(inputId = 'PCy',
                                       label = 'Select which PC to plot on the y-axis',
                                       value = 2, 
                                       min = 1,
                                       max = 69),
                          hr('PCA biplot'), 
                          plotOutput('pca_plot'),
                          hr('PCA projection'),
                          plotOutput('bee')))), 
      
      
      
      
      # DE tab ------------------------------------------------------------------
     
      tabPanel('DE',
               sidebarLayout(
                 sidebarPanel(
                   fileInput(
                     inputId = 'diff_ex', 
                     label = 'Please load a differential expression file',
                     accept = 'DESeq2.csv',
                     placeholder = 'DE_results'
                   ),
                   radioButtons(
                     inputId = 'volc_x', 
                     label = 'Select a variable to plot on the x-axis',
                     choices = de_columns, 
                     selected = 'log2FoldChange'
                   ),
                   radioButtons(
                     inputId = 'volc_y', 
                     label = 'Select a variable to plot on the y-axis',
                     choices = de_columns, 
                     selected = 'padj'
                   ),
                   colourInput(
                     inputId = 'de_c1',
                     label = 'Base point color', 
                     value = '#000000', 
                     closeOnClick = T
                   ),
                   colourInput(
                     inputId = 'de_c2',
                     label = 'Highlight point color', 
                     value = '#f900e6', 
                     closeOnClick = T
                   ),
                   sliderInput(
                     'de_slider', 
                     'Select the magnitude of the p-adjusted coloring:',
                     min = -35, max = 0, value = -4
                   ),
                   submitButton(
                     text = 'Submit', 
                     icon = icon('heart')
                   )
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel(
                       'DE Filtered Table',
                       DTOutput('diff_ex_table'), 
                       br()
                     ), 
                     tabPanel(
                       'DE Volcano Plot',
                       hr('Volcano Plot'),
                       plotOutput('diff_ex_volc')))))),
               
       
      
      
      
      # GSEA -----------------------------------------------------------------
      

      
  
  tabPanel(tags$code('FGSEA'), 'Gene Set Enrichment Analysis',
           sidebarLayout(
             sidebarPanel(
               fileInput(
                 inputId = 'DESeq2_results',
                 label = 'Load DESeq2 results',
                 accept = '.csv',
                 placeholder = 'GSE64810_HD_DESeq2.csv'
               ),
               radioButtons(
                 inputId = 'pos_or_neg',
                 label = 'Choose which sign for the NES results 
                 you would like to visualize on the scatter plot, 
                 and download as a CSV (+, -, or All)',
                 choices = c('Negative', 'Positive', 'All'),
                 selected = 'All'
               ),
               sliderInput(
                 'fgsea_slider',
                 'Select the magnitude of the p-adjusted results:',
                 min = -35, max = 1, value = 1
               ),
               
               numericInput(inputId = 'num_paths',
                            label = 'Select number of enriched pathways to visualize',
                            value = 10, 
                            min = 1,
                            max = 50),
               
               submitButton(
                 text = 'Submit',
                 icon = icon('heart')
               )
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel(
                   'FGSEA Bar Plot',
                   plotOutput('fgsea_bar'),
                   br()
                 ),
                 tabPanel(
                   'FGSEA Scatter Plot',
                   plotOutput('fgsea_scatter')),
                 tabPanel(
                   'FGSEA Results table',
                   DTOutput('fgsea_table'),
                   downloadButton("download_csv", "Download as CSV")
                 )))))))
  )

  




# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  # load in the metadata---------------
  load_file <- reactive ({
    req(input$file)
    HD_data <- data.table::fread(input$file$datapath) %>% 
      dplyr::rename('Age of Death' = 'age_of_death',
                    'PMI' = 'pmi',
                    'mRNA-seq Reads' = 'mrna-seq_reads',
                    'Age of Onset' = 'age_of_onset',
                    'Vonsattel Grade' = 'vonsattel_grade')
    
    return(HD_data)  
  })
  
  # Summary table functions----------------------------------
  Summary <- function (HD_data) {
    # get column names
    column_names <- c('Age of Death',
                      'PMI',
                      'RIN',
                      'mRNA-seq Reads',
                      'Age of Onset',
                      'Vonsattel Grade')
    
    
    # get averages
    ## avg age of death
    avg_age_of_death <- mean(HD_data$`Age of Death`) #
    ## avg pmi
    avg_pmi <- mean(as.numeric(HD_data$PMI), na.rm = TRUE)
    ## avg RIN value
    avg_rin <- avg_rin <- mean(HD_data$RIN) #
    ## avg mRNA-seq reads
    avg_mRNA_seq_reads <- mean(HD_data$`mRNA-seq Reads`) #
    ## avg age of onset
    avg_AoO <- mean(HD_data$`Age of Onset`[51:nrow(HD_data)]) #
    ## avg Vonsattel grade
    avg_v <- mean(HD_data$`Vonsattel Grade`[51:nrow(HD_data)])
    
    
    classes <- c("double", "double", "double", "double", 
                 "double", "double")
    
    
    avg_vector <- c(avg_age_of_death, avg_pmi,
                    avg_rin, avg_mRNA_seq_reads, 
                    avg_AoO, avg_v)
    
    summarized <- tibble::tibble(column_names, classes, avg_vector) %>% 
      dplyr::rename('Column Name' = 'column_names',
                    'Class' = 'classes',
                    'Average Value' = 'avg_vector')
    return(summarized)
  }
  
  ## Filtered Data Table-----------------------------
  
  filtered_data_table <- function (HD_data) {
    data_table <- HD_data
    
    return (data_table)
  }
  
  
  ## density function----------------------
  plot_histogram <- function (HD_data, x_axis, color_1, color_2) {
    # defining ggplot object
    hist <- ggplot(data = HD_data, 
                   aes(x = !!sym(x_axis), 
                       fill = Diagnosis)) +
      geom_density(alpha = 0.80) +
      scale_fill_manual(values = c(color_1, color_2))
    
    return (hist)
  }
  
  
  ## violin plot function----------------------
  plot_violin <- function (HD_data, y_axis, color_1, color_2) {
    # defining the ggplot object
    violin <- ggplot(data = HD_data,
                     aes(x = Diagnosis,
                         y = !!sym(y_axis),
                         fill = Diagnosis)) +
      geom_violin(alpha = 0.80) +
      scale_fill_manual(values = c(color_1, color_2))
    
    return (violin)
  }
  
  
  # Counts matrix table functions---------------------------
  load_counts_mat <- reactive ({
    req(input$file1)
    # loading in the counts data
    counts_matrix <- data.table::fread(input$file1$datapath)
    # rownames are not supported w/ tibbles so it will become first column
    # renaming it as 'gene'
    colnames(counts_matrix)[1] <- 'gene'
    return (counts_matrix)
    
  })
  
  
  ## filtering based on variance function------------------
  variance_filtering <- function (counts_matrix, slider) {
    # we want everything minus the first column
    # remember, we set this column as the gene names
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-gene)
    
    # adding a column for variance
    counts_matrix <- counts_matrix %>% 
      dplyr::mutate(variance = rowVars(as.matrix(counts_matrix)), .before = C_0002)
    
    
    # sorting rows based on decreasing order of the variance column
    counts_matrix <- counts_matrix %>% 
      dplyr::arrange(desc(variance))
    
    # making sure the slider value is between 1 and 100
    slider_value <- slider / 100
    # making sure they are whole numbers
    counts_matrix <- counts_matrix %>% 
      dplyr::slice(1:ceiling(nrow(.) * slider_value))
    
    # selecting the 1st and 71st columns
    # storing them in a new data frame called `counts_variance`
    counts_variance <- counts_matrix %>% 
      dplyr::select(1:ncol(.))
    
    # number of columns from counts matrix
    counts_columns <- ncol(counts_matrix)
    
    # number of rows from counts_matrix
    counts_rows <- nrow(counts_matrix)
    
    # number of columns from counts variance data frame
    counts_variance_columns <- ncol(dplyr::select(counts_variance, -variance))
    
    # number of rows from counts variance data frame
    counts_variance_rows <- nrow(counts_variance)
    
    not_pass <- (28087 - counts_variance_rows)
    pass <- (counts_variance_rows / 28087 * 100)
    p_n_pass <- (100 - (counts_variance_rows / 28087 * 100))
    
    # information about the counts data to be returned to user
    # after they make adjustments with the slider
    sample_info <- c('Number of samples', 
                     'Total number of genes',
                     'Number of genes passing threshold',
                     'Number of gene not passing threshold',
                     'Percent of genes passing threshold',
                     'Percent of genes not passing threshold')
    
    variables <- c(counts_columns, 
                   28087, 
                   counts_variance_rows, 
                   not_pass,
                   pass,
                   p_n_pass)
    # putting together in a data frame for the user
    variance_summary <- data_frame(sample_info, variables)
    
    return(variance_summary)  
    
  }
  
  ## filtering for no zeros---------------------------------
  filtering_zero <- function(counts_matrix, slider_2) {
    # remove gene column, same as before
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-gene)
    
    # new column `frequency` with count of non-zero elements in each row
    counts_matrix$frequency <- rowSums(counts_matrix != 0)
    # counts_matrix <- dplyr::mutate(frequency = rowSums(counts_matrix != 0))
    # slider_2 <- (slider_2 / 100) * max(counts_matrix$frequency)
    # selecting rows from the `counts_matrix` where the value of the `frequency` column
    # is greater than the `slider_2` threshold
    # we then create a new data frame called `zeros`
    # this new df contains only the values of the first column of the selected rows from before
    zeros <- counts_matrix %>%
      dplyr::filter(frequency < slider_2) %>%
      dplyr::select(1)
    
    # number of samples
    N_samples <- ncol(counts_matrix)
    # number of genes
    N_genes <- nrow(counts_matrix)
    # which genes pass
    genes_passing <- nrow(zeros)
    # percent of genes passing
    percent_passing <- ((genes_passing / 28087) * 100)
    # number of genes not passing
    genes_not_passing <- (28087 - genes_passing)
    # number of genes not passing
    percent_not_passing <- (100 - percent_passing)
    
    sample_info <- c('Number of samples', 
                     'Total number of genes',
                     'Number of genes passing threshold',
                     'Number of genes not passing threshold',
                     'Percent of genes passing threshold',
                     'Percent of genes not passing threshold')
    
    variables <- c(N_samples,
                   N_genes,
                   genes_passing,
                   genes_not_passing,
                   percent_passing,
                   percent_not_passing)
    
    zeros_summary <- tibble(sample_info, variables)
    
    return(zeros_summary)  
  }
  
  ## diagnostic scatter plots-----------------------------------
  scat_plot_med_vs_var <- function(counts_matrix, slider, color1, color2) {
    # remove gene column, same as before
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-gene)
    
    # adding a column for variance
    counts_matrix <- counts_matrix %>% 
      dplyr::mutate(variance = rowVars(as.matrix(counts_matrix)), .before = C_0002)
    
    # sorting rows based on decreasing order of the variance column
    counts_matrix <- counts_matrix %>% 
      dplyr::arrange(desc(variance))
    
    # making sure the slider value is between 1 and 100
    slider_value <- ((slider / 100) * max(log10(counts_matrix$variance)))
    
    # new column for median
    counts_matrix$median <- apply(counts_matrix, 1, median)
    
    # assigning groups for plotting
    group <- ifelse(counts_matrix$variance < 10^slider_value, TRUE, FALSE)
    
    plot <- ggplot(data = counts_matrix, aes(x = log10(median), y = log10(variance), group = group, color = group))+
      geom_point() +
      scale_color_manual(name = glue('Variance < {slider}%'), values=c(color2, color1)) + 
      theme(legend.position = 'bottom') +
      xlab(label = 'Median Count (log10 scale)') +
      ylab(label = 'log10(Variance)') +
      ggtitle(label = 'Median Count (log10 scale) vs log10(Variance)')
    
    return(plot)
  } 
  
  scat_plot_2 <- function(counts_matrix, slider, color1, color2) {
    # remove gene column, same as before
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-gene)
    # new frequency column
    counts_matrix$frequency <- rowSums(counts_matrix != 0)
    # new median column
    counts_matrix$median <- apply(counts_matrix, 1, median)
    
    # percentage
    slider_value <- ((slider / 100) * max(counts_matrix$frequency))
    # assigning groups for plotting
    group <- ifelse(counts_matrix$frequency < slider_value, TRUE, FALSE)
    
    plot <- ggplot(data = counts_matrix, aes(x = log10(median), y = frequency, group = group, color = group))+
      geom_point() + 
      scale_color_manual(name=glue('Frequency of Non-Zeros < {slider}%'), values=c(color2, color1)) + 
      theme(legend.position = "bottom") +
      xlab(label = 'Median Count (log10 scale)') +
      ylab(label = 'Frequency') +
      ggtitle(label = 'log10(Median Count) vs. Frequency of Non-Zeros')
    
    return(plot)
  }
  
  ## heatmap------------------------------------
  counts_heatmap <- function (counts_matrix, slider_1) {
    
    # remove gene column
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-gene)
    
    # add the variance column
    counts_matrix <- counts_matrix %>% 
      dplyr::mutate(variance = rowVars(counts_matrix), .before = 'C_0002')
    
    counts_matrix <- counts_matrix %>% dplyr::mutate_all(function(x) log2(x+1))
    
    
    
    slider_value_1 <- ((slider_1 / 100) * max(counts_matrix$variance))
    
    # Filter rows where variance is less than or equal to slider_value_1
    counts_matrix <- counts_matrix[counts_matrix$variance <= slider_value_1, ]
    
    # Remove the variance column
    
    counts_matrix <- counts_matrix[, -1]
    
    
    
    
    
    colors <- colorRampPalette(rev(brewer.pal(n = 11, name = 'PiYG')))(50)
    # colors <- rev(colors)
    
    hm <- pheatmap(counts_matrix,
                   scale = 'row',
                   color = colors,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   # clustering_distance_rows = 'manhattan',
                   # clustering_distance_cols = 'correlation',
                   treeheight_row = 0,
                   show_rownames = FALSE)

    return (hm)
    
  }
  
  ## PCA--------------------------
  pca <- function (counts_matrix, X, Y) {

    # remove gene column
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-gene) %>% 
      dplyr::mutate_all(function(x) log2(x + 1))

    # Create metadata with the "Neuro" column
    metadata <- tibble('sample' = colnames(counts_matrix), 
                       'Neuro' = if_else(str_starts(sample, "C_"), "Neurologically Normal", "Huntington's Disease"))
    
    # Performing PCA
    pca_results <- prcomp(scale(t(counts_matrix)), center = FALSE, scale. = FALSE)
    
    # Extract the user desired principal components
    PCx <- pca_results$x[, X]
    PCy <- pca_results$x[, Y]
    
    # Combine the principal component scores with the metadata
    plot_df <- cbind(data.frame(PCx = PCx, PCy = PCy))
    plot_df <- rownames_to_column(plot_df, var = 'sample')
    
    # match IDs
    plot_df <- dplyr::inner_join(plot_df, metadata, by = 'sample')
    
    # Create the plot
    p <- ggplot(plot_df, aes(x = PCx, y = PCy, color = Neuro)) +
      geom_point() +
      scale_color_manual(values = c("#fa0000", "#240032")) +
      stat_ellipse(aes(x = PCx, y = PCy, color = Neuro)) +
      scale_x_continuous(name = paste0("PC", X, " (", round(100*pca_results$sdev[X]^2/sum(pca_results$sdev^2), 1), "%)")) +
      scale_y_continuous(name = paste0("PC", Y, " (", round(100*pca_results$sdev[Y]^2/sum(pca_results$sdev^2), 1), "%)")) +
      labs(title = paste0("PC", X, " vs. PC", Y, "; Grouped by Neurological Function"), color = 'Neurological Function') +
      geom_vline(xintercept = 0, linetype = 3) +
      geom_hline(yintercept = 0, linetype = 3)
    
    return(p)
  }
  
  ## PCA bee------------------------
  pca_bee <- function (counts_matrix) {
    
    # remove gene column
    counts_matrix <- counts_matrix %>% 
      dplyr::select(-1) %>% 
      dplyr::mutate_all(function(x) log2(x + 1))
    
    
    # Performing PCA
    pca_results <- prcomp(scale(t(counts_matrix)), center = TRUE, scale. = TRUE)
    
    plot <- pca_results$x %>% as_tibble() %>% 
      pivot_longer(everything(), names_to = 'PC', values_to = 'Projection') %>% 
      # slice(1:1000) %>% 
      mutate(PC = fct_relevel(PC, str_c('PC', 1:69))) %>% 
      ggplot(aes(x = PC, y = Projection)) +
      geom_beeswarm() +
      labs(title = 'PCA Projection Plot') +
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
    
    return (plot)
    
  }
  
  
  # Differential expression functions------------------------
    DE_volcano <- function (de_file, x_name, y_name, slider, color1, color2) {
      # dropping NAs
      de_file <- de_file %>% drop_na()
      #   filter(padj != 0)
      
      # assigning groups for plotting
      group <- ifelse(de_file[[y_name]]<1*(10^slider), TRUE, FALSE)
      # ggplot object
      vp <- ggplot(de_file, aes(x=de_file[[x_name]], y=-log10(de_file[[y_name]]), group=group, color=group)) + 
        # setting point geometry
        geom_point() + 
        # setting colors for volcano plot
        scale_color_manual(name=glue('padj < 1e{slider}'), values=c(color1, color2)) + 
        # labelling axes
        labs(x = x_name, y = y_name) +
        theme(plot.title = element_blank(), legend.position = 'bottom')
        # return the volcano plot
        return(vp)

    }
  
  ## DE table----------------------------------------
    de_table <- function (de_file, slider) {
      table <- de_file %>% 
        dplyr::filter(padj < 1*(10^slider))
    
      return (table)
    }
  
  
  
  # FGSEA -----------------------------------------------
 
  run_gsea <- function(labeled_results, gmt, slider, sign) {
    # I think we can add the slider after `arrange`
    GENES <- labeled_results %>%
      arrange(desc(log2FoldChange))
    
    # pathways for `fgsea`
    pathways_fgsea <- fgsea::gmtPathways(gmt)
    
    gsea_list <- setNames(GENES$log2FoldChange, GENES$symbol)
    
    fgsea_results <- fgsea(pathways_fgsea, na.omit(gsea_list), minSize = 15, maxSize = 500)
    
    fgsea_results <- fgsea_results %>% as_tibble() %>%
      arrange(padj) %>%
      dplyr::filter(padj <= 1*(10^slider))
    # making the results more reader-friendly for plotting
    fgsea_results$pathway <- base::gsub('_', ' ', fgsea_results$pathway)
    
    return (fgsea_results)
    
  }
  
  
  ## plot fgsea results------------------------------
  plot_fgsea <- function (fgsea_results, num_paths) {
    fgseaResTidy <- fgsea_results %>% dplyr::select(-leadingEdge, -ES, -size)
    top_pos <- fgseaResTidy %>% top_n(n = num_paths, wt = NES) %>% dplyr::mutate(sign = "pos")
    top_neg <- fgseaResTidy %>% top_n(n = -num_paths, wt = NES) %>% dplyr::mutate(sign = "neg")
    top_nes <- rbind(top_pos, top_neg)
    # Format the pathway columns to make axis labels look nice
    hall_NES <- top_nes %>%
      dplyr::arrange(NES) %>%
      # Format the pathway columns to make axis labels look nice
      dplyr::mutate(pathway = gsub("_"," ", pathway), 
                    pathway = stringr::str_wrap(pathway, width = 80))
    
    # Make the bar plot of the top pathways
    plot <- hall_NES %>% ggplot2::ggplot() + geom_col(aes(NES, pathway, fill=sign)) + 
      # Set the columns to the correct order with highest scores at the top and 
      # lowest scores at the bottom
      scale_y_discrete(limits = pull(hall_NES, pathway)) +
      # Set the colors for the respective NES groups
      scale_fill_manual(values=c("#bf1616","#2d43fa")) +
      # Add labels and formatting 
      ggtitle("fgsea results for Hallmark MSigDB gene sets") +
      xlab("Normalized Enrichment Score (NES)") + 
      theme(axis.text.y=element_text(size = 6), legend.position="none",
            axis.title.y=element_blank(), axis.title.x=element_text(size = 8),
            title=element_text(size = 8.5), axis.ticks = element_blank())
    
    return(plot)
    
  }
  
  ## fgsea scatter
  fgsea_scatter <- function(fgsea_results, sign) {
    
    res <- fgsea_results %>% 
      mutate(sign = ifelse(NES > 0, 'Pos', 'Neg')) %>% 
      as_tibble()
    
    if(sign == 'Positive') {
      p_results <- filter(res, sign == 'Pos')
      p <- ggplot(data = p_results) +
        geom_point(mapping = aes(x = NES, y = -log10(padj), color = 'Pos')) +
        scale_color_manual(values = "#2d43fa", labels = "Pos") +
        ggtitle("Positives") +
        labs(color = 'Pos')
      pp <- plot(p)
      return(pp)
    }
    else if(sign == 'Negative') {
      n_results <- filter(res, sign == 'Neg')
      n <- ggplot(data = n_results) +
        geom_point(mapping = aes(x = NES, y = -log10(padj), color = 'Neg')) +
        scale_color_manual(values = "#bf1616", labels = "Neg") +
        ggtitle("Negatives") +
        labs(color = 'Neg')
      nn <- plot(n)
      return(nn)
    }
    else if(sign == 'All') {
      a <- ggplot(data = res) +
        geom_point(mapping = aes(x = NES, y = -log10(padj), color = sign)) +
        scale_color_manual(values = c("#bf1616", "#2d43fa", "#9700cf"),
                           labels = c("Neg", "Pos", "All")) +
        labs(color = "Sign") +
        ggtitle("All")
      aa <- plot(a)
      return(aa)
    }
    
    
  }
  
  
  
  # Outputs -----------------------------------------------------------------
  ## Summary-----------------------------
  ### summary table output--------------------------
  output$summary_table <- renderTable({
    # req(input$file)
    table <- load_file()
    return(Summary(table))
    
  })
  
  
  ### filtered table output---------------------------
  output$filtered_data_table <- renderDT({
    # req(input$file)
    f_table <- load_file()
    return(filtered_data_table(f_table))
  }, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
  
  
  
  
  isolate({
    observeEvent(load_file(), {
      updateRadioButtons(session, 'x_axis', choices = meta_columns, selected = "Age of Death")
    })
  })
  
  
  
  ### histogram output-------------------------
  output$histogram_plot <- renderPlot({
    HD_data <- load_file()
    
    plot_histogram(HD_data, input$x_axis, input$fill_color1, input$fill_color2)
    
  })
  observeEvent(input$histogram_submit, {
    
  })
  
  
  
  ### violin plot output-------------------------
  output$violin_plot <- renderPlot({
    HD_data <- load_file()
    
    plot_violin(HD_data, input$y_axis, input$fill_color3, input$fill_color4)
  })
  observeEvent(input$violin_submit, {
    
  })
  
  ## Counts tab-----------------------------------
  ### variance filter table output-------------------
  output$variance_table <- renderTable({
    c_mat <- load_counts_mat()
    
    res <- variance_filtering(c_mat, input$slider_1) %>% 
      dplyr::rename('Sample Information' = 'sample_info',
                    'Results' = 'variables')
    return(res)
    
  })
  
  
  ### zero filter table output----------------------
  output$filter_zero <- renderTable({
    c_mat <- load_counts_mat()
    
    res <- filtering_zero(c_mat, input$slider_2) %>% 
      dplyr::rename('Sample Information' = 'sample_info',
                    'Results' = 'variables')
    
    return(res)
    
  })
  
  ### scatterplot outputs---------------------------
  output$scat1 <- renderPlot({
    
    c_mat <- load_counts_mat()
    res <- scat_plot_med_vs_var(c_mat, input$slider_1, input$c1, input$c2)
    
    return(res)
    
  })
  
  output$scat2 <- renderPlot({
    
    c_mat <- load_counts_mat()
    res <- scat_plot_2(c_mat, input$slider_2, input$c1, input$c2)
    
    return(res)
    
  })
  
  
  ### heatmap--------------------------------------
  output$heatmap_vis <- renderPlot({
    
    c_mat <- load_counts_mat()
    HEAT <- counts_heatmap(c_mat, input$slider_1)
    
    return (HEAT)
  })
  
  
  ### pca--------------------------------------------
  output$pca_plot <- renderPlot({
    
    c_mat <- load_counts_mat()
    PCA <- pca(c_mat, input$PCx, input$PCy)
    
    return (PCA)
  })
  
  output$bee <- renderPlot ({
    
    c_mat <- load_counts_mat()
    bee <- pca_bee(c_mat)
    
    return (bee)
  })
  
  ## Differential Expression---------------
  
  ### volcano plot-------------------------
  
  # load in the DESeq2 data
  load_diff_ex <- reactive ({
    req(input$diff_ex)
    # loading in the counts data
    deseq2 <- utils::read.csv(input$diff_ex$datapath, sep = '')
    
    return (deseq2)
    
  })
  
  # DE table here
  output$diff_ex_table <- renderDT({
    
    DE_table <- load_diff_ex()
    
    res <- de_table(DE_table, input$de_slider)
    
    return (res)

  }, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
  
  
  # volcano here
  output$diff_ex_volc <- renderPlot({
    
    deseq2_res <- load_diff_ex()
    vp <- DE_volcano(deseq2_res, input$volc_x, input$volc_y, input$de_slider, input$de_c1, input$de_c2)
    
    return (vp)
    
  })
  
  ## GSEA--------------------------------------
  ### load in DESeq2 data--------------------
  DS2_data <- reactive({
    req(input$DESeq2_results)
    
    deseq2 <- utils::read.csv(input$DESeq2_results$datapath, sep = '')
    
    return (deseq2)
  })
  
  ### GSEA table-------------------------
  output$fgsea_table <- renderDT({
    
    table <- run_gsea(DS2_data(), gmt, input$fgsea_slider, input$pos_or_neg) %>% 
      dplyr::select(-leadingEdge) %>% 
      mutate(sign = ifelse(NES > 0, 'Pos', 'Neg')) %>% 
      as_tibble()
    
    if (input$pos_or_neg == 'Positive') {
      table <- table %>% 
        dplyr::filter(sign == 'Pos')
      return(table)
    }
      
    else if (input$pos_or_neg == 'Negative') {
      table <- table %>% 
        dplyr::filter(sign == 'Neg')
      
      return(table)
    }
    
    
    
    
    return (table)
    
  }, options = list(scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
  
  ### GSEA bar-----------------
  output$fgsea_bar <- renderPlot({
    
    table <- run_gsea(DS2_data(), gmt, input$fgsea_slider)
    bp <- plot_fgsea(table, input$num_paths)
    return (bp)
  })
  
  
  ### GSEA scatter-------------------------
  output$fgsea_scatter <- renderPlot({
    
    table <- run_gsea(DS2_data(), gmt, input$fgsea_slider)
    
    sp <- fgsea_scatter(table, input$pos_or_neg) 
    
    return(sp)
  })
  
  table <- reactive({
    run_gsea(DS2_data(), gmt, input$fgsea_slider) %>% 
      dplyr::select(-leadingEdge) %>% 
      mutate(sign = ifelse(NES > 0, 'Pos', 'Neg')) %>% 
      as_tibble()
  })
  
  
  ### download table results as CSV
  output$download_csv <- downloadHandler(
    filename = function() {
      paste("fgsea_table_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(table(), file, row.names = FALSE)
    }
  )
  
  
  
}
# Run the application 
shinyApp(ui = ui, server = server)


