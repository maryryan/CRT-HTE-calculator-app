#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(plotly)
library(DT)
library(patchwork)
library(latex2exp)
library(RColorBrewer)
library(shinybrowser)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # detect window size #
  tags$head(tags$script('
                                  var dimension = [0, 0];
                                  $(document).on("shiny:connected", function(e) {
                                      dimension[0] = window.innerWidth;
                                      dimension[1] = window.innerHeight;
                                      Shiny.onInputChange("dimension", dimension);
                                  });
                                  $(window).resize(function(e) {
                                      dimension[0] = window.innerWidth;
                                      dimension[1] = window.innerHeight;
                                      Shiny.onInputChange("dimension", dimension);
                                  });
                              ')),
  # Application title #
  #titlePanel("Heterogeneous CRT Calculator"),
  tags$title("CRT HTE Calculator"),
  div(h2("CRT HTE Calculator", style="margin:0;"),
      h4("Power and sample size for effect modification in CRTs", style="margin:0;")),
  #),
  
  # Input sidebar
  sidebarLayout(
    #### type of trial ####
    sidebarPanel(#style="overflow-y:scroll; max-height: 800px; position:relative;",
      radioButtons(inputId="trial",
                   label="Design Type",
                   choiceNames=c("Parallel two-level CRT", "Parallel three-level CRT",
                                 "Multiple-period parallel CRT",
                                 "Two-period cluster cross-over",
                                 "Multiple-period cluster cross-over",
                                 "Stepped-wedge",
                                 "Individually randomized group treatment",
                                 "Parallel two-level CRT (sample size and ICC heterogeneous by arm)",
                                 "Upload custom design"),
                   choiceValues = c("parallel", "three_level","parallel_m","crossover_2","crossover_m","SWD", "irgt", "het_two","upload")),#can also add help alt/hover text later using helpText()
      
      ## helptext for designs ##
      conditionalPanel(
        condition = "input.trial == 'parallel'",
        helpText("A parallel two-level CRT randomizes clusters of individuals (i.e., students within schools) to either receive the control or intervention condition.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'three_level'",
        helpText("A parallel three-level CRT randomizes clusters of individuals, who are nested within subclusters (i.e., patients within physicians within clinics), to either receive the control or intervention condition.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'parallel_m'",
        helpText("A multiple-period parallel CRT randomizes clusters of individuals to either receive the control or intervention condition over J periods.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'crossover_2'",
        helpText("A two-period cross-over CRT randomizes clusters of individuals to one of two treatment sequences, either receiving the control first followed by the intervention, or receiving the intervention condition first followed by the control.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'crossover_m'",
        helpText("A multiple-period cross-over CRT randomizes clusters of individuals to one of two treatment sequences alternated over J periods: either the receipt of the control condition followed by the intervention, or the receipt of the intervention condition followed by the control.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'SWD'",
        helpText("A stepped-wedge CRT randomizes clusters of individuals to the time at which the cluster initiates the intervention. All clusters begin on the control condition and, as time progresses, batches of clusters are transitioned to the intervention.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'het_two'",
        helpText("A parallel two-level CRT randomizes clusters of individuals to either receive the control or intervention condition, with allowance for sample size (number of clusters, cluster size) and ICCs to vary by study arm.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'irgt'",
        helpText("A parallelt two-level CRT randomizes individuals to either receive treatment or intervention condition, but treatment conditions are administered in a clustered or grouped way (e.g., group instruction).", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      conditionalPanel(
        condition = "input.trial == 'upload'",
        helpText("Designs not fitting the description of any of the above choices can by accommodated by uploading the design as a CSV file.")
      ),
      
      ## upload own design ##
      conditionalPanel(
        condition = "input.trial == 'upload'",
        fileInput("file1", "Upload a design matrix:", accept=c('text/plain', '.csv')),
        helpText("The file must be a comma separated .csv file consisting of 0s (control) and 1s (treatment), with a column for each time period and a row for each treatment sequence. There may not be any mising cluster-periods. Do not include row or column names.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      
      ## clustered control for IRGT ##
      # when treatment arm only, this will trigger m0=1 #
      conditionalPanel(
        condition = "input.trial == 'irgt'",
        radioButtons(inputId="control_cluster",
                     label="Clustering prevalence:",
                     choiceNames = c("Both arms",
                                     "Treatment arm only"),
                     choiceValues = c("cluster","indiv")),
        helpText("In individually randomized group treatment trials, both the treatment and control arms may be clustered, or the control arm could have no clustering whatsoever if the standard of care is not administered in groups or by a common provider.", style="margin-top:-0.5em; margin-bottom:1em;")),
      
      ## cross-sectional vs closed cohort ##
      conditionalPanel(
        condition = "input.trial == 'parallel_m' || input.trial == 'crossover_2' || input.trial == 'crossover_m' || input.trial == 'SWD' || input.trial == 'upload'",
        radioButtons(inputId="cohort",
                     label="Sampling scheme:",
                     choiceNames = c("Cross-sectional",
                                     "Closed-cohort"),
                     choiceValues = c("cross","closed")),
        helpText("A cross-sectional sampling scheme assumes individuals within each cluster only contribute data to a single time period. A closed-cohort sampling scheme assumes indivdiuals with each cluster are followed longitudinally across all periods.", style="margin-top:-0.5em; margin-bottom:1em;")
      ),
      ## Plot display options ##
      radioButtons(inputId="plot_display",
                   label="Plot display",
                   choiceNames = c("Cluster size vs Power",
                                   "Number of clusters vs Power",
                                   "Cluster size vs Number of clusters"),#,
                   #"Number of clusters vs Cluster size (fixed power)"),
                   choiceValues=c("m_v_power", "n_v_power", "fixed_power")),# "n_v_m"),
      helpText("Hover over the plot lines to obtain precise design parameter information", style="margin-top:-0.5em; margin-bottom:1em;"),
      
      #### Parallel designs ####
      conditionalPanel(
        condition = "input.trial == 'parallel'",
        
        ## sample size ##
        # numeric n #
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power'",
          numericInput(inputId="n",
                       label="Total number of clusters (n)",
                       10, min=1, max=9999999999)
        ),
        # m range #
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power' || input.plot_display == 'fixed_power'",
          sliderInput(inputId="m_range",
                      label="Plot cluster size range",
                      min=1, max=3000,
                      value=c(5,100))
        ),
        #numeric m, n range #
        conditionalPanel(
          condition = "input.plot_display == 'n_v_power'",
          numericInput(inputId="m",
                       label="Cluster size (m)",
                       20, min=1, max=9999999999),
          sliderInput(inputId="n_range",
                      label="Plot number of clusters range",
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # numeric power #
        conditionalPanel(
          condition = "input.plot_display == 'fixed_power'",
          numericInput(inputId="power",
                       label="Power",
                       0.9, min=0, max=1)
        ),
        ## ICCs ##
        tags$u(h3("ICC options")),
        numericInput(inputId="oicc_est",
                     label="Assumed outcome ICC",
                     0.1,min=0, max=1, step=0.001),
        numericInput(inputId="cicc_est",
                     label="Assumed covariate ICC",
                     0.1,min=0, max=1, step=0.001),
        
        radioButtons(inputId="icc_sensitivity",
                     label="ICC sensitivity analyses",
                     choiceNames = c("Only display results for assumed ICCs",
                                     "Display results for ICC ranges"),
                     choiceValues = c("est_only", "sensitivity"),
                     inline=T),
        conditionalPanel(
          condition = "input.icc_sensitivity == 'sensitivity'",
          
          radioButtons(inputId="icc_display",
                       label="ICC display",
                       choiceNames = c("Outcome ICC constant within-plot", "Covariate ICC constant within-plot"),
                       choiceValues=c("oICC_constant", "cICC_constant"),
                       inline=T),
          
          ## play with font bolding/italics/color for section titles to make sections pop ##
          # h4("Outcome ICC"),
          sliderInput(inputId="oicc_range",
                      label="Outcome ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          # numericInput(inputId="oicc_min",
          #              label="Minimum",
          #              0.05, min=0, max=1, step=0.001),
          # numericInput(inputId="oicc_max",
          #              label="Maximum ",
          #              0.2,min=0, max=1, step=0.001),
          # h4("Covariate ICC"),
          # numericInput(inputId="cicc_min",
          #              label="Minimum",
          #              0.05, min=0, max=1, step=0.001),
          # numericInput(inputId="cicc_max",
          #              label="Maximum ",
          #              0.2,min=0, max=1, step=0.001)
          sliderInput(inputId="cicc_range",
                      label="Covariate ICC range",
                      min=0, max=1,
                      value=c(0,0.9))
        )#end ICC sensitivity conditional
      ),# end of parallel
      
      # tags$u(h3("Outcome and variable options")),
      # ## outcomes ##
      # radioButtons(inputId="outcome",
      #              label="Outcome type",
      #              choiceNames = c("Continuous", "Binary"),
      #              choiceValues=c("continuous", "binary"),
      #              inline=T),
      # # effect size and sd for continuous #
      # conditionalPanel(
      #   condition = "input.outcome == 'continuous'",
      #   # numericInput(inputId="mean_diff_ATE",
      #   #              label="Assumed ATE",
      #   #              1, min=0, max=999999),
      #   numericInput(inputId="sd_outcome",
      #                label="Outcome standard deviation",
      #                1, min=0, max=999999)
      # ),#end conditional
      # # effect size for binary #
      # conditionalPanel(
      #   condition = "input.outcome == 'binary'",
      #   numericInput(inputId="prop_control",
      #                # outcome proportion/event prop of those on control #
      #                label=withMathJax("Mean outcome under control (\\(p_0\\))"),
      #                0.5, min=0, max=1,step=0.001),
      #   numericInput(inputId="prop_trt",
      #                # outcome proportion/event prop of those on intervention #
      #                label=withMathJax("Mean outcome under intervention (\\(p_1\\))"),
      #                0.5, min=0, max=1, step=0.001)
      # ),#end conditional
      # 
      # # covariate #
      # radioButtons(inputId="covar",
      #              label="Covariate type",
      #              choiceNames = c("Continuous", "Binary"),
      #              choiceValues=c("continuous", "binary"),
      #              inline=T),
      # # can keep this as specifying HTE effect size since HTE for binary outcome doesn't translate directly into diff of proportions
      # numericInput(inputId="mean_diff_HTE",
      #              label="Assumed HTE",
      #              1, min=0, max=999999),
      # # covariate sd for continuous #
      # conditionalPanel(
      #   condition = "input.covar == 'continuous'",
      #   numericInput(inputId="sd_covar",
      #                label="Covariate standard deviation",
      #                1, min=0, max=999999)
      # ),#end conditional
      # # covariate sd for binary #
      # conditionalPanel(
      #   condition = "input.covar == 'binary'",
      #   numericInput(inputId="prop_covar",
      #                label="Covariate proportion",
      #                0.5, min=0, max=1,step=0.001)
      # )#end conditional
      
      # tags$u(h3("Other design options")),
      # #treatment allocation #
      # numericInput(inputId="w",
      #              label="Intervention allocation",
      #              0.5, min=0.000001, max=1, step=0.001),
      # # sig level #
      # numericInput(inputId="sig",
      #              label="Significance level",
      #              0.05, min=0.000001, max=1, step=0.001)
      #),
      #END PARALLEL DESIGN
      # normal or t #
      # radioButtons(inputId="approx",
      #              label="Normal approximation",
      #              choiceNames = c("Normal approximation", "T-distribution"),
      #              choiceNames = c("normal", "t"))
      
      #### THREE-LEVEL DESIGN ####
      conditionalPanel(
        condition = "input.trial == 'three_level'",
        # radioButtons(inputId="plot_display_three",
        #              label="Plot display",
        #              choiceNames = c("Individuals vs Power (fixed clusters and subclusters)",
        #                              #"Subclusters vs Power (fixed clusters and individuals)",
        #                              "Clusters vs Power (fixed subclusters and individuals)",
        #                              #"Subclusters vs Clusters (fixed power and individuals)",
        #                              "Individuals vs Clusters (fixed power and subclusters)"#,
        #                              #"Individuals vs Subclusters (fixed power and clusters)"
        #                              ),#,
        #              #"Number of clusters vs Cluster size (fixed power)"),
        #              choiceValues=c("m_v_power", #"ns_v_power",
        #                             "nc_v_power", #"ns_v_nc",
        #                             "m_v_nc"#,"m_v_ns"
        #                             )),# "n_v_m"),
        ## sample size ##
        # ns numeric #
        numericInput(inputId="ns",
                     label=withMathJax("Number of subclusters per cluster (\\(n_s\\))"),
                     10, min=1, max=9999999999),
        #nc numeric #
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power'",
          numericInput(inputId="nc",
                       label=withMathJax("Total number of clusters (\\(n_c\\))"),
                       10, min=1, max=9999999999)
        ),
        
        
        # m numeric and nc range #
        conditionalPanel(
          condition = "input.plot_display == 'nc_v_power'",
          numericInput(inputId="m_three",
                       label="Individuals per subcluster (m)",
                       20, min=1, max=9999999999),
          sliderInput(inputId="nc_range",
                      label="Plot number of clusters range",
                      min=1, max=3000,
                      value=c(5,100))
        ),
        
        # m range #
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power' || input.plot_display == 'm_v_nc'",
          sliderInput(inputId="m_three_range",
                      label="Plot individuals (per subcluster) range",
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # fixing power #
        # conditionalPanel(
        #   condition = "input.plot_display_three == 'ns_v_nc'",
        #   numericInput(inputId="power_three",
        #                label="Power",
        #                0.9, min=0, max=1),
        #   numericInput(inputId="m3_three",
        #                label="Individuals (m)",
        #                100, min=1, max=9999999999),
        #   sliderInput(inputId="ns_range_p",
        #               label="Plot subcluster size range",
        #               min=1, max=3000,
        #               value=c(5,100))
        # ),
        conditionalPanel(
          condition = "input.plot_display =='m_v_nc'",
          numericInput(inputId="power_three",
                       label="Power",
                       0.9, min=0, max=1)
        ),
        # conditionalPanel(
        #   condition = "input.plot_display_three == 'm_v_ns'",
        #   numericInput(inputId="power3_three",
        #                label="Power",
        #                0.9, min=0, max=1),
        #   numericInput(inputId="nc3",
        #                label="Number of clusters (nc)",
        #                100, min=1, max=9999999999),
        #   sliderInput(inputId="m2_range_p",
        #               label="Plot individuals range",
        #               min=1, max=3000,
        #               value=c(5,100))
        # ),
        radioButtons(inputId="randomization_three",
                     label="Randomization level",
                     choiceNames=c("Randomized at cluster level", "Randomized at subcluster level"),
                     choiceValues = c("cluster", "subcluster"),
                     inline=T),
        tags$u(h3("ICC options")),
        tags$i(h4("Outcome")),
        numericInput(inputId="oicc_wsub_est_three",
                     label=withMathJax("Assumed within-subcluster ICC (\\(\\alpha_0\\))"),
                     0.1,min=0, max=1, step=0.001),
        numericInput(inputId="oicc_ratio_est_three",
                     label=withMathJax("Assumed CAC (\\(\\alpha_1/\\alpha_0\\))"),
                     0.5,min=0, max=1, step=0.001),
        helpText("The ratio of the between-subcluster outcome ICC over the within-subcluster outcome ICC.", style="margin-top:-0.5em; margin-bottom:1em;"),
        
        tags$i(h4("Covariate")),
        numericInput(inputId="cicc_wsub_est_three",
                     label=withMathJax("Assumed within-subcluster ICC (\\(\\rho_0\\))"),
                     0.1,min=0, max=1, step=0.001),
        numericInput(inputId="cicc_ratio_est_three",
                     label=withMathJax("Assumed CAC (\\(\\rho_1/\\rho_0\\))"),
                     0.5,min=0, max=1, step=0.001),
        helpText("The ratio of the between-subcluster covariate ICC over the within-subcluster covariate ICC.", style="margin-top:-0.5em; margin-bottom:1em;"),
        
        radioButtons(inputId="icc_sensitivity_three",
                     label="ICC sensitivity analyses",
                     choiceNames = c("Only display results for assumed ICCs",
                                     "Display results for ICC ranges"),
                     choiceValues = c("est_only", "sensitivity"),
                     inline=T),
        conditionalPanel(
          condition = "input.icc_sensitivity_three == 'sensitivity'",
          
          radioButtons(inputId="icc_display_three",
                       label="ICC display",
                       choiceNames = c("Outcome ICCs constant within-plot", "Covariate ICCs constant within-plot"),
                       choiceValues=c("oICC_constant", "cICC_constant"),
                       inline=T),
          ## play with font bolding/italics/color for section titles to make sections pop ##
          # h4("Outcome ICC"),
          tags$i(h4("Outcome")),
          sliderInput(inputId="oicc_wsub_range_three",
                      label="Within-subcluster ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          sliderInput(inputId="oicc_ratio_range_three",
                      label="CAC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          
          
          tags$i(h4("Covariate")),
          sliderInput(inputId="cicc_wsub_range_three",
                      label="Within-subcluster ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          sliderInput(inputId="cicc_ratio_range_three",
                      label="CAC range",
                      min=0, max=1,
                      value=c(0,0.9))
        )#end ICC sensitivity conditional
        
        
        
      ),# end three-level conditional
      
      # tags$u(h3("Other design options")),
      
      
      #### IRGT ####
      conditionalPanel(
        condition = "input.trial == 'irgt'",
        
        ## sample size ##
        # treatment clusters #
        tags$u(h3("Sample Size/Power options")),
        tags$i(h4("Treatment arm")),
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power' || input.plot_display == 'n_v_power' || (input.plot_display == 'fixed_power' & input.control_cluster == 'cluster')",
          numericInput(inputId="n1_fix",
                       label=withMathJax("Assumed number of clusters in treatment arm (\\(n_1\\))"),
                       10, min=1, max=9999999999),
          
        ),
        conditionalPanel(
          condition = "input.plot_display == 'n_v_power'",
          sliderInput(inputId="n1_slide",
                      label=withMathJax("Number of clusters range in treatment arm (\\(n_1\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # treatment cluster size #
        conditionalPanel(
          condition = "(input.plot_display == 'm_v_power' & input.control_cluster == 'cluster') || input.plot_display == 'n_v_power' || (input.plot_display == 'fixed_power' && input.control_cluster == 'cluster')",
          numericInput(inputId="m1_fix",
                       label=withMathJax("Assumed cluster size in treatment arm (\\(m_1\\))"),
                       10, min=1, max=9999999999),
          
        ),
        conditionalPanel(
          condition = "input.plot_display != 'n_v_power'",
          sliderInput(inputId="m1_slide",
                      label=withMathJax("Cluster size range in treatment arm (\\(m_1\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # control clusters #
        tags$i(h4("Control arm")),
        numericInput(inputId="n0_fix",
                     label=withMathJax("Assumed number of clusters in control arm (\\(n_0\\))"),
                     10, min=1, max=9999999999),
        conditionalPanel(
          condition = "input.plot_display == 'n_v_power'",
          sliderInput(inputId="n0_slide",
                      label=withMathJax("Number of clusters range in control arm (\\(n_0\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # control cluster size #
        # only show control cluster size when clustering in control arm #
        conditionalPanel(
          condition = "input.control_cluster == 'cluster'",
          numericInput(inputId="m0_fix",
                       label=withMathJax("Assumed cluster size in control arm (\\(m_0\\))"),
                       10, min=1, max=9999999999),
          conditionalPanel(
            condition = "input.plot_display != 'n_v_power'",
            sliderInput(inputId="m0_slide",
                        label=withMathJax("Cluster size range in control arm (\\(m_0\\))"),
                        min=1, max=3000,
                        value=c(5,100))
          )
        ),
        conditionalPanel(
          condition = "input.control_cluster == 'indiv'",
          helpText("When clustering only occurs in the treatment arm, the cluster size of the control arm is set to 1.")
        ),
        conditionalPanel(
          condition = "input.plot_display == 'fixed_power'",
          numericInput(inputId="power_irgt",
                       label="Power",
                       0.9, min=0, max=1)
        ),
        
        
        
        
        
        tags$u(h3("ICC options")),
        #tags$i(h4("Outcome")),
        numericInput(inputId="oicc_trt_est_irgt",
                     label=withMathJax("Assumed treatment-arm outcome ICC (\\(\\rho_{y|x,1}\\))"),
                     0.5,min=0, max=1, step=0.001),
        numericInput(inputId="oicc_ctrl_est_irgt",
                     label=withMathJax("Assumed control-arm outcome ICC (\\(\\rho_{y|x,0}\\))"),
                     0.1,min=0, max=1, step=0.001),
        
        
        radioButtons(inputId="icc_sensitivity_irgt",
                     label="ICC sensitivity analyses",
                     choiceNames = c("Only display results for assumed ICCs",
                                     "Display results for ICC ranges"),
                     choiceValues = c("est_only", "sensitivity"),
                     inline=T),
        conditionalPanel(
          condition = "input.icc_sensitivity_irgt == 'sensitivity'",
          
          radioButtons(inputId="icc_display_irgt",
                       label="ICC display",
                       choiceNames = c("Treatment-arm ICC constant within-plot", "Control-arm ICC constant within-plot"),
                       choiceValues=c("oICC1_constant", "oICC0_constant")),
          
          sliderInput(inputId="oicc_trt_range_irgt",
                      label="Treatment-arm outcome ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          sliderInput(inputId="oicc_ctrl_range_irgt",
                      label="Control-arm outcome ICC range",
                      min=0, max=1,
                      value=c(0,0.9))
          
        )#end ICC sensitivity conditional
        
        
        
      ),# end irgt conditional
      
      # tags$u(h3("Other design options")),
      
      
      #### HET CRT ####
      conditionalPanel(
        condition = "input.trial == 'het_two'",
        
        ## sample size ##
        tags$u(h3("Sample Size/Power options")),
        tags$i(h4("Treatment arm")),
        # treatment clusters #
        conditionalPanel(
          condition = "input.plot_display != 'fixed_power'",
          numericInput(inputId="n1_fix_het",
                       label=withMathJax("Assumed number of clusters in treatment arm (\\(n_1\\))"),
                       10, min=1, max=9999999999)
        ),
        
        
        conditionalPanel(
          condition = "input.plot_display == 'n_v_power'",
          sliderInput(inputId="n1_slide_het",
                      label=withMathJax("Number of clusters range in treatment arm (\\(n_1\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # treatment cluster size #
        numericInput(inputId="m1_fix_het",
                     label=withMathJax("Assumed cluster size in treatment arm (\\(m_1\\))"),
                     10, min=1, max=9999999999),
        conditionalPanel(
          condition = "input.plot_display != 'n_v_power'",
          sliderInput(inputId="m1_slide_het",
                      label=withMathJax("Cluster size range in treatment arm (\\(m_1\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # control clusters #
        tags$i(h4("Control arm")),
        conditionalPanel(
          condition = "input.plot_display != 'fixed_power'",
          numericInput(inputId="n0_fix_het",
                       label=withMathJax("Assumed number of clusters in control arm (\\(n_0\\))"),
                       10, min=1, max=9999999999)
        ),
        conditionalPanel(
          condition = "input.plot_display == 'n_v_power'",
          sliderInput(inputId="n0_slide_het",
                      label=withMathJax("Number of clusters range in control arm (\\(n_0\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # control cluster size #
        numericInput(inputId="m0_fix_het",
                     label=withMathJax("Assumed cluster size in control arm (\\(m_0\\))"),
                     10, min=1, max=9999999999),
        conditionalPanel(
          condition = "input.plot_display != 'n_v_power'",
          sliderInput(inputId="m0_slide_het",
                      label=withMathJax("Cluster size range in control arm (\\(m_0\\))"),
                      min=1, max=3000,
                      value=c(5,100))
        ),
        conditionalPanel(
          condition = "input.plot_display == 'fixed_power'",
          numericInput(inputId="power_het",
                       label="Power",
                       0.9, min=0, max=1)
        ),
        
        
        
        
        
        tags$u(h3("ICC options")),
        #tags$i(h4("Outcome")),
        numericInput(inputId="oicc_trt_est_het",
                     label=withMathJax("Assumed treatment-arm outcome ICC (\\(\\rho_{y|x,1}\\))"),
                     0.5,min=0, max=1, step=0.001),
        numericInput(inputId="oicc_ctrl_est_het",
                     label=withMathJax("Assumed control-arm outcome ICC (\\(\\rho_{y|x,0}\\))"),
                     0.1,min=0, max=1, step=0.001),
        numericInput(inputId="cicc_est_het",
                     label=withMathJax("Assumed covariate ICC (\\(\\rho_x\\))"),
                     0.5,min=0, max=1, step=0.001),
        
        
        radioButtons(inputId="icc_sensitivity_het",
                     label="ICC sensitivity analyses",
                     choiceNames = c("Only display results for assumed ICCs",
                                     "Display results for ICC ranges"),
                     choiceValues = c("est_only", "sensitivity"),
                     inline=T),
        conditionalPanel(
          condition = "input.icc_sensitivity_het == 'sensitivity'",
          sliderInput(inputId="oicc_trt_range_het",
                      label="Treatment-arm outcome ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          sliderInput(inputId="oicc_ctrl_range_het",
                      label="Control-arm outcome ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          sliderInput(inputId="cicc_range_het",
                      label="Covariate ICC range",
                      min=0, max=1,
                      value=c(0,0.9)),
          radioButtons(inputId="icc_constant",
                       label="ICC to stay constant across plots",
                       choiceNames = c("Treatment-arm outcome ICC",
                                       "Control-arm outcome ICC",
                                       "Covariate ICC"),
                       choiceValues = c("oicc1","oicc0","cicc")),
          radioButtons(inputId="icc_constant_within",
                       label="ICC to stay constant within plot",
                       choiceNames = c(#"Treatment-arm outcome ICC",
                         "Control-arm outcome ICC",
                         "Covariate ICC"),
                       choiceValues = c(#"oicc1",
                         "oicc0","cicc"))
          
        )#end ICC sensitivity conditional
        
        
        
      ),# end irgt conditional
      
      # tags$u(h3("Other design options")),
      
      
      #### crossovers & SW-CRT & upload DESIGN ####
      conditionalPanel(
        condition = "input.trial == 'parallel_m' || input.trial == 'crossover_2' || input.trial == 'crossover_m' || input.trial == 'SWD'||input.trial=='upload'",
        
        ## sample size ##
        # always numeric periods #
        conditionalPanel(
          condition = "input.trial == 'SWD'",
          numericInput(inputId="J_1",
                       label="Number of sequences (i.e., 'steps')",
                       2, min=2, max=100)
        ),
        
        # specify number of periods of multi-period crossover #
        conditionalPanel(
          condition = "input.trial == 'parallel_m' || input.trial == 'crossover_m'",
          numericInput(inputId="J",
                       label="Number of periods",
                       2, min=2, max=100)
        ),
        
        
        # numeric n #
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power'",
          numericInput(inputId="ns_swd",
                       label="Number of clusters (per sequence)",
                       1, min=1, max=9999999999)
        ),
        # m range #
        conditionalPanel(
          condition = "input.plot_display == 'm_v_power' || input.plot_display == 'fixed_power'",
          sliderInput(inputId="m_swd_range",
                      label="Plot cluster size (per period) range",
                      min=1, max=3000,
                      value=c(5,100))
        ),
        #numeric m, n range #
        conditionalPanel(
          condition = "input.plot_display == 'n_v_power'",
          numericInput(inputId="m_swd",
                       label="Cluster size (per period)",
                       20, min=1, max=9999999999),
          sliderInput(inputId="ns_swd_range",
                      label="Plot number of clusters (per sequence) range",
                      min=1, max=3000,
                      value=c(5,100))
        ),
        # numeric power #
        conditionalPanel(
          condition = "input.plot_display == 'fixed_power'",
          numericInput(inputId="power_swd",
                       label="Power",
                       0.9, min=0, max=1)
        ),
        
        
        
        tags$u(h3("ICC options")),
        conditionalPanel(
          condition = "input.cohort == 'cross'",
          tags$i(h4("Outcome")),
          numericInput(inputId="oicc_wperiod_est_swd",
                       label=withMathJax("Assumed within-period ICC (\\(\\alpha_0\\))"),
                       0.1,min=0, max=1, step=0.001),
          numericInput(inputId="oicc_ratio_est_swd",
                       label=withMathJax("Assumed CAC (\\(\\alpha_1/\\alpha_0\\))"),
                       0.5,min=0, max=1, step=0.001),
          helpText("The ratio of the between-period outcome ICC over the within-period outcome ICC.", style="margin-top:-0.5em; margin-bottom:1em;"),
          
          tags$i(h4("Covariate")),
          numericInput(inputId="cicc_wperiod_est_swd",
                       label=withMathJax("Assumed within-period ICC (\\(\\rho_0\\))"),
                       0.1,min=0, max=1, step=0.001),
          numericInput(inputId="cicc_ratio_est_swd",
                       label=withMathJax("Assumed CAC (\\(\\rho_1/\\rho_0\\))"),
                       0.5,min=0, max=1, step=0.001),
          helpText("The ratio of the between-period covariate ICC over the within-period covariate ICC.", style="margin-top:-0.5em; margin-bottom:1em;"),
        ),#end cross conditional
        
        conditionalPanel(
          condition = "input.cohort == 'closed'",
          tags$i(h4("Outcome")),
          numericInput(inputId="oicc_wperiod_est_swd_cc",
                       label=withMathJax("Assumed within-period ICC (\\(\\alpha_0\\))"),
                       0.1,min=0, max=1, step=0.001),
          numericInput(inputId="oicc_ratio_est_swd_cc",
                       label=withMathJax("Assumed CAC (\\(\\alpha_1/\\alpha_0\\))"),
                       0.5,min=0, max=1, step=0.001),
          helpText("The ratio of the between-period outcome ICC over the within-period outcome ICC.", style="margin-top:-0.5em; margin-bottom:1em;"),
          numericInput(inputId="oicc_windiv_est_swd_cc",
                       label=withMathJax("Assumed within-individual ICC (\\(\\alpha_2\\))"),
                       0.1,min=0, max=1, step=0.001),
          
          tags$i(h4("Covariate")),
          numericInput(inputId="cicc_wperiod_est_swd_cc",
                       label=withMathJax("Assumed ICC (\\(\\rho_0\\))"),
                       0.1,min=0, max=1, step=0.001)
        ),#end closed conditional
        
        radioButtons(inputId="icc_sensitivity_swd",
                     label="ICC sensitivity analyses",
                     choiceNames = c("Only display results for assumed ICCs",
                                     "Display results for ICC ranges"),
                     choiceValues = c("est_only", "sensitivity"),
                     inline=T),
        conditionalPanel(
          condition = "input.icc_sensitivity_swd == 'sensitivity'",
          
          radioButtons(inputId="icc_display_swd",
                       label="ICC display",
                       choiceNames = c("Outcome ICCs constant within-plot", "Covariate ICCs constant within-plot"),
                       choiceValues=c("oICC_constant", "cICC_constant"),
                       inline=T),
          ## play with font bolding/italics/color for section titles to make sections pop ##
          # h4("Outcome ICC"),
          conditionalPanel(
            condition="input.cohort == 'cross'",
            tags$i(h4("Outcome")),
            sliderInput(inputId="oicc_wperiod_range_swd",
                        label="Within-period ICC range",
                        min=0, max=1,
                        value=c(0,0.9)),
            sliderInput(inputId="oicc_ratio_range_swd",
                        label="CAC range",
                        min=0, max=1,
                        value=c(0.01,0.9)),
            
            
            tags$i(h4("Covariate")),
            sliderInput(inputId="cicc_wperiod_range_swd",
                        label="Within-period ICC range",
                        min=0, max=1,
                        value=c(0,0.9)),
            sliderInput(inputId="cicc_ratio_range_swd",
                        label="CAC range",
                        min=0, max=1,
                        value=c(0.01,0.9))
          ),
          
          conditionalPanel(
            condition="input.cohort == 'closed'",
            
            tags$i(h4("Outcome")),
            sliderInput(inputId="oicc_wperiod_range_swd_cc",
                        label="Within-period ICC range",
                        min=0, max=1,
                        value=c(0,0.9)),
            sliderInput(inputId="oicc_ratio_range_swd_cc",
                        label="CAC range",
                        min=0, max=1,
                        value=c(0.01,0.9)),
            sliderInput(inputId="oicc_windiv_range_swd_cc",
                        label="Within-individual ICC range",
                        min=0, max=1,
                        value=c(0.01,0.9)),
            
            
            tags$i(h4("Covariate")),
            sliderInput(inputId="cicc_wperiod_range_swd_cc",
                        label="ICC range",
                        min=0, max=1,
                        value=c(0,0.9))
          )
          
        )#end ICC sensitivity conditional
        
        
      ),# end swd conditional
      #### UNIVERSAL OPTIONS ####
      ## outcomes ##
      # irgt can only do continuous outcome - have done updateRadioButtons() for this server-side #
      tags$u(h3("Outcome and variable options")),
      radioButtons(inputId="outcome",
                   label="Outcome type",
                   choiceNames = c("Continuous", "Binary"),
                   choiceValues=c("continuous", "binary"),
                   inline=T),
      # effect size and sd for continuous #
      conditionalPanel(
        condition = "input.outcome == 'continuous' & input.trial != 'irgt' & input.trial != 'het_two'",
        # numericInput(inputId="mean_diff_ATE",
        #              label="Assumed ATE",
        #              1, min=0, max=999999),
        numericInput(inputId="sd_outcome",
                     label="Outcome standard deviation",
                     1, min=0, max=999999)
      ),#end conditional
      conditionalPanel(
        condition = "input.trial == 'irgt' || input.trial == 'het_two'",
        numericInput(inputId="sd_outcome1",
                     label="Treatment-arm outcome standard deviation",
                     1, min=0, max=999999),
        numericInput(inputId="sd_outcome0",
                     label="Control-arm outcome standard deviation",
                     1, min=0, max=999999)
      ),#end conditional
      # effect size for binary #
      conditionalPanel(
        condition = "input.outcome == 'binary'",
        numericInput(inputId="prop_control",
                     # outcome proportion/event prop of those on control #
                     label=withMathJax("Mean outcome under control (\\(p_0\\))"),
                     0.5, min=0, max=1,step=0.001),
        numericInput(inputId="prop_trt",
                     # outcome proportion/event prop of those on intervention #
                     label=withMathJax("Mean outcome under intervention (\\(p_1\\))"),
                     0.5, min=0, max=1, step=0.001)
      ),#end conditional
      
      # covariate #
      radioButtons(inputId="covar",
                   label="Covariate type",
                   choiceNames = c("Continuous", "Binary"),
                   choiceValues=c("continuous", "binary"),
                   inline=T),
      # can keep this as specifying HTE effect size since HTE for binary outcome doesn't translate directly into diff of proportions
      numericInput(inputId="mean_diff_HTE",
                   label="Assumed HTE",
                   1, min=0, max=999999),
      helpText("Specify the target effect size for the effect modification, e.g., the difference in difference estimate for a binary effect modifier.", style="margin-top:-0.5em; margin-bottom:1em;"),
      # covariate sd for continuous #
      conditionalPanel(
        condition = "input.covar == 'continuous'",
        numericInput(inputId="sd_covar",
                     label="Covariate standard deviation",
                     1, min=0, max=999999)
      ),#end conditional
      # covariate sd for binary #
      conditionalPanel(
        condition = "input.covar == 'binary'",
        numericInput(inputId="prop_covar",
                     label="Covariate proportion",
                     0.5, min=0, max=1,step=0.001),
        helpText("e.g., prevalence of the binary covariate.", style="margin-top:-0.5em; margin-bottom:1em;"),
      ),#end conditional
      
      #treatment allocation #
      conditionalPanel(
        condition = "input.trial == 'parallel' || input.trial == 'three_level'",
        numericInput(inputId="w",
                     label="Intervention allocation",
                     0.5, min=0.000001, max=1, step=0.001),
        helpText("Enter the proportion of total participants in the intervention arm.", style="margin-top:-0.5em; margin-bottom:1em;"),
      ),
      # sig level #
      numericInput(inputId="sig",
                   label="Significance level",
                   0.05, min=0.000001, max=1, step=0.001),
      ## loading message when the app is calculating ##
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 8%;
               left: 33.5%;
               width: 50%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 150%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),
      conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                       tags$div("Loading...",id="loadmessage"))
      
    ),#end sidebar panel
    
    
    #### Output panel ####
    mainPanel(
      style="position:fixed;margin-left:32vw;",
      
      tabsetPanel(
        tabPanel("Power",#br(),
                 plotlyOutput("powerPlot_hte")), 
        tabPanel("Design Matrix",#br(),
                 tableOutput("design_matrix")),
        tabPanel("References and Resources", 
                 h3("References and Resources"),
                 p(em("This application has been written by Mary M. Ryan (University of Wisconsin - Madison USA), Fan Li (Yale School of Public Health USA), and with input from Monica Taljaard (University of Ottawa - Ottawa Hospital Research Institute CA). Please email mary.ryan@wisc.edu if you need to report errors or would like to submit comments or feedback."), style = "font-size:13pt;"),
                 p(tags$ul(
                   tags$li(strong("Parallel two-level design"),": Yang S., Li F., Starks M., Hernandez A.F., Mentz R.J., Choudhury K.R. Sample size requirements for detecting treatment effect heterogeneity in cluster randomized trials. Stat Med. 2020 July 16;39(28):4218-4237. doi: 10.1002/sim.8721. Epub 2020 Aug 21. PubMed PMID: 32823372", style = "font-size:13pt;"),
                   tags$ul(
                     tags$li(strong("Heterogeneous ICCs and variances"),":  Tong G., Taljaard M., Li F. Sample size considerations for assessing treatment effect heterogeneity in randomized trials with heterogeneous intracluster correlations and variances. Stat Med. 2023; 42(19): 3392–3412. doi: 10.1002/sim.9811. PubMed PMID: 37316956", style = "font-size:13pt;"),
                   ),
                   tags$li(strong("Parallel three-level design"),": Li F., Chen X., Tian Z., Esserman D., Heagerty P.J., Wang R. Designing three-level cluster randomized trials to assess treatment effect heterogeneity . Biostatistics. 2022 July 21. doi: 10.1093/biostatistics/kxac026. PubMed PMID: 35861621", style = "font-size:13pt;"),
                   tags$li(strong("Two-period and multiple-period cross-over designs"),": Li F., Chen X., Tian Z., Wang R., Heagerty P.J. Planning stepped wedge cluster randomized trials to detect treatment effect heterogeneity. Stat Med. 2024; 43(5): 890-911. doi: 10.1002/sim.9990. PubMed PMID: 38115805", style = "font-size:13pt;"),
                   tags$li(strong("Stepped-wedge design"),": Li F., Chen X., Tian Z., Wang R., Heagerty P.J. Planning stepped wedge cluster randomized trials to detect treatment effect heterogeneity. Stat Med. 2024; 43(5): 890-911. doi: 10.1002/sim.9990. PubMed PMID: 38115805", style = "font-size:13pt;"),
                   tags$li(strong("Individually randomized group treatment design"),":  Tong G., Taljaard M., Li F. Sample size considerations for assessing treatment effect heterogeneity in randomized trials with heterogeneous intracluster correlations and variances. Stat Med. 2023; 42(19): 3392–3412. doi: 10.1002/sim.9811. PubMed PMID: 37316956", style = "font-size:13pt;"),
                   tags$li(strong("How to estimate ICCs"),":  Ouyang Y., Hemming K., Li F., Taljaard M.  Estimating intra-cluster correlation coefficients for planning longitudinal cluster randomized trials: A tutorial. International Journal of Epidemiology. 2023; 52(5): 1634–1647. doi: 10.1093/ije/dyad062. PubMed PMID: 37196320, PubMed Central ID: PMC10555741", style = "font-size:13pt;")
                 )),
                 br(),
                 p("To calculate ATE power for CRTs, please visit The Shiny CRT Calculator developed by Karla Hemming (University of Birmingham UK) and Jesica Kasza (Monash University Australia):",a("https://clusterrcts.shinyapps.io/rshinyapp/", href="https://clusterrcts.shinyapps.io/rshinyapp/"), style = "font-size:13pt;")
        )
      )
      
      
      
      # fluidRow(
      #   splitLayout(cellWidths = c("50%", "50%"),
      #               plotlyOutput("powerPlot_hte_n"), plotlyOutput("powerPlot_hte_m"))
      # ),
      # fluidRow(
      #   HTML("<br>")
      # ),
      # fluidRow(
      #   splitLayout(cellWidths = c("50%", "50%"),
      #               plotlyOutput("powerPlot_hte_m"), plotlyOutput("powerPlot_hte_m"))
      # )
    )
  )
))
