#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)
library(ggplot2)
library(ggpubr)
library(scales)
library(viridis)
library(gt)
library(tidyverse)
library(dplyr)
library(rmarkdown)
library(markdown)
#library(huxtable)

model = function(t, x, model_par){
    
    #Effectiveness
    pi_alpha = model_par[1]
    pi_beta = model_par[2]
    #Relative value
    phi_alpha = 1 - pi_alpha
    phi_beta = 1 - pi_beta
    
    #R0 N = intervenable, O= other
    R0_N = model_par[3]
    R0_O = model_par[4]
    
    
    
    I = x[1]
    I_N = x[2]
    
    S = x[3]
    S_N = x[4]
    
    dxdt = numeric(length(x))
    dxdt[1]  =  S * ( R0_N * (I + I_N*phi_alpha) + R0_O * (I + I_N)) -  I
    dxdt[2]  =  S_N * (phi_beta * R0_N * (I + I_N*phi_alpha) + 
                           R0_O * (I + I_N)) -  I_N
    
    dxdt[3] =   I - S * ( R0_N * (I + I_N*phi_alpha) + R0_O * (I + I_N)) 
    dxdt[4] =   I_N - S_N * (phi_beta * R0_N * (I + I_N*phi_alpha) + 
                                 R0_O * (I + I_N)) 
    
    return(list(dxdt))
    
}

simulate = function(par){
    
    par_baseline_condition = par[1]
    par_intervention_compliance = par[2]
    R0 = par[3]
    R0_ratio = par[4]
    efficacy_alpha = par[5]
    efficacy_beta= par[6]
    coverage=par[7]
    
    rho_vec = c(1-par_intervention_compliance, par_intervention_compliance)
    baseline_adherence_vec = c(1-par_baseline_condition, par_baseline_condition)
    
    prev= 0.06
    x0 = c(rep(prev,2)*(rho_vec*coverage + baseline_adherence_vec*(1-coverage)),
           rep(1-prev, 2)*(rho_vec*coverage + baseline_adherence_vec*(1-coverage)))
    
    R0_N= R0*R0_ratio
    R0_O = R0*(1-R0_ratio)
    model_par = c(efficacy_alpha,efficacy_beta, R0_N,R0_O)
    out_coverage = ode(x0, times = seq(0,100), model, model_par,method="vode")
    steady_state = tail(out_coverage[,2:5],1)
    prevalence  = sum(steady_state[1:2])
    
    return(prevalence)
}

fill_matrix=function(par0,factor1,factor2,res,intervention_effectiveness_actual){
  flag = par0[5]==0
  
    matrix = matrix(NA, res*res,3)
    
    for (i in 1:res){
        for (j in 1:res){
            par=par0
            
            if (factor1=="conditions"){par[1] = seq(0,1,length.out = res)[i]
            factor1_temp=par[1]}
            if (factor1=="compliance"){par[2] = seq(0,1,length.out = res)[i]
            factor1_temp=par[2]}
            if (factor1=="R0"){par[3] = seq(1,1.5,length.out = res)[i]
            factor1_temp=par[3]}
            if (factor1=="completeness"){par[4] = seq(0,1,length.out = res)[i]
            factor1_temp=par[4]}
            if (factor1=="efficacy"){par[5] = seq(0,1,length.out = res)[i]*(1-flag)
            par[6] = seq(0,1,length.out = res)[i]*flag
            factor1_temp=seq(0,1,length.out = res)[i]}
            if (factor1=="coverage"){par[7] = seq(0,1,length.out = res)[i]
            factor1_temp=par[7]}
            
            if (factor2=="conditions"){par[1] = seq(0,1,length.out = res)[j]
            factor2_temp=par[1]}
            if (factor2=="compliance"){par[2] = seq(0,1,length.out = res)[j]
            factor2_temp=par[2]}
            if (factor2=="R0"){par[3] = seq(1,1.5,length.out = res)[j]
            factor2_temp=par[3]}
            if (factor2=="completeness"){par[4] = seq(0,1,length.out = res)[j]
            factor2_temp=par[4]}
            if (factor2=="efficacy"){par[5] = seq(0,1,length.out = res)[j]*(1-flag)
            par[6] = seq(0,1,length.out = res)[j]*flag
            factor2_temp=seq(0,1,length.out = res)[j]}
            if (factor2=="coverage"){par[7] = seq(0,1,length.out = res)[j]
            factor2_temp=par[7]}
            
            if (par[1]>par[2]){
                change_intervention_effectiveness = NA
            }else{
                prev_control = simulate(c(par[1:6],0))
                if (prev_control<1E-4){
                    change_intervention_effectiveness = NA
                }
                prev_intervention = simulate(par)
                relative_risk_counterfactual = prev_intervention/prev_control
                intervention_effectiveness_counterfactual = 1- relative_risk_counterfactual
                
                change_intervention_effectiveness = 100*(intervention_effectiveness_counterfactual-intervention_effectiveness_actual)
            }
            
            matrix[j+(i-1)*res,]=c(factor2_temp, factor1_temp, change_intervention_effectiveness)
        }
    }
    matrix = as.data.frame(matrix)
    colnames(matrix)=c("x","y","value")
    return(matrix)
}


# Define UI for application that draws a histogram
ui <- fluidPage(
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style-page.css")),
    navbarPage(HTML("Modeling Intervention Effectiveness"),
    # Application title
    #titlePanel("Change in Intervention Effectiveness"),
    tabPanel("Intervention Effectiveness",
            # Sidebar with a slider input for number of bins 
            sidebarLayout(
                sidebarPanel(id = "module_sidebar",
                    sliderInput("par_baseline_condition",
                                "Baseline WASH conditions. The proportion of population with existing WASH infrastructure similar to the intervention.",
                                min = 0,
                                max = 100,
                                value = 25, post = "%"), 
                    sliderInput("par_intervention_compliance",
                                "Intervention compliance. The proportion of the intervention group that receives and uses WASH infrastructure.",
                                min = 0,
                                max = 100,
                                value = 75, post = "%"), 
                    #Display error if par_intervention_compliance < par_baseline_condition
                    #"Intervention compliance must be greater than baseline conditions.
                    sliderInput("R0", 
                                 "Basic reproduction number. The number of expected cases generated by one case in an otherwise susceptible population.", 
                                 min = 1, 
                                 max = 2, 
                                value = 1.25), 
                    sliderInput("R0_ratio", 
                                 "Intervenable fraction. Proportion of overall transmission controlled if the intervention had perfect efficacy for those using it.", 
                                min = 0,
                                max = 100,
                                value = 75, post = "%"), 
                    radioButtons("flag", 
                                 "Pathogen control mechanism. Whether the intervention acts to reduce  shedding into the environment (e.g., latrines) or transmission from the environment to people (e.g., handwashing).",
                                 choices = list("Shedding" = 0, "Transmission" = 1),selected = 1), 
                    sliderInput("efficacy",
                                "Intervention efficacy. The percent of shedding/transmission prevented by the intervention.",
                                min = 0,
                                max = 100,
                                value = 75, post = "%"), 
                    
                    sliderInput("coverage",
                                "Community coverage. The percent of total population within the intervention area that is targeted to receive the intervention.",
                                min = 0,
                                max = 100,
                                value = 11, post = "%"), 
                    radioButtons("make_plot", 
                                 "Run sensitivity analysis?",
                                 choices = list("Yes", "No"),selected = "No")
                    
                ),
        
                # Show a plot of the generated distribution
                mainPanel(
                   plotOutput("dp1"), 
                   br(),
                   gt_output("dp1_table"),
                   br(), 
                   br(),
                   conditionalPanel(
                       condition = "input.make_plot == 'Yes'",
                       plotOutput("main_matrix", height = "800px")
                   )
                )
            )
        ), # tab
    tabPanel("About",
             mainPanel(fluidRow(column(12, includeMarkdown(rmarkdown::render("about.Rmd")))))
    ) # second tab
)
)
# Define server logic required to draw a histogram
server <- function(input, output) {

    
    data_one <- reactive({
        
        par = c(as.numeric(input$par_baseline_condition)/100, 
                as.numeric(input$par_intervention_compliance)/100,
                as.numeric(input$R0), 
                as.numeric(input$R0_ratio)/100, 
                (1-as.numeric(input$flag))*as.numeric(input$efficacy)/100,  
                as.numeric(input$flag)*as.numeric(input$efficacy)/100, 
                as.numeric(input$coverage)/100)
        
        
        
        
        #Baseline control
        prev_control = simulate(c(par[1:6],0))
        #Baseline intervention
        prev_intervention = simulate(par)
        #Baseline max intervention
        prev_intervention_max = simulate(c(par[1],1,par[3:4],(1-as.numeric(input$flag)),as.numeric(input$flag),par[7]))
        
        data_plot1=as.data.frame(cbind(c("Control","Control","Intervention","Intervention"),
                                       100*c(prev_intervention_max,prev_control-prev_intervention_max,
                                             prev_intervention_max,prev_intervention-prev_intervention_max),
                                       c("Non-intervenable","Intervenable","Non-intervenable","Intervenable")))
        colnames(data_plot1)=c("Arm","Prevalence","Transmission")
        data_plot1$Prevalence = as.numeric(data_plot1$Prevalence )
        
        
        
        return(data_plot1)
        
    })
    
    data_one_table <- reactive({
        
        par = c(as.numeric(input$par_baseline_condition)/100, 
                as.numeric(input$par_intervention_compliance)/100,
                as.numeric(input$R0), 
                as.numeric(input$R0_ratio)/100, 
                (1-as.numeric(input$flag))*as.numeric(input$efficacy)/100,  
                as.numeric(input$flag)*as.numeric(input$efficacy)/100, 
                as.numeric(input$coverage)/100)
        
        #Baseline control
        prev_control = simulate(c(par[1:6],0))
        #Baseline intervention
        prev_intervention = simulate(par)
        #Baseline max intervention
        prev_intervention_max = simulate(c(par[1],1,par[3:4],(1-as.numeric(input$flag)),as.numeric(input$flag),par[7]))
        
        relative_risk_actual = prev_intervention/prev_control
        intervention_effectiveness_actual = 1- relative_risk_actual
        
        #Display: "Prevalence in control arm"
        a <- round(prev_control*100,digits=1)
        #Display: "Prevalence in intervention arm"
        b <- round(prev_intervention*100,digits=1)
        #Dispaly: "Non-intervenable prevalence"
        c <- round(prev_intervention_max * 100,digits=1)
        #Display: "Intervention effectiveness. What percent of disease is prevented by the intervention?"
        d <- round(intervention_effectiveness_actual * 100,digits=1) #Display %

        
        ftable <- data.frame(prev_control = a, prev_intervention = b, prev_min = c, actual_effect = d)

        return(ftable)

    })
    
    output$dp1 <- renderPlot({
        
        validate(
            need(input$par_intervention_compliance > input$par_baseline_condition, 'Intervention compliance must be greater than baseline conditions.'), 
            need(input$R0 >= 1 & input$R0 <= 2, "Basic reproduction number needs to be between 1 and 2."), 
            need(input$R0_ratio >= 0 & input$R0_ratio <= 100, "Intervenable fraction needs to be between 0 and 1")
        )
        
        data_plot1 <- data_one()
        
        #print(data_plot1)
        
        #Display this plot
        ggplot(data_plot1,aes(x = Arm, y = Prevalence, fill = Transmission)) +
            geom_col() +
            theme_classic() +
            #ylab("Prevalence (%)")+
            labs(x = "", 
                 y = "Prevalence (%)", 
                 caption = "Simulated infection prevalence and intervention effectiveness in the randomized controlled trial.\nNon-intervenable prevalence is the infection prevalence that would remain with perfect intervention compliance and efficacy.") +
            scale_fill_manual(values = c("grey75", "grey25")) +
          theme(axis.text.x = element_text(size = 20), 
                axis.text.y = element_text(size = 20), 
                axis.title.x = element_text(size = 20), 
                axis.title.y = element_text(size = 20)) 
        
    })
    
    output$dp1_table <- render_gt({
        
        validate(
          need(input$par_intervention_compliance > input$par_baseline_condition, 'Intervention compliance must be greater than baseline conditions.'), 
          need(input$R0 >= 1 & input$R0 <= 2, "Basic reproduction number needs to be between 1 and 2."), 
          need(input$R0_ratio >= 0 & input$R0_ratio <= 100, "Intervenable fraction needs to be between 0 and 1")
        )
        
        data_point_list <- data_one_table()
        
        data_point_list %>% 
            gt() %>% 
            cols_label(
                prev_control = "Prevalence in control arm",
                prev_intervention = "Prevalence in intervention arm",
                prev_min = "Non-intervenable prevalence",
                actual_effect = "Intervention effectiveness"
            ) %>% cols_align(
                align = c("center"),
                columns = everything()
            ) %>% 
            tab_style(
                style = cell_borders(
                    sides = c("all"),
                    color = "#000000",
                    weight = px(3),
                    style = "solid"
                ),
                locations = cells_body(
                    columns = everything(),
                    rows = everything()
                )
            ) %>% tab_style(
                locations = cells_column_labels(columns = everything()),
                style     = list(
                    #Give a thick border below
                    cell_borders(sides = c("all"),
                                 color = "#000000",
                                 weight = px(3)),
                    #Make text bold
                    cell_text(weight = "bold")
                )
            ) %>% fmt_number(
                columns = c("actual_effect"),
                decimals = 1,
                pattern = "{x}%"
                
            )
        
    })
    
    
    output$main_matrix <- renderPlot({
        
        breaks <- input$make_plot
        if (breaks == "Yes") {
        
                validate(
                  need(input$par_intervention_compliance > input$par_baseline_condition, 'Intervention compliance must be greater than baseline conditions.'), 
                  need(input$R0 >= 1 & input$R0 <= 2, "Basic reproduction number needs to be between 1 and 2."), 
                  need(input$R0_ratio >= 0 & input$R0_ratio <= 100, "Intervenable fraction needs to be between 0 and 1")
                )
                
                
                par = c(as.numeric(input$par_baseline_condition)/100, 
                        as.numeric(input$par_intervention_compliance)/100,
                        as.numeric(input$R0), 
                        as.numeric(input$R0_ratio)/100, 
                        (1-as.numeric(input$flag))*as.numeric(input$efficacy)/100,  
                        as.numeric(input$flag)*as.numeric(input$efficacy)/100, 
                        as.numeric(input$coverage)/100)
                
                #Baseline control
                prev_control = simulate(c(par[1:6],0))
                #Baseline intervention
                prev_intervention = simulate(par)
                #Baseline max intervention
                prev_intervention_max = simulate(c(par[1],1,par[3:4],(1-as.numeric(input$flag)),as.numeric(input$flag),par[7]))
                
                relative_risk_actual = prev_intervention/prev_control
                intervention_effectiveness_actual = 1- relative_risk_actual
                
        
                #conditions, compliance, R0, R0_1/R0, efficacy, coverage
                col_limits = 100*c(0-intervention_effectiveness_actual,1-intervention_effectiveness_actual)
                
                #Resolution
                res = 25
                
                withProgress(message = 'Making plot...', value = 0, {
                
                  incProgress(1/15, detail = "Generating Plot 1")
                  
                #conditions, compliance
                matrix1= fill_matrix(par,"conditions","compliance",res,intervention_effectiveness_actual)
                p1= ggplot()+
                  geom_raster(data=matrix1, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[2],y=par[1],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  xlab("Compliance")+ylab("")+#ylab("Baseline conditions")+
                  theme_classic()
                
                incProgress(1/15, detail = "Generating Plot 2")
                
                #conditions, R0
                matrix2= fill_matrix(par,"conditions","R0",res,intervention_effectiveness_actual)
                p2= ggplot()+
                  geom_raster(data=matrix2, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[3],y=par[1],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  xlab(expression(R[0]))+ylab("")+#+ylab("Baseline conditions")+
                  theme_classic()
                
                incProgress(1/15, detail = "Generating Plot 3")
                
                #conditions, completeness
                matrix3= fill_matrix(par,"conditions","completeness",res,intervention_effectiveness_actual)
                p3= ggplot()+
                  geom_raster(data=matrix3, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[4],y=par[1],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  xlab("Intervenable fraction")+ylab("")+#+ylab("Baseline conditions")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 4")
                
                #conditions, efficacy
                matrix4= fill_matrix(par,"conditions","efficacy",res,intervention_effectiveness_actual)
                p4= ggplot()+
                  geom_raster(data=matrix4, aes(x=x,y=y,fill=value))+
                  annotate("point",x=c(par[5],par[6])[par[5]+1],y=par[1],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  xlab("Efficacy")+ylab("")+#+ylab("Baseline conditions")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 5")
                
                #conditions, coverage
                matrix5= fill_matrix(par,"conditions","coverage",res,intervention_effectiveness_actual)
                p5= ggplot()+
                  geom_raster(data=matrix5, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[7],y=par[1],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  xlab("Coverage fraction")+ylab("Baseline conditions")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 6")
                
                #compliance, R0
                matrix6= fill_matrix(par,"compliance","R0",res,intervention_effectiveness_actual)
                p6= ggplot()+
                  geom_raster(data=matrix6, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[3],y=par[2],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab(expression(R[0]))+ylab("Compliance")+
                  xlab("")+ylab("")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 7")
                
                #compliance, completeness
                matrix7= fill_matrix(par,"compliance","completeness",res,intervention_effectiveness_actual)
                p7= ggplot()+
                  geom_raster(data=matrix7, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[4],y=par[2],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Intervenable fraction")+ylab("Compliance")+
                  xlab("")+ylab("")+
                  theme_classic()
                
                incProgress(1/15, detail = "Generating Plot 8")

                #compliance, efficacy
                matrix8= fill_matrix(par,"compliance","efficacy",res,intervention_effectiveness_actual)
                p8= ggplot()+
                  geom_raster(data=matrix8, aes(x=x,y=y,fill=value))+
                  annotate("point",x=c(par[5],par[6])[par[5]+1],y=par[2],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Efficacy")+ylab("Compliance")+
                  xlab("")+ylab("")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 9")
                
                #compliance, coverage
                matrix9= fill_matrix(par,"compliance","coverage",res,intervention_effectiveness_actual)
                p9= ggplot()+
                  geom_raster(data=matrix9, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[7],y=par[2],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Coverage fraction")
                  xlab("")+ylab("Compliance")+
                  theme_classic()
                
                incProgress(1/15, detail = "Generating Plot 10")
                
                #R0, completeness
                matrix10= fill_matrix(par,"R0","completeness",res,intervention_effectiveness_actual)
                p10= ggplot()+
                  geom_raster(data=matrix10, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[4],y=par[3],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Intervenable fraction")+ylab(expression(R[0]))+
                  xlab("")+ylab("")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 11")
                
                #R0, efficacy
                matrix11= fill_matrix(par,"R0","efficacy",res,intervention_effectiveness_actual)
                p11= ggplot()+
                  geom_raster(data=matrix11, aes(x=x,y=y,fill=value))+
                  annotate("point",x=c(par[5],par[6])[par[5]+1],y=par[3],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Efficacy")+ylab(expression(R[0]))+
                  xlab("")+ylab("")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 12")
                
                #R0, coverage
                matrix12= fill_matrix(par,"R0","coverage",res,intervention_effectiveness_actual)
                p12= ggplot()+
                  geom_raster(data=matrix12, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[7],y=par[3],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Coverage fraction")
                  xlab("")+ylab(expression(R[0]))+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 13")
                
                #completeness, efficacy
                matrix13= fill_matrix(par,"completeness","efficacy",res,intervention_effectiveness_actual)
                p13= ggplot()+
                  geom_raster(data=matrix13, aes(x=x,y=y,fill=value))+
                  annotate("point",x=c(par[5],par[6])[par[5]+1],y=par[4],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Efficacy")+ylab("Intervenable fraction")+
                  xlab("")+ylab("")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 14")
                
                #completeness, coverage
                matrix14= fill_matrix(par,"completeness","coverage",res,intervention_effectiveness_actual)
                p14= ggplot()+
                  geom_raster(data=matrix14, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[7],y=par[4],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  # xlab("Coverage fraction")
                  xlab("")+ylab("Intervenable fraction")+
                  theme_classic()

                incProgress(1/15, detail = "Generating Plot 15")
                
                #efficacy, coverage
                matrix15= fill_matrix(par,"efficacy","coverage",res,intervention_effectiveness_actual)
                p15= ggplot()+
                  geom_raster(data=matrix15, aes(x=x,y=y,fill=value))+
                  annotate("point",x=par[7],y=c(par[5],par[6])[par[5]+1],pch=19,col="white",size=3)+
                  scale_fill_gradientn(colors = viridis(10),values=rescale(c(-1,-0.5,-0.25,0,0.25,0.5,1)),limits=col_limits)+
                  guides(fill=guide_colorbar(title="Change in\nintervention\neffectiveness\n(percentage points)        "))+
                  xlab("")+ylab("Efficacy")+
                  theme_classic()

                })
                
                #Can we put a loading bar to indicate % of graphs calculated?
                
                p=ggarrange(p15,NULL,NULL,NULL,NULL,
                            p14,p13,NULL,NULL,NULL,
                            p12,p11,p10,NULL,NULL,
                            p9,p8,p7,p6,NULL,
                            p5,p4,p3,p2,p1,
                            ncol=5,nrow=5,common.legend = TRUE,legend="top")
                p <- annotate_figure(p, bottom = "Intervention effectiveness as a function of WASH intervention factors.\nThe heatmaps denote how intervention effectiveness depends on each pair of WASH factors,\ncompared to the original scenario indicated by the white points.")
                return(p)
        }
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
