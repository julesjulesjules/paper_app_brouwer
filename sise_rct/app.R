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


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Change in Intervention Effectiveness"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("par_baseline_condition",
                        "Baseline conditions. What fraction of population already has WASH infrastructure similar to the intervention?",
                        min = 0,
                        max = 100,
                        value = 25, post = "%"), 
            sliderInput("par_intervention_compliance",
                        "Intervention compliance. What fraction of the intervention group receives and uses the already has WASH infrastructure similar to the intervention?",
                        min = 0,
                        max = 100,
                        value = 75, post = "%"), 
            #Display error if par_intervention_compliance < par_baseline_condition
            #"Intervention compliance must be greater than baseline conditions.
            numericInput("R0", 
                         "Basic reproduction number. Controls the burden of disease. (1 to 1.5)", 
                         value = 1.25, 
                         min = 1, 
                         max = 1.5), 
            numericInput("R0_ratio", 
                         "Intervenable fraction. How much of the transmission could be controlled if the intervention was perfectly effective? (0 to 1)", 
                         value = 0.75, 
                         min = 0, 
                         max = 1), 
            radioButtons("flag", 
                         "Does the intervention reduce shedding into the environment or transmission from the environment to people?",
                         choices = list("Shedding" = 0, "Transmission" = 1),selected = 1), 
            sliderInput("efficacy",
                        "Intervention efficacy. What percent of shedding/transmission is prevented by the intervention?",
                        min = 0,
                        max = 100,
                        value = 75, post = "%"), 
            sliderInput("efficacy_beta",
                        "Efficacy at reducing infection. What percent of potential transmission is prevented by the intervention?",
                        min = 0,
                        max = 100,
                        value = 75, post = "%"), 
            sliderInput("coverage",
                        "What percent of the population is included in the intervention?",
                        min = 0,
                        max = 100,
                        value = 11, post = "%")
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("dp1"), 
           br(),
           gt_output("dp1_table")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    
    data_one <- reactive({
        
        par = c(as.numeric(input$par_baseline_condition)/100, 
                as.numeric(input$par_intervention_compliance)/100,
                as.numeric(input$R0), 
                as.numeric(input$R0_ratio), 
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
                as.numeric(input$R0_ratio), 
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
        #Display: "Intervention effectiveness. What percent of disease is prevented by the intervention?"
        c <- round(intervention_effectiveness_actual * 100,digits=1) #Display %
        
        ftable <- data.frame(prev_control = a, prev_intervention = b, actual_effect = c)
        #colnames(ftable) <- c("prev_control", "prev_intervention", "actual_effect")                  
        
        return(ftable)

    })
    
    output$dp1 <- renderPlot({
        
        validate(
            need(input$par_intervention_compliance > input$par_baseline_condition, 'Intervention compliance must be greater than baseline conditions.')
        )
        
        data_plot1 <- data_one()
        
        #print(data_plot1)
        
        #Display this plot
        ggplot(data_plot1,aes(x=Arm,y=Prevalence,fill=Transmission))+
            geom_col()+theme_classic()+ylab("Prevalence (%)")+
            scale_fill_manual(values=c("grey75","grey25"))
        
    })
    
    output$dp1_table <- render_gt({
        
        data_point_list <- data_one_table()
        
        data_point_list %>% 
            gt() %>% 
            cols_label(
                prev_control = "Prevalence in control arm",
                prev_intervention = "Prevalence in intervention arm",
                actual_effect = "Intervention effectiveness. What percent of disease is prevented by the intervention?"
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
    
}

# Run the application 
shinyApp(ui = ui, server = server)
