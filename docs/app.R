library(shiny)
library(bslib)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reactable)

# Function to calculate combined test performance metrics
calculate_metrics <- function(sens_a, spec_a, sens_b, spec_b, prevalence, sample_size) {
  # Calculate proportions
  base_metrics <- calculate_base_metrics(sens_a, spec_a, sens_b, spec_b, prevalence)
  
  # Add counts based on sample size
  lapply(base_metrics, function(df) {
    df$count <- round(df$value * sample_size)
    df$percentage <- round(df$value * 100, 1)
    df
  })
}

# Base metrics calculation
calculate_base_metrics <- function(sens_a, spec_a, sens_b, spec_b, prevalence) {
  # Single test calculations
  test_a <- calc_proportions(sens_a, spec_a, prevalence)
  test_b <- calc_proportions(sens_b, spec_b, prevalence)
  
  # Parallel AND combination
  p_and_sens <- sens_a * sens_b
  p_and_spec <- spec_a + spec_b - (spec_a * spec_b)
  
  # Parallel OR combination
  p_or_sens <- sens_a + sens_b - (sens_a * sens_b)
  p_or_spec <- spec_a * spec_b
  
  # Serial AND combination
  s_and_sens <- sens_a * sens_b
  s_and_spec <- spec_a + (1 - spec_a) * spec_b
  
  # Serial OR combination
  s_or_sens <- sens_a + (1 - sens_a) * sens_b
  s_or_spec <- spec_a * spec_b
  
  list(
    test_a = calc_proportions(sens_a, spec_a, prevalence),
    test_b = calc_proportions(sens_b, spec_b, prevalence),
    p_and = calc_proportions(p_and_sens, p_and_spec, prevalence),
    p_or = calc_proportions(p_or_sens, p_or_spec, prevalence),
    s_and = calc_proportions(s_and_sens, s_and_spec, prevalence),
    s_or = calc_proportions(s_or_sens, s_or_spec, prevalence)
  )
}

# Calculate proportions for confusion matrix
calc_proportions <- function(sens, spec, prev) {
  tp <- sens * prev
  fn <- prev - tp
  tn <- spec * (1 - prev)
  fp <- (1 - prev) - tn
  
  data.frame(
    category = c("True Positive", "False Negative", "True Negative", "False Positive"),
    value = c(tp, fn, tn, fp)
  )
}

# Function to create mosaic plot
create_mosaic_plot <- function(data, title) {
  ggplot(data, aes(x = "", y = value, fill = category)) +
    geom_bar(stat = "identity", width = 1) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "True Positive" = "#2ECC71",
        "False Negative" = "#E74C3C",
        "True Negative" = "#3498DB",
        "False Positive" = "#F1C40F"
      )
    ) +
    labs(
      title = title,
      x = "",
      y = "Proportion",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = "bottom"
    )
}

# Function to create summary table
create_summary_table <- function(data) {
  reactable(
    data.frame(
      Category = data$category,
      Count = data$count,
      Percentage = sprintf("%.1f%%", data$percentage)
    ),
    compact = TRUE,
    striped = TRUE,
    defaultPageSize = 4
  )
}

ui <- page_sidebar(
  title = "Diagnostic Testing Visualization",
  sidebar = sidebar(
    title = "Test Parameters",
    
    # Disease parameters
    card(
      "Disease Parameters",
      numericInput("prevalence", "Disease Prevalence:", 0.1, min = 0, max = 1, step = 0.01),
      numericInput("sample_size", "Sample Size:", 1000, min = 100, max = 1000000, step = 1)
    ),
    
    # Test A inputs
    card(
      "Test A",
      textInput("name_a", "Test Name:", "Test A"),
      selectInput("timepoint_a", "Time Point:", 
                  choices = c("Baseline", "Follow-up 1", "Not used"),
                  selected = "Baseline"),
      numericInput("sens_a", "Sensitivity:", 0.8, min = 0, max = 1, step = 0.01),
      numericInput("spec_a", "Specificity:", 0.9, min = 0, max = 1, step = 0.01)
    ),
    
    # Test B inputs
    card(
      "Test B",
      textInput("name_b", "Test Name:", "Test B"),
      selectInput("timepoint_b", "Time Point:", 
                  choices = c("Baseline", "Follow-up 1", "Not used"),
                  selected = "Not used"),
      numericInput("sens_b", "Sensitivity:", 0.7, min = 0, max = 1, step = 0.01),
      numericInput("spec_b", "Specificity:", 0.85, min = 0, max = 1, step = 0.01)
    )
  ),
  
  # Add the tags$style() here for custom styles
  tags$style(HTML("
    .container {
      width: 95% !important; /* Expands the app width */
    }
    .card {
      margin: 15px; /* Adds space around cards */
    }
  ")),
  
  # Dynamic testing section
  uiOutput("strategy_plots"),
  
  theme = bs_theme(version = 5, bootswatch = "minty"), # Optional theme customization
  fluid = TRUE
)

server <- function(input, output) {
  
  # Reactive expression to determine testing strategy
  testing_strategy <- reactive({
    timepoints <- c(input$timepoint_a, input$timepoint_b)
    
    if ("Not used" %in% timepoints) {
      if (input$timepoint_a == "Not used") return("Test B only")
      if (input$timepoint_b == "Not used") return("Test A only")
    }
    
    if (input$timepoint_a == input$timepoint_b) {
      return("Parallel")
    }
    
    # If we get here, tests are at different timepoints
    return("Serial")
  })
  
  # Reactive expression for metrics calculation
  metrics <- reactive({
    calculate_metrics(
      input$sens_a, input$spec_a,
      input$sens_b, input$spec_b,
      input$prevalence,
      input$sample_size
    )
  })
  
  # Dynamic UI based on strategy selection
  output$strategy_plots <- renderUI({
    strategy <- testing_strategy()
    
    if (strategy %in% c("Test A only", "Test B only")) {
      test_name <- if(strategy == "Test A only") input$name_a else input$name_b
      timepoint <- if(strategy == "Test A only") input$timepoint_a else input$timepoint_b
      
      card(
        card_header(paste(test_name, "Performance at", timepoint)),
        layout_columns(
          col_widths = 12,
          card(
            full_screen = TRUE,
            plotOutput("single_plot"),
            reactableOutput("single_summary")
          )
        )
      )
    } else {
      strategy_type <- if(strategy == "Parallel") "Parallel" else "Serial"
      card(
        card_header(paste(strategy_type, "Testing Strategies")),
        layout_columns(
          col_widths = c(6, 6),
          card(
            full_screen = TRUE,
            plotOutput("and_plot"),
            reactableOutput("and_summary")
          ),
          card(
            full_screen = TRUE,
            plotOutput("or_plot"),
            reactableOutput("or_summary")
          )
        )
      )
    }
  })
  
  # Generate single test plot and summary
  output$single_plot <- renderPlot({
    req(testing_strategy() %in% c("Test A only", "Test B only"))
    
    if (testing_strategy() == "Test A only") {
      create_mosaic_plot(metrics()$test_a, 
                         paste(input$name_a, "Performance at", input$timepoint_a))
    } else {
      create_mosaic_plot(metrics()$test_b, 
                         paste(input$name_b, "Performance at", input$timepoint_b))
    }
  })
  
  output$single_summary <- renderReactable({
    req(testing_strategy() %in% c("Test A only", "Test B only"))
    
    if (testing_strategy() == "Test A only") {
      create_summary_table(metrics()$test_a)
    } else {
      create_summary_table(metrics()$test_b)
    }
  })
  
  # Generate AND plot and summary
  output$and_plot <- renderPlot({
    req(testing_strategy() %in% c("Parallel", "Serial"))
    
    if (testing_strategy() == "Parallel") {
      create_mosaic_plot(
        metrics()$p_and,
        sprintf("Parallel Testing at %s: AND Strategy\n(%s AND %s must be positive)", 
                input$timepoint_a, input$name_a, input$name_b)
      )
    } else {
      create_mosaic_plot(
        metrics()$s_and,
        sprintf("Serial Testing (%s then %s): AND Strategy\n(%s AND %s must be positive)", 
                input$timepoint_a, input$timepoint_b, input$name_a, input$name_b)
      )
    }
  })
  
  output$and_summary <- renderReactable({
    req(testing_strategy() %in% c("Parallel", "Serial"))
    
    if (testing_strategy() == "Parallel") {
      create_summary_table(metrics()$p_and)
    } else {
      create_summary_table(metrics()$s_and)
    }
  })
  
  # Generate OR plot and summary
  output$or_plot <- renderPlot({
    req(testing_strategy() %in% c("Parallel", "Serial"))
    
    if (testing_strategy() == "Parallel") {
      create_mosaic_plot(
        metrics()$p_or,
        sprintf("Parallel Testing at %s: OR Strategy\n(Either %s OR %s positive)", 
                input$timepoint_a, input$name_a, input$name_b)
      )
    } else {
      create_mosaic_plot(
        metrics()$s_or,
        sprintf("Serial Testing (%s then %s): OR Strategy\n(Either %s OR %s positive)", 
                input$timepoint_a, input$timepoint_b, input$name_a, input$name_b)
      )
    }
  })
  
  output$or_summary <- renderReactable({
    req(testing_strategy() %in% c("Parallel", "Serial"))
    
    if (testing_strategy() == "Parallel") {
      create_summary_table(metrics()$p_or)
    } else {
      create_summary_table(metrics()$s_or)
    }
  })
}

shinyApp(ui, server)