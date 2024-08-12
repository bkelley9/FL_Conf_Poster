library(shiny)
library(shinycssloaders)
library(tidyverse)
library(viridis)
library(akima)
library(plotly)
library(gridlayout)
library(bslib)
library(reactable)
library(climaemet)
library(ggpubr)
library(bslib)

wqdata <- read.csv("wqapp.csv") %>%
  mutate(Date.Time = as.POSIXct(Date.Time),
         Date = as.Date(Date)) %>%
  filter(Date > "2014-07-05", !(Date > "2015-04-01" & Date < "2015-06-11")) %>%
  select(Date.Time, Date, Time, Depth, Temp, Chla, DO, DOsat, PC, SpCond, Turbidity)

metdata <- read.csv("metapp.csv") %>%
  mutate(Date.Time = as.POSIXct(Date.Time),
         Time = format(Date.Time, "%H:%M:%S"),
         Date = as.Date(Date)) %>%
  filter(Date > "2014-07-05", !(Date > "2015-04-01" & Date < "2015-06-11"), !(Date > "2023-04-01" & Date < "2023-05-03")) %>%
  select(Date.Time, Date, Time, AirTemp, Rhumid, Pressure, SolarRad, WindSpeed, WindDir)

ui <- grid_page(
  layout = c(
    "sidebar heatPlot heatPlot heatPlot",
    "table   table    metplot  windrose  ",
    "table   table    metplot  windrose  "
  ),
  row_sizes = c(
    "1.7fr",
    "0.78fr",
    "1fr"
  ),
  col_sizes = c(
    "325px",
    "0.44fr",
    "1.2fr",
    "1.36fr"
  ),
  gap_size = "1rem",
  grid_card(
    area = "sidebar",
    card_header("Settings"),
    card_body(
      # selectInput(
      #   inputId = "year_select",
      #   label = "Select Year",
      #   choices = 2014:2023
      # ),
      selectInput("heat_param", "WQ Parameter",
                  choices = c("Temperature (\u00B0C)" = "Temp",
                              "Chlorophyll (\u00B5g/L)" = "Chla",
                              "Dissolved Oxygen (mg/L)" = "DO",
                              "Dissolved Oxygen Saturation (%)" = "DOsat",
                              "Phycocyanin (RFU)" = "PC",
                              "Specific Conductivity (\u00B5S/cm)" = "SpCond",
                              "Turbidity (FNU)" = "Turbidity")),
      selectInput("met_param", "Met Parameter",
                  choices = c("Air Temperature (\u00B0C)" = "AirTemp",
                              "Relative Humidity (%)" = "Rhumid",
                              "Pressure (mbar)" = "Pressure",
                              "Solar Radiation (W/m2)" = "SolarRad",
                              "Wind Speed (m/s)" = "WindSpeed",
                              "Wind Direction (\u00B0 from North)" = "WindDir")),
      sliderInput("date_range", "Date Range",
                  min = as.Date("2023-05-03"),
                  max = as.Date("2023-11-01"),
                  value = c(as.Date("2023-05-03"), as.Date("2023-11-01")))
    )
  ),
  grid_card(
    area = "table",
    card_header("Table"),
    full_screen = TRUE,
    card_body(
      tabsetPanel(
        nav_panel(
          title = "WQ Data",
          withSpinner(reactableOutput(outputId = "wq_table", width = "100%"), type = 5)
        ),
        nav_panel(
          title = "Met Data",
          withSpinner(reactableOutput(outputId = "met_table", width = "100%"), type = 4)
        )
      )
    )
  ),
  grid_card(
    area = "windrose",
    full_screen = TRUE,
    card_header("Wind Rose"),
    card_body(
      withSpinner(plotOutput(outputId = "windrose_plot", width = "100%"))
    )
  ),
  grid_card(
    area = "heatPlot",
    card_header("Water Quality Plot"),
    card_body(
      card(
        full_screen = TRUE,
        card_body(
          withSpinner(plotOutput(outputId = "heatmap_plot")))
      )
    )
  ),
  grid_card(
    area = "metplot",
    full_screen = TRUE,
    card_header("Meterological Plot"),
    card_body(
      withSpinner(plotOutput(outputId = "met_plot", width = "100%"))
    )
  )
)

server <- function(input, output, session){
  
  wq_filtered <- reactive({
    data <- wqdata %>% filter(year(Date) == 2023)
    data
  })
  
  met_filtered <- reactive({
    data <- metdata %>% filter(year(Date) == 2023)
    data
  })
  
  output$dynamic_date <- renderUI({
    sliderInput("date_range", "Date Range",
                min = min(wq_filtered()$Date, na.rm = T),
                max = max(wq_filtered()$Date, na.rm = T),
                value = c(min(wq_filtered()$Date, na.rm = T), max(wq_filtered()$Date, na.rm = T)))
  })
  
  # observeEvent(input$year_select, {
  #   freezeReactiveValue(input, "date_range")
  #   updateSliderInput(inputId = "date_range", 
  #                     min = min(wq_filtered()$Date, na.rm = T),
  #                     max = max(wq_filtered()$Date, na.rm = T),
  #                     value = c(min(wq_filtered()$Date, na.rm = T), max(wq_filtered()$Date, na.rm = T)))
  # })
  output$heatmap_plot <- renderPlot({
    browser()
    variable <- input$heat_param
    
    data1 <- wq_filtered() %>%
      filter(between(Date, input$date_range[1], input$date_range[2])) %>%
      select(Date.Time, Date, Time, Depth, !!sym(variable)) %>%
      mutate(depth_bin = cut(Depth, breaks = c(0, 1.25, 2.25, seq(3.75, 47.25, by = 1.5)),
                             labels = c(1, seq(1.5, 46.5, by = 1.5))) %>%
               as.character() %>%
               as.numeric(),
             Time = ifelse(Time < "11:00:00", "00:00:00", "12:00:00"),
             Date.Time = as.POSIXct(paste(Date, Time))
      ) %>%
      select(Date.Time, depth_bin, variable = variable) %>%
      filter(!is.na(variable), !is.na(depth_bin))
    
    if (variable == "SpCond"){
      data1 <- data1 %>%
        with(
          akima::interp(Date.Time, depth_bin, variable, duplicate = "mean",
                        xo = seq(min(Date.Time), max(Date.Time), by = "12 hours"),
                        yo = c(1, seq(1.5, 45, by = 1.5)),
                        linear = F,
                        remove = F))
    } else {
      data1 <- data1 %>%
        with(
          akima::interp(Date.Time, depth_bin, variable, duplicate = "mean",
                        xo = seq(min(Date.Time), max(Date.Time), by = "12 hours"),
                        yo = c(1, seq(1.5, 45, by = 1.5)),
                        linear = T,
                        remove = F))
    }
    
    base_plot <- data.frame(expand.grid(Date.Time = data1$x, Depth = data1$y), fill = c(data1$z)) %>%
      mutate(Date.Time = as.POSIXct(Date.Time) %>% format("%m-%d %H:%M:%S") %>% factor(),
             Depth = Depth %>% factor() %>% fct_rev()
      ) %>%
      filter(Depth != "46.5") %>% 
      ggplot(aes(x = Date.Time, y = Depth, fill = fill)) + 
      geom_tile() + 
      scale_x_discrete(breaks = c("04-01 00:00:00", "05-01 00:00:00","06-01 00:00:00", "07-01 00:00:00", "08-01 00:00:00", "09-01 00:00:00", "10-01 00:00:00"),
                       labels  = c("April", "May", "Jun", "July", "Aug", "Sept", "Oct")) +
      scale_y_discrete(limits = levels(data1$Depth),
                       breaks = c(1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45)) +
      scale_fill_viridis(option = "turbo") +
      theme(text = element_text(size = 20)) +
      labs(x = "Date",
           y = "Depth (m)",
           fill = input$heat_param)
    
    print(base_plot)
  })
  
  output$met_plot <- renderPlot({
    
    data <- met_filtered() %>%
      filter(between(Date, input$date_range[1], input$date_range[2]))
    
    variable <- input$met_param
    
    plot <- ggplot(data, aes(x = Date.Time, y = !!sym(variable))) +
          geom_line() +
          geom_smooth() +
          theme_pubclean() +
          theme(text = element_text(size = 20)) + 
          labs(y = variable,
               x = "Date")
    print(plot)
  })
  
  output$windrose_plot <- renderPlot({
    
    data <- met_filtered() %>%
      filter(between(Date, input$date_range[1], input$date_range[2]))
    
    plot <- ggwindrose(data$WindSpeed, data$WindDir,
                       speed_cuts = c(0, 2, 4, 6, 8, 10, 12, 14),
                       col_pal = "Viridis") +
      theme(text = element_text(size = 20))
    
    print(plot)
  })
  
  output$wq_table <- renderReactable({
    data <- wq_filtered() %>%
      filter(year(Date) == 2023)
    
    table <- reactable(data,
                       columns = list(
                         Date.Time = colDef(show = F),
                         Temp = colDef(name = "Temperature"),
                         Chla = colDef(name = "Chl-a", format = colFormat(digits = 3)),
                         DO = colDef(name = "Dissolved Oxygen"),
                         DOsat = colDef(name = "Dissolved Oxygen Saturation"),
                         SpCond = colDef(name = "Specific Conductivity"),
                         PC = colDef(name = "Phycocyanin", format = colFormat(digits = 3))
                         
                       ))
    
    return(table)
    
  })
  
  output$met_table <- renderReactable({
    data <- met_filtered() %>%
      filter(year(Date) == 2023)
    
    table <- reactable(data,
                       columns = list(
                         Date.Time = colDef(show = F),
                         AirTemp = colDef(name = "Air Temperature"),
                         Rhumid = colDef(name = "Relative Humidity"),
                         SolarRad = colDef(name = "Solar Radiation"),
                         WindSpeed = colDef(name = "Wind Speed"),
                         WindDir = colDef(name = "Wind Direction")
                         
                       ))
    
    return(table)
  })

}

shinyApp(ui, server)