library(tidyverse)
library(shiny)
library(dplyr)

#Define the user interface
ui= fluidPage(
  
  titlePanel("Star Wars movie data"),
  
  selectInput("color", label="Skin color", choices= unique(starwars$skin_color)),
  
  submitButton("Get results"),

  p("Home world of characters with the chosen skin color"),
  
  verbatimTextOutput("tablehomeworld"),
  
  h2("Data"),
  dataTableOutput("selecteddata"),
  
  plotOutput("myplot")

)

#Define the server calculations
server= function(input, output, session) {

  #Avoid code duplication and repeating computations. 
  #Since several output components use the same data, define a reactive expression (mydata below) 
  #The result is a dataset but it's accessed like a function, i.e. mydata()
  mydata= reactive({
    filter(starwars, skin_color == input$color)
  })
  
  output$tablehomeworld= renderPrint({
    table(mydata()$homeworld)
  })
  
  output$selecteddata= renderDataTable({
    mydata()
  })
  
  output$myplot= renderPlot({
    mytitle= paste("Characters with skin color ", input$color, collapse="")
    ggplot(mydata(), aes(mass, height)) + geom_point() + labs(x="Mass",y="Height",title=mytitle)
  })
  
}

#run App
shinyApp(ui, server)