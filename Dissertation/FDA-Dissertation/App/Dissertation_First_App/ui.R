#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("Functional Data Analysis Walk-Through"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            numericInput("Constant",
                        "Constant Term:",
                        min = -50,
                        max = 50,
                        value = 2),
            numericInput("x_coeff",
                         "Coefficient of X Term:",
                         min = -50,
                         max = 50,
                         value = 2),
            numericInput("x2_coeff",
                         "Coefficient of $$X^2$$ Term:",
                         min = -50,
                         max = 50,
                         value = 2)
        ),

        # Show a plot of the generated distribution
        mainPanel(

        )
    )
)
