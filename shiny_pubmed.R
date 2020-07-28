library(shiny)
library(shinythemes)
library(ggplot2)
library(RISmed)
library(stringr)
library(tokenizers)
library(purrr)
library(DT)
library(dplyr)
library(tidyr)

values <- reactiveValues()
values$fetch <- NULL
values$res <- NULL
values$Text <- NULL
values$D1 <- NULL
values$D2 <- NULL

app <- shinyUI(fluidPage(theme=shinytheme("united"),
                         headerPanel("PubMed Search"),
                         sidebarLayout(
                           sidebarPanel(
                             helpText("Type a word below and search PubMed to find documents that contain that word in the text. You can even type multiple words. You can search authors, topics, any acronym, etc. Please note, that if you do a broad search, only 1000 publications will be returned by the program"),
                             textInput("text", label = h3("Keyord(s)"), value = "Type some text"),
                             helpText("You can specify the start and end dates of your search, use the format YYYY/MM/DD"),
                             textInput("date1", label = h3("From"),value="1990/01/01"),
                             textInput("date2", label = h3("To"),  value = paste0(str_split(paste(Sys.Date()),"-")[[1]][1],"/",str_split(paste(Sys.Date()),"-")[[1]][2],"/",str_split(paste(Sys.Date()),"-")[[1]][3])),
                             helpText("If you want to return PMID's for abstracts containing specific text, insert your search parameters below"),
                             textInput("search_text", label = h3("Search string as regex expression. Set to 'NULL' to return all"), value = "NULL"),
                             checkboxInput("search_for_significance", "Return sentences indicating significance?", TRUE),
                             helpText("Now select the output you'd like to see. You can see a barplot of articles per year, a table of the top 20 authors or a table of top 20 journals"),
                             
                             
                             actionButton("goButton","PLOT"),
                             actionButton("authButton","AUTHORS"),
                             actionButton("authButton1","JOURNALS"),
                             actionButton("authButton2","SENTENCES")),
                           
                           mainPanel(
                             tabsetPanel(id = "inTabset",
                                         tabPanel("PLOT", 
                                                  br(),
                                                  textOutput("querry_plot"), 
                                                  br(),
                                                  plotOutput("distPlot")),
                                         tabPanel("AUTHORS", 
                                                  br(),
                                                  textOutput("querry_auth"),
                                                  br(),
                                                  DTOutput("authList")),
                                         tabPanel("JOURNALS", 
                                                  br(),
                                                  textOutput("querry_auth1"),
                                                  br(),
                                                  DTOutput("journalList")),
                                         tabPanel("SENTENCES", 
                                                  br(),
                                                  textOutput("querry_sentences"),
                                                  br(),
                                                  DTOutput("abstract_of_interest"))
                             )
                           )
                         )))

server <- shinyServer(function(input, output, session) {
  
  observeEvent(input$goButton, {
    updateTabsetPanel(session, "inTabset",
                      selected = "PLOT")
    output$distPlot <- renderPlot({
      withProgress(message = 'Making plot', value = 0, {
        
        if(is.null(isolate(values$fetch))) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        } else if (isolate(input$text) != values$Text | isolate(input$date1) != values$D1 | isolate(input$date2) != values$D2) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        }
        
        count<-table(YearPubmed(values$fetch))
        count<-as.data.frame(count)
        
        if(nrow(count) == 0) {
          showNotification("Please reduce the search. PubMed finds no articles with this combination of words", type = "error")
        } else {
          names(count)<-c("Year", "Counts")
          num <- data.frame(Year=count$Year, Counts=cumsum(count$Counts)) 
          num$g <- "g"
          names(num) <- c("Year", "Counts", "g")
          
          incProgress(1/3, detail = paste("Creating graph"))
          q <- ggplot(aes(x=Year, y=Counts), data=count) +
            geom_bar(stat="identity")
          q <- q + geom_line(aes(x=Year, y=Counts, group=g), data=num) +
            ggtitle(paste("PubMed articles containing '", isolate(input$text), "' ", "= ", max(num$Counts), sep="")) +
            ylab("Number of articles") +
            xlab(paste("Year n Query date: ", Sys.time(), sep="")) +
            labs(colour="") + theme_bw()
          q
        }
      })
    })
    output$querry_plot <- renderText({
      if(!is.null(values$res)) {
        paste0("Pudmed Search string: ",slot(values$res,"querytranslation"),"\n","Publication counts: ", QueryCount(values$res)) 
      } else {
        paste(NULL)
      }
    })
  })
  
  observeEvent(input$authButton, {
    updateTabsetPanel(session, "inTabset",
                      selected = "AUTHORS")
    output$authList<-renderDT({
      withProgress(message = 'Generating table', value = 0, {
        if(is.null(isolate(values$fetch))) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        } else if (isolate(input$text) != values$Text | isolate(input$date1) != values$D1 | isolate(input$date2) != values$D2) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        }
        
        AuthorList<-Author(values$fetch)
        Last<-sapply(AuthorList, function(x)paste(x$ForeName,x$LastName))
        auths<-as.data.frame(sort(table(unlist(Last)), dec=TRUE))
        auths <- auths[!(auths$Var1 == "NA NA"),]
        colnames(auths)<-"Count"
        incProgress(1/3, detail = paste("Generating table"))
        auths <- cbind(Author = rownames(auths), auths)
        rownames(auths) <- NULL
        auths<-head(auths, 20)
        auths
      })
    })
    output$querry_auth <- renderText({
      if(!is.null(values$res)) {
        paste0("Pudmed Search string: ",slot(values$res,"querytranslation")) 
      } else {
        paste(NULL)
      }
    })
  })
  
  observeEvent(input$authButton1, {
    updateTabsetPanel(session, "inTabset",
                      selected = "JOURNALS")
    output$journalList<-renderDT({
      withProgress(message = 'Generating table', value = 0, {
        
        if(is.null(isolate(values$fetch))) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        } else if (isolate(input$text) != values$Text | isolate(input$date1) != values$D1 | isolate(input$date2) != values$D2) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        }
        
        AuthorList<-Title(values$fetch)
        data <- data.frame(journal_names = names(table(AuthorList)),Count = table(AuthorList))
        incProgress(1/3, detail = paste("Creating table"))
        data <- data[order(data$Count.Freq,decreasing = TRUE),]
        auths<-head(data[,c(1,3)], 20)
        auths
      })
    })
    output$querry_auth1 <- renderText({
      if(!is.null(values$res)) {
        paste0("Pudmed Search string: ",slot(values$res,"querytranslation")) 
      } else {
        paste(NULL)
      }
    })
  })
  
  observeEvent(input$authButton2, {
    updateTabsetPanel(session, "inTabset",
                      selected = "SENTENCES")
    output$abstract_of_interest <- renderDT({
      withProgress(message = 'Getting abstracts', value = 0, {
        
        search_for_significant <- as.character(isolate(input$search_for_significance))
        extra_search_patterns <- as.character(isolate(input$search_text))
        
        if(is.null(isolate(values$fetch))) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        } else if (isolate(input$text) != values$Text | isolate(input$date1) != values$D1 | isolate(input$date2) != values$D2) {
          incProgress(1/3, detail = paste("Getting data"))
          values$Text <- isolate(input$text)
          values$D1 <- isolate(input$date1)
          values$D2 <- isolate(input$date2)
          values$res <- EUtilsSummary(isolate(input$text), type="esearch", db="pubmed", datetype='pdat', mindate=isolate(input$date1), maxdate=isolate(input$date2), retmax=1000)
          values$fetch <- EUtilsGet(values$res, type="efetch", db="pubmed")      
          incProgress(1/3, detail = paste("Preparing data"))
        }
        
        articles <- data.frame('Abstract'=AbstractText(values$fetch))
        abstracts <- paste((articles$Abstract))
        pubmedID <- paste(RISmed::ArticleId(values$fetch))
        name_publication <- paste(RISmed::ArticleTitle(values$fetch))
        
        #Tokenize the sentences
        list_sentences <- tokenize_sentences(abstracts)
        
        #Create tibble with data
        data <- mapply(cbind, "PMID" = pubmedID,"Abstract_sentences" = list_sentences, "Title" = name_publication, SIMPLIFY=F)
        data <- data %>% 
          map(~ as_tibble(.x)) %>% 
          bind_rows() %>% 
          drop_na()
        
        search_pattern_significant <- regex("([p]( |)(=|<)( |)\\d)|significant|( )tend")
        if(search_for_significant == "TRUE") data <- data[str_detect(data$Abstract_sentences, search_pattern_significant),]
        if(extra_search_patterns != "NULL") data <- data[str_detect(x, extra_search_patterns),]
        datatable(tibble(PMID = data$PMID,Title = data$Title,`Abstract sentences` = data$Abstract_sentences), options = list( pageLength = 100, columnDefs = list(list(targets = c(0,1,2), searchable = FALSE))))
      })
    })
    output$querry_sentences <- renderText({
      if(!is.null(values$res)) {
        paste0("Pudmed Search string: ",slot(values$res,"querytranslation")) 
      } else {
        paste(NULL)
      }
    })
  })
  
})

shinyApp(ui = app, server = server)

#___________#

# d1 <- "2010/01/01"
# word1 <- "remote virtual adherence retention dropout e-consent recruitment"
# d2 <- "2019/01/01"
