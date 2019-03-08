  MASTER.list <- list(bayes_jeffreys, bayes_uniform, coit,
                      chao_huwang, easterling, lindstrom_madden, lindstrom_madden_AC,
                      madansky, mann_grubbs, normal_approximation,
                      myhre_rennie1, myhre_rennie2, nishime,
                      rice_moore)
  names(MASTER.list) <- c("Bayes (Jeffrey's prior)", "Bayes (Uniform prior)", "Coit",
                          "Chao-Huwang", "Easterling", "Lindstrom-Madden (Clopper-Pearson)", "Lindstrom-Madden (Agresti-Coull)",
                          "Madansky", "Mann-Grubbs", "MLE Normal Approximation",
                          "Myhre-Rennie (modified ML)", "Myhre-Rennie (reliability invariant)", "Nishime",
                          "Rice-Moore")


  MonteCarlo <- getShinyOption("MonteCarlo", default=1000)
  backup.method <- getShinyOption("backup.method", default=lindstrom_madden_AC)
  use.backup <- getShinyOption("use.backup", default=TRUE)


  ## input is results from main() function
  ## output is all the summary stats and histograms of each method; which method is "good"
  summary.LCB <- function(results, alpha, delta, mode){
    p.mat <- results$"p.mat"
    if(is.matrix(p.mat)==F) p.mat <- t(p.mat)
    big.list <- results$"big.list"
    num.p <- nrow(p.mat)
    M <- ncol(big.list[[1]])
    M.names <- colnames(big.list[[1]])
    p.names <- NULL

    LCB.avg.mat <- array(NA, dim=c(num.p, M), dimnames=list(p.names, M.names))
    cov0.mat <- array(NA, dim=c(num.p, M), dimnames=list(p.names, M.names))
    cov.minus.mat <- array(NA, dim=c(num.p, M), dimnames=list(p.names, M.names))
    cov.plus.mat <- array(NA, dim=c(num.p, M), dimnames=list(p.names, M.names))
    rmse.mat <- array(NA, dim=c(num.p, M), dimnames=list(p.names, M.names))

    delta.fix <- delta
    delta.hold <- 0 #place-holding delta, because may increase delta
    who.made.it <- NULL
    ranked.methods <- NULL #keep track of methods as delta increases from 0
    ranked.deltas <- NULL #corresponding delta values for the ranked methods
    out.list <- vector("list", length=M)
    out.list.ind <- 1
    out <- list(NULL)
    prepare.to.break <- FALSE


    ## this newer "best" is going to be the RANKED version.
    ## The user's delta value is the largest tolerable delta value to pass a method
    ## Any methods that pass for this largest value will be ranked according to their
    ## best delta value that they pass.
    ## In this function, just rank all selected methods. In the "renderprint" function later,
    ## that's when we can cut-off all the methods that are above the input delta
    if(mode=="best"){
      # calculate all LCB avg's and delta coverages

      while(is.null(out[[1]])){
        if( delta.hold>=delta ){
          delta.fix <- 1
          prepare.to.break <- TRUE
        }

        while(delta.hold<delta.fix){
          delta.hold <- delta.hold+.0001
          for(j in 1:num.p){
            Rs <- prod(p.mat[j,])
            LCB.mat <- big.list[[j]]
            LCB.avg.mat[j,] <- colMeans(LCB.mat)
            cov0.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs)) #regular coverage
            cov.minus.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs-delta.hold)) #minus delta coverage
            cov.plus.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs+delta.hold)) #plus delta coverage
          }

          #method is "good" if satisifies both coverage conditions across all p-vectors
          first <- colSums((cov0.mat>=(1-alpha)) & (cov.minus.mat<(1-alpha)))
          second <- colSums((cov0.mat<(1-alpha)) & (cov.plus.mat>=(1-alpha)))
          made.it <- (first+second)>=num.p
          sum.made.it <- sum(made.it)

          if(sum.made.it>0){
            for(s in which(made.it)){
              this.method <- M.names[s]
              if(!(this.method %in% ranked.methods)){
                summary <- rbind(colMeans(LCB.avg.mat), colMeans(cov0.mat),
                                 colMeans(cov.minus.mat), colMeans(cov.plus.mat))
                summary <- rbind(summary, rep(delta.hold, ncol(summary)))
                rownames(summary) <- c("Average LCB", "Coverage", "Delta Coverage (minus)",
                                       "Delta Coverage (plus)", "Best Delta Value")
                #out.list[[out.list.ind]] <- cbind(summary[c(1,2,5), s, drop=FALSE])
                out.list[[out.list.ind]] <- cbind(summary[c(5), s, drop=FALSE])
                out.list.ind <- out.list.ind + 1

                ranked.methods <- c(ranked.methods, this.method)
                ranked.deltas <- c(ranked.deltas, delta.hold)

              }
            }
          }
        if(prepare.to.break){
          if(!is.null(out.list[[1]])) break
        }

        }

        out <- out.list

      }

    }



    ## the alternative mode to rank the methods
    ## to pass, the method must have coverage at least 1-alpha for all R grid points
    ## then calculate the RMSE and average across all R grid points
    ## the method with the lowest average RMSE is ranked better
    if(mode=="rmse"){
      for(j in 1:num.p){
        Rs <- prod(p.mat[j,])
        LCB.mat <- big.list[[j]]
        LCB.avg.mat[j,] <- colMeans(LCB.mat)
        cov0.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs)) #regular coverage
        rmse.mat[j,] <- apply(X=LCB.mat, MARGIN=2, FUN=rmse.LCB, R=Rs)
      }

      avg.rmse <- colMeans(rmse.mat) #average the RMSE's across the R grid points


      # check which methods have coverage greater than or equal 1-alpha across all R grid points
      # this is a vector of T and F corresponding to methods which made it
      isit <- apply(X=cov0.mat, MARGIN=2, FUN=">=", 1-alpha)
      if(is.null(dim(isit))) isit <- t(isit)
      made.it <- colSums(isit) == num.p

      names.made.it <- M.names[made.it]
      rmse.made.it <- avg.rmse[made.it]
      to.display <- names.made.it[order(rmse.made.it)]
      sorted.rmse.made.it <- sort(rmse.made.it)


      listy <- vector("list", length=length(made.it))
      for(j in 1:sum(made.it)){
        flannel <- as.matrix(sorted.rmse.made.it[j])
        rownames(flannel) <- "Average RMSE"
        colnames(flannel) <- to.display[j]
        listy[[j]] <- flannel
      }

      out <- listy
      #out <- list(listy, sort(rmse.made.it))
    }


    ### Ranking the methods the same as Delta Coverage but with only criterion (3)
    ## that is, only consider methods that standard coverage at least 1-alpha
    ## then rank those methods according to the delta criterion
    if(mode=="exact"){
      # calculate all LCB avg's and delta coverages

      while(is.null(out[[1]])){
        if( delta.hold>=delta ){
          delta.fix <- 1
          prepare.to.break <- TRUE
        }

        while(delta.hold<delta.fix){
          delta.hold <- delta.hold+.0001
          for(j in 1:num.p){
            Rs <- prod(p.mat[j,])
            LCB.mat <- big.list[[j]]
            LCB.avg.mat[j,] <- colMeans(LCB.mat)
            cov0.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs)) #regular coverage
            cov.minus.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs-delta.hold)) #minus delta coverage
            #cov.plus.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs+delta.hold)) #plus delta coverage
          }

          #method is "good" if satisifies both coverage conditions across all p-vectors
          first <- colSums((cov0.mat>=(1-alpha)) & (cov.minus.mat<(1-alpha)))
          #second <- colSums((cov0.mat<(1-alpha)) & (cov.plus.mat>=(1-alpha)))
          made.it <- (first)>=num.p
          sum.made.it <- sum(made.it)

          if(sum.made.it>0){
            for(s in which(made.it)){
              this.method <- M.names[s]
              if(!(this.method %in% ranked.methods)){
                summary <- rbind(colMeans(LCB.avg.mat), colMeans(cov0.mat),
                                 colMeans(cov.minus.mat))
                summary <- rbind(summary, rep(delta.hold, ncol(summary)))
                rownames(summary) <- c("Average LCB", "Coverage", "Delta Coverage (minus)",
                                       "Best Delta Value")
                #out.list[[out.list.ind]] <- cbind(summary[c(1,2,5), s, drop=FALSE])
                out.list[[out.list.ind]] <- cbind(summary[c(4), s, drop=FALSE])
                out.list.ind <- out.list.ind + 1

                ranked.methods <- c(ranked.methods, this.method)
                ranked.deltas <- c(ranked.deltas, delta.hold)

              }
            }
          }
          if(prepare.to.break){
            if(!is.null(out.list[[1]])) break
          }

        }

        out <- out.list

      }

    }



    ## detailed output of all summary statistics;
    ## delta coverage is calculated based on input delta value
    if(mode=="detailed"){
      for(j in 1:num.p){
        Rs <- prod(p.mat[j,])
        LCB.mat <- big.list[[j]]
        LCB.avg.mat[j,] <- colMeans(LCB.mat)
        cov0.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs)) #regular coverage
        cov.minus.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs-delta)) #minus delta coverage
        cov.plus.mat[j,] <- colMeans(apply(X=LCB.mat, MARGIN=2, FUN="<", Rs+delta)) #plus delta coverage
        #rmse.mat[j,] <- apply(X=LCB.mat, MARGIN=2, FUN=rmse.LCB)
      }

      listy <- vector("list", length=num.p)
      listy.pnames <- NULL
      for(j in 1:num.p){
        foobar <- rbind(LCB.avg.mat[j,], cov0.mat[j,],
                        cov.minus.mat[j,], cov.plus.mat[j,],
                        rep(delta, ncol(LCB.mat)) )
        rownames(foobar) <- c("Average LCB", "Coverage", "Delta Coverage (minus)",
                              "Delta Coverage (plus)", "Delta Value")
        listy[[j]] <- foobar
        listy.pnames[j] <- paste(" p = (", paste(round(p.mat[j,], 3), collapse=", "), ") ", sep="")
      }
      names(listy) <- listy.pnames
      out <- listy
    }

    return(out)
  }




  ui <- navbarPage(title="serieslcb",
                   tabPanel(title="Simulations", fluidPage(

                     fluidRow(
                       column(width=5,
                              textInput(inputId = "Rs.int",
                                        label = "System Reliability Interval (lower, upper)",
                                        value="0.90, 0.95"),

                              textInput(inputId = "n",
                                        label = "Component Sample Sizes (comma delimited)",
                                        value = "35, 100, 60"),

                              textInput(inputId = "nsim",
                                        label = "Number of Simulations",
                                        value="1000"),

                              #textInput(inputId = "montecarlo",
                              #          label = "Monte Carlo Number of Samples",
                              #          value="1000"),

                              textInput(inputId = "seed",
                                        label = "Seed Number",
                                        value = "17"),

                              textInput(inputId = "confidence",
                                        label = "Confidence Level",
                                        value = "0.90"),

                              textInput(inputId = "delta",
                                        label = "Delta Level",
                                        value = "0.01")
                       ),

                       column(width=5,
                              checkboxGroupInput(inputId = "checkbox",
                                                 label = "Methods",
                                                 choices= c("Bayes (Jeffrey's prior)" = "Bayes (Jeffrey's prior)",
                                                            "Bayes (Uniform prior)" = "Bayes (Uniform prior)",
                                                            "Coit" = "Coit",
                                                            "Chao-Huwang" = "Chao-Huwang",
                                                            "Easterling" = "Easterling",
                                                            "Lindstrom-Madden (Clopper-Pearson)" = "Lindstrom-Madden (Clopper-Pearson)",
                                                            "Lindstrom-Madden (Agresti-Coull)" = "Lindstrom-Madden (Agresti-Coull)",
                                                            "Madansky" = "Madansky",
                                                            "Mann-Grubbs" = "Mann-Grubbs",
                                                            "MLE Normal Approximation" = "MLE Normal Approximation",
                                                            "Nishime" = "Nishime",
                                                            "Myhre-Rennie (modified ML)" = "Myhre-Rennie (modified ML)",
                                                            "Myhre-Rennie (reliability invariant)" = "Myhre-Rennie (reliability invariant)",
                                                            "Rice-Moore" = "Rice-Moore"))
                       )
                     ),

                     img(src='deltacoveragedefinitions.png', align = "right", width=450),


                     radioButtons(inputId = "radio", label = "Display Options",
                                  choices = c("Best Delta Coverage Methods (criteria (*) and (**))" = "best",
                                              "Best Delta Coverage Methods (criterion (*) only)" = "exact",
                                              "Detailed Output" = "detailed") ),

                     actionButton(inputId = "submit",
                                  label = "Submit"),

                     downloadButton(outputId = "histogramplot", label = "Download Histograms"),

                     verbatimTextOutput(outputId = "submitted.text"),
                     verbatimTextOutput(outputId = "summary.text"),
                     verbatimTextOutput(outputId = "delta.message"),

                     verbatimTextOutput(outputId = "test") #to test and make sure it's working

                     #mainPanel(
                     #   bsModal(id = "histogrammodal", title = "LCB Histograms",
                     #           trigger = "histogrambutton", size = "large",
                     #           plotOutput(outputId = "histogramplot"),
                     #           downloadButton(outputId = "downloadplot", label="Download"))
                     # )

                   )),
                   tabPanel(title="LCB Calculator", fluidPage(

                     fluidRow(
                       column(width=5,

                              textInput(inputId = "s.calc",
                                        label = "Component Number of Successes (comma delimited)",
                                        value = "35, 97, 59"),
                              textInput(inputId = "n.calc",
                                        label = "Component Sample Sizes (comma delimited)",
                                        value = "35, 100, 60"),
                              #textInput(inputId = "montecarlo.calc",
                              #          label = "Monte Carlo Number of Samples",
                              #          value="1000"),
                              textInput(inputId = "confidence.calc",
                                        label = "Confidence Level",
                                        value = "0.90")


                       ),

                       column(width=5,
                              checkboxGroupInput(inputId = "checkbox.calc",
                                                 label = "Methods",
                                                 choices= c("Bayes (Jeffrey's prior)" = "Bayes (Jeffrey's prior)",
                                                            "Bayes (Uniform prior)" = "Bayes (Uniform prior)",
                                                            "Coit" = "Coit",
                                                            "Chao-Huwang" = "Chao-Huwang",
                                                            "Easterling" = "Easterling",
                                                            "Lindstrom-Madden (Clopper-Pearson)" = "Lindstrom-Madden (Clopper-Pearson)",
                                                            "Lindstrom-Madden (Agresti-Coull)" = "Lindstrom-Madden (Agresti-Coull)",
                                                            "Madansky" = "Madansky",
                                                            "Mann-Grubbs" = "Mann-Grubbs",
                                                            "MLE Normal Approximation" = "MLE Normal Approximation",
                                                            "Nishime" = "Nishime",
                                                            "Myhre-Rennie (modified ML)" = "Myhre-Rennie (modified ML)",
                                                            "Myhre-Rennie (reliability invariant)" = "Myhre-Rennie (reliability invariant)",
                                                            "Rice-Moore" = "Rice-Moore"))
                       )
                     ),
                     actionButton(inputId = "submit.calc",
                                  label = "Submit"),

                     verbatimTextOutput(outputId = "submitted.calc.text"),

                     verbatimTextOutput(outputId = "calc.text")

                   ))
  )


  server <- function(input, output) {

    button <- reactiveValues(submitted = 0, submitted.calc = 0)

    user <- reactiveValues(Rs.int=NULL,
                           n=NULL,
                           nsim=NULL,
                           #montecarlo=NULL,
                           seed=NULL,
                           alpha=NULL,
                           delta=NULL,
                           checkbox=NULL,
                           state=NULL,
                           radio=NULL)

    hold <- reactiveValues(name="abcdefg",
                           my.list=NULL,
                           LCB=NULL)

    ## priority = 3 (highest priority)
    observeEvent(eventExpr = input$"submit", priority = 3, handlerExpr = {
      user$"Rs.int" <- as.numeric(unlist(strsplit(input$"Rs.int", ",")))
      user$"n" <- as.numeric(unlist(strsplit(input$"n", ",")))
      user$"nsim" <- as.numeric(unlist(input$"nsim"))
     # user$"montecarlo" <- as.numeric(unlist(input$"montecarlo"))
      user$"seed" <- as.numeric(unlist(input$"seed"))
      user$"alpha" <- 1 - as.numeric(unlist(input$"confidence"))
      user$"delta" <- as.numeric(unlist(input$"delta"))
      user$"checkbox" <- names(MASTER.list) %in% input$"checkbox" #make vector of T and F
      user$"radio" <- input$"radio"

      state.vec <- NULL
      # generated state = 1 if good, = 0 if invalid input(s)
      Rs.low <- min(user$"Rs.int")
      Rs.up <- max(user$"Rs.int")
      n <- user$"n"
      nsim <- user$"nsim"
      alpha <- user$"alpha"
      #delta <- user$"delta" #don't think we technically need a check for a "valid" delta
      out <- 1 #this means the state is "good"
      if( (Rs.low<0) | (Rs.up>1) ) out <- 0
      if( (sum(n<=0) > 0) ) out <- 0
      if(nsim<1) out <- 0
      if( (alpha<0) | (alpha>1) ) out <- 0
      state.vec[1] <- out

      # states$"submitted" = 1 if one or more methods are selected; = 0 otherwise
      state.vec[2] <- ifelse(sum(user$"checkbox") > 0, yes = 1, no = 0)

      user$"state" <- state.vec
    })

    #priority = 2 (middle priority, calculated after user is updated)
    observeEvent(eventExpr = input$"submit", priority = 2, handlerExpr = {
      req( prod(user$state)==1 ) #require valid inputs and at least one checked box

      thing <- vector(mode="list", length=5)
      thing[[1]] <- user$"Rs.int"
      thing[[2]] <- user$"n"
      thing[[3]] <- user$"nsim"
      thing[[4]] <- user$"seed"
      thing[[5]] <- user$"alpha"
      new.name <- paste(unlist(thing), collapse="/")

      #if need to refresh so all LCB's are re-calculated; e.g. data sims change, alpha changes...
      if(new.name != hold$"name"){
        #hold$"checkbox" <- rep(0, length(MASTER.list))
        hold$"name" <- new.name

        set.seed(seed = user$"seed")
        p.mat <- pm(Rs.int=user$"Rs.int", m=length(user$"n"))
        out <- vector(mode="list", length=nrow(p.mat))
        for(i in 1:nrow(p.mat)){
          out[[i]] <- matrix(rbinom(n=length(user$"n")*user$"nsim", size=rep(user$"n", user$"nsim"),
                                    prob=rep(as.numeric(p.mat[i,]), user$"nsim")),
                             byrow=TRUE, nrow=user$"nsim")
        }
        hold$"my.list" <- list(p.mat=p.mat, data.list=out)

        mat <- matrix(rep(NA, user$"nsim"*length(MASTER.list)), nrow=user$"nsim")
        colnames(mat) <- names(MASTER.list)
        out <- vector(mode = "list", length = nrow(p.mat))
        for(i in 1:length(out)){
          out[[i]] <- mat
        }
        hold$"LCB" <- out
      }

      #hold$"checkbox" <- hold$"checkbox" + user$"checkbox"
    })

    #priority = 1 (lowest priority)
    observeEvent(eventExpr = input$"submit", priority = 1, handlerExpr = {
      button$"submitted" = input$"submit"
    })

    # main function; calculates LCB's and the summary
    # invalidates every time "submit" is pressed, even if no inputs change
    # note reactives() are executed AFTER observes()
    handle <- eventReactive(eventExpr = input$"submit", valueExpr = {
      req( prod(user$"state")==1 )

      my.list <- hold$"my.list"
      data.list <- my.list$"data.list"
      num.p <- length(data.list) #should be the number of rows of p.mat

      #todo <- hold$"checkbox"==1 #if =0, don't calculate, if >1, already calculated it.
      todo <- user$"checkbox" & is.na(colSums(hold$"LCB"[[1]]))

      prog <- 0
      total <- sum(todo) * user$"nsim" * num.p
      withProgress(message = "Computing...", value = 0, {
        count <- 0
        for(z in 1:num.p){
          for(m in which(todo)){
            foo <- NULL
            for(b in 1:nrow(data.list[[z]])){
              count <- count + 1
              foo[b] <- MASTER.list[[m]](s=data.list[[z]][b,], n=user$"n", alpha=user$"alpha", MonteCarlo=MonteCarlo, use.backup=use.backup, backup.method=backup.method)
              incProgress(1/total, detail = paste(100*round(count/total, 2), "%", sep=""))
            }
            hold$"LCB"[[z]][,m] <- foo
          }
        }
      })

      LCB.list <- vector(mode="list", length=num.p)
      for(i in 1:num.p){
        LCB.list[[i]] <- hold$"LCB"[[i]][, user$"checkbox", drop=FALSE]
      }
      results <- list(p.mat=my.list$"p.mat", big.list=LCB.list)

      summary <- summary.LCB(results=results, alpha=user$"alpha",
                             delta=user$"delta", mode=user$"radio")

      out <- list(summary=summary, results=results)
      return(out)
    })


    ## begin subsection for the LCB Calculator tab...
    user.calc <- reactiveValues(s=NULL,
                                n=NULL,
                                alpha=NULL,
                                #montecarlo=NULL,
                                checkbox=NULL,
                                state=NULL)

    ## priority = 3 (highest priority)
    observeEvent(eventExpr = input$"submit.calc", priority = 3, handlerExpr = {
      user.calc$"s" <- as.numeric(unlist(strsplit(input$"s.calc", ",")))
      user.calc$"n" <- as.numeric(unlist(strsplit(input$"n.calc", ",")))
      user.calc$"alpha" <- 1 - as.numeric(unlist(input$"confidence.calc"))
      #user.calc$"montecarlo" <- as.numeric(unlist(input$"montecarlo.calc"))
      user.calc$"checkbox" <- names(MASTER.list) %in% input$"checkbox.calc" #make vector of T and F

      state.vec <- NULL
      # generated state = 1 if good, = 0 if invalid input(s)
      s <- user.calc$"s"
      n <- user.calc$"n"
      alpha <- user.calc$"alpha"
      out <- 1 #this means the state is "good"
      if(length(s) != length(n) ) out <- 0
      if( sum(s>n) > 0 ) out <- 0
      if( (alpha<0) | (alpha>1) ) out <- 0
      state.vec[1] <- out

      # states$"submitted" = 1 if one or more methods are selected; = 0 otherwise
      state.vec[2] <- ifelse(sum(user.calc$"checkbox") > 0, yes = 1, no = 0)

      user.calc$"state" <- state.vec
    })

    #priority = 1 (lowest priority)
    observeEvent(eventExpr = input$"submit.calc", priority = 1, handlerExpr = {
      button$"submitted.calc" = input$"submit.calc"
    })

    ## Main function to calculate LCB's for input data
    calculator <- eventReactive(eventExpr = input$"submit.calc", valueExpr = {
      req( prod(user.calc$"state")==1 )

      these.names <- input$"checkbox.calc"

      out <- NULL
      for(j in 1:length(these.names)){
        ind <- which(these.names[j]==names(MASTER.list))
        out[j] <- MASTER.list[[ind]](s=user.calc$"s", n=user.calc$"n",
                                     alpha=user.calc$"alpha", MonteCarlo=MonteCarlo)
      }
      out <- t(as.matrix(out))
      colnames(out) <- these.names
      rownames(out) <- paste(100*(1-user.calc$alpha), "%", " Lower Confidence Bound", sep="")

      out
    })


    ## begin output section...
    output$"submitted.text" <- renderPrint({
      if(button$"submitted" == 0) print <- "No methods have been submitted."
      if(button$"submitted" > 0){
        A <- user$state[1]
        B <- user$state[2]
        if( A==0 ) print <- "Error: Invalid input(s)."
        if( (A==1) & (B==0) ) print <- "Error: No methods have been selected."
        if( (A==1) & (B==1) ) print <- "Submission successful."
      }
      cat(print)
    })


    ## the actual numerical summary output
    output$"summary.text" <- renderPrint({
      #req( handle() )
      summary <- handle()$"summary"
      delta.vec <- do.call(cbind, summary)[1,]

      if(is.null(delta.vec)){
        out <- "No methods passed."
        cat(out)
      }
      if(length(delta.vec)>0){
        good.ind <- !sapply(summary, is.null)
        out <- summary[good.ind]
        print(out, digits=4)
      }

    })

    ## any message about the numerical summary output
    output$"delta.message" <- renderPrint({
      #req( handle() )
      summary <- handle()$"summary"
      delta.vec <- do.call(cbind, summary)[1,]

      #if(length(delta.vec)==0) out <- "Consider a larger value of delta."
      #if(length(delta.vec)>0) out <- "Displaying all methods (ranked according to their best delta value) which passed for a delta smaller than the input delta value. Both criteria (*) and (**) were used."
      if(user$"radio" == "best"){
        out <- "Displaying all methods (ranked according to their best delta value) which passed for a delta smaller than the input delta value. Both criteria (*) and (**) were used."
        if(delta.vec[1]>user$"delta") out <- "Displaying the best method(s). Note that delta was increased by .0001 until at least one method passed. Both criteria (*) and (**) were used."
      }

      if(user$"radio" == "exact"){
        out <- "Displaying all methods (ranked according to their best delta value) which passed for a delta smaller than the input delta value. Only criterion (*) was used."
        if(delta.vec[1]>user$"delta") out <- "Displaying the best method(s). Note that delta was increased by .0001 until at least one method passed. Only criterion (*) was used."
      }

      if(user$"radio" == "detailed"){
        out <- "Results for each p-vector combination are displayed. Delta coveraged are calculated using the input delta value."
      }



      cat(out)
    })


    ## function to make prettier breakpoints in the histogram function
    #LCB.breaks <- function(){
    #
    #}


    ## function to calculate histograms based on handle()$"results"
    ## reacts to "histogram" button
    LCB.histogram <- function(){
      p.mat <- handle()$"results"$"p.mat"
      big.list <- handle()$"results"$"big.list" # length(big.list) = nrow(p.mat)
      num.methods <- ncol(big.list[[1]])
      method.names <- colnames(big.list[[1]])

      #create dataframe to plot a table legend for the p-vector combinations
      listy.pnames <- NULL
      for(z in 1:nrow(p.mat)){
        listy.pnames[z] <- paste("(", paste(round(p.mat[z,], 3), collapse=", "), ") ", sep="")
      }
      df <- data.frame(1:length(listy.pnames), listy.pnames, row.names=NULL)
      colnames(df) <- c("Combo ID", "p-vector")
      gplots::textplot(df, show.rownames = FALSE) # in "gplots" package...

      ## collect all histogram data and determine appropriate xlim and ylim
      ## this is to make all histograms looks "pretty" and comparable to one another
      x.min <- NULL
      x.max <- NULL
      y.max <- NULL
      for(j in 1:num.methods){
        for(i in 1:nrow(p.mat)){
          h <- hist(big.list[[i]][,j], plot=FALSE, breaks=seq(.34, 1, .02))
          h$counts <- h$counts / sum(h$counts)
          y.max <- max(y.max, max(h$counts))
          x.min <- min(c(x.min, big.list[[i]][,j]))
          x.max <- max(c(x.max, big.list[[i]][,j]))
        }
      }
      nice.seq <- seq(0, 1, .02)
      x.min <- max(c(0, nice.seq[max(which(nice.seq < x.min))] - .06))
      x.max <- min(c(1, nice.seq[min(which(nice.seq > x.max))] + .06))
      y.max <- y.max + .06

      for(j in 1:num.methods){
        counter <- 0
        for(i in 1:nrow(p.mat)){
          counter <- counter + 1
          R.true <- prod(p.mat[i,])
          ## calculate relative frequencies... not sure how to make more pretty in general though
          h <- hist(big.list[[i]][,j], plot=FALSE, breaks=seq(x.min, x.max, .02))
          h$counts <- h$counts / sum(h$counts)
          plot(h, freq=TRUE, ylab="Relative Frequency", xlab="LCB",
               xlim=c(x.min, x.max), ylim=c(0, y.max),
               main=paste(method.names[j], "\np-vector Combination", counter))
          abline(v=R.true, col="red")

          #hist(big.list[[i]][,j], xlab="LCB", main=paste(method.names[j], "\np-vector Combination", counter))
          #points(c(R.true, R.true), c(0, user$"nsim"), type='l', col="red", lwd=2)
        }
      }

    }

    output$"histogramplot" <- downloadHandler(
      filename = function(){
        paste("LCB Histograms ", Sys.time(), ".pdf", sep="")
      },
      content = function(file){
        pdf(file)
        LCB.histogram()
        dev.off()
      })


    output$"submitted.calc.text" <- renderPrint({
      if(button$"submitted.calc" == 0) print <- "No calculations have been submitted."
      if(button$"submitted.calc" > 0){
        A <- user.calc$state[1]
        B <- user.calc$state[2]
        if( A==0 ) print <- "Error: Invalid input(s)."
        if( (A==1) & (B==0) ) print <- "Error: No methods have been selected."
        if( (A==1) & (B==1) ) print <- "Submission successful."
      }
      cat(print)
    })

    output$"calc.text" <- renderPrint({
      print(calculator(), digits=4)
    })
  }


  shinyApp(ui = ui, server = server)


