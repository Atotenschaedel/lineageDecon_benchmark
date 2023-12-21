## function for prosprediction plots
library(dplyr)
library(tidyr)
library(tidyr)
library(reshape2)
library(wesanderson)
library(stringi)

# first letter capitalized
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

getLookUP <- function(data.tool, sim.lineages, json.data, level) {
  json_data <- NULL
  json_data <- data.frame(lineage = sim.lineages, ancestor = sim.lineages, unaliased = NA)
  
  for (l in sim.lineages) {
    if (l != "others"){
      ancestor <- l
      descendant <- NULL
      
      curr_nodes <- l
      # get children (level 1 or grandchildren level 2 ect...
      for (i in 1:level){
        next_node <- NULL
        
        for (node in curr_nodes) {
          children <- json.data[[node]]$children
          descendant <- c(descendant, children)
          next_node <- c(next_node, children)
        }
        curr_nodes <- next_node
      }
    }
    
    if(length(descendant) == 0) {}
    else {
      temp <- data.frame(lineage = unlist(descendant, use.names=FALSE), ancestor=ancestor, unaliased=json.data[[ancestor]]$unaliased)
      json_data <- rbind(json_data, temp)
      
    }
  }
  
  json_clean <- json_data %>% group_by(lineage) %>% filter(n() == 1)
  json_unclear <- json_data %>% group_by(lineage) %>% filter(n() > 1)
  
  if (nrow(json_unclear) > 0)  {
    for (lin in unique(json_unclear$lineage)) {
      z <- filter(json_unclear, lineage==lin) %>%
        rowwise() %>% mutate(n = nchar(unaliased)) %>%
        mutate(n = coalesce(n, 0)) %>%
        unique() %>% ungroup() %>% 
        slice_min(n=1, n) %>% select(-c("n"))
      json_clean <- rbind(json_clean, z)
    }
  }

  json_data <- json_clean %>% select(-c("unaliased"))
  return (json_data)
}

densityPlot <- function(data, x.value, group, Colors) {
  ggplot(data, aes_string(color=group)) +
    geom_density(aes_string(x=x.value)) +
    scale_colour_manual(values=Colors)
}

tool_predictionScatterplot <- function(data, sim.data.long, x.value, mode){
  # change measures and legend based on variable  
  if (mode=="strict") {
    measures <- c("jaccard_index", "TPR", "FPR", "relAb_RMSE")
    label.legend <- c(FPR="False Positive Rate", jaccard_index="Jaccard Index", 
                TPR="True Positive Rate", relAb_RMSE="Relative Abundance - Mean Square Error")
    label.title <- "Prediction Error - Strict Mode"
  }
  
  else if (mode == "adjusted") {
    measures <- c("jaccard_index_adj_child", "TPR_adj_child", "FPR_adj_child", "relAb_RMSE_adj_child")
    label.legend <- c(FPR_adj_child="False Positive Rate", jaccard_index_adj_child="Jaccard Index",
                      TPR_adj_child="True Positive Rate", relAb_RMSE_adj_child="Relative Abundance - Mean Square Error")
    label.title <- "Prediction Error - Adjusted Mode"
  }
  
  # data filter
  data.filter <- data %>%
    group_by(tool_name, paste(x.value)) %>%
    gather(key="measure", value="value", measures)
  
  # plot 
  p <- ggplot(data.filter, aes_string(x=x.value, y="value", color="tool_name")) +
    geom_point()
  
  if (x.value == "timepoint") {
    label.title <- paste(label.title, " - per Timepoint")
    label.x <- "Timepoint"
    p <- p + geom_smooth(method="loess", se=TRUE, fullrange=TRUE, level=0.95, aes(fill=tool_name))
  }
  
  else if (x.value == "uniformity_wg_per"){
    label.title <- paste(label.title, " - genomic coverage")
    label.x <- "genomic coverage (in %)"
    p <- p + geom_smooth(method=lm, se=TRUE, fullrange=TRUE, level=0.95, aes(fill=tool_name))
  }
  
  
  # finish plot
  p + facet_wrap(~measure, labeller=as_labeller(label.legend)) + 
    scale_y_continuous(limits = c(0,1)) +
  #scale_y_continuous(trans='log10') +
  #geom_line(data= sim.data.long, aes_string(x=x.value, y="rel_abun100", color="lineage"), linetype="dotted") +
  #scale_y_continuous(limits = c(0,1), sec.axis = sec_axis(~.*100, name="relative Abundance (in %)")) +
    labs(
      title=label.title,
      y="",
      x=label.x) +
    scale_color_manual(values=myColors) +
    theme_minimal() + guides(fill=FALSE)
}

tool_predictionLineplot <- function(data, sim.data.long, mode) {
  scale=100
  
  if (mode == "strict") {
    d <- data %>% 
      group_by(tool_name, timepoint) %>%
      summarise(jaccard_index = mean(jaccard_index), RMSE = mean(relAb_RMSE), .groups = 'keep') %>%
      gather(key="measure", value = "value", c("jaccard_index", "RMSE"))
    label.title <- "Prediction Error - Identification - strict mode"
  } 
  
  else if (mode == "adjusted") {
    d <- data %>% 
      group_by(tool_name, timepoint) %>%
      summarise(jaccard_index = mean(jaccard_index_adj_child), RMSE = mean(relAb_RMSE_adj_child), .groups = 'keep')  %>%
      gather(key="measure", value = "value", c("jaccard_index", "RMSE"))
    label.title <- "Prediction Error - Identification - adjusted mode"
  }
  
  ggplot(sim.data.long, aes(x=timepoint, y=rel_abun, color=lineage)) +
    geom_line(linetype="dotted") +
    geom_point(data = d %>% filter(tool_name=="simulation", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "simulation")) +
    geom_point(data = d %>% filter(tool_name=="vaquero", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "vaquero")) +
    geom_point(data = d %>% filter(tool_name=="VLQ", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "VLQ")) +
    geom_point(data = d %>% filter(tool_name=="lollipop", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "lollipop")) +
    geom_point(data = d %>% filter(tool_name=="freyja", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "freyja")) +
    
    geom_line(data = d %>% filter(tool_name=="simulation", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "simulation")) +
    geom_line(data = d %>% filter(tool_name=="vaquero", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "vaquero")) +
    geom_line(data = d %>% filter(tool_name=="VLQ", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "VLQ")) +
    geom_line(data = d %>% filter(tool_name=="lollipop", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "lollipop")) +
    geom_line(data = d %>% filter(tool_name=="freyja", measure == "jaccard_index"), aes(x=timepoint, y=value*scale, colour = "freyja")) +
    
    geom_point(data = d %>% filter(tool_name=="simulation", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "simulation")) +
    geom_point(data = d %>% filter(tool_name=="vaquero", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "vaquero")) +
    geom_point(data = d %>% filter(tool_name=="VLQ", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "VLQ")) +
    geom_point(data = d %>% filter(tool_name=="lollipop", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "lollipop")) +
    geom_point(data = d %>% filter(tool_name=="freyja", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "freyja")) +
    
    geom_line(data = d %>% filter(tool_name=="simulation", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "simulation")) +
    geom_line(data = d %>% filter(tool_name=="vaquero", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "vaquero")) +
    geom_line(data = d %>% filter(tool_name=="VLQ", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "VLQ")) +
    geom_line(data = d %>% filter(tool_name=="lollipop", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "lollipop")) +
    geom_line(data = d %>% filter(tool_name=="freyja", measure == "RMSE"), aes(x=timepoint, y=value*scale, colour = "freyja")) +
    
    scale_color_manual(name="tool name", breaks=c("simulation", "vaquero", "VLQ", "lollipop", "freyja"), values=c("simulation"="grey" , "vaquero"="blue", "VLQ"="orange", "lollipop"="green", "freyja"="red")) +
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name="Jaccard Index")) +
    facet_wrap(~factor(measure, levels=c('RMSE', 'jaccard_index')), labeller=as_labeller(c(jaccard_index="Jaccard Index", RMSE="Relative Abundance - Mean Square Error"))) +
    labs(
      title=label.title,
      y="Relative Abundace (in %)",
      x="Timepoint") +
    theme_minimal() 
}

tool_twoDimensionPlot <- function(data.metrics, colScale, title) {
  data.filter <- data.metrics %>%
    group_by(tool_name) %>%
    summarise(jaccard_index = median(jaccard_index),
              relAb_RMSE = median(RMSE), .groups="drop")
  
  ggplot(data.filter, aes(x=jaccard_index, y=relAb_RMSE)) +
    geom_point(aes(colour=tool_name), size = 3, shape=4, stroke=2) +
    facet_grid(~time-course) + 
    labs(
      title=title,
      y="Root Mean Squared Error - relative Abundance",
      x="Jaccard Index") +
    colScale +
    theme_minimal()
  
}

tool_timepPlot <- function(data.long, sim.data.long, tool, colScale, mode) {
  caption <- "Errorbars are standard error over all experiments."
  
  data.tool <- data.long %>% filter(tool_name == tool) %>% select(-c("tool_name", "uniformity_wg_per")) %>% mutate(series = tool)
  data.sim <- sim.data.long %>% mutate(series = "simulation")
  
  # remove signals only observed in one timepoint
  sim.lineages <- filter(data.sim, series=="simulation") %>% 
    group_by(lineage) %>% 
    summarise(count = n_distinct(timepoint)) %>% 
    filter(count > 1) %>% pull(lineage)
  
  data.tool <- filter(data.tool, lineage %in% (data.tool %>% 
                                                 group_by(lineage) %>% 
                                                 summarise(count = n_distinct(timepoint)) %>% 
                                                 filter(count > 1) %>% pull(lineage)))
  
  zero_count <- length(filter(data.sim, series=="simulation") %>% 
                         group_by(lineage) %>% summarise(count = n_distinct(timepoint)) %>% 
                         filter(count == 1) %>% pull(lineage))
  if(zero_count != 0){
    caption <- paste(caption,
                     (paste0(zero_count," lineages only simulated in 1 Timepoint and removed from visualisation.")))
  } else{}
  
  if (length(union(unique(data.tool$lineage), sim.lineages)) > 36) {
    
    # get top 36 lineage by relative abundance
    top_lineages <- data.tool %>% 
      group_by(lineage) %>% 
      summarise(rel_abun = mean(rel_abun)) %>% 
      arrange(desc(rel_abun)) %>%
      slice_head(n=36) %>% 
      pull(lineage)
    
    caption <- paste(caption, 
                     paste0("Showing top ", length(top_lineages)," lineage of ", 
                            length(unique(data.tool$lineage)),
                            " lineages (top lineage according to mean of abundance over all timepoints and predicted at more then 1 timepoint.)"))
    
    data.tool <- data.tool %>% filter(lineage %in% top_lineages)
    
    # check for false negatives
    if (length(sim.lineages) - length(intersect(top_lineages, sim.lineages)) > 0){
      caption <- paste(caption, 
                       length(sim.lineages) - length(intersect(top_lineages, sim.lineages)), " false negative signals.")
      data.long <- rbind(data.tool, filter(data.sim, lineage %in%  intersect(top_lineages, sim.lineages)))
    }
    else {data.long <- rbind(data.tool, data.sim)}
    
  } else {
    data.long <- rbind(data.tool, data.sim)
  }
  
  
  
  d <- data.long %>% group_by(timepoint, lineage, series) %>% summarise_at(vars(rel_abun), list(sd=sd, mean=mean))
  d <- merge(data.long, d, by=c("timepoint", "lineage", "series"))
  
  d$error_min <- d$mean - d$sd
  d$error_max <- d$mean + d$sd
  
  d <- d %>% mutate("error_min" = case_when(error_min < 0 ~ 0, TRUE ~ error_min), 
                    "error_max" = case_when(error_max > 1 ~ 1, TRUE ~ error_max))
  
  transparency <- 0.3
  p <- ggplot(filter(d, series==tool), aes(x=timepoint, y=mean, color=series)) + 
    #geom_line(data=d %>% filter(series=="simulation", lineage %in% sim.lineages), aes(x=timepoint, y=mean)) +
    geom_point(data=d %>% filter(series=="simulation", lineage %in% sim.lineages), aes(x=timepoint, y=mean), size=0.8) +
    geom_smooth(data=d %>% filter(series=="simulation"), se=FALSE, size=0.8, alpha=transparency) +
    #geom_line() + 
    geom_point() +
    geom_errorbar(aes(ymin=error_min, ymax=error_max), width=.2) +
    facet_wrap(~lineage, ncol=6) +
    labs(
      title=paste0("Simulated Lineage Dynamic - ",tool, " - ", mode),
      y="Relative Abundance (in %)",
      x="Timepoint",
      caption = caption
    ) +
    colScale +
    theme_minimal()
  
  print(caption)
  
  return(p)
}

calculateConfusionMatrix <- function(data.long, sim.data.long, grouping=TRUE) {
  data.metrics <- NULL
  
  for (exp in unique(data.long$experiment)) {
    for (rep in unique(data.long$replicate)) {
      for (tp in unique(data.long$timepoint)) {
        for (tool in unique(data.long$tool_name)) {
          data.temp <- full_join(
            filter(data.long, experiment==exp, replicate==rep, tool_name==tool, timepoint==tp), 
            filter(sim.data.long, experiment==exp, replicate==rep, timepoint==tp), 
            by=c("timepoint", "lineage", "experiment", "replicate"), 
            suffix = c("","_sim")) %>%
            mutate("P" = case_when(!is.na(rel_abun_sim) ~ 1),
                   "PP" = case_when(!is.na(rel_abun) ~ 1),
                   "TP" = case_when(!is.na(rel_abun) & !is.na(rel_abun_sim) ~ 1),
                   "FP" = case_when(!is.na(rel_abun) & is.na(rel_abun_sim) ~ 1),
                   "FN" = case_when(is.na(rel_abun) & !is.na(rel_abun_sim) ~ 1)) %>%
            mutate_at(c("rel_abun","rel_abun_sim", "P", "PP", "TP", "FP", "FN"), ~replace_na(.,0)) %>%
            mutate("Error" = case_when(P == 1 ~ (rel_abun - rel_abun_sim)))
            
            if (grouping){
              data.temp <- data.temp %>%
                select(-c("rel_abun", "rel_abun_sim")) %>%
                group_by(across(all_of(c("timepoint", "experiment", "replicate")))) %>% summarise(
                  P = sum(P),
                  PP = sum(PP),
                  TP = sum(TP),
                  FP = sum(FP),
                  FN = sum(FN),
                  RMSE = sqrt(mean(Error^2, na.rm=TRUE)), .groups = "drop") %>%
                mutate(tool_name = tool)
            }
            else {}
          data.metrics <- rbind(data.metrics, data.temp)
        }
      }
    }
  }
  return(data.metrics)
}

getAdjustedData <- function(data.long, sim.lineages, json.data, level) {
  
  json_data <- getLookUP(data.long, sim.lineages, json.data, level)

  d <- data.long %>%
    left_join(., y=json_data, by="lineage", multiple = "all") %>%
    mutate(lineage = ifelse(!is.na(ancestor), ancestor, lineage)) %>% 
    select(-c("ancestor")) %>% group_by(timepoint, experiment, replicate, tool_name, lineage) %>%
    summarise(rel_abun=sum(rel_abun),
              uniformity_wg_per=mean(uniformity_wg_per), .groups="drop")

  return(d)
}

calculateMetrics <- function(data.metrics){
  data.metrics$TPR <- data.metrics$TP / (data.metrics$TP + data.metrics$FN)
  data.metrics$FNR <- data.metrics$FN / data.metrics$P
  data.metrics$PPV <- data.metrics$TP / data.metrics$PP #precision
  data.metrics$FDR <- data.metrics$FP / data.metrics$PP
  data.metrics$jaccard_index <- data.metrics$TP / (data.metrics$TP + data.metrics$FP + data.metrics$FN)
  data.metrics$F1_score <- (data.metrics$TP*2) / (data.metrics$TP*2 + data.metrics$FP + data.metrics$FN)
  return(data.metrics)
}

vennDiagramReplicates <- function(data.long) {
  names <- NULL
  for (i in (1:length(unique(data.long$replicate)))) {
    names <- c(names, paste("Replicate", unique(data.long$replicate)[i]))
  }
  rp <- data.long %>% select(c("lineage", "replicate")) %>% unstack(lineage ~ replicate)
  
  ggVennDiagram(rp, label_size= 3, label_alpha = 0, category.names = c("R1", "R2", "R3"), set_size = 4) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#3A9AB2") +
    theme(legend.position = "none", plot.title=element_text(hjust=0.5))
}

lineageIdentificationMetrics <- function(data.metrics, colScale){
  data.temp <- data.metrics %>% 
    select(tool_name, PP, TP, FP, FN, TPR, FNR, PPV, FDR, jaccard_index, F1_score) %>% 
    group_by(tool_name) %>% 
    summarise(
      PP = median(PP, na.rm = TRUE),
      TP = median(TP, na.rm = TRUE),
      FP = median(FP, na.rm = TRUE),
      FN = median(FN, na.rm = TRUE),
      TPR = median(TPR, na.rm = TRUE),
      FNR = median(FNR, na.rm = TRUE),
      PPV = median(PPV, na.rm = TRUE),
      FDR = median(FDR, na.rm = TRUE),
      jaccard_index = median(jaccard_index, na.rm = TRUE),
      F1_score =  median(F1_score, na.rm = TRUE)) %>%
    gather(key="metric", value="value", -c("tool_name"))
  
  label.legend <- c(PP="Predicted Pos.", TP="True Pos.", FP = "False Pos.", FN = "False Neg.",
                    TPR = "True Pos. Rate", FNR="False Neg. Rate", PPV = "Pos. Predictive Value", FDR = "False Discovery Rate",
                    jaccard_index = "Jaccard Index", F1_score = "F1 score")
  
  ggarrange(
    ggplot(filter(data.temp, metric %in% c("PP", "TP", "FP", "FN")), aes(x=tool_name, y=value, color=tool_name, text=paste0("Tool: ", tool_name, 
                                                                                                                    "\n", "Result: ", round(value, digits = 2)))) +
      geom_point(aes(shape=metric, color=tool_name), size=4, stroke=1) +
      theme_minimal() + ylab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15)) +
      geom_vline(aes(xintercept="simulation"),linetype="dashed")  + 
      facet_wrap(~metric, labeller=as_labeller(label.legend)) +
      theme(
        legend.position = "none",
        text = element_text(size=12)) +
      labs(x = "") +
      colScale,
    
    ggplot(filter(data.temp, metric %in% c("TPR", "FNR", "PPV", "jaccard_index")), aes(x=tool_name, y=value, color=tool_name, text=paste0("Tool: ", tool_name, 
                                                                                                                    "\n", "Result: ", round(value, digits = 2)))) +
      geom_point(aes(shape=metric, color=tool_name), size=4, stroke=1) +
      theme_minimal() + ylab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15)) +
      geom_vline(aes(xintercept="simulation"),linetype="dashed")  + 
      facet_wrap(~metric, labeller=as_labeller(label.legend)) +
      theme(
        legend.position = "none",
        text = element_text(size=12)) +
      labs(x = "") +
      colScale,
    
    ggplot(filter(data.temp, metric %in% c("FDR", "F1_score")), aes(x=tool_name, y=value, color=tool_name, text=paste0("Tool: ", tool_name, 
                                                                                                      "\n", "Result: ", round(value, digits = 2)))) +
      geom_point(aes(shape=metric, color=tool_name), size=4, stroke=1) +
      theme_minimal() + ylab("") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15)) +
      geom_vline(aes(xintercept="simulation"),linetype="dashed")  + 
      facet_wrap(~metric, labeller=as_labeller(label.legend)) +
      theme(
        legend.position = "none",
        text = element_text(size=12)) +
      labs(x = "") +
      colScale
  )
}
  
tool_predictionHeatmap <- function(data.frames, sim.lineages, tool) {
  
  series <- list("Strict Mode", "Adj. 1 level", "Adj. 2 level", "Adj. 3 level", "Adj. 4 level")
  
  for (i in 1:length(data.frames)){
    s <- filter(data.frames[[i]], tool_name=="simulation") %>% group_by(lineage) %>% summarise(rel_abun = mean(rel_abun))
    t <- filter(data.frames[[i]], tool_name==tool) %>% group_by(lineage) %>% summarise(rel_abun = mean(rel_abun))
    
    comb.temp <- merge(x=s, y=t, by="lineage", all=TRUE, suffixes =  c(".simul", ".tool")) %>% 
      mutate("count" =  case_when(!is.na(rel_abun.simul) & !is.na(rel_abun.tool) ~ "Match",
                                  is.na(rel_abun.simul) & !is.na(rel_abun.tool) ~ "FP",
                                  !is.na(rel_abun.simul) & is.na(rel_abun.tool) ~ "FN"),
             rel_abun = rel_abun.tool,
             "series" = series[[i]]) %>% select("series", "lineage", "rel_abun", "count")
    
    if(i == 1){ comb.data <- comb.temp }
    else { comb.data <- rbind(comb.data, comb.temp) }
  }
  
  comb.data$count <- factor(comb.data$count, levels=c("Match", "FP", "FN"))
  comb.data$series <- factor(comb.data$series, levels=c("Strict Mode", "Adj. 1 level", "Adj. 2 level", "Adj. 3 level", "Adj. 4 level"))
  comb.data$lineage <- factor(comb.data$lineage, levels=unique(c(sim.lineages, unique(comb.data$lineage))))
  
  ggplot(comb.data, aes(x=lineage, y=count, fill=rel_abun)) +
    geom_tile() + facet_grid(series ~ .) +
    scale_fill_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous")) +
    labs(title = paste("Lineage identification -", tool)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          plot.title=element_text(hjust=0.5))
}

abundanceEstimation <- function(data.long, sim.data.long, title){
  
  plot_d <- NULL
  
  for (exp in unique(data.long$experiment)) {
    for (rep in unique(data.long$replicate)) {
      for (tp in unique(data.long$timepoint)) {
        for (tool in unique(data.long$tool_name)) {
          data.temp <- full_join(
            filter(data.long, experiment==exp, replicate==rep, tool_name==tool, timepoint==tp), 
            filter(sim.data.long, experiment==exp, replicate==rep, timepoint==tp), 
            by=c("timepoint", "lineage", "experiment", "replicate"), 
            suffix = c("","_sim")) %>%
            mutate("P" = case_when(!is.na(rel_abun_sim) ~ 1)) %>%
            mutate_at(c("rel_abun","rel_abun_sim", "P"), ~replace_na(.,0)) %>%
            mutate("Error" = case_when(P == 1 ~ (rel_abun - rel_abun_sim))) %>%
            mutate("tool_name" = tool) %>%
            filter(P == 1) %>%
            select(c("experiment", "replicate", "timepoint","tool_name","rel_abun","rel_abun_sim"))
          
          plot_d <- rbind(plot_d, data.temp)
        }
      }
    }
  }
  
  plot_d <- plot_d %>% mutate(timecourse = case_when(.$experiment %in% c("Ex09", "Ex10", "Ex11", "Ex12", "Ex13", "Ex14", "Ex15") ~ "WideQual",
                                                     .$experiment %in% c("Ex16", "Ex17", "Ex18", "Ex19", "Ex20", "Ex21", "Ex22") ~ "realTimecourse",
                                                     .$experiment %in% c("Ex23", "Ex24", "Ex25", "Ex26", "Ex27", "Ex28", "Ex29") ~ "nearQual")) %>%
    left_join(., unique(select(data.long, c("experiment", "replicate", "timepoint", "uniformity_wg_per"))), by=c("experiment", "replicate", "timepoint")) 
  
  label.legend <- c(nearQual = "Omicron time-course", realTimecourse = "Vienna time-course", WideQual = "VOC time-course",
                    freyja="Freyja", lollipop = "Lollipop", vaquero = "VaQuERo", vlq = "VLQ") 
  
  
  ggplot(filter(plot_d, tool_name != "simulation"), aes(x=rel_abun_sim, y=rel_abun)) +
    geom_jitter(alpha=0.1, aes(color=uniformity_wg_per)) +
    geom_abline(intercept = 0, slope = 1, linetype="longdash", alpha=0.2) +
    scale_color_gradient(name = "genomic\ncoverage\n(in %)\n", 
                         low = "red", high = "blue") +
    facet_grid(tool_name ~ timecourse, labeller=as_labeller(label.legend)) +
    labs(title = "tool") +
    theme(legend.position = "none") +
    labs(
      title=title,
      y="Predicted Relative Abundance (in %)",
      x="Simulated Relative Abundance (in %)",
    ) +
    theme_bw()

}

errorEstimation <- function(data.long, sim.data.long, title, colors){
  
  data.relAb.comb <- NULL
  
  for (exp in unique(data.long$experiment)) {
    for (rep in unique(data.long$replicate)) {
      for (tp in unique(data.long$timepoint)) {
        for (tool in unique(data.long$tool_name)) {
          data.temp <- full_join(
            filter(data.long, experiment==exp, replicate==rep, tool_name==tool, timepoint==tp), 
            filter(sim.data.long, experiment==exp, replicate==rep, timepoint==tp), 
            by=c("timepoint", "lineage", "experiment", "replicate"), 
            suffix = c("","_sim")) %>%
            mutate("P" = case_when(!is.na(rel_abun_sim) ~ 1)) %>%
            mutate_at(c("rel_abun","rel_abun_sim", "P"), ~replace_na(.,0)) %>%
            mutate("Error" = rel_abun - rel_abun_sim) %>%
            mutate("tool_name" = tool) %>%
            filter(rel_abun != 0.00) %>% 
            filter(rel_abun_sim != 0.00) %>%
            select(c("timepoint", "experiment", "replicate", "lineage", "rel_abun_sim", "tool_name", "rel_abun", "Error"))
          
          data.relAb.comb <- rbind(data.relAb.comb, data.temp)
        }
      }
    }
  }
  
  ggplot(data.relAb.comb, aes(x=rel_abun_sim, y=Error, fill=tool_name)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_point(shape=21, size=3, stroke=0.1) +
    stat_smooth(method="lm", aes(color=tool_name)) +
    #geom_line(stat="smooth", method="lm", se=TRUE, alpha=0.4) +
    xlim(0, 1) +
    ylim(-1, 1) +
    facet_wrap(~tool_name) +
    labs(
      x="True Abundance (%)",
      y="Error",
      title = title
    ) +
    scale_color_manual(values=colors) +
    scale_fill_manual(values=colors) +
    theme_bw() +
    theme(
      legend.position = "none",
      text = element_text(size=12)
    )
}  

tool_metricsPer <- function(data.metrics, x_value, colScale) {
  ggplot(select(data.metrics, -c("P", "PP", "TP", "FP", "FN")) %>% 
           gather(key="metric", value = "value", -c(x_value, "tool_name")), 
         aes_string(x=x_value, y="value", color="tool_name")) +
    geom_step(alpha=0.8, stroke=2) + 
    colScale +
    facet_wrap(~metric, ncol=2) +
    theme_minimal()
}

getDistance <- function(var1, var2, json.data){
  
  var1_nuc <- c(json.data[[var1]]$nucSubstitutions, json.data[[var1]]$nucDeletions)
  var2_nuc <- c(json.data[[var2]]$nucSubstitutions, json.data[[var2]]$nucDeletions)
  
  if (length(var1_nuc) == 0 || length(var2_nuc) == 0){
    return(Inf)
  }
  
  else {
    
    if (var1 == "B" || var2 == "B"){
      var1_nuc <- c(var1_nuc, "SARS_COV_2")
      var2_nuc <- c(var2_nuc, "SARS_COV_2")
    }
    
    var1_nuc <- gsub("-", "del", var1_nuc)
    var2_nuc <- gsub("-", "del", var2_nuc)
    
    var1_nuc <- as.character(stri_remove_empty(var1_nuc, na_empty = FALSE))
    var2_nuc <- as.character(stri_remove_empty(var2_nuc, na_empty = FALSE))
    
    percent_similarityPerIndex <- length(intersect(var1_nuc, var2_nuc)) / (length(union(var1_nuc, var2_nuc)))
    return (percent_similarityPerIndex)
  }
}

CheckDistanceMatrix <- function(lineages, update=FALSE) {
  if (!(file.exists("reference/distanceMatrix.csv")) | update) {
    distanceMatrix <- matrix()
  }
  
  else {
    distanceMatrix <- read.csv("reference/distanceMatrix.csv")
  }

  for (i in 1:length(lineages)){
    
    #print(paste("Calculating lineage ", i, "of", length(lineages)))
    
    lineage <- lineages[i]
    if(lineage != "others"){
      if (lineage %in% colnames(distanceMatrix)) {}
      else {
        if (length(colnames(distanceMatrix)) == 0){
          
          colnames(distanceMatrix) <- lineage
          rownames(distanceMatrix) <- lineage
        }
        else {
          distanceData <- NULL
          
          for (j in 1:length(colnames(distanceMatrix))){
            distanceData <- append(distanceData, getDistance(lineage, colnames(distanceMatrix)[j], json.data))
          }
          
          distanceMatrix <- cbind(distanceMatrix, matrix(data=distanceData, 
                                                         ncol=1, nrow=nrow(distanceMatrix),
                                                         dimnames=list(rownames(distanceData),
                                                                       lineage)
                                                         )
                                  )
          
          
          distanceMatrix <- rbind(distanceMatrix, matrix(data=NA, 
                                                         ncol=ncol(distanceMatrix), nrow=1,
                                                         dimnames = list(lineage, 
                                                                         colnames(distanceMatrix))
                                                         )
                                  )
          
        }
        index <- match(lineage, rownames(distanceMatrix))
        distanceMatrix[index, index]  <- getDistance(lineage, lineage, json.data)
        distanceMatrix <- round(distanceMatrix, 4)
        write.table(distanceMatrix, file="reference/distanceMatrix.csv", sep=",")
      }
    }
  }
  distanceMatrix[lower.tri(distanceMatrix)] <- t(distanceMatrix)[lower.tri(distanceMatrix)]
  distanceMatrix <- as.data.frame(distanceMatrix)
  return(distanceMatrix)
}

getminDistance <- function(lin, tp, sim.data.long) {
  
  Distance <- sim.data.long %>% 
    filter(timepoint == tp) %>% select("lineage") %>% filter(lineage !="others") %>% 
    unique() %>% rowwise() %>%
    mutate(similarity= distanceMatrix[lin, lineage])
  
  minDistance <- max(Distance$similarity)
  return(minDistance)
}

