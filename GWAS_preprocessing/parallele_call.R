# Functions that contact the NCBI should not be called in parallel

library(tidyverse)
library(RCurl)
library(curl)
library(furrr)
library(magrittr)
library(future)


links <- read_csv('links.csv')

## OMDB API for downloading images
times <- tibble(
  plan = character(),
  number = numeric(),
  iteration = numeric(),
  cores = numeric(),
  time = numeric()
)

imgAPI <- 'http://img.omdbapi.com/?apikey=[your-own-API-key]' #get your Own API key and replace it here
rootDir <- 'images/'
get_img <- function(id){
  url <- paste0(imgAPI,str_pad(id, 7, pad = "0"))
  if(url.exists(url)){
    curl_download(url = url ,destfile = paste0(rootDir,'tt',id,'.png'))
  }
}

for(n_movies in c(500, 1000, 2500)){
  for(cores in c(3, 5, 7)){
    plan(multiprocess, workers = cores)
    for(i in seq(1, 5)){
      start_time <- Sys.time()
      
      links %>%
        top_n(n_movies) %>%
        extract2(2) %>%
        future_map(get_img)
      end_time <- Sys.time()
      
      times <- times %>% add_row(
        plan = 'multiprocess',
        number = n_movies,
        iteration = i,
        cores = cores,
        time = as.numeric(end_time - start_time))
    }
  }
  
  for(i in seq(1, 5)){
    start_time <- Sys.time()
    
    links %>%
      top_n(n_movies) %>%
      extract2(2) %>%
      map(get_img)
    end_time <- Sys.time()
    
    times <- times %>% add_row(
      plan = 'sequential',
      number = n_movies,
      iteration = i,
      cores = 1,
      time = as.numeric(end_time - start_time))
  }
}



times %>%          
  mutate(time = ifelse(plan == 'sequential' & number == 2500, time * 60, time)) %>%
  group_by(plan, number, cores) %>%
  summarise(avg.time = mean(time)) %>%
  ggplot(aes(x = as.factor(number), y = avg.time, fill = as.factor(cores))) +
  geom_bar(position = 'dodge', stat = 'identity', alpha = 0.8) +
  labs(title = 'Performance of using parallel processing for API connection', 
       fill = '# of Cores',
       x = '# of Movies Poster Downloaded',
       y = 'Avg. Time to download all movies (seconds)') +
  theme(legend.position = 'bottom',
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white",
                                        colour = "white"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line( size=.25, color="lightgray" ),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line( size=.25, color="lightgray" )) +
  scale_fill_brewer(palette = 'BrBG')

ggsave('PerformanceComparison.jpg')
