library(tidyverse)
library(lubridate)

specimens <- read.csv("epidemiology.csv") %>% 
  separate(week, into = c("year", "week"), sep = 4) %>% 
  mutate(week = as.integer(week)) %>% 
  mutate(year = as.integer(year))

# dates <- read.csv("presentation_dates.csv") %>% 
#   mutate(week = str_pad(X2022_week, 2, pad = "0")) %>% 
#   mutate(week = as.integer(paste0("2022", week))) %>% 
#   mutate(count = 1) %>% 
#   group_by(week) %>% 
#   summarise(sum(count)) %>% 
#   rename(cases = 'sum(count)')


dates <- read.csv("presentation_dates.csv") %>%
  select(-date) %>% 
  rename(week = X2022_week) %>% 
  mutate(year = 2022) %>% 
  #mutate(count = 1) %>% 
  #group_by(year, week) %>% 
  #summarise(sum(count)) %>% 
  #rename(cases = 'sum(count)') %>% 
  mutate(year = as.integer(year))

numbers <- seq(1, 21)
randomized_numbers <- sample(numbers)
dates$random_id <- randomized_numbers

epidemiology <- full_join(specimens, dates, by = c("year", "week")) %>% 
  mutate(x = make_date(year)) %>% 
  mutate(date = x + lubridate::weeks(week-1)) %>% 
  mutate(transplant = factor(transplant, levels = c("Transplant", "No transplant"), 
                                 labels = c("Transplant", "No transplant")))


plot <- ggplot(epidemiology, aes(x = date)) +
  geom_line(aes(y = number_of_specimens)) +
  geom_point(aes(y = random_id*5 + 100, color = transplant), size = 2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        #legend.title = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_date(date_labels = "%b\n%Y", date_breaks = "2 months") +
  scale_color_discrete( na.translate = F, type = c("#D95F02", "#1F78B4")) +
  ylab("Weekly adenovirus cases") +
  labs(color = "Unknown\nhepatitis cases")

ggsave(filename = "Figures/epidemiology.pdf", plot = plot, device = "pdf", 
       units = "in", width = 7, height = 4, dpi = 300)



            