library(tidyverse)
library(lubridate)

# Read in data for adenovirus cases over time
specimens <- read.csv("epidemiology.csv") %>% 
  separate(week, into = c("year", "week"), sep = 4) %>% 
  mutate(week = as.integer(week)) %>% 
  mutate(year = as.integer(year))

# Read in dates of presentation for cases
dates <- read.csv("presentation_dates.csv") %>%
  select(-date) %>% 
  rename(week = X2022_week) %>% 
  mutate(year = 2022) %>% 
  mutate(year = as.integer(year))

# Random heights to display the presentation dates to avoid crowding
numbers <- seq(1, 21)
randomized_numbers <- sample(numbers)
dates$random_id <- randomized_numbers

# Join together
epidemiology <- full_join(specimens, dates, by = c("year", "week")) %>% 
  mutate(x = make_date(year)) %>% 
  mutate(date = x + lubridate::weeks(week-1)) %>%
  mutate(transplant = factor(transplant, levels = c("Transplant", "No transplant"), 
                                 labels = c("Transplant", "No transplant")))

# Plot
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



            