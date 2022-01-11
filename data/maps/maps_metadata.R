library(sp)
library(rgdal)
library(OpenStreetMap)
library(raster)
library(here)
library(dplyr)
library(ggplot2)
library(ggmap)
library(gganimate)
library(lubridate)

select <- dplyr::select

# Get lat and lon from metadata
meta <- read.csv(here("uploads", "metadata_fixed.csv")) %>%
    select(sample = specimen.collector.sample.id, lab, folder, lat, lon, 
        date = sample.collection.date)

coco <- read.csv(here("data", "coco.csv"))
coco <- left_join(coco, meta, by = c("sample", "lab"))
coco <- coco[!is.na(coco$lat), ]


# Prepare for the results
res_df <- coco %>% 
    select(sample, folder, lab, date, location, lon, lat) %>%
    distinct() %>%
    mutate(prop = rep(NA, n()), 
        lo = rep(NA, n()), hi = rep(NA, n()),
        status = rep(NA, n()),
        totol_succ = rep(NA, n()),
        totol_fail = rep(NA, n())
        )

for (i in seq_len(nrow(res_df))) {
    thiscoco <- coco[coco$sample == res_df$sample[i] & coco$lab == res_df$lab[i] & coco$folder == res_df$folder[i], ]
    success <- round(as.numeric(thiscoco$count))
    failure <- round(as.numeric(thiscoco$coverage - success))
    res_df$totol_fail[i] <- sum(failure, na.rm = TRUE)
    res_df$totol_succ[i] <- sum(success, na.rm = TRUE)

    if (sum(success, na.rm = TRUE) < 3 | 
            sum(failure, na.rm = TRUE) < 10 |
            any(success > failure, na.rm = TRUE)) {
        res_df$status[i] <- "Insufficient data"
        next
    }

    cocomod <- glm(cbind(success, failure) ~ 1, family = binomial)
    res_df$prop[i] <- exp(cocomod$coef) / (1+exp(cocomod$coef))
    ci <- confint(cocomod)
    res_df$lo[i] <- exp(ci[1]) / (1+exp(ci[1]))
    res_df$hi[i] <- exp(ci[2]) / (1+exp(ci[2]))
}

res_df <- res_df[!is.na(res_df$prop),]

my_bbox <- make_bbox(lon = range(as.numeric(res_df$lon)), 
    lat = range(as.numeric(res_df$lat)), f = 0.1)
ontario3 <- get_map(location = my_bbox, source = "stamen")

animap <- ggmap(ontario3) +
    geom_point(data = res_df, size = 5,
        aes(x = as.numeric(lon), 
            y = as.numeric(lat),
            group = paste0(sample, lab),
            colour = as.numeric(prop))) +
    scale_colour_viridis_c(option = "B") +
    transition_time(ymd(res_df$date))
animap2 <- animate(animap)

print(animap2)
