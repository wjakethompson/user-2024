library(tidyverse)
library(modelr)
library(measr)
library(mirt)
library(spotifyr)
library(taylor)
library(here)
library(fs)

attributes <- 3
items_per_att <- 7
single_att_items <- 2
sample_size <- 500

set.seed(121389)

# Functions --------------------------------------------------------------------
logit <- function(x) log(x / (1 - x))
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Generate Q-matrix ------------------------------------------------------------
possible_items <- create_profiles(attributes = attributes) |> 
  rowwise() |> 
  filter(between(sum(c_across(everything())), 1, 2)) |> 
  ungroup()

qmatrix <- map_dfr(
  paste0("att", seq_len(attributes)),
  \(x, items, s_wt = 0.3, total, single) {
    bind_rows(
      items |> 
        rowwise() |> 
        filter(!!sym(x) == 1, sum(c_across(everything())) == 1) |> 
        ungroup() |> 
        slice_sample(n = single_att_items, replace = TRUE),
      items |> 
        rowwise() |> 
        filter(!!sym(x) == 1) |> 
        ungroup() |> 
        mutate(weight = c(s_wt, rep((1 - s_wt) / (n() - 1), n() - 1))) |> 
        slice_sample(n = (items_per_att - single_att_items), weight_by = weight,
                     replace = TRUE) |> 
        select(-weight)
    )
  },
  items = possible_items, total = items_per_att, single = single_att_items
) |> 
  slice_sample(n = (items_per_att * attributes), replace = FALSE)

# Generate item parameters -----------------------------------------------------
all_params <- qmatrix |> 
  model_matrix(~.^2) |> 
  rowid_to_column(var = "item_id") |> 
  pivot_longer(cols = -item_id, names_to = "parameter", values_to = "needed") |> 
  filter(needed == 1)

intercepts <- all_params |> 
  filter(parameter == "(Intercept)") |> 
  mutate(prob = runif(n = n(), min = 0.05, max = 0.45),
         value = logit(prob))

maineffects <- all_params |> 
  filter(parameter != "(Intercept)", !str_detect(parameter, ":")) |> 
  add_count(item_id, name = "num_att") |> 
  left_join(select(intercepts, item_id, prob, value), by = join_by(item_id)) |> 
  rename(int_prob = prob, int_value = value) |> 
  mutate(prob = map2_dbl(num_att, int_prob,
                         \(att, prob) {
                           cur_min <- ifelse(att == 1,
                                             max(prob, 0.60), max(prob, 0.2))
                           cur_max <- ifelse(att == 1,
                                             0.95, .7)
                           runif(1, cur_min, cur_max)
                         }),
         value = logit(prob) - int_value) |> 
  select(-c(num_att:int_value))

interactions <- all_params |> 
  filter(str_detect(parameter, ":")) |> 
  left_join(select(intercepts, item_id, value), by = join_by(item_id)) |> 
  rename(int_value = value) |> 
  left_join(select(maineffects, item_id, prob, value), by = join_by(item_id)) |>
  rename(mef_prob = prob, mef_value = value) |> 
  nest(params = c(int_value, mef_prob, mef_value)) |> 
  mutate(ints = map(params,
                    \(x) {
                      int_prob = runif(n = 1, min = max(max(x$mef_prob, 0.6)),
                                       max = 0.95)
                      tibble(prob = int_prob,
                             value = logit(int_prob) - unique(x$int_value) -
                               sum(x$mef_value))
                    })) |> 
  select(-params) |> 
  unnest(ints)

item_params <- bind_rows(intercepts, maineffects, interactions) |> 
  select(-prob, - needed) |> 
  arrange(item_id)

# Generate profiles ------------------------------------------------------------
profiles <- create_profiles(attributes = attributes) |> 
  slice_sample(n = sample_size, replace = TRUE) |> 
  rowid_to_column(var = "resp_id")

# Generate data ----------------------------------------------------------------
prof_profiles <- profiles |> 
  select(-resp_id) |> 
  model_matrix(~.^2) |> 
  rowid_to_column(var = "resp_id")

example_dat <- expand_grid(resp_id = seq_len(nrow(profiles)),
                           item_id = seq_len(nrow(qmatrix))) |> 
  left_join(prof_profiles, by = join_by(resp_id)) |> 
  pivot_longer(-c(resp_id, item_id), names_to = "parameter",
               values_to = "proficient") |> 
  left_join(item_params, by = join_by(item_id, parameter)) |> 
  mutate(final_value = value * proficient) |> 
  summarize(log_odds = sum(final_value, na.rm = TRUE),
            .by = c(resp_id, item_id)) |> 
  mutate(prob = inv_logit(log_odds),
         rand = runif(n = n()),
         score = case_when(rand < prob ~ 1L,
                           rand >= prob ~ 0L)) |> 
  select(resp_id, item_id, score) |> 
  pivot_wider(names_from = item_id, values_from = score)

# Estimate models --------------------------------------------------------------
## DCM
initial_dcm <- measr_dcm(data = example_dat, qmatrix = qmatrix,
                         resp_id = "resp_id",
                         type = "lcdm", backend = "cmdstanr",
                         chains = 4, parallel_chains = 4, seed = 121389,
                         iter_warmup = 1500, iter_sampling = 500,
                         file = here("fits", "data", "initial-lcdm"))

## IRT
mirt_type_2pl <- rep("2PL", nrow(qmatrix))
model_spec_2pl <- mirt.model(
  paste(
    paste0("F1 = 1-", attributes * items_per_att),
    paste0("PRIOR = (1-", attributes * items_per_att, ", a1, lnorm, 0, 1),"),
    paste0("        ", "(1-", attributes * items_per_att, ", d, norm, 0, 2)"),
    sep = "\n"
  )
)
initial_irt <- mirt(data = select(example_dat, -resp_id),
                    model = model_spec_2pl,
                    itemtype = mirt_type_2pl,
                    technical = list(set.seed = 121389))

# Identify albums --------------------------------------------------------------
irt_theta <- fscores(initial_irt) %>%
  as_tibble(.name_repair = ~"theta") %>%
  rowid_to_column(var = "resp_id")

dcm_profile <- predict(initial_dcm, summary = FALSE) |> 
  pluck("attribute_probabilities") |> 
  mutate(across(where(posterior::is_rvar), posterior::E)) |> 
  mutate(across(where(is.double), \(x) case_when(x < 0.5 ~ 0L, x >= 0.5 ~ 1L)))

profiles <- dcm_profile |> 
  mutate(resp_id = as.integer(as.character(resp_id))) |> 
  left_join(irt_theta, join_by(resp_id)) |> 
  nest(resp = c(resp_id, theta)) |> 
  mutate(min_theta = map_dbl(resp, \(x) min(x$theta)),
         max_theta = map_dbl(resp, \(x) max(x$theta)),
         range = max_theta - min_theta)


resp_albums <- tribble(
  ~album,                          ~era,             ~songwriting, ~production, ~vocals,
  "Taylor Swift",                  "Debut",               "xmark",     "xmark", "xmark",
  "Fearless",                      "Fearless",            "check",     "xmark", "xmark",
  "Fearless (Taylor's Version)",   "Fearless (TV)",       "check",     "xmark", "xmark",
  "Speak Now",                     "Speak Now",           "check",     "xmark", "xmark",
  "Speak Now (Taylor's Version)",  "Speak Now (TV)",      "check",     "xmark", "xmark",
  "Red",                           "Red",                 "check",     "check", "check",
  "Red (Taylor's Version)",        "Red (TV)",            "check",     "check", "check",
  "1989",                          "1989",                "xmark",     "check", "check",
  "1989 (Taylor's Version)",       "1989 (TV)",           "xmark",     "check", "check",
  "reputation",                    "reputation",          "check",     "xmark", "check",
  "Lover",                         "Lover",               "xmark",     "xmark", "check",
  "folklore",                      "folklore",            "check",     "check", "check",
  "evermore",                      "evermore",            "check",     "check", "check",
  "Midnights",                     "Midnights",           "check",     "xmark", "check",
  "THE TORTURED POETS DEPARTMENT", "TTPD",                "check",     "xmark", "check"
) |> 
  mutate(across(-c(album, era),
                \(x) case_when(x == "xmark" ~ 0L, x == "check" ~ 1L))) |> 
  mutate(target_theta = case_when(era == "Fearless" ~ -2.500,
                                  era == "Debut" ~ -2.080,
                                  era == "Red" ~ -1.670,
                                  era == "Speak Now" ~ -1.400,
                                  era == "1989" ~ -1.000,
                                  era == "reputation" ~ -0.4170,
                                  era == "evermore" ~ -0.200,
                                  era == "Fearless (TV)" ~ 0.317,
                                  era == "Lover" ~ 0.817,
                                  era == "Midnights" ~ 1.350,
                                  era == "TTPD" ~ 1.450,
                                  era == "folklore" ~ 1.500,
                                  era == "Speak Now (TV)" ~ 1.750,
                                  era == "1989 (TV)" ~ 2.000 ,
                                  era == "Red (TV)" ~ 2.500)) |> 
  left_join(profiles, join_by(songwriting == att1, production == att2,
                              vocals == att3)) |> 
  mutate(selection = map2(resp, target_theta,
                          \(r, t) {
                            r |> 
                              slice_min(n = 1, order_by = abs(theta - t),
                                        with_ties = FALSE)
                          })) |> 
  unnest(selection) |> 
  select(resp_id, album, era, songwriting, production, vocals,
         target_theta, theta)

# Taylor-fy the data -----------------------------------------------------------
access_token <- get_spotify_access_token()
playlists <- c("Hot Hits USA"     = "37i9dQZF1DX0kbJZpiYdZl",
               "Today's Top Hits" = "37i9dQZF1DXcBWIGoYBM5M",
               "Top Hits of 2024" = "774kUuKDzLa8ieaSmi8IfS",
               "Top Hits of 2023" = "6unJBM7ZGitZYFJKkO0e4P",
               "Top Hits of 2022" = "56r5qRUv3jSxADdmBkhcz7",
               "Top Hits of 2021" = "5GhQiRkGuqzpWZSE7OU4Se",
               "Top Hits of 2020" = "2fmTTbBkXi8pewbUvG3CeZ",
               "Top Hits of 2019" = "37i9dQZF1DWVRSukIED0e9",
               "Top Hits of 2018" = "5l771HfZDqlBsDFQzO0431",
               "Top Hits of 2017" = "0cNzbPLBSPIFUrJkqxkmUw",
               "Top Hits of 2016" = "0MvxEHiGaSFsKXry6Vjvb1",
               "Top Hits of 2015" = "4xYPZBwLtndo3MqjM0kMz1",
               "Top Hits of 2014" = "5oN9E0LneOjbrYkGrZmUc9",
               "Top Hits of 2013" = "1pQPyEwAz8Vs5ebAZY6y4G",
               "Top Hits of 2012" = "0rZUzeYjHNvHQhgMkq5C9W",
               "Top Hits of 2011" = "45yWuxLgXRJD3HnhOVZ3xz",
               "Top Hits of 2010" = "07lujzK1sKPxqzZqBq5ANi",
               "All Out 2000s"    = "37i9dQZF1DX4o1oenSJRJd",
               "All Out 2010s"    = "37i9dQZF1DX5Ejj0EkURtP")

albums <- map_dfr(playlists,
                  \(x) {
                    tracks <- get_playlist(playlist_id = x)
                    
                    tibble(album = tracks$tracks$items$track.album.name,
                           artist = tracks$tracks$items$track.album.artists,
                           type = tracks$tracks$items$track.album.album_type) |> 
                      filter(type == "album") |> 
                      select(-type) |> 
                      unnest(artist) |> 
                      select(album, artist = name)
                  }) |> 
  summarize(artist = paste(unique(artist), collapse = ","),
            .by = album) |> 
  filter(!str_detect(artist, "Taylor Swift")) |> 
  anti_join(resp_albums, join_by(album == era)) |> 
  distinct()

taylor_data <- example_dat |> 
  mutate(album = sample(albums$album, size = n(), replace = FALSE),
         .after = resp_id) |> 
  left_join(select(resp_albums, resp_id, ts_alb = album), join_by(resp_id)) |> 
  mutate(album = case_when(is.na(ts_alb) ~ album,
                           !is.na(ts_alb) ~ ts_alb)) |> 
  select(-resp_id, -ts_alb) |> 
  mutate(album = factor(album, levels = c(resp_albums$album, sample(albums$album)))) |> 
  arrange(album) |> 
  mutate(album = as.character(album)) |> 
  write_rds(here("data", "taylor-data.rds"))

taylor_qmatrix <- qmatrix |> 
  rename(songwriting = att1, production = att2, vocals = att3) |> 
  write_rds(here("data", "taylor-qmatrix.rds"))

# Taylor models ----------------------------------------------------------------
taylor_dcm <- measr_dcm(data = taylor_data, qmatrix = taylor_qmatrix,
                        resp_id = "album",
                        type = "lcdm", backend = "cmdstanr",
                        chains = 4, parallel_chains = 4, seed = 121389,
                        iter_warmup = 1000, iter_sampling = 500,
                        file = here("fits", "data", "taylor-lcdm"))

# Create final results ---------------------------------------------------------
popularity <- tribble(
  ~album_name,                           ~album_uri,
  "Taylor Swift",                        "7mzrIsaAjnXihW3InKjlC3",
  "Fearless",                            "43OpbkiiIxJO8ktIB777Nn",
  "Fearless (Taylor's Version)",         "4hDok0OAJd57SGIT8xuWJH",
  "Speak Now",                           "5EpMjweRD573ASl7uNiHym",
  "Speak Now (Taylor's Version)",        "5AEDGbliTTfjOB8TSm1sxt",
  "Red",                                 "1KlU96Hw9nlvqpBPlSqcTV",
  "Red (Taylor's Version)",              "6kZ42qRrzov54LcAk4onW9",
  "1989",                                "34OkZVpuzBa9y40DCy0LPR",
  "1989 (Taylor's Version)",             "64LU4c1nfjz1t4VnGhagcg",
  "reputation",                          "6DEjYFkNZh67HP7R9PSZvv",
  "Lover",                               "1NAmidJlEaVgA3MpcPFYGq",
  "folklore",                            "1pzvBxYgT6OVwJLtHkrdQK",
  "evermore",                            "6AORtDjduMM3bupSWzbTSG",
  "Midnights",                           "1fnJ7k0bllNfL1kVdNVW1A",
  "THE TORTURED POETS DEPARTMENT",       "5H7ixXZfsNMGbIE5OBSpcb"
) |> 
  mutate(popularity = map_int(album_uri,
                              function(.x) {
                                album <- get_album(.x)
                                album$popularity
                              }))

album_meta <- here("figure", "taylor", "eras") |> 
  dir_ls() |> 
  as_tibble() |> 
  rename(img = value) |> 
  mutate(img = str_replace_all(img, "^.*(?=figure)", ""),
         img = path_file(img),
         era = path_ext_remove(img),
         era = str_replace(era, "[0-9]{2}-", ""),
         era = str_replace_all(era, "-", " "),
         era = case_when(era %in% c("reputation", "folklore", "evermore") ~
                           era,
                         era == "ttpd" ~ "TTPD",
                         str_detect(era, "tv") ~ str_replace(str_to_title(era),
                                                             "Tv", "(TV)"),
                         TRUE ~ str_to_title(era)),
         album = case_when(str_detect(era, "\\(TV\\)") ~
                             str_replace(era, "TV", "Taylor's Version"),
                           era == "Debut" ~ "Taylor Swift",
                           era == "TTPD" ~ "THE TORTURED POETS DEPARTMENT",
                           TRUE ~ era)) |> 
  relocate(img, .after = last_col()) |> 
  left_join(select(popularity, -album_uri),
            join_by(album == album_name)) |> 
  left_join(select(taylor_albums, -ep),
            join_by(album == album_name)) |>
  relocate(album_release, .after = album)

predict(taylor_dcm, summary = FALSE) |> 
  pluck("attribute_probabilities") |> 
  mutate(across(where(posterior::is_rvar), posterior::E)) |> 
  semi_join(resp_albums, join_by(album)) |> 
  mutate(across(where(is.double),
                \(x) case_when(x < 0.5 ~ 0L, x >= 0.5 ~ 1L))) |> 
  left_join(select(resp_albums, album, era, songwriting:vocals, theta),
            join_by(album, songwriting, production, vocals)) |> 
  relocate(era, .after = album) |> 
  left_join(album_meta, join_by(era, album)) |> 
  select(album, era, album_release, img,
         songwriting:vocals, popularity:user_score, theta) |> 
  arrange(album_release) |> 
  write_rds(here("data", "taylor-results.rds"))
