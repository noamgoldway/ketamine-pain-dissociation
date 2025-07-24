# Load required libraries
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(corrplot)
library(pheatmap)
library(readr)
library(lme4)
library(emmeans)
library(afex)
library(Hmisc)   
library(ggpattern)
library(patchwork)
library(ggstatsplot)
library(ggside)
library(pheatmap)
library(gridExtra)
library(readxl)
library(broom)
library(R.matlab)
library(ComplexHeatmap)
library(circlize)




one_session_mri_list <- c("0435", "0614", "0642", "0839", "0863", "0873", "0991", "1008", "1137", "1175", "1405",
                      "1842", "2060", "2494", "2711", "3571", "3911", "4900", "5005", "5071", "6684", "7826",
                      "7932", "8295", "8298", "8446", "8754", "9018", "9235", "9364", "9477", "9501", "9616")

both_sessions_list <- c("0435", "0614", "0642", "0839", "0863", "0991", "1008", "1137", "1175", "1405",
                        "1842", "2060", "2494", "2711", "3440", "3571", "3911", "5005", "5071", "5407",
                        "7826", "7856", "7932", "8295", "8298", "8446", "8754", "9018", "9235", "9364",
                        "9501", "9616", "9653")


demog <- read.csv(
  "/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/demog.csv"
)

# ⚠️ Note: Gender is coded as 1 = Male, 2 = Female

# Define list of participants with both sessions
both_sessions_list <- c("0435", "0614", "0642", "0839", "0863", "0991", "1008", "1137", "1175", "1405",
                        "1842", "2060", "2494", "2711", "3440", "3571", "3911", "5005", "5071", "5407",
                        "7826", "7856", "7932", "8295", "8298", "8446", "8754", "9018", "9235", "9364",
                        "9501", "9616", "9653")
both_numeric <- as.numeric(both_sessions_list)

# Filter groups
demog_both <- demog %>% filter(subji %in% both_numeric)

# Function to compute summary stats
get_demog_summary <- function(df) {
  list(
    N = nrow(df),
    mean_age = round(mean(df$age, na.rm = TRUE), 2),
    sd_age = round(sd(df$age, na.rm = TRUE), 2),
    gender_counts = table(df$gender)
  )
}

# Compute summaries
summary_both <- get_demog_summary(demog_both)
summary_one <- get_demog_summary(demog)

# Print results
cat("Demographic Summary for Participants with One Session Only:\n")
cat("N =", summary_one$N, "\n")
cat("Mean Age =", summary_one$mean_age, "\n")
cat("SD Age =", summary_one$sd_age, "\n")
cat("Gender Distribution (1 = Male, 2 = Female):\n")
print(summary_one$gender_counts)


cat("Demographic Summary for Participants with Both Sessions:\n")
cat("N =", summary_both$N, "\n")
cat("Mean Age =", summary_both$mean_age, "\n")
cat("SD Age =", summary_both$sd_age, "\n")
cat("Gender Distribution (1 = Male, 2 = Female):\n")
print(summary_both$gender_counts)


########## Pain calibration##########

library(readxl)

# Load the dataset
Pain_calibration <- read_excel("/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/pain_calibration.xlsx")

# Paired t-test for nonpainful temperature (VAS = 2)
t_test_vas2 <- t.test(Pain_calibration$vas2_p, Pain_calibration$vas2_k, paired = TRUE)

# Paired t-test for painful temperature (VAS = 8)
t_test_vas8 <- t.test(Pain_calibration$vas8_p, Pain_calibration$vas8_k, paired = TRUE)

# Summary statistics
mean_vas2_p <- mean(Pain_calibration$vas2_p)
sd_vas2_p <- sd(Pain_calibration$vas2_p)
mean_vas2_k <- mean(Pain_calibration$vas2_k)
sd_vas2_k <- sd(Pain_calibration$vas2_k)

mean_vas8_p <- mean(Pain_calibration$vas8_p)
sd_vas8_p <- sd(Pain_calibration$vas8_p)
mean_vas8_k <- mean(Pain_calibration$vas8_k)
sd_vas8_k <- sd(Pain_calibration$vas8_k)

# Print results
cat("Nonpainful temp; placebo =", round(mean_vas2_p, 1), "±", round(sd_vas2_p, 2),
    "°C, ketamine =", round(mean_vas2_k, 1), "±", round(sd_vas2_k, 2),
    "°C,", "t(", t_test_vas2$parameter, ") =", round(t_test_vas2$statistic, 2),
    ", p =", round(t_test_vas2$p.value, 2), "\n")

cat("Painful temp; placebo =", round(mean_vas8_p, 1), "±", round(sd_vas8_p, 2),
    "°C, ketamine =", round(mean_vas8_k, 1), "±", round(sd_vas8_k, 2),
    "°C,", "t(", t_test_vas8$parameter, ") =", round(t_test_vas8$statistic, 2),
    ", p =", round(t_test_vas8$p.value, 2), "\n")



##########Effect of ketamine on pain perception ##########

Pain_ratings <- read.csv("/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/pain_ratings.csv", header = TRUE, stringsAsFactors = FALSE)

# Select relevant columns
Pain_ratings_wide <- dplyr::select(Pain_ratings, subji, rating.response_heat_high_1, rating.response_heat_high_2, rating.response_heat_low_1, rating.response_heat_low_2)

# Rename columns
colnames(Pain_ratings_wide) <- gsub("rating.response_heat_", "", colnames(Pain_ratings_wide))



# Create the Pain_data data frame
Pain_data <- Pain_ratings_wide %>%
  pivot_longer(
    cols = starts_with("high") | starts_with("low"),
    names_to = c("intensity", "session"),
    names_pattern = "(\\w+)_(\\d+)",
    values_to = "rating"
  )

# Convert subji and session to factors
Pain_data$subji <- as.factor(Pain_data$subji)
Pain_data$session <- as.factor(Pain_data$session)

# Map session levels to meaningful labels
Pain_data$session <- ifelse(Pain_data$session == '1', 'placebo', 'ketamine')

Pain_data <- Pain_data %>%
  mutate(intensity = case_when(
    intensity == "high" ~ "Pain",
    intensity == "low" ~ "No pain",
    TRUE ~ intensity
  ))

# Calculate summary statistics
RatingBySessionIntensity <- Pain_data %>% 
  group_by(session, intensity) %>% 
  dplyr::summarise(N = n(), 
                   meanRating = mean(rating, na.rm = TRUE), 
                   sdRating = sd(rating, na.rm = TRUE),
                   seRating = sdRating / sqrt(N))

# Define the order of session levels
session_order <- c("placebo", "ketamine")

# Reorder session levels in the data frames
RatingBySessionIntensity$session <- factor(RatingBySessionIntensity$session, levels = session_order)
Pain_data$session <- factor(Pain_data$session, levels = session_order)

# Create the rating plot
rating_plot <- ggplot(RatingBySessionIntensity, aes(x = intensity, y = meanRating, fill = session, color = intensity)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), size = 1) +
  scale_fill_manual(values = c("placebo" = "#ADD8E6", "ketamine" = "#FFB6C1")) +
  geom_errorbar(aes(ymin = meanRating - seRating, ymax = meanRating + seRating), 
                position = position_dodge(width = 0.9), width = 0.2) +
  labs(y = "Subjective pain rating", x = "", title = "") +
  scale_x_discrete(limits = c("Pain", "No pain")) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# Add data points to the rating plot
rating_plot + 
  geom_point(data = Pain_data, aes(x = intensity, y = rating, color = session), 
             position = position_dodge(width = 0.9), size = 2, color = "darkgray", na.rm = TRUE)



# Base plot with bars
rating_plot <- ggplot(RatingBySessionIntensity, aes(x = intensity, y = meanRating, fill = session)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 1.2) +
  scale_fill_manual(values = c("placebo" = "#ADD8E6", "ketamine" = "#FFB6C1")) +
  labs(y = "Subjective pain rating", x = "", title = "") +
  scale_x_discrete(limits = c("Pain", "No pain")) +
  theme_classic(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add points, then re-add error bars on top
rating_plot +
  geom_point(data = Pain_data, 
             aes(x = intensity, y = rating, fill = session), 
             position = position_dodge(width = 0.9), 
             shape = 21,
             size = 3,
             stroke = 0.3,
             color = "black",
             na.rm = TRUE) +
  geom_errorbar(data = RatingBySessionIntensity,
                aes(ymin = meanRating - seRating, ymax = meanRating + seRating),
                position = position_dodge(width = 0.9),
                width = 0.2,
                size = 0.7)  # Optional: adjust thickness for visibility




###this is providing 'isSingular' reducing complexity 
Pain_ratings.m <- mixed(rating~intensity*session+(intensity+session+1|subji),data=Pain_data,control=lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=1e6)),method="KR")

Pain_ratings.m <- mixed(rating~intensity*session+(intensity+1|subji),data=Pain_data,control=lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=1e6)),method="KR")

nice(Pain_ratings.m)
Pain_ratings.m.post<-emmeans(Pain_ratings.m, consec ~ session|intensity)
test(Pain_ratings.m.post, by = NULL, adjust = "fdr")




##########Effect of ketamine on dissociative state##########

# Read the CADSS data
CADSS_wide <- read.csv('/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/CADSS.csv')


CADSS<-CADSS_wide %>%pivot_longer(
  !subji,
  names_to = c("scale",  "sub_scale","session","time_point"),
  names_sep="_",
  values_to = "rating"
)

# Map sub_scale labels
CADSS$sub_scale <- ifelse(CADSS$sub_scale == 'Amnesiasum', 'Amnesia',
                          ifelse(CADSS$sub_scale == 'Depersonalizationsum', 'Depersonalization',
                                 ifelse(CADSS$sub_scale == 'Derealisationsum', 'Derealisation', NA)))

# Convert to factors
CADSS$sub_scale <- as.factor(CADSS$sub_scale)
CADSS$session <- ifelse(CADSS$session == '2', 'Placebo', 'Ketamine')
CADSS$session <- as.factor(CADSS$session)
levels(CADSS$session) <- c("Placebo", "Ketamine")

# Map time_point labels
CADSS$time_point <- ifelse(CADSS$time_point == 'baseline', 'Pre-infusion',
                           ifelse(CADSS$time_point == 'post', 'Post-bolus',
                                  ifelse(CADSS$time_point == "end", 'End of infusion', NA)))

# Convert to factors
CADSS$time_point <- as.factor(CADSS$time_point)
CADSS <- na.omit(CADSS)

# Calculate standard errors
RatingBySessionSubScaleTimePoint <- CADSS %>% 
  group_by(session, sub_scale, time_point) %>%
  dplyr::summarise(N = n(), 
                   meanRating = mean(rating, na.rm = TRUE), 
                   sdRating = sd(rating, na.rm = TRUE),
                   seRating = sdRating / sqrt(N))

CADSS_rating <- ggplot(CADSS, aes(x = time_point, y = rating)) +
  geom_point(aes(color = session, group = subji), size = 2) +
  geom_line(aes(color = session, group = subji), size = 0.5) +
  geom_line(data = RatingBySessionSubScaleTimePoint, aes(x = time_point, y = meanRating, group = session), size = 1, color = "black") +
  geom_errorbar(data = RatingBySessionSubScaleTimePoint, aes(x=time_point, y=meanRating, ymin = meanRating - seRating, ymax = meanRating + seRating), width=0.2) +
  
  facet_grid(session ~ sub_scale, labeller = label_wrap_gen(15)) + # Use facet_grid() instead of facet_wrap()
  labs(y = "Rating", x = "Time point", title = "") +
  scale_x_discrete(limits = c("Pre-infusion", "Post-bolus", "End of infusion")) +
  scale_color_manual(values = c("lightgray", "darkgray")) +
  theme_minimal(base_size = 24) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(), # Remove the background of facet labels
        strip.text = element_text(size = 24)) # Increase the size of facet labels

CADSS_rating



#more complex models gives singular values
CADSS_ratings.m <- mixed(rating~time_point*session*sub_scale+(1+session|subji),data=CADSS,control=lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=1e6)),method="KR")
nice(CADSS_ratings.m)

CADSS_ratings.m.session<-emmeans(CADSS_ratings.m, consec ~ session|time_point, adjust = "fdr") #all subslces are significant only in the ketamine session.
test(CADSS_ratings.m.session, by = NULL, adjust = "fdr")


CADSS_ratings.m.post<-emmeans(CADSS_ratings.m, consec ~ time_point|session|sub_scale, adjust = "fdr") #all subslces are significant only in the ketamine session.
test(CADSS_ratings.m.post, by = NULL, adjust = "fdr")


##########Correlation between pain and dissociation ratings ##########


# Prepare data for correlation analysis
CADSS_wide$Amnesia_k_min_p <- CADSS_wide$CADSS_Amnesiasum_2_post - CADSS_wide$CADSS_Amnesiasum_1_post
CADSS_wide$Depersonalization_k_min_p <- CADSS_wide$CADSS_Depersonalizationsum_2_post - CADSS_wide$CADSS_Depersonalizationsum_1_post
CADSS_wide$Derealisation_k_min_p <- CADSS_wide$CADSS_Derealisationsum_2_post - CADSS_wide$CADSS_Derealisationsum_1_post
CADSS_wide$sum <- CADSS_wide$CADSS_Amnesiasum_2_post + CADSS_wide$CADSS_Depersonalizationsum_2_post + CADSS_wide$CADSS_Derealisationsum_2_post
CADSS_for_corr <- dplyr::select(CADSS_wide, subji, CADSS_Amnesiasum_2_post, CADSS_Depersonalizationsum_2_post, CADSS_Derealisationsum_2_post, sum)
Pain_ratings_wide$high_k_min_p <- Pain_ratings_wide$high_2 - Pain_ratings_wide$high_1
Pain_ratings_wide$high_low_k_min_p <- (Pain_ratings_wide$high_2 - Pain_ratings_wide$low_2) - (Pain_ratings_wide$high_1 - Pain_ratings_wide$low_1)
Pain_for_corr <- dplyr::select(Pain_ratings_wide, subji, high_k_min_p, high_low_k_min_p, high_2)
Pain_CADSS_for_corr <- merge(Pain_for_corr, CADSS_for_corr, by = "subji")


# Rename columns
colnames(Pain_CADSS_for_corr) <- c("subji", "Pain high ketamine-placebo", "Pain high-low ketamine-placebo", "Pain high ketamine", "CADSS amnesia", "CADSS depersonalization", "CADSS derealisation", "CADSS sum")
Pain_CADSS_for_corr_nosubji <- Pain_CADSS_for_corr[, !(names(Pain_CADSS_for_corr) %in% "subji")]

# Perform correlation analysis
Cor_pain_data <- rcorr(as.matrix(Pain_CADSS_for_corr_nosubji), type = "spearman")



# Define the color palette
col <- colorRampPalette(c("lightblue", "white", "red"))(100)

col_fun <- colorRamp2(c(-1, 0, 1), c("lightblue", "white", "red"))

# Draw heatmap with legend at the bottom
Heatmap(Cor_pain_data$r,
        name = "Correlation",
        col = col_fun,
        show_column_names = FALSE,
        show_row_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        heatmap_legend_param = list(
          at = c(-1, -0.5, 0, 0.5, 1),
          title_position = "topcenter",
          legend_direction = "horizontal"
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Cor_pain_data$r[i, j]), x, y, gp = gpar(fontsize = 10))
        }) %>%
  draw(heatmap_legend_side = "bottom")



pheatmap(Cor_pain_data$P, 
         color = col, 
         # symm = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         show_colnames = TRUE,
         legend = FALSE)



##########fMRI results##########
##########Pain response: univariate analysis##########
##########Effects of Ketamine on pain related brain activationn##########


ROI_data_wide <- read.csv("/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/roi_beta_values_by_condition.csv")



# Pivot long the ROI_data_wide data frame

ROI_data <- ROI_data_wide %>%
  pivot_longer(
    cols = -c(subji, Condition),
    names_to = "ROI",
    values_to = "beta"
  ) %>%
  mutate(
    subji = as.factor(subji),
    ROI = as.factor(ROI),
    intensity = if_else(str_detect(Condition, "high"), "pain", "no-pain"),
    session = if_else(str_detect(Condition, "ketamine"), "ketamine", "placebo")
  )

ROI_data$beta_scale <- scale(ROI_data$beta)


ROI_activation_m <- mixed(
  beta_scale ~ intensity * session * ROI + (1 | subji),
  data = ROI_data,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
  method = "KR"
)

nice(ROI_activation_m)


# Ensure data is in correct format
ROI_data <- ROI_data %>%
  mutate(across(c(intensity, session, ROI), as.factor))

# Pivot wider to get separate columns for high and low pain
ROI_diff_data <- ROI_data %>%
  filter(intensity %in% c("pain", "no-pain")) %>%
  pivot_wider(
    id_cols = c(subji, session, ROI),
    names_from = intensity,
    values_from = beta
  ) %>%
  mutate(
    diff_beta = pain - `no-pain`  # Correct way to compute pain - no-pain
  )



pain_diff_model <- mixed(
  diff_beta ~ session * ROI + (1 | subji),
  data = ROI_diff_data,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
  method = "KR"
)

nice(pain_diff_model)


emm <- emmeans(pain_diff_model, ~ session | ROI)
contrast(emm, method = "pairwise", adjust = "fdr")


roi_labels <- c(
  "ant_insula_L" = "Ant. Insula (L)",
  "ant_insula_R" = "Ant. Insula (R)",
  "dlPFC_R" = "DLPFC (R)",
  "dACC" = "Dorsal ACC",
  "s1_R" = "S1 (R)",
  "s2_L" = "S2 (L)",
  "s2_R" = "S2 (R)"
)

ROI_diff_data <- ROI_diff_data %>%
  mutate(
    ROI_clean = recode(ROI, !!!roi_labels),
    ROI_clean = factor(ROI_clean, levels = unname(roi_labels)),
    session = factor(session, levels = c("placebo", "ketamine"))
  )

emms <- emmeans(pain_diff_model, ~ session | ROI)
summary_emms_df <- contrast(emms, method = "pairwise", adjust = "fdr") %>%
  as.data.frame() %>%
  mutate(
    ROI = recode(ROI, "dosal_stiatum_L" = "dorsal_stiatum_L"),
    ROI_clean = recode(ROI, !!!roi_labels),
    ROI_clean = factor(ROI_clean, levels = unname(roi_labels)),
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    session = factor("ketamine", levels = c("placebo", "ketamine"))
  )

ROI_summary <- ROI_diff_data %>%
  dplyr::group_by(session, ROI_clean) %>%
  dplyr::summarise(
    mean_diff = mean(diff_beta, na.rm = TRUE),
    se_diff = sd(diff_beta, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

asterisk_positions <- ROI_summary %>%
  group_by(ROI_clean) %>%
  summarise(y_pos = max(mean_diff + se_diff, na.rm = TRUE) + 0.1)

summary_emms_df <- summary_emms_df %>%
  left_join(asterisk_positions, by = "ROI_clean")

ggplot(ROI_summary, aes(x = session, y = mean_diff, fill = session)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6),
           width = 0.5, color = "black") +
  geom_errorbar(
    aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff),
    width = 0.2, position = position_dodge(width = 0.6)
  ) +
  geom_text(
    data = summary_emms_df,
    aes(x = 1.5, y = 0.9, label = sig),
    size = 12, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("placebo" = "#ADD8E6", "ketamine" = "#FFB6C1")) +
  facet_wrap(~ ROI_clean, ncol = 4, scales = "fixed") +
  labs(
    # title = "Effect of Session on Pain-Evoked Activation",
    x = "Session",
    y = "Mean Beta Difference (High Pain – Low Pain)",
  ) +
  theme_classic(base_size = 20) +
  theme(
    # plot.title = element_text(hjust = 0.5, size = 0, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


##########Correspondence between pain-related brain activation and pain ratings.##########




pain_long <- Pain_ratings_wide %>%
  mutate(
    pain_diff_placebo  = high_1 - low_1,
    pain_diff_ketamine = high_2 - low_2,
    pain_diff_diff     = pain_diff_ketamine - pain_diff_placebo
  ) %>%
  dplyr::select(
    subji,
    placebo_rating = high_1,
    ketamine_rating = high_2,
    rating_diff = high_k_min_p,
    pain_diff_placebo,
    pain_diff_ketamine,
    pain_diff_diff
  )

# Merge pain and ROI data by subject and session
pain_long_diff <- pain_long %>%
  dplyr::select(subji, pain_diff_placebo, pain_diff_ketamine) %>%
  pivot_longer(
    cols = starts_with("pain_diff"),
    names_to = "session",
    names_pattern = "pain_diff_(.*)",
    values_to = "pain_diff_rating"
  )

roi_diff_long <- ROI_diff_data %>%
  dplyr:: rename(pain_diff_activation = diff_beta) %>%
  dplyr:: select(subji, session, ROI, pain_diff_activation)

pain_roi_diff <- left_join(
  roi_diff_long %>% mutate(subji = as.character(subji)),
  pain_long_diff %>% mutate(subji = as.character(subji)),
  by = c("subji", "session")
)


# Define ROI labels
roi_labels <- c(
  "ant_insula_L" = "Ant. Insula (L)",
  "ant_insula_R" = "Ant. Insula (R)",
  "dlPFC_R"      = "DLPFC (R)",
  "dACC"         = "Dorsal ACC",
  "s1_R"         = "S1 (R)",
  "s2_L"         = "S2 (L)",
  "s2_R"         = "S2 (R)"
)

# Step 1: Compute correlation + p-value
cor_results <- pain_roi_diff %>%
  group_by(ROI, session) %>%
  summarise(
    p = cor.test(pain_diff_activation, pain_diff_rating, method = "spearman")$p.value,
    r = cor.test(pain_diff_activation, pain_diff_rating, method = "spearman")$estimate,
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    x_pos = 3,
    y_pos = case_when(
      session == "placebo"  ~ 1.9,
      session == "ketamine" ~ 1.6
    ),
    color = case_when(
      session == "placebo"  ~ "#136580",
      session == "ketamine" ~ "#b82339"
    )
  )

ggplot(pain_roi_diff, aes(x = pain_diff_rating, y = pain_diff_activation, color = session)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +  # Thicker line, subtle CI
  geom_text(
    data = cor_results %>% filter(sig != ""),
    aes(x = x_pos, y = y_pos, label = sig, color = session),
    inherit.aes = FALSE,
    size = 10
  ) +
  facet_wrap(~ ROI, scales = "fixed", ncol = 7, labeller = as_labeller(roi_labels)) +
  scale_color_manual(values = c("placebo" = "#136580", "ketamine" = "#b82339")) +
  labs(
    x = "Difference in pain rating (high - low)",
    y = "Difference in ROI activation \n(high - low)",
    color = "Session"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )


# Correspondence between pain-related brain activation and dissociation ratings.#####


# Step 1: Merge data
# Assuming 'pain_roi_diff' has a column named 'subj' that corresponds to 'subji' in CADSS_wide
merged_data <- merge(pain_roi_diff, CADSS_wide, by = "subji")
library(dplyr)
library(broom)

# Filter to ketamine session only
ketamine_data <- merged_data %>%
  filter(session == "ketamine")

# List of ROIs and CADSS predictors
roi_list <- unique(ketamine_data$ROI)
cadss_vars <- c("Depersonalization_k_min_p", "Derealisation_k_min_p", "Amnesia_k_min_p")


# Initialize a more informative results data frame
model_results <- data.frame()

for (roi in roi_list) {
  data_roi <- ketamine_data %>% filter(ROI == roi)
  
  # Model 1: pain_diff_rating only
  model1 <- lm(pain_diff_activation ~ pain_diff_rating, data = data_roi)
  
  # Model 2: pain_diff_rating + 3 dissociation scales
  model2 <- lm(pain_diff_activation ~ pain_diff_rating + 
                 Depersonalization_k_min_p + 
                 Derealisation_k_min_p + 
                 Amnesia_k_min_p, data = data_roi)
  
  # ANOVA model comparison
  model_comparison <- anova(model1, model2)
  p_anova <- model_comparison$`Pr(>F)`[2]
  f_stat <- model_comparison$F[2]
  df1 <- model_comparison$Df[2]
  df2 <- model_comparison$Res.Df[2]
  
  # Extract adjusted R² and AIC
  adj_r2_model1 <- summary(model1)$adj.r.squared
  adj_r2_model2 <- summary(model2)$adj.r.squared
  aic_model1 <- AIC(model1)
  aic_model2 <- AIC(model2)
  
  # Extract p-values for Model 2 predictors
  coef_summary <- summary(model2)$coefficients
  pain_p <- coef_summary["pain_diff_rating", "Pr(>|t|)"]
  dep_p <- coef_summary["Depersonalization_k_min_p", "Pr(>|t|)"]
  der_p <- coef_summary["Derealisation_k_min_p", "Pr(>|t|)"]
  amn_p <- coef_summary["Amnesia_k_min_p", "Pr(>|t|)"]
  
  # Append to results
  model_results <- rbind(model_results, data.frame(
    ROI = roi,
    pain_p = pain_p,
    dep_p = dep_p,
    der_p = der_p,
    amn_p = amn_p,
    p_anova = p_anova,
    f_stat = f_stat,
    df1 = df1,
    df2 = df2,
    adj_r2_model1 = adj_r2_model1,
    adj_r2_model2 = adj_r2_model2,
    aic_model1 = aic_model1,
    aic_model2 = aic_model2
  ))
}

# Optionally, add significance stars
model_results <- model_results %>%
  mutate(
    pain_sig = case_when(
      pain_p < 0.001 ~ "***",
      pain_p < 0.01  ~ "**",
      pain_p < 0.05  ~ "*",
      TRUE           ~ ""
    ),
    dissoc_sig = case_when(
      p_anova < 0.001 ~ "***",
      p_anova < 0.01  ~ "**",
      p_anova < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

model_results


######## Pain response: multivariate analysis########

NPS_response<-read.csv("/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/NPS.csv")

NPS_long <- NPS_response %>%
  pivot_longer(
    !subji,
    names_to = c("signature", "session", "intensity"),
    names_sep = "_",
    values_to = "value"
  )


# Perform mixed-effects model analysis
NPS_m <- mixed(value ~ intensity * session + (intensity + session + 1 | subji),
               data = NPS_long,
               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
               method = "KR")
nice(NPS_m)

# Calculate estimated marginal means (EMMs)
NPS_emm_session_intensity <- emmeans(NPS_m, consec ~ session | intensity, adjust = "fdr")
NPS_emm_intensity_session <- emmeans(NPS_m, consec ~ intensity | session, adjust = "fdr")

NPS_emm_session_intensity
NPS_emm_intensity_session
# Summarize the data for plotting
NPS_summary <- NPS_long %>%
  group_by(session, intensity) %>%
  dplyr::summarize(mean_value = mean(value),
                   se = sd(value) / sqrt(n()))

# Define color palette for the plot



NPS_plot <- ggplot(NPS_summary, aes(x = session, y = mean_value,color = intensity, fill = interaction(session, intensity))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), size = 1) +
  scale_fill_manual(values = fill_colors) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se), 
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.75) +
  labs(y = "Pain Signature", x = "") +
  scale_x_discrete(limits = c("placebo", "ketamine")) +
  theme_classic(base_size = 30) +
  #facet_wrap(~ ROI) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("") # add title


NPS_plot<-NPS_plot + 
  geom_point(data = NPS_long, aes(x = session, y = value), 
             position = position_dodge(width = 0.9), size = 2)
NPS_plot


#### NPS correlations #####

# Merge NPS response data with pain ratings data
Pain_NPS_for_corr <- merge(NPS_response, Pain_ratings_wide, by = "subji")

# Calculate difference scores
Pain_NPS_for_corr$high_p_min_k <- Pain_NPS_for_corr$high_1 - Pain_NPS_for_corr$high_2



nps_corr <- tibble(
  session = c("placebo", "ketamine"),
  x = c("NPS_placebo_high", "NPS_ketamine_high"),
  y = c("high_1", "high_2")
) %>%
  rowwise() %>%
  mutate(
    cor_test = list(cor.test(Pain_NPS_for_corr[[x]], Pain_NPS_for_corr[[y]], method = "spearman")),
    r = cor_test$estimate,
    p = cor_test$p.value,
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    r_label = paste0("ρ = ", round(r, 2), sig),
    x_pos = 3,
    y_pos = if_else(session == "placebo", 1.6, 1.6),
  )
# Plot placebo
p1 <- ggplot(Pain_NPS_for_corr, aes(x = NPS_placebo_high, y = high_1)) +
  geom_point(color = "#383637", alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#383637") +
  geom_text(data = nps_corr %>% filter(session == "placebo"),
            aes(x = x_pos, y = y_pos, label = r_label),
            inherit.aes = FALSE, size = 12, color = "#383637") + 
  labs(x = "Pain Signature", y = "Subjective Rating") +
  ylim(0, NA) +  # start y-axis at 0
  theme_minimal(base_size = 24) +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot ketamine
p2 <- ggplot(Pain_NPS_for_corr, aes(x = NPS_ketamine_high, y = high_2)) +
  geom_point(color = "#383637", alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "#383637") +
  geom_text(data = nps_corr %>% filter(session == "ketamine"),
            aes(x = x_pos, y = y_pos, label = r_label),
            inherit.aes = FALSE, size = 12, color = "#383637") +
  labs(x = "Pain Signature", y = "Subjective Rating") +
  ylim(0, NA) +  # start y-axis at 0
  theme_minimal(base_size = 24) +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine plots
p1 + p2


# Merge NPS and CADSS with pain ratings
merged_nps <- merge(Pain_NPS_for_corr,CADSS_wide, by = "subji") 
  #  mutate(
  #   pain_diff_rating = high_2 - low_2  # or just use high_2 if preferred
  # )

# MODEL 1: subjective pain ~ NPS only
model1 <- lm(NPS_ketamine_high ~ high_2, data = merged_nps)

# MODEL 2: subjective pain ~ NPS + dissociation subscales
model2 <- lm(NPS_ketamine_high ~ high_2 +
               Depersonalization_k_min_p +
               Derealisation_k_min_p +
               Amnesia_k_min_p,
             data = merged_nps)

# ANOVA comparison
model_comparison <- anova(model1, model2)
p_anova <- model_comparison$`Pr(>F)`[2]
f_stat <- model_comparison$F[2]
df1 <- model_comparison$Df[2]
df2 <- model_comparison$Res.Df[2]

# Model summaries
adj_r2_model1 <- summary(model1)$adj.r.squared
adj_r2_model2 <- summary(model2)$adj.r.squared
aic_model1 <- AIC(model1)
aic_model2 <- AIC(model2)

# p-values from full model
coef_summary <- summary(model2)$coefficients
Pain_p   <- coef_summary["high_2", "Pr(>|t|)"]
dep_p   <- coef_summary["Depersonalization_k_min_p", "Pr(>|t|)"]
der_p   <- coef_summary["Derealisation_k_min_p", "Pr(>|t|)"]
amn_p   <- coef_summary["Amnesia_k_min_p", "Pr(>|t|)"]

# Compile results
model_results <- data.frame(
  outcome = "NPS",
  Pain_p = Pain_p,
  dep_p = dep_p,
  der_p = der_p,
  amn_p = amn_p,
  p_anova = p_anova,
  f_stat = f_stat,
  df1 = df1,
  df2 = df2,
  adj_r2_model1 = adj_r2_model1,
  adj_r2_model2 = adj_r2_model2,
  aic_model1 = aic_model1,
  aic_model2 = aic_model2
) %>%
  mutate(
    dissoc_sig = case_when(
      p_anova < 0.001 ~ "***",
      p_anova < 0.01  ~ "**",
      p_anova < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

model_results



#######Neural correlates of dissociative states#######
file_path <- "/Users/noamgoldway/Library/CloudStorage/Box-Box/Goldway, Noam/tlvphd/manuscript-pain/data/git/within_connectivity.xlsx"
all_connectivity <- read_excel(path = file_path)

all_connectivity <- all_connectivity %>%
  rename_with(~ "Subject", all_of("...1"))


# Retrieve the first row to get the session names (PLC, KET)
session_names <- all_connectivity[1, -1]

# Drop the first row as it's going to be used for column names
all_connectivity <- all_connectivity[-1, ]

# Create a vector of the current network names (placeholders for now)
network_names <- colnames(all_connectivity)[-1]


# Ensure the Subject column is numeric
all_connectivity <- all_connectivity %>%
  mutate(Subject = as.numeric(Subject))

# Get valid subject IDs that appear in both ROI and NPS data
# valid_ids <- intersect(unique(ROI_data_wide$subji), unique(NPS_response$subji))

# Filter all_connectivity to keep only valid subjects
# all_connectivity <- all_connectivity %>%
#   filter(Subject %in% valid_ids)
# Replace placeholders with actual network names. Assuming you know the network names, you'd replace 'network1', 'network2', etc., with the actual names.
# Here's an example where we assume there are five networks, and they repeat twice (PLC then KET):
actual_network_names <- rep(c("DefaultMode", "Salience", "Frontoparietal", "SensoryMotor", "DorsalAttention"), each = 2)

# Rename the columns by pasting the actual network names with the session names
colnames(all_connectivity)[-1] <- paste(actual_network_names, session_names, sep = "_")
colnames(all_connectivity)[1] <- "subji"

# Pivot longer
all_connectivity_long <- all_connectivity %>%
  pivot_longer(cols = -subji, names_to = c("network", "session"), names_sep = "_", values_to = "connectivity")



all_connectivity_long$connectivity<-as.numeric(all_connectivity_long$connectivity)


connectivty.m <- mixed(connectivity~session*network+(1|subji),data=all_connectivity_long,control=lmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=1e6)),method="KR")
nice(connectivty.m)

# Split data by network and run t-tests
t_tests_results <- all_connectivity_long %>%
  group_by(network) %>%
  do(tidy(t.test(connectivity ~ session, data = .)))

t_tests_results <- t_tests_results %>%
  mutate(p.value.adjusted = p.adjust(p.value, method = "BH"))

t_tests_results



# Adjusting the factor levels so PLC appears on the left and KET on the right
all_connectivity_long$session <- factor(all_connectivity_long$session, levels = c("PLC", "KET"))


# Update session labels
all_connectivity_long$session <- recode(
  all_connectivity_long$session,
  "PLC" = "Placebo",
  "KET" = "Ketamine"
)

# Plot
all_networks_plot <- ggplot(all_connectivity_long, aes(x = session, y = connectivity, group = subji)) +
  geom_point(aes(color = session), position = position_dodge(width = 0.2)) +
  geom_line(aes(group = subji, color = session), alpha = 0.5, position = position_dodge(width = 0.2)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = 1), color = "red", shape = 18, position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(group = 1), position = position_dodge(width = 0.2)) +
  labs(y = "Within-Network Connectivity", x = "Condition") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_manual(values = c("Placebo" = "lightgray", "Ketamine" = "darkgray")) +
  facet_wrap(~network, scales = "free_y")

# Display
all_networks_plot



dmn_data <- all_connectivity_long %>%
  filter(network == "DefaultMode") 
# Calculate summary stats

summary_stats <- dmn_data %>%
  dplyr::mutate(session = forcats::fct_drop(session)) %>%
  dplyr::group_by(session) %>%
  dplyr::summarise(
    mean = mean(connectivity, na.rm = TRUE),
    se = sd(connectivity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Ensure summary_stats has session as factor to match plotting order
summary_stats$session <- factor(summary_stats$session, levels = c("Placebo", "Ketamine"))

# Plot
ggplot(dmn_data, aes(x = session, y = connectivity)) +
  # Individual subject lines
  geom_line(aes(group = subji), color = "grey70", alpha = 0.6) +
  geom_point(color = "grey40", size = 2) +
  
  # Group mean and SE
  geom_errorbar(data = summary_stats,
                aes(x = session, ymin = mean - se, ymax = mean + se),
                width = 0.2, size = 0.8, inherit.aes = FALSE) +
  geom_line(data = summary_stats,
            aes(x = session, y = mean, group = 1),
            color = "black", size = 1.2, inherit.aes = FALSE) +
  geom_point(data = summary_stats,
             aes(x = session, y = mean),
             color = "black", size = 3, shape = 18, inherit.aes = FALSE) +
  
  # Significance marker
  annotate("text", x = 1.5, y = max(dmn_data$connectivity, na.rm = TRUE) * 0.98,
           label = "***", size = 10, fontface = "bold") +
  
  labs(x = "", y = "Within Network Connectivity") +
  theme_minimal(base_size = 24) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# List of network names
# Define only the networks that actually exist in your data
networks <- c("DefaultMode", "Salience", "Frontoparietal", "SensoryMotor", "DorsalAttention")

for (network in networks) {
  # Construct expected column names
  plc_col <- paste(network, "PLC", sep = "_")
  ket_col <- paste(network, "KET", sep = "_")
  diff_col <- paste(network, "diff", sep = "_")
  
  # Check if both columns exist in the data
  if (plc_col %in% colnames(all_connectivity) && ket_col %in% colnames(all_connectivity)) {
    # Ensure numeric conversion
    all_connectivity[[plc_col]] <- as.numeric(all_connectivity[[plc_col]])
    all_connectivity[[ket_col]] <- as.numeric(all_connectivity[[ket_col]])
    
    # Calculate and store the difference
    all_connectivity[[diff_col]] <- all_connectivity[[plc_col]] - all_connectivity[[ket_col]]
  } else {
    message(paste("Skipping", network, "- columns not found"))
  }
}



names(all_connectivity) <- gsub(" ", "_", names(all_connectivity))




# Selecting subji and DMN_diff from all_connectivity
DMN_diff_subset <- all_connectivity %>%
  dplyr::select(subji, DefaultMode_diff)

# Merging DMN_diff_subset with Pain_NPS_CADSS by subji
#Pain_NPS_CADSS <- left_join(Pain_NPS_CADSS, DMN_diff_subset, by = "subji")
# Convert subji column to character in Pain_CADSS_for_corr
Pain_CADSS_for_corr <- Pain_CADSS_for_corr %>%
  mutate(subji = as.character(subji))

# Convert subji column to character in DMN_diff_subset
DMN_diff_subset <- DMN_diff_subset %>%
  mutate(subji = as.character(subji))

# Perform the left join
Pain_NPS_CADSS <- left_join(Pain_CADSS_for_corr, DMN_diff_subset, by = "subji")





Pain_NPS_CADSS$DMN_diff<-as.numeric(Pain_NPS_CADSS$DefaultMode_diff)
Pain_NPS_CADSS$CADSS_sum<-as.numeric(Pain_NPS_CADSS$"CADSS sum")
# Compute Spearman correlation
dmn_cadss_test <- suppressWarnings(
  cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS sum`, method = "spearman")
)

# Extract R and p-value
r_value <- round(dmn_cadss_test$estimate, 2)
p_value <- dmn_cadss_test$p.value

# Create significance label
sig_label <- ifelse(p_value < 0.001, "***",
                    ifelse(p_value < 0.01, "**",
                           ifelse(p_value < 0.05, "*", "")))

r_text <- paste0("r = ", r_value, sig_label)

# Plot with consistent aesthetics
ggplot(Pain_NPS_CADSS, aes(x = DefaultMode_diff, y = CADSS_sum)) +
  geom_point(color = "#383637", alpha = 0.6,size=3) +
  geom_smooth(method = "lm", se = TRUE, color = "#383637", size = 2) +
  geom_text(aes(x = 0.03, y = min(CADSS_sum, na.rm = TRUE)*0.95, label = r_text),
            inherit.aes = FALSE, size = 10, color = "#383637") +
  labs(x = "DMN within-network connectivity \n (placebo - ketamine)",
       y = "CADSS sum") +
  ylim(0, NA) +
  theme_minimal(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS amnesia`, method = "spearman")
cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS depersonalization`, method = "spearman")
cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS derealisation`, method = "spearman")


