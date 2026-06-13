library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(dplyr)
library(tidyr)
library(httr)
library(jsonlite)
library(ggplot2)

load("3_data_analysis/5_MTBLS749/annotation_combined2.rda")

length(unique(annotation_combined2$variable_id))


annotation_combined2 %>%
  group_by(variable_id) %>%
  summarise(any_true = any(super_class_match)) %>%
  pull(any_true) %>%
  sum()

annotation_combined2 %>%
  group_by(variable_id) %>%
  summarise(any_true = any(class_match)) %>%
  pull(any_true) %>%
  sum()

annotation_combined2 %>%
  group_by(variable_id) %>%
  summarise(any_true = any(sub_class_match)) %>%
  pull(any_true) %>%
  sum()

n_total <- 14

class_data <- data.frame(
  level = factor(
    c("Kingdom", "Superclass", "Class", "Subclass"),
    levels = c("Subclass", "Class", "Superclass", "Kingdom")  # top to bottom
  ),
  matched = c(14, 9, 4, 2),
  total   = n_total
) |>
  mutate(pct = matched / total * 100,
         label_bar = sprintf("%d/%d", matched, total),
         label_pct = sprintf("%.1f%%", pct))

p <- ggplot(class_data,
            aes(x = pct, y = level)) +
  geom_col(fill = "#0F6E56", width = 0.6, alpha = 0.8) +
  geom_text(aes(label = ifelse(
    pct == round(pct),
    paste0(round(pct), "%"),
    paste0(round(pct, 1), "%")
  )),
  hjust = -0.15, size = 3, colour = "grey20") +
  scale_x_continuous(limits = c(0, 115),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25%", "50%", "75%", "100%"),
                     expand = c(0, 0)) +
  labs(x = "Consistency rate (%)",
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.line.y  = element_line(linewidth = 0.4),
        axis.ticks.y = element_blank(),
        axis.line.x  = element_line(linewidth = 0.4),
        axis.ticks.x = element_line(linewidth = 0.4),
        axis.title   = element_text(size = 10),
        plot.margin  = margin(8, 12, 8, 8))

ggsave("4_manuscript/Figures/Figure_5/top1_remaining_class_accuracy.pdf", plot = p, width = 14, height = 6)

# p <- ggplot(class_data,
#             aes(x = pct, y = level)) +
#   geom_col(fill = "#0F6E56", width = 0.6, alpha = 0.8) +
#   geom_text(aes(label = label_bar),
#             hjust = -0.15, size = 3, colour = "grey20") +
#   scale_x_continuous(limits    = c(0, 115),
#                      breaks    = c(0, 25, 50, 75, 100),
#                      labels    = c("0", "25%", "50%", "75%", "100%"),
#                      expand    = c(0, 0)) +
#   labs(x = "Consistency rate (%)",
#        y = NULL,
#        tag = "B") +
#   theme_classic(base_size = 11) +
#   theme(axis.line.y  = element_line(linewidth = 0.4),
#         axis.ticks.y = element_blank(),
#         axis.line.x  = element_line(linewidth = 0.4),
#         axis.ticks.x = element_line(linewidth = 0.4),
#         axis.title   = element_text(size = 10),
#         plot.tag     = element_text(face = "bold", size = 12),
#         plot.margin  = margin(8, 12, 8, 8))



# p <- ggplot(class_data,
#                   aes(x = pct, y = level)) +
#   geom_col(fill = "#0F6E56", width = 0.6, alpha = 0.8) +
#   geom_text(aes(label = label_bar),
#             hjust = -0.15, size = 3, colour = "grey20") +
#   geom_vline(xintercept = 100, colour = "grey70",
#              linewidth = 0.4, linetype = "dashed") +
#   scale_x_continuous(limits    = c(0, 115),
#                      breaks    = c(0, 25, 50, 75, 100),
#                      labels    = c("0", "25%", "50%", "75%", "100%"),
#                      expand    = c(0, 0)) +
#   labs(x = "Consistency rate (%)",
#        y = NULL,
#        tag = "B") +
#   theme_classic(base_size = 11) +
#   theme(axis.line.y  = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.x  = element_line(linewidth = 0.4),
#         axis.ticks.x = element_line(linewidth = 0.4),
#         axis.title   = element_text(size = 10),
#         plot.tag     = element_text(face = "bold", size = 12),
#         plot.margin  = margin(8, 12, 8, 8))

# data <- data.frame(
#   level = c("Kingdom", "Super Class", "Class", "Sub Class"),
#   matched = c(14, 9, 4, 2),
#   total = 14
# )
# 
# data$rate <- data$matched / data$total
# data$label <- paste0(scales::percent(data$rate, accuracy = 1), " (", data$matched, "/", data$total, ")")
# data$level <- factor(data$level, levels = c("Kingdom", "Super Class", "Class", "Sub Class"))
# 
# p <- 
# ggplot(data, aes(x = level, y = rate)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   geom_text(aes(label = label), vjust = -0.5) +
#   scale_y_continuous(labels = scales::percent, limits = c(0, 1.15)) +
#   labs(y = "Match Rate", title = "Feature Match Rate by Classification Level") +
#   theme_minimal()
# 
# ggsave("4_manuscript/Figures/Figure_5/class_match_percentage.pdf", p, width = 6, height = 4)

# ── 1. Function definitions ──────────────────────────────────────────
get_embedding <- function(text, api_key) {
  if (is.na(text) || text == "" || is.null(text)) return(NULL)
  
  res <- POST(
    url = "https://api.openai.com/v1/embeddings",
    add_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ),
    body = toJSON(list(
      model = "text-embedding-ada-002",
      input = text
    ), auto_unbox = TRUE),
    encode = "raw"
  )
  
  if (status_code(res) != 200) {
    message("HTTP Error ", status_code(res), ": ", content(res, "text"))
    return(NULL)
  }
  
  parsed <- fromJSON(content(res, "text", encoding = "UTF-8"))
  parsed$data$embedding[[1]]
}

cosine_similarity <- function(vec1, vec2) {
  if (is.null(vec1) || is.null(vec2)) return(NA_real_)
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

# ── 2. Set API Key ───────────────────────────────────────────────────
api_key <- Sys.getenv("OPENAI_API_KEY")

# ── 3. Combine description text ──────────────────────────────────────
annotation_combined2 <- annotation_combined2 %>%
  mutate(
    MS2_combined_desc = trimws(paste(
      ifelse(is.na(MS2_PUBCHEM_Description), "", MS2_PUBCHEM_Description),
      ifelse(is.na(ms2_description),         "", ms2_description)
    )),
    MS1_combined_desc = trimws(paste(
      ifelse(is.na(MS1_PUBCHEM_Description), "", MS1_PUBCHEM_Description),
      ifelse(is.na(ms1_description),         "", ms1_description)
    )),
    MS2_combined_desc = na_if(MS2_combined_desc, ""),
    MS1_combined_desc = na_if(MS1_combined_desc, "")
  )

# ── 4. Collect all unique texts (merged and deduplicated from 3 sources) ─
all_texts <- unique(na.omit(c(
  # HMDB
  annotation_combined2$ms2_description,
  annotation_combined2$ms1_description,
  # PubChem
  annotation_combined2$MS2_PUBCHEM_Description,
  annotation_combined2$MS1_PUBCHEM_Description,
  # Combined
  annotation_combined2$MS2_combined_desc,
  annotation_combined2$MS1_combined_desc
)))

message("Total unique texts to embed: ", length(all_texts))

# ── 5. Request embeddings (global cache) ─────────────────────────────
embedding_cache <- list()
for (txt in all_texts) {
  Sys.sleep(0.1)
  embedding_cache[[txt]] <- get_embedding(txt, api_key)
  message("Done: ", substr(txt, 1, 40))
}

# ── 6. Compute three similarity scores ───────────────────────────────
annotation_combined2 <- annotation_combined2 %>%
  rowwise() %>%
  mutate(
    similarity_hmdb = cosine_similarity(
      embedding_cache[[ms2_description]],
      embedding_cache[[ms1_description]]
    ),
    similarity_pubchem = cosine_similarity(
      embedding_cache[[MS2_PUBCHEM_Description]],
      embedding_cache[[MS1_PUBCHEM_Description]]
    ),
    similarity_combined = cosine_similarity(
      embedding_cache[[MS2_combined_desc]],
      embedding_cache[[MS1_combined_desc]]
    )
  ) %>%
  ungroup()


# Reshape to long format
df_long <- annotation_combined2 %>%
  dplyr::select(similarity_combined, similarity_hmdb, similarity_pubchem) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Type",
    values_to = "Similarity_Score"
  ) %>%
  mutate(Type = factor(Type,
                       levels = c("similarity_combined", "similarity_hmdb", "similarity_pubchem"),
                       labels = c("Combined", "HMDB", "PubChem")
  ))

# Draw violin plot
p <-
ggplot(df_long, aes(x = Type, y = Similarity_Score, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c(
    "Combined" = "#D4781A",
    "HMDB"     = "#2E8B57",
    "PubChem"  = "#2C6EA6"
  )) +
  scale_y_continuous(limits = c(0.6, 1.1), breaks = seq(0.6, 1.1, 0.1)) +
  labs(
    x = NULL,
    y = "Similarity Score",
    title = "Distribution of Biotext Similarity",
    fill = "Source"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12)
  )

ggsave("4_manuscript/Figures/Figure_5/biotext_similarity.pdf", p, width = 6, height = 4, dpi = 300)



# library(dplyr)
# library(httr)
# library(jsonlite)
# library(ggplot2)
# 
# # ── 1. Get embedding function ────────────────────────────────────────
# get_embedding <- function(text, api_key) {
#   if (is.na(text) || text == "" || is.null(text)) return(NULL)
#   
#   res <- POST(
#     url = "https://api.openai.com/v1/embeddings",
#     add_headers(
#       Authorization = paste("Bearer", api_key),
#       `Content-Type` = "application/json"
#     ),
#     body = toJSON(list(
#       model = "text-embedding-ada-002",
#       input = text
#     ), auto_unbox = TRUE),
#     encode = "raw"
#   )
#   
#   if (status_code(res) != 200) {
#     message("HTTP Error ", status_code(res), ": ", content(res, "text"))
#     return(NULL)
#   }
#   
#   parsed <- fromJSON(content(res, "text", encoding = "UTF-8"))
#   parsed$data$embedding[[1]]
# }
# 
# # ── 2. Cosine similarity function ───────────────────────────────────
# cosine_similarity <- function(vec1, vec2) {
#   if (is.null(vec1) || is.null(vec2)) return(NA_real_)
#   sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
# }
# 
# # ── 3. Set API Key ───────────────────────────────────────────────────
# api_key <- Sys.getenv("OPENAI_API_KEY")
# 
# # ── 4. Combine PubChem and HMDB description text ────────────────────
# annotation_combined2 <- annotation_combined2 %>%
#   mutate(
#     MS2_combined_desc = trimws(paste(
#       ifelse(is.na(MS2_PUBCHEM_Description), "", MS2_PUBCHEM_Description),
#       ifelse(is.na(ms2_description),         "", ms2_description)
#     )),
#     MS1_combined_desc = trimws(paste(
#       ifelse(is.na(MS1_PUBCHEM_Description), "", MS1_PUBCHEM_Description),
#       ifelse(is.na(ms1_description),         "", ms1_description)
#     )),
#     # Convert empty strings to NA after combining
#     MS2_combined_desc = na_if(MS2_combined_desc, ""),
#     MS1_combined_desc = na_if(MS1_combined_desc, "")
#   )
# 
# # ── 5. Get embeddings for all unique texts (with cache) ──────────────
# all_texts_combined <- unique(na.omit(c(
#   annotation_combined2$MS2_combined_desc,
#   annotation_combined2$MS1_combined_desc
# )))
# 
# message("Total unique texts to embed: ", length(all_texts_combined))
# 
# embedding_cache <- list()
# for (txt in all_texts_combined) {
#   Sys.sleep(0.1)
#   embedding_cache[[txt]] <- get_embedding(txt, api_key)
#   message("Done: ", substr(txt, 1, 40))
# }
# 
# # ── 6. Compute cosine similarity per row ─────────────────────────────
# annotation_combined2 <- annotation_combined2 %>%
#   rowwise() %>%
#   mutate(
#     similarity = cosine_similarity(
#       embedding_cache[[MS2_combined_desc]],
#       embedding_cache[[MS1_combined_desc]]
#     )
#   ) %>%
#   ungroup()
