# -------------------------------------------------------------------
# SETUP
# -------------------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(here)
source(here('analyses/prior_item.R'))



extract_study <- function() {
  subj_dirs <- dir('data/peers-data/', pattern = 'LTP')
  study <- data.frame()
  for (subj in subj_dirs) {
    for (session in 0:19) {
      file <- paste0('data/peers-data/',subj,'/session_',session,'/session.log')
      if (file.exists(file)) {
        tmp <- read.delim(file, header=F, fill=T, col.names=paste0('V', 1:12))
        tmp <- tmp %>% 
          filter(V3 == 'FR_PRES') %>% 
          select(V1, V4, V5, V6, V7, V8, V9)
        names(tmp) <- c('timestamp', 'listid','word','wordid', 'task','resp','rt')
        tmp$subject <- subj
        tmp$session <- session
        tmp$resp <- as.numeric(tmp$resp)
        study <- bind_rows(study, tmp)
      }
    }
  }
  study <- study %>% 
    group_by(subject, listid, session) %>% 
    mutate(trial = 1:length(word),
           timestamp = timestamp-min(timestamp),
           duration = c(timestamp[2:length(timestamp)]-timestamp[1:(length(timestamp)-1)],NA)) %>% 
    select(-timestamp)
  
  study <- study %>% 
    mutate(task = ifelse(task == -1, 'notask','task'),
           study_resp = ifelse(task == 'notask', NA, resp),
           study_rt = ifelse(task == 'notask', NA, rt)) %>% 
    select(-resp, -rt)
  return(study)
}

extract_test <- function() {
  subj_dirs <- dir('data/peers-data/', pattern = 'LTP')
  test <- data.frame()
  for (subj in subj_dirs) {
    for (session in 0:19) {
      for (listid in 0:15) {
        file <- paste0('data/peers-data/',subj,'/session_',session,'/',listid,'.par')
        if (file.exists(file)) {
          tmp <- read.delim(file, header=F, fill=T, col.names=paste0('V', 1:3))
          tmp$rt <- tmp$V1-c(0,tmp$V1[-length(tmp$V1)])
          names(tmp) <- c('rttime','wordid','word', 'rt')
          if (nrow(tmp) > 0) {
            tmp$subject <- subj
            tmp$session <- session
            tmp$listid <- listid
            test <- bind_rows(test, tmp)
          }
        }
      }
    }
  }
  test <- filter(test, wordid > 0)
  return(test)
}

study <- extract_study()
test <- extract_test()

write.csv(study, 'data/peers-cleaned/study.csv', row.names=F)
write.csv(test, 'data/peers-cleaned/test.csv', row.names=F)

# get study position
study <- study %>% 
  ungroup() %>% 
  mutate(listid = as.numeric(listid)) %>% 
  group_by(subject, session, listid) %>% 
  mutate(study_position = 1:length(subject))

# get test position
test <- test %>% 
  group_by(subject, session, listid) %>% 
  mutate(test_position = 1:length(subject))

# remove duplicate responses
test <- test %>% 
  group_by(subject, session, listid, word) %>% 
  mutate(rep = 1:length(word)) %>% 
  filter(rep == 1) %>% 
  select(-rep)

dat <- left_join(study, test)
dat$acc = ifelse(!is.na(dat$rt), 1, 0)

# get frequency info
subtlex <- read.csv('misc/prior-item-effects/materials/SUBTLEX_all_words.csv')
subtlex <- select(subtlex, Word, Lg10WF)
subtlex$Word <- toupper(subtlex$Word)
dat <- left_join(dat, subtlex, by=c('word'='Word'))

# get info about preceding freq
dat <- dat %>% 
  group_by(subject, session, listid) %>% 
  prior_item_analysis('Lg10WF','freq', max_lag = 4)

dat <- dat %>% 
  mutate(freq_prioritem1_bins = cut_number(as.numeric(freq_prioritem1), 20),
         freq_prioritem2_bins = cut_number(as.numeric(freq_prioritem2), 20),
         freq_prioritem3_bins = cut_number(as.numeric(freq_prioritem3), 20),
         freq_prioritem4_bins = cut_number(as.numeric(freq_prioritem4), 20)) %>% 
  group_by(freq_prioritem1_bins) %>% 
  mutate(ave_freq_prioritem1 = mean(as.numeric(freq_prioritem1), na.rm=T)) %>% 
  group_by(freq_prioritem2_bins) %>% 
  mutate(ave_freq_prioritem2 = mean(as.numeric(freq_prioritem2), na.rm=T)) %>% 
  group_by(freq_prioritem3_bins) %>% 
  mutate(ave_freq_prioritem3 = mean(as.numeric(freq_prioritem3), na.rm=T))%>% 
  group_by(freq_prioritem4_bins) %>% 
  mutate(ave_freq_prioritem4 = mean(as.numeric(freq_prioritem4), na.rm=T))

write.csv(dat, 'data/peers-cleaned/alldata.csv', row.names=F)
