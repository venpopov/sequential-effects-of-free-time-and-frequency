####
# The data reported here are from the Penn Electrophysiology of Encoding and Retrieval dat (PEERS)
# Experiment 1. 7 sessions, each session contains 16 lists of 16 words presented one at a time
# each dat list was followed by a free recall test
# check this paper for description http://memory.psych.upenn.edu/files/pubs/HealEtal14.pdf

# E1 - sessions 0 to 6
# E2 - sessions 7 to 12
# E3 - session 13 to 19

# From http://memory.psych.upenn.edu/files/pubs/KuhnEtal18.pdf :
# We analyze final free-recall data collected as part of the Penn
# Electrophysiology of Encoding and Retrieval Study (PEERS).
# Subjects recruited to PEERS took part in three subsequently administered multisession experiments, comprising a total of 20
# experimental sessions (of the 171 participants who completed the
#                        seven sessions of Experiment 1, 158 also completed six sessions of
#                        Experiment 2, and 151 completed the remaining sessions of Experiment 3). Subjects in these experiments studied and then freely
# recalled lists of 16 common words under various conditions (immediate free recall in Experiments 1 and 3; immediate, delayed,
#                                                             and continual-distractor free recall [CDFR] in Experiment 2).
# During half of the experimental sessions, subjects were also given
# a final-free-recall test, and the present paper is the first to report on
# these data. Experiment 3 differed from Experiment 1 in that a
# subset of subjects were asked to verbalize all words that came to
# mind during recall


#############################################################################
# SETUP
#############################################################################

# rm(list=ls())
library(tidyverse)
library(ggtext)
library(here)
library(lme4)
library(patchwork)
library(brms)
source(here('analyses/prior_item.R'))

mean_se2 <- function (x) {
  x <- stats::na.omit(x)
  se <- 1.96 * sqrt(stats::var(x)/length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}

theme_set(theme_classic(base_size=9))

#############################################################################
# DATA
#############################################################################
dat1 <- read.csv('data/peers-cleaned/alldata_part1.csv')
dat2 <- read.csv('data/peers-cleaned/alldata_part2.csv')
dat <- rbind(dat1,dat2)


dat <- dat %>% 
  group_by(subject, session, listid) %>% 
  arrange(subject, session, listid, test_position) %>% 
  mutate(resplag = study_position-c(NA,study_position[-length(study_position)])) %>% 
  mutate(resplag = ifelse(!is.na(test_position), resplag, NA)) %>% 
  ungroup()

dat <- dat %>% 
  mutate(exp = case_when(
    session <= 6 ~ 1,
    session > 6 & session <= 12 ~ 2,
    TRUE ~ 3
  ))

dat$sp_cat <- ifelse(dat$study_position <= 8, 'first half','second half')

dat <- dat %>% 
  mutate_at(vars(starts_with("ave_")), list(~round(., 1)))


exp1 <- dat %>% filter(exp == 1) %>%
  group_by(subject, session, listid) %>% 
  mutate(list_acc = mean(acc, na.rm=T)) %>% 
  ungroup() %>% 
  # filter(duration <= 4200 | is.na(duration)) %>%
  # filter(is.na(duration)) %>% 
  # mutate(duration = ifelse(duration>4200,4200,duration)) %>% 
  mutate(post = duration-3000)

exp1 <- exp1 %>% 
  group_by(subject, listid, session) %>% 
  arrange(study_position) %>% 
  prior_item_analysis('post','pre')
names(exp1)[names(exp1) == 'pre_prioritem'] <- 'pre1'

exp1 <- exp1 %>% 
  mutate(post = as.numeric(post),
         pre1 = as.numeric(pre1),
         post_bin = ceiling(post/100),
         pre1_bin = ceiling(pre1/100))

exp1 <- exp1 %>% 
  arrange(subject, session, listid, study_position) %>% 
  group_by(subject, session,listid) %>% 
  arrange(study_position) %>% 
  mutate(lagm1acc = acc,
         lagm2acc = c(NA, acc[-length(acc)]),
         lagm3acc = c(NA,NA, acc[-((length(acc)-1):length(acc))]),
         lagm4acc = c(NA,NA,NA, acc[-((length(acc)-2):length(acc))]),
         lag1acc = c(acc[-1],NA),
         lag2acc = c(acc[-(1:2)],NA,NA),
         lag3acc = c(acc[-(1:3)],NA,NA,NA),
         lag4acc = c(acc[-(1:4)],NA,NA,NA,NA))

exp1 <- exp1 %>% 
  arrange(subject, session, listid, study_position)

exp1 <- exp1 %>% 
  group_by(subject, listid, session) %>% 
  arrange(desc(study_position)) %>% 
  prior_item_analysis('Lg10WF','postfreq')

exp1 <- exp1 %>% 
  arrange(subject, session, listid, study_position) %>% 
  ungroup() %>% 
  mutate(postfreq_prioritem = as.numeric(postfreq_prioritem),
         postfreq_prioritem_q = ecdf(postfreq_prioritem)(postfreq_prioritem),
         postfreq_prioritem_bin = ceiling(postfreq_prioritem_q*20)) %>% 
  group_by(postfreq_prioritem_bin) %>% 
  mutate(postfreq_ave = mean(postfreq_prioritem, na.rm=T))
# -------------------------------------------------------------------
# Modeling prep
# -------------------------------------------------------------------
study <- exp1 %>% 
  select(subject, listid, word, wordid, task, session, trial, duration, 
         study_position, acc, Lg10WF, freq_prioritem1, ave_freq_prioritem1, 
         postfreq_prioritem1, postfreq_ave, post, pre1, post_bin, pre1_bin) %>% 
  mutate(procedure = 'study') %>% 
  arrange(subject, session, listid, study_position) %>% 
  group_by(subject, session, listid) %>% 
  mutate(onset = cumsum(ifelse(!is.na(pre1),pre1/1000,0)+3)-3)
test <- study %>% 
  mutate(procedure = 'test') %>% 
  group_by(subject, session, listid) %>% 
  mutate(onset = max(onset)+1.3 + 4*(0:(length(listid)-1)))
df1 <- bind_rows(study, test) %>% 
  arrange(subject, session, listid, procedure, study_position)

df1 <- df1 %>% 
  rename(list = listid,
         stim = word)

#remove lists with extreme ISIs
df1 <- df1 %>% 
  group_by(subject, session, list) %>% 
  mutate(postmax = max(post_bin,na.rm=T)) %>% 
  filter(postmax < 14)
  
write.csv(df1, 'data/peers_for_modeling.csv', row.names=F)

# -------------------------------------------------------------------
# SECTION
# -------------------------------------------------------------------


# sac_fit <- read.csv('output/peers_results_fit.csv')
# sac_fit$listid <- as.numeric(sac_fit$list)-1
# sac_fit$subject <- toupper(sac_fit$subject)
# 
# dat1 <- exp1 %>% 
#   filter(exp == 1) %>% 
#   left_join(select(sac_fit, subject, session, listid, wordid, acc_pred)) %>% 
#   gather(data, acc, acc, acc_pred) %>% 
#   mutate(data = recode(data, acc='Observed', acc_pred='SAC Model fit'))



# -------------------------------------------------------------------
# FIGURES
# -------------------------------------------------------------------


# main effect of prior freq
(f1 <- exp1 %>%
   # filter(subject == 'LTP079') %>% 
   ggplot(aes(ave_freq_prioritem1, acc)) +
   stat_summary(geom='point') +
   geom_smooth(method='lm', se=F) +
   scale_x_continuous('Mean log(freq) of preceding item during study', breaks=c(1,1.5,2,2.5,3,3.5,4)) +
   scale_y_continuous('P(recall)'))


# does the preceding isi matter overall?
exp1 %>% 
  filter(pre1_bin <= 13, post_bin <= 13) %>% 
  group_by(subject, pre1_bin) %>% 
  summarise(acc = mean(acc, na.rm=T)) %>% 
  filter(!is.na(pre1_bin)) %>% 
  ungroup() %>% 
  mutate(overall_mean = mean(acc, na.rm=T)) %>% 
  group_by(subject) %>% 
  mutate(acc = acc-mean(acc)+overall_mean,
         pre1_bin = as.factor(pre1_bin)) %>%
  ggplot(aes(pre1_bin, acc, group=1)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_discrete(name='**Pre**-ISI duration (ms.)', labels=c('800-900', '900-1000', '1000-1100', '1100-1200', '1200-1250')) +
  scale_y_continuous('Recall probability') +
  coord_cartesian(ylim=c(0.53,0.555)) +
  theme(axis.title.x = element_markdown(),
        axis.text.x = element_text(angle = 45,hjust=1))

ggsave('figures/peers_preisi_acc.png', w=3.25, h=2.75)


exp1 %>% 
  filter(acc==1, pre1_bin <= 13) %>% 
  group_by(subject, pre1_bin) %>% 
  summarise(rt = mean(rt, na.rm=T)) %>% 
  filter(!is.na(pre1_bin)) %>% 
  ungroup() %>% 
  mutate(overall_mean = mean(rt, na.rm=T)) %>% 
  group_by(subject) %>% 
  mutate(rt = rt-mean(rt)+overall_mean) %>% 
  ggplot(aes(as.factor(pre1_bin), rt, group=1)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_discrete(name='PRE* ISI duration (ms.)', labels=c('800-900', '900-1000', '1000-1100', '1100-1200', '1200-1250')) +
  scale_y_continuous('Response time (ms.)')


ggsave('figures/peers_isi_rt.png', w=4, h=3.5)

# effect of post_isi?
exp1 %>% 
  filter(post_bin <= 13,pre1_bin <= 13) %>% 
  group_by(subject, post_bin) %>% 
  summarise(acc = mean(acc, na.rm=T)) %>% 
  filter(!is.na(post_bin)) %>% 
  ungroup() %>% 
  mutate(overall_mean = mean(acc, na.rm=T)) %>% 
  group_by(subject) %>% 
  mutate(acc = acc-mean(acc)+overall_mean,
         post_bin = as.factor(post_bin)) %>%
  ggplot(aes(post_bin, acc, group=1)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_discrete(name='**Post**-ISI duration (ms.)', labels=c('800-900', '900-1000', '1000-1100', '1100-1200', '1200-1250')) +
  scale_y_continuous('Recall probability')  +
  coord_cartesian(ylim=c(0.53,0.555))  +
  theme(axis.title.x = element_markdown(),
        axis.text.x = element_text(angle = 45,hjust=1))
ggsave('figures/peers_postisi_acc.png', w=3.25, h=2.75)

# main effect of prior freq
(f1 <- exp1 %>% 
    filter(!is.na(pre1_bin), pre1_bin <= 13,post_bin <= 13) %>% 
    ggplot(aes(ave_freq_prioritem1, acc)) +
    stat_summary(geom='point') +
    geom_smooth(method='lm', se=F) +
    scale_x_continuous('Mean log(freq) of preceding item during study', breaks=c(1,1.5,2,2.5,3,3.5,4)) +
    scale_y_continuous('P(recall)'))

# main effect of prior freq by pre_isi
(f2 <- exp1 %>% 
    filter(!is.na(pre1_bin), pre1_bin <= 13,post_bin <= 13) %>% 
    mutate(overall_mean = mean(acc, na.rm=T)) %>% 
    # group_by(pre1_bin) %>% 
    # mutate(acc = acc-mean(acc)+overall_mean) %>% 
    ggplot(aes(ave_freq_prioritem1, acc, color=as.factor(pre1_bin))) +
    stat_summary(geom='point') +
    geom_smooth(method='lm', se=F) +
    scale_x_continuous('Mean log(freq) of preceding item during study', breaks=c(1,1.5,2,2.5,3,3.5,4)) +
    scale_y_continuous('P(recall)') +
    scale_color_discrete(name='Preceding ISI duration', labels=c('800-900 ms.', '900-1000 ms.', '1000-1100 ms.', '1100-1200 ms.', '1200-1250ms.')) +
    coord_cartesian(ylim=c(0.515,0.565)))

ggsave('figures/peers_isi_by_priorfreq.png', w=4.5, h=3)


# interaction of prior freq by pre_isi (main effect of pre_isi removed)
(f3 <- exp1 %>% 
    filter(!is.na(pre1_bin), pre1_bin <= 13,post_bin <= 13) %>% 
    mutate(overall_mean = mean(acc, na.rm=T)) %>% 
    group_by(pre1_bin) %>% 
    mutate(acc = acc-mean(acc)) %>% 
    ggplot(aes(ave_freq_prioritem1, acc, color=as.factor(pre1_bin))) +
    # stat_summary(geom='point') +
    geom_smooth(method='lm', se=F) +
    scale_x_continuous('Mean log(freq) of preceding item during study', breaks=c(1,1.5,2,2.5,3,3.5,4)) +
    scale_y_continuous('P(recall)') +
    scale_color_discrete(name='Preceding ISI duration', labels=c('800-900 ms.', '900-1000 ms.', '1000-1100 ms.', '1100-1200 ms.','1200-1250 ms.')))

ggsave('figures/peers_isi_by_priorfreq_centered.png', w=4.5, h=3)

f1+theme(legend.position='none') + f3 + theme(legend.position=c(0.05,1), legend.justification = c(0.05,1))
ggsave('figures/peers_isi_by_freqprior_both.png', w=7.5, h=3.5)

#### LAG PLOTS #########################

exp1 %>% 
  filter(post_bin <= 12, !is.na(post_bin), study_position >= 4, study_position <= 12) %>% 
  gather(lag,acc,lagm1acc:lag4acc) %>% 
  mutate(lag = case_when(
    lag == 'lagm4acc' ~ -4,
    lag == 'lagm3acc' ~ -3,
    lag == 'lagm2acc' ~ -2,
    lag == 'lagm1acc' ~ -1,
    lag == 'lag1acc' ~ 1,
    lag == 'lag2acc' ~ 2,
    lag == 'lag3acc' ~ 3,
    lag == 'lag4acc' ~ 4
  )) %>% 
  ggplot(aes(lag, study_position, color=as.factor(post_bin))) +
  stat_summary(geom="pointrange") +
  stat_summary(geom="line")
    

# -------------------------------------------------------------------
# models
# -------------------------------------------------------------------


mldata1 <- exp1 %>% 
  filter(!is.na(freq_prioritem1), 
         !is.na(pre1_bin),
         !is.na(postfreq_prioritem), 
         !is.na(post_bin), 
         study_position <= 12, 
         pre1_bin <= 12, 
         post_bin <= 12) %>% 
  mutate(Lg10WF = as.numeric(Lg10WF),
         pre1 = pre1/1000,
         post = post/1000,
         pre1 = pre1-min(pre1),
         post = post-min(post))


gml_fpre <- glmer(acc ~ Lg10WF + postfreq_prioritem + pre1 + post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_fpost <- glmer(acc ~ Lg10WF + freq_prioritem1 + pre1 + post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_isipre <- glmer(acc ~ Lg10WF + freq_prioritem1 + postfreq_prioritem + post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_isipost <- glmer(acc ~ Lg10WF + freq_prioritem1 + postfreq_prioritem + pre1 + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_main <- glmer(acc ~ Lg10WF + freq_prioritem1 + postfreq_prioritem + pre1 + post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_inter <- glmer(acc ~ Lg10WF + freq_prioritem1 * pre1 + postfreq_prioritem + post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_interpost <- glmer(acc ~ Lg10WF + freq_prioritem1 * pre1 + postfreq_prioritem * post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_interpostlg10wf <- glmer(acc ~ Lg10WF*post + freq_prioritem1 * pre1 + postfreq_prioritem + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)
gml_interprelg10wf <- glmer(acc ~ Lg10WF*pre1 + freq_prioritem1 * pre1 + postfreq_prioritem + post + (1|subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)

anova(gml_fpre, gml_main)
anova(gml_fpost, gml_main)
anova(gml_isipre, gml_main)
anova(gml_isipost, gml_main)
anova(gml_main, gml_inter)
anova(gml_inter, gml_interpost)
anova(gml_inter, gml_interpostlg10wf)
anova(gml_inter, gml_interprelg10wf)
summary(gml_main)
summary(gml_inter)
summary(gml_interpost)


cont=matrix(c(0,0,1,0,-1,0,0),ncol=7)
rownames(cont)="freq_prioritem1 - freq_postitem"
library(multcomp)
summary(glht(gml_inter,linfct=cont))



gml_main_slopes <- glmer(acc ~ Lg10WF + freq_prioritem1 + postfreq_prioritem + pre1 + post + (Lg10WF + freq_prioritem1 + postfreq_prioritem + pre1||subject) +(1|word), data=mldata1, family='binomial', nAGQ=0)



#### SUBJECTS WHO SHOW A POSITIVE preFreq do they show a null postFreq and vice versa? #########################

fits <- mldata1 %>% 
  nest(-subject) %>% 
  mutate(fit = map(data, function(x) glm(acc ~Lg10WF + freq_prioritem1+postfreq_prioritem+pre1+post, data=x, family='binomial')),
         coefs = map(fit, function(x) as.data.frame(t(x$coefficients)))) %>% 
  unnest(coefs)
