### CNV analysis ####

library("readxl")
library("writexl")
library(tidyr)
library(naniar)

## open data -------------------------------------------------------------------------
CNV <- read.xlsx("./input/CNV analysis.xlsx")

## plots Chromosomes

CNV %>%
  group_by(location, Chr1) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr1),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr1.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr2) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr2),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr2.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr3) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr3),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr3.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr4) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr4),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr4.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr5) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr5),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr5.pdf", path= "./plots/")
  
CNV %>%
  group_by(location, Chr6) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr6),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr6.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr7) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr7),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr7.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr8) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr8),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr8.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr9) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr9),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr9.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr10) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr10),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr10.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr11) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr11),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr11.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr12) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr12),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr12.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr13) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr13),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr13.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr14) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr14),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr14.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr15) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr15),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr15.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr16) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr16),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr16.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr17) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr17),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr17.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr18) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr18),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr18.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr19) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr19),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr19.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr20) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr20),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr20.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr21) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr21),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr21.pdf", path= "./plots/")

CNV %>%
  group_by(location, Chr22) %>%
  summarise(n()) %>%
  rename("count" = "n()") %>%
  ggplot(aes(y=count, x=location, fill=location)) +
  geom_bar(
    aes(fill = Chr22),
    stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle=90))
ggsave("Chr22.pdf", path= "./plots/")