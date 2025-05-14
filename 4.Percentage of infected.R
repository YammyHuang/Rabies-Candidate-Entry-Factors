library(dplyr)
library(ggplot2)
RAAV_RV.INTE<-readRDS("H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/adult.inte.RDS")
Idents(RAAV_RV.INTE) <- 'sample'
Idents(RAAV_RV.INTE) <- factor(Idents(RAAV_RV.INTE), 
                               levels = c("RV_infected", "rAAV_infected", "rAAV_MOCK", "RV_MOCK"))
nMaintype <- RAAV_RV.INTE@meta.data$Maintype
nsample <- RAAV_RV.INTE@meta.data$sample
nMaintype

custom_order <- c('Pvalb','Sst','Vip/Lamp5',"Astro","Microglia",'Oligo','L2/3IT','L4/5IT','L6IT','L5NP','L6CT','L5ET')
nMaintype <- factor(nMaintype, levels = custom_order)
# 转换为data.frame对象
data <- data.frame(nMaintype, nsample)

levels(data)

# 计算每个nMaintype中nsample1和nsample2的数量
summary_data <- data %>%
  group_by(nMaintype, nsample) %>%
  summarise(count = n()) %>%
  ungroup()
summary_data_3 <- summary_data %>%
  filter(nsample %in% c("rAAV_infected", "RV_infected"))
# 计算百分比
percent_data <- summary_data %>%
  group_by(nMaintype) %>%
  mutate(percentage = count / sum(count))

percent_data$nMaintype
library(dplyr)

percent_data_3 <- summary_data_3 %>%
  group_by(nMaintype) %>%
  mutate(percentage = count / sum(count))

#不同亚群占比
p1<- ggplot(percent_data, aes(x = nMaintype, y = percentage, fill = nsample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f8766d","#f8766d","#00bfc4", '#eade5e')) +
  ylab("Percentage") +
  theme_bw() +
  labs(
    title = "Percentage of RV rAAV infected cells",  # 设置标题
    x = "Cluster Type",  # 设置x轴标签
    y = "Percentage of virus infected cells"
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),  # 设置标题大小为20，居中对齐
    axis.title.x = element_text(size = 14),  # 设置x轴标签大小为14
    axis.title.y = element_text(size = 14),  # 设置y轴标签大小为14
    axis.text.x = element_text(angle = 45, hjust = 1,size = 14,color ='black'),  # 横坐标倾斜45度
    legend.text = element_text(size = 14))

# 过滤掉 "rAAV_MOCK" 和 "RV_MOCK"
filtered_data <- data %>%
  filter(nsample %in% c("RV_infected", "rAAV_infected")) %>%  # 只保留 RV_infected 和 rAAV_infected
  filter(!nMaintype %in% c("Pvalb", "Sst", "Vip/Lamp5"))  # 去掉这几类细胞


# 重新绘制柱状图
p2<- ggplot(filtered_data, aes(x = nMaintype, fill = nsample)) +
  geom_bar() +
  scale_fill_manual(values = c("#00bfc4", '#eade5e')) +  # 仅保留 RV_infected 和 rAAV_infected 的颜色
  theme_bw() +
  labs(
    title = "Number of Virus infected cells",
    x = " ",
    y = "Cell Numbers"
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color ='black'),
    legend.text = element_text(size = 14)
  )

print(p2)
#Virus Percentage===RV/rAAV
# 计算每个nMaintype中nsample1和nsample2的数量
summary_data <- filtered_data%>%
  group_by(nMaintype, nsample) %>%
  summarise(count = n()) %>%
  ungroup()
summary_data_3 <- summary_data %>%
  filter(nsample %in% c("rAAV_infected", "RV_infected"))

percent_data_3 <- summary_data_3 %>%
  group_by(nMaintype) %>%
  mutate(percentage = count / sum(count))

p3 <- ggplot(percent_data_3, aes(x = nMaintype, y = percentage, fill = nsample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#00bfc4", '#eade5e')) +
  ylab("Percentage") +
  theme_bw() +
  labs(
    title = "Percentage of RV rAAV infected cells",  # 设置标题
    x = "Cluster Type",  # 设置x轴标签
    y = "Percentage of virus infected cells"
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),  # 设置标题大小为20，居中对齐
    axis.title.x = element_text(size = 14),  # 设置x轴标签大小为14
    axis.title.y = element_text(size = 14),  # 设置y轴标签大小为14
    axis.text.x = element_text(angle = 45, hjust = 1,size = 14,color ='black')  # 横坐标倾斜45度
  )

p3


p1
p2
p3
ggsave("H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/p6.Percentage of virus infected cells.pdf", 
       plot = p1, width =7, height = 8, dpi = 300)
ggsave("H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/p7.Percentage of virus infected cells.pdf", 
       plot = p2, width =6, height = 8, dpi = 300)
ggsave("H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/p8.Percentage of virus infected cells.pdf", 
       plot = p3, width =6, height = 8, dpi = 300)
