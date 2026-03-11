test<- readRDS("pairwise_sharing_m.s.RDS")
row_dist <- dist(test)  
col_dist <- dist(t(test))  
row_clust <- hclust(row_dist)
col_clust <- hclust(col_dist)  

# 按聚类顺序重新排列矩阵
reordered_test <- test[row_clust$order, col_clust$order]

# 设置上三角为 NA
lower_tri <- reordered_test
diag(lower_tri) <- 0
lower_tri[upper.tri(lower_tri)] <- NA

# 自定义颜色
heatmap_colors <- colorRampPalette(c("white", "firebrick3"))(50)

# 绘制热图
library(pheatmap)

# 设置颜色断点
breaks <- seq(0, 1, length.out = 100)  # 数据范围为 0 到 1
heatmap_colors <- c(
  colorRampPalette(c("white", "red"))(100)   # 白色到红色
)

# 生成热图
pdf("heatmap.pdf", width = 12, height = 10)
pheatmap(lower_tri,
         cluster_rows = FALSE,  # 关闭行聚类
         cluster_cols = FALSE,  # 关闭列聚类
         color = heatmap_colors,  # 设置颜色
         breaks = breaks,  # 使用自定义颜色断点
         na_col = "transparent",  # 设置 NA 为透明
         border_color = NA,  # 去掉网格线
         main = "Lower Triangle Heatmap",
         labels_row = rownames(lower_tri),
         angle_col = 45)  # 设置标题
dev.off()