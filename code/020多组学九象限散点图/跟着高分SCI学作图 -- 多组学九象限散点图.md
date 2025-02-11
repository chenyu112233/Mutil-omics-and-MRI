# 跟着高分SCI学作图 -- 堆积柱状图+折线+热图注释

> 已经付费加群的小伙伴无需二次付费，等待师兄后续更新即可！

![封面](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905224626692.png)

从这个系列开始，师兄就带着大家从各大顶级期刊中的Figuer入手，从仿照别人的作图风格到最后实现自己游刃有余的套用在自己的分析数据上！这一系列绝对是高质量！还不赶紧**点赞+在看**，学起来！

插入公众号

> 本期分享的是期刊：**Nucleic Acids Research（IF = 19.16）**上面一篇文章中的一个**多组学九象限散点图**！
>
> 本系列所有代码和示例数据将会和生信常用图形系列绘图放在一起，公众号右下角添加师兄微信，**付费149元**，即可加入生信绘图交流群。**群内不仅提供生信常用图形系列的代码，还会提供本系列后续所有Figure的实例数据和代码，我会在文章更新后第一时间上传。**
>
> 当然了！如果你还想白嫖，师兄的文章中代码已经写的很清楚了！但是师兄还是希望你点个赞再走呗！
>
> 以上就是本期的全部内容啦！**欢迎点赞，点在看！**师兄会尽快更新哦！制作不易，你的打赏将成为师兄继续更新的十足动力！
>
> **优惠方式：点赞+在看，并转发这两个系列任意一篇文章至朋友圈，集赞30个，即可享受￥119入群！**

![参考文献](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905224824538.png)

话不多说，直接上图！

### 读图

![原图](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905225018809.png)

> 九象限图经常用于多组学比较，两种组学的差异分析结果对比，这张图应该算是我见过的颜值最高的一个，难点主要表现在以下几个方面：
>
> - 散点分组，如何根据FC值对散点进行分组？
> - 最难的地方是，ggplot2如何设置中心为原点的坐标轴？
> 
> 图形整体难度不大，但是样式值得大家学习和应用！



### 效果展示

![效果展示](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905230026057.png)



### 数据构建+分组

```R
######### 多组学散点图（九象限散点图） #############
# 创建示例数据(具体到你们自己的数据就是两种组学差异分析的FC值)
library(MASS)

covariance <- matrix(c(1,0.8,0.8,1),nrow=2,byrow=TRUE)
data <- mvrnorm(n=1000, mu=c(0,0), covariance)
data <- as.data.frame(data)
colnames(data) <- c("mRNA_FC", "RPF_FC")

head(data)
# mRNA_FC     RPF_FC
# [1,] -1.0069132 -1.1906012
# [2,]  1.5264603  0.6807646
# [3,] -0.7447173 -0.6654957
# [4,]  0.9627548  1.2683072
# [5,]  0.8246237  0.6616568
# [6,]  0.9687847  1.2512485

# 分组处理：设置阈值，以正负log10(2)为阈值：
group <- ifelse((abs(data[,1]) > log10(2) & abs(data[,2]) > log10(2))|(abs(data[,1]) < -log10(2) & abs(data[,2]) < -log10(2)),
                "mRNA+RPF_both", ifelse(
                  (data[,1] > log10(2) & data[,2] < log10(2) & data[,2] > -log10(2))|(data[,1] < -log10(2) & data[,2] < log10(2) & data[,2] > -log10(2)),
                "mRNA_only", ifelse(
                  (data[,2] > log10(2) & data[,1] < log10(2) & data[,1] > -log10(2))|(data[,2] < -log10(2) & data[,1] < log10(2) & data[,1] > -log10(2)),
                "RPF_only", NA)))

data$group <- group
```



### 绘图

```R
# 绘图：
library(ggplot2)

# 定义绘制坐标轴函数：
draw_axis_line <- function(length_x, length_y, 
                           tick_step = NULL, lab_step = NULL){
  axis_x_begin <- -1*length_x
  axis_x_end <- length_x
  
  axis_y_begin  <- -1*length_y
  axis_y_end    <- length_y
  
  if (missing(tick_step))
    tick_step <- 1
  
  if (is.null(lab_step))
    lab_step <- 2
  
  # axis ticks data
  tick_x_frame <- data.frame(ticks = seq(axis_x_begin, axis_x_end, 
                                         by = tick_step))
  
  tick_y_frame <-  data.frame(ticks = seq(axis_y_begin, axis_y_end, 
                                          by = tick_step))
  
  # axis labels data
  lab_x_frame <- subset(data.frame(lab = seq(axis_x_begin, axis_x_end, 
                                             by = lab_step), zero = 0), 
                        lab != 0)
  
  lab_y_frame <- subset(data.frame(lab = seq(axis_y_begin, axis_y_end,
                                             by = lab_step),zero = 0), 
                        lab != 0)
  
  # set tick length
  tick_x_length = 0.05
  tick_y_length = 0.05
  
  # set zero point
  
  data <- data.frame(x = 0, y = 0)
  p <- ggplot(data = data) +
    
    # draw axis line
    geom_segment(y = 0, yend = 0, 
                 x = axis_x_begin, 
                 xend = axis_x_end,
                 size = 0.5) + 
    geom_segment(x = 0, xend = 0, 
                 y = axis_y_begin, 
                 yend = axis_y_end,
                 size = 0.5) +
    # x ticks
    geom_segment(data = tick_x_frame, 
                 aes(x = ticks, xend = ticks, 
                     y = 0, yend = 0 - tick_x_length)) +
    # y ticks
    geom_segment(data = tick_y_frame, 
                 aes(x = 0, xend = 0 - tick_y_length, 
                     y = ticks, yend = ticks)) + 
    
    # labels
    geom_text(data=lab_x_frame, aes(x=lab, y=zero, label=lab), vjust = 1.5) +
    geom_text(data=lab_y_frame, aes(x=zero, y=lab, label=lab), hjust= 1.5) +
    theme_minimal()+
    theme(panel.grid = element_blank(),axis.text = element_blank())
  return(p)
}

p <- draw_axis_line(4, 4)

p1 <- p + geom_point(data=data, aes(mRNA_FC, RPF_FC, color = group))+
  scale_color_manual(values = c("mRNA+RPF_both" = "#dd8653", 
                                "mRNA_only" = "#59a5d7", 
                                "RPF_only" = "#aa65a4", 
                                "#878787"),
                     breaks = c("mRNA+RPF_both","mRNA_only","RPF_only"))+
  xlab("mRNA:FC(P42/E15.5)")+
  ylab("RPF:FC(P42/E15.5)")+
  theme(legend.position = "bottom")+
  annotate("text", label = "bolditalic(Brain)", parse = TRUE, 
           x = -2, y = 2, size = 4, colour = "black")+
  guides(color = guide_legend(title = "", ncol = 1, byrow = TRUE))

p1

ggsave("plot.pdf", plot = p1, height = 6, width = 6)
```

![plot1](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905230259100.png)

### 多图拼接

```R
# 多图拼接：这里不再多画了，直接都用p1
p_list <- list(p1=p1, p2=p1, p3=p1)

library(cowplot)

p_all <- plot_grid(plotlist = p_list, ncol = 3)

ggsave("plot2.pdf", plot = p_all, height = 4, width = 9)
```

![效果展示](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905230026057.png)

### 绘图群附加福利

凡是**已经加群**的小伙伴，你们在看文献的时候如果看到好看的Figure，可以发到群里！师兄会及时关注的，**如果被师兄选中，就会在推文中更新！**

#### 今天的图也是来自粉丝推荐

![选中案例](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220905230648635.png)

#### **Figure的要求如下（被选中的前提）：**

- 首先肯定是要符合大众审美的，在图形样式上要过关。
- 新颖独特，有与常见图形不一样的地方。
- 有一定难度，太简单的大家都会做，没什么挑战性哈！

#### 往期选中案例

![往期选中案例](https://picgo-1312459003.cos.ap-shanghai.myqcloud.com/img/image-20220828213019162.png)
