# 目前已经做完了基础代码1~5的全部功能，接下来还需要实现1.按要求输出数据 2.作图 3.bootstrapping approach不知道是干嘛的




# 加载数据
#data <- read.table("F:/R-file/engcov.txt") # 请确保将数据文件放在工作目录，取消注释这行代码
data <- read.table("D:/rstudio/engcov.txt") 
data <- data[1:150, ]  # 只使用前150行

# 参数设置
meanlog <- 3.152
sdlog <- 0.451
lockdown_day <- 84  # 英国封锁的日期

# 1. 初始化感染到死亡时间的概率分布
duration_probs <- dlnorm(1:80, meanlog, sdlog)
duration_probs <- duration_probs / sum(duration_probs)  # 标准化

# 2. 初始感染时间猜测
# death 生成一个向量，显示的是每个死人是哪一天去世的，共29442个死人
death <- rep(data$julian, data$nhs)  # 将死亡人数展开成一个个个体
n <- 29442 # 总死亡人数

# 计算初始t0, t0是预计每个患者的初始感染时间
infection_to_death <- sample(1:80, n, replace=TRUE, prob=duration_probs)
t0 <- death - infection_to_death

# 3. 模拟及优化过程的核心函数 （这里的参数设置有问题，没有严格按照pdf要求来搞）
deconv_optimized_plot <- function(t, deaths, n.rep=100, bs=FALSE, t0=NULL) {
  
  #按照pdf中要求，t是有死亡记录的天数，deaths是每天死了多少人，t0如果为空使用前面2行代码再随机生成一个,bs是后面要用的检验
  colors <- rainbow(n.rep)
  #这些代码是之前作图需要的变量
  P_history <- numeric(n.rep)  # 保存每次迭代后的 P 值
  inft <- matrix(0, nrow=310, ncol=n.rep)  # 保存每日感染人数的历史
  # 创建初始图像，绘制实际的死亡人数曲线
  if (bs == FALSE){
    di <- tabulate(death, nbins=310)
    plot(data$julian, data$nhs, type='h', col="blue", lwd=2, xlab="Day", ylab="Deaths",
       main="Real and Simulated Deaths Over Iterations")
  #legend("topright", legend=c("Actual Deaths"), col=c("red"), lty=1, lwd=2)
  }
   if (is.null(t0)) {
    t0 <- rep(t, deaths)  # 初始猜测感染时间
  }
  
  # 执行多次迭代，重复4、5步骤100次
  for (rep in 1:n.rep) {
    # 随机打乱 t0 的顺序
    indices <- sample(seq_along(t0))  # 随机化索引顺序，生成长度为t0，顺序为随机乱序的序号向量
    
    # 再生成一个新的感染到死亡时间向量，再算出模拟的死亡时间simulated_deaths
    infection_to_death <- sample(1:80, n, replace=TRUE, prob=duration_probs)
    simulated_deaths <- t0 + infection_to_death
    # 确保模拟的死亡天数在 1 到 310 之间
    t0 <- pmin(pmax(t0, 1), 310)
    
    # 计算初始的P值
    #dsi是模拟的每天死了多少人，di是实际上每天死了多少人
    dsi <- tabulate(simulated_deaths, nbins=310)
    
    # 添加bootstrap模块
    if (bs){
      #di = rpois(length(di), lambda = di)
      #di <- smooth.spline(1:length(di), di)$y # 平滑处理
      death_bs <- rpois(n, death)
      
      # 筛选出在区间 [a, b] 内的随机数
      filtered_values <- death_bs[death_bs >= 1 & death_bs <= 310]
      while(length(filtered_values) < n) {
        new_values <- rpois(n, death)
        filtered_values <- c(filtered_values, new_values[new_values >= 1 & new_values <= 310])
      }
      death_bs <- filtered_values[1:n]  # 保留前 n 个符合条件的随机数
      
      # 查看结果
      di <- tabulate(death_bs, nbins=310)
      
      if(rep == 1){
      plot(1:310, di, type='h', col=rgb(0.4,0.1,0,0.1), lwd=2, xlab="Day", ylab="Deaths",
           main="Real and Simulated Deaths Over Iterations////")
      }
    }
    
    p <- sum((di - dsi)^2 / pmax(1, dsi),na.rm=TRUE)
     
    #对p的计算作了更改，加了一个for循环，严格按照要求来 （这个功能先舍弃了）
    # p <- 0
    # for (i in 1:length(dsi)){
    #   p <- p + (di[i] - dsi[i])^2 / max(1, dsi[i])
    # }
    
    #声明一个长度为n的随机步长向量，用循环中每一个词加上该向量中的步长，避免循环中调用sample语句
    if (rep > 75) {
      step_size <- sample(c(-2, -1, 1, 2), length(t0), replace=TRUE)
    } else if (rep > 50) {
      step_size <- sample(c(-4, -2, -1, 1, 2, 4), length(t0), replace=TRUE)
    } else {
      step_size <- sample(c(-8, -4, -2, -1, 1, 2, 4, 8), length(t0), replace=TRUE)
    }
 
    # 遍历 t0 中的每一个元素并应用步长
    for (i in indices) {
      stepp <- step_size[i]
      former_death_day <- simulated_deaths[i]
      current_death_day <-  pmin(pmax(t0[i] + stepp, 1), 310) + infection_to_death[i]
      # 使用更简单的更新方式
      dsi[former_death_day] <- dsi[former_death_day] - 1
      dsi[current_death_day] <- dsi[current_death_day] + 1

      new_p <- sum((di - dsi)^2 / pmax(1, dsi),na.rm=TRUE)

      if (new_p < p) {
        t0[i] <- pmin(pmax(t0[i] + stepp, 1), 310)
        p <- new_p
      } else {
        dsi[former_death_day] <- dsi[former_death_day] + 1
        dsi[current_death_day] <- dsi[current_death_day] - 1
      }

      
      # # 更新 t0[i] 并确保值在合理范围内
      # proposed_t0 <- t0  # 创建 t0 副本
      # proposed_t0[i] <- pmin(pmax(t0[i] + step_size[i], 1), 310)  # 更新单个 t0[i]，确保不越界
      # 
      # #为了避免在循环中调用tabulate语句，在原先的dsi中使当前的第i天的人数减一，再使移动后的那一天的人数加一，作为新的dsi
      # former_death_day <- simulated_deaths[i]  # 更新前的模拟死亡日期
      # current_death_day <- proposed_t0[i] + infection_to_death[i]  # 更新后的死亡日期
      # 
      # proposed_dsi <- dsi  # 创建 dsi 副本
      # proposed_dsi[former_death_day] <- proposed_dsi[former_death_day] - 1  # 减少原死亡日期的人数
      # proposed_dsi[current_death_day] <- proposed_dsi[current_death_day] + 1  # 增加新死亡日期的人数
      # 
      # # 计算新的 P 值
      # new_p <- sum((di - proposed_dsi)^2 / pmax(1, proposed_dsi))
      # 
      # # 如果新的 P 值更小，则更新 t0[i] 和 P
      # if (new_p < p) {
      #   t0[i] <- proposed_t0[i]  # 接受新的 t0 值
      #   p <- new_p  # 更新 P 值
      #   dsi <- proposed_dsi  # 更新死亡人数分布
      # }
      # # 否则保留原值，什么也不做
    }
    
    P_history[rep] <- p
    inft[, rep] <- tabulate(t0, nbins=310)
    cat("Iteration:", rep, "p value:", p, "\n")
    # 每次迭代后用不同颜色叠加模拟的死亡数据
    if (bs){
      lines(1:310, dsi[1:310], col="grey", lwd=1)
    }
    else{
    lines(1:211, dsi[1:211], col=colors[rep], lwd=1)  # 使用不同颜色绘制模拟死亡数据
    }
    
  }
  print(P_history)
  print(t0)
  
    # 返回结果
  return(list(P=P_history, inft=inft, t0=t0))
  
}

# 运行函数
# 算时间
start_time = Sys.time()
result_optimized <- deconv_optimized_plot(data$julian, data$nhs, bs=FALSE)
end_time = Sys.time()
print(end_time-start_time)
# 提取实际死亡数据和最终模拟死亡数据
#deaths <- data$nhs[1:310]  # 实际死亡数据
final_simulated_deaths <- result_optimized$inft[, ncol(result_optimized$inft)]  # 模拟的最终感染率
# 绘制图像，实际死亡数据和模拟数据

# plot(1:310, final_simulated_deaths[1:310],type='l' ,col="red", lwd=2, xlab="Day", ylab="Deaths/Incidence",
#      main="Actual Deaths vs Simulated Incidence")  #模拟的感染曲线
# 
# lines(data$julian, data$nhs, col="blue", lwd=2 ) # 真实死亡曲线
# # lines(result_optimized$dsi,col="pink")  #模拟的死亡曲线
# 
# # 标出封锁的第一天（第84天）
# abline(v=84, col="black", lty=2,lwd=2)
# 
# # 添加图例
# legend("topright", legend=c("Simulated Incidence", "Actual Deaths", "Lockdown Start"),
#        col=c("red", "blue", "black"), lty=c(1,1,2), lwd=2)






###############################################
# p和t0收敛后进行泊松分布模拟（将已收敛的t0代入）
bootstrap_result <- deconv_optimized_plot(data$julian, data$nhs, bs=TRUE, t0=result_optimized$t0)
final_simulated_deaths_bs <- bootstrap_result$inft[, ncol(bootstrap_result$inft)]

# 画图
plot(1:310, final_simulated_deaths[1:310],type='l' ,col="red", lwd=2, xlab="Day", ylab="Number",
     main="Actual Deaths vs Simulated Incidence")  #模拟的感染曲线

lines(data$julian, data$nhs, col="blue", lwd=2 ) # 真实死亡曲线
lines(1:310, final_simulated_deaths_bs[1:310],type='l' ,col="green", lwd=2)




# 绘制图像，实际死亡数据和模拟数据
# plot(1:310, result[1:310], type='l', col="red", lwd=2, 
#      xlab="Day", ylab="Deaths/Incidence", main="Actual Deaths vs Simulated Incidence")

# 标出封锁的第一天（第84天）
abline(v=84, col="black", lty=2, lwd=2)

# 添加图例
# legend("topright", 
#        legend=c("Simulated Incidence", "Actual Deaths", "Lockdown Start", "95% Confidence Interval"),
#        col=c("red", "blue", "black", NA), lty=c(1, 1, 2, NA), lwd=c(2, 2, 2, NA),
#        fill=c(NA, NA, NA, rgb(0.1, 0.1, 0.9, 0.15)), border=NA, bty="n")


# infection_rate_se <- sqrt(final_simulated_deaths_bs)
# 
# #  infection_rate ± 1.96 * 标准误差
# infection_rate_lower <- final_simulated_deaths_bs - 1.96 * infection_rate_se
# infection_rate_upper <- final_simulated_deaths_bs + 1.96 * infection_rate_se
# 
# # 确保下限不小于 0 （感染人数不能为负数）
# infection_rate_lower[infection_rate_lower < 0] <- 0

# 使用 polygon 添加置信区间带
# polygon(c(1:length(final_simulated_deaths_bs), rev(1:length(final_simulated_deaths_bs))),  # X 轴坐标
#         c(infection_rate_lower, rev(infection_rate_upper)),  # Y 轴下限和上限
#         col=rgb(0.7, 0.2, 0.4, 0.5), border=NA)  # 设置颜色和透明度


# lines(result_optimized$dsi,col="pink")  #模拟的死亡曲线

# 添加图例
legend("topright", legend=c("Simulated Incidence", "Actual Deaths", "Lockdown Start", "95% Confidence Interval"),
       col=c("red", "blue", "black", rgb(0.7, 0.2, 0.4, 0.2)), lty=c(1,1,2,1),
       lwd=c(2,2,2,5),
       border=NA,
       y.intersp=1.2, # 设置图例的垂直间距
       box.lwd=0)




