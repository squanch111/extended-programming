# 目前已经做完了基础代码1~5的全部功能，接下来还需要实现1.按要求输出数据 2.作图 3.bootstrapping approach不知道是干嘛的

# 加载数据
# data <- read.table("F:/R-file/engcov.txt") # 请确保将数据文件放在工作目录，取消注释这行代码
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
deconv_optimized_plot <- function(t, deaths, n.rep=100, steps=c(8,4,2,1), bs=False, t0=NULL) {
  #按照pdf中要求，t是有死亡记录的天数，deaths是每天死了多少人，t0如果为空使用前面2行代码再随机生成一个,bs是后面要用的检验
  if (is.null(t0)) {
    t0 <- rep(t, deaths)  # 初始猜测感染时间
  }
  
  #这些代码是之前作图需要的变量
  # P_history <- numeric(length(steps))  # 存储每次更新后的P值
  # inft <- matrix(0, 310, length(steps))  # 存储每次更新后的感染率
  
  # 执行多次迭代，重复4、5步骤100次
  for (rep in 1:n.rep) {
    # 随机打乱 t0 的顺序
    indices <- sample(seq_along(t0))  # 随机化索引顺序，生成长度为t0，顺序为随机乱序的序号向量
    
    # 再生成一个新的感染到死亡时间向量，再算出模拟的死亡时间simulated_deaths
    infection_to_death <- sample(1:80, n, replace=TRUE, prob=duration_probs)
    simulated_deaths <- t0 + infection_to_death
    
    # 计算初始的P值
    #dsi是模拟的每天死了多少人，di是实际上每天死了多少人
    dsi <- tabulate(simulated_deaths, nbins=310)
    di <- tabulate(death, nbins=310)
    p <- sum((di - dsi)^2 / pmax(1, dsi))
    
    #对p的计算作了更改，加了一个for循环，严格按照要求来 （这个功能先舍弃了）
    # p <- 0
    # for (i in 1:length(dsi)){
    #   p <- p + (di[i] - dsi[i])^2 / max(1, dsi[i])
    # }
    
    #声明一个长度为n的随机步长向量，用循环中每一个词加上该向量中的步长，避免循环中调用sample语句
    step_size <- sample(c(-4, -2, -1, 1, 2, 4), n, replace = TRUE)
    
    # 遍历 t0 中的每一个元素并应用步长
    for (i in indices) {
      
      # 更新 t0[i] 并确保值在合理范围内
      proposed_t0 <- t0  # 创建 t0 副本
      proposed_t0[i] <- pmin(pmax(t0[i] + step_size[i], 1), 310)  # 更新单个 t0[i]，确保不越界
      
      #为了避免在循环中调用tabulate语句，在原先的dsi中使当前的第i天的人数减一，再使移动后的那一天的人数加一，作为新的dsi
      former_death_day <- simulated_deaths[i]  # 更新前的模拟死亡日期
      current_death_day <- proposed_t0[i] + infection_to_death[i]  # 更新后的死亡日期
      
      proposed_dsi <- dsi  # 创建 dsi 副本
      proposed_dsi[former_death_day] <- proposed_dsi[former_death_day] - 1  # 减少原死亡日期的人数
      proposed_dsi[current_death_day] <- proposed_dsi[current_death_day] + 1  # 增加新死亡日期的人数
      
      # 计算新的 P 值
      new_p <- sum((di - proposed_dsi)^2 / pmax(1, proposed_dsi))
      
      # 如果新的 P 值更小，则更新 t0[i] 和 P
      if (new_p < p) {
        t0[i] <- proposed_t0[i]  # 接受新的 t0 值
        p <- new_p  # 更新 P 值
        dsi <- proposed_dsi  # 更新死亡人数分布
      }
      # 否则保留原值，什么也不做
    }
    
    # 每次迭代后可以输出 P 值和当前 t0 的状态用于检查
    cat("Iteration:", rep, "p value:", p, "\n")
  }
  
  #准备输出，前面的有些我删去了，麻烦再对照pdf中的要求重新声明需要输出的数据，并准备绘图
  #return(list(P=P_history, inft=inft, t0=t0))
  
}
# 运行函数
result_optimized <- deconv_optimized_plot(data$julian, data$nhs)




