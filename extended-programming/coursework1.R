## reprocessing 
setwd("/Users/huijiezhang/Desktop/extended-programming") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("
#####
library(stringr)
split_punct = function(words,punct){
  punct_escaped <- gsub("([\\.^$|?*+(){}\\[\\]])", "\\\\\\1", punct) 
  words <- gsub(paste0("([", punct_escaped, "])"), " \\1", words)
  split_words <- unlist(strsplit(words, " "))
  split_words <- split_words[split_words != ""]
  return(split_words)
}
a2 <- split_punct(a, ",.?!;:-_")
##改成小写
lower_a = tolower(a2)
dictionary = unique(lower_a)

index = match(lower_a,dictionary)

frequency = tabulate(index)

## take the first 1000 words
sort_index = order(frequency,decreasing = TRUE)
sort_freq = frequency[sort_index]
sort_words = dictionary[sort_index]
threshold = sort_freq[1000]
## 全是小写格式
b = sort_words[sort_freq>=threshold]
# b的索引
index_vector = match(lower_a,b)

freq = tabulate(index_vector)

############# #############
mlag = 4
n = length(index_vector)
M = matrix(NA,nrow = n-mlag,ncol = mlag+1)
for (i in 1:(mlag + 1)) {
  M[, i] <- index_vector[i:(n - mlag + i - 1)]
}
############## M 已生成################
# 判断是否是字母
is_word <- function(word) {
  grepl("^[a-z]+$", word)
}
# 判断是否是标点符号
is_punctuation <- function(word) {
  grepl("^[^a-z]+$", word)
}

#生成函数，避免连续选择相同类型单词#####
text_simulator <- function(b, M, mlag, nw) {
  generator <- vector("character", nw)
  
  # 选择第一个单词索引，并从 b 中提取相应的单词
  valid_first_words <- na.omit(M[, 1])
  
  # 检查 valid_first_words 是否为空
  if (length(valid_first_words) > 0) {
    first_word_index <- sample(valid_first_words, 1)
  } else {
    first_word_index <- sample(seq_along(b), 1)  # 从 b 中随机选择索引
  }
  
  generator[1] <- b[first_word_index]
  
  for (i in 2:nw) {
    word_found <- FALSE
    
    # 内层循环，逐渐减少延迟（从 mlag 到 1）
    for (j in mlag:1) {
      if (i > j) {
        # 匹配当前序列生成的单词索引
        current_seq <- match(generator[(i - j):(i - 1)], b)
        
        # 使用 drop=FALSE 保证结果为矩阵，并找到匹配的行
        match_rows <- which(rowSums(M[, 1:j, drop = FALSE] == current_seq) == j)
        
        # 检查 match_rows 是否为空
        if (length(match_rows) > 0) {
          valid_next_words <- na.omit(M[match_rows, j + 1])
          
          # 避免选择与上一个单词相同的单词
          valid_next_words <- valid_next_words[b[valid_next_words] != generator[i - 1]]
          
          # 避免连续生成标点符号
          if (is_punctuation(generator[i - 1])) {
            valid_next_words <- valid_next_words[is_word(b[valid_next_words])]
          } 
          # 检查 valid_next_words 是否为空
          if (length(valid_next_words) > 0) {
            next_word_index <- sample(valid_next_words, 1)
            generator[i] <- b[next_word_index]
            word_found <- TRUE
            break
          }
        }
      }
    }
    
    # 如果没有找到符合条件的单词，则从 b 中随机选择一个单词
    if (!word_found) {
      # 检查 b 中是否有有效单词
      valid_b <- na.omit(b)
      if (length(valid_b) > 0) {
        generator[i] <- sample(valid_b, 1)
      } else {
        stop("Error: no valid words found in 'b'.")
      }
    }
  }
  # 格式化输出并去除标点符号前的空格
  formatted_output <- generator[1]
  
  for (i in 2:nw) {
    formatted_word <- generator[i]
    
    if (is_punctuation(formatted_word)) {
      # 标点符号前不加空格，直接拼接
      formatted_output <- paste0(formatted_output, formatted_word)
    } else {
      # 普通单词前加空格
      formatted_output <- paste(formatted_output, formatted_word)
    }
  }
  # 使用 cat 打印生成的文本
  cat(formatted_output,"\n")
  
}


#Q10#########
# 寻找大写词库
indicies = match(b,lower_a)
b2 = a2[indicies]
# 判断是否是大写字母开头的单词
is_capitalized <- function(word) {
  grepl("^[A-Z]", word)
}
 mark_capitalized = function(b2){
  
   # 忽视单词大小写，统计在全文中总出现次数
   lowercase_b <- tolower(b)
   total_count <- table(lowercase_b)
   
   # 统计大写字母开头的单词次数
   capitalized_words <- b[is_capitalized(b)]
   capitalized_count <- table(tolower(capitalized_words))
   
   # 创建一个新的向量，用来存储修改后的词库
   new_b <- b2
   
   # 标记那些大写出现次数占比 50% 或以上的单词
   for (word in names(total_count)) {
     total_freq <- total_count[word]  # 该单词的总出现次数
     capitalized_freq <- capitalized_count[word]  # 该单词首字母大写的次数
     
     # 处理 NA 值，假设如果该单词没有大写形式，则 capitalized_freq 为 0
     if (is.na(capitalized_freq)) {
       capitalized_freq <- 0
     }
     # 处理 total_freq 的 NA 情况
     if (!is.na(total_freq) && total_freq > 0 && capitalized_freq / total_freq > 0.5) {
       # 将词库 b2 中符合条件的单词转换为大写
       new_b[lowercase_b == word] <- str_to_title(word)
     }
   }

   return(new_b)  # 返回更新后的词库 b
 }
   new_b = mark_capitalized(b2)
   # 调用函数生成文本
 text_simulator(new_b,M,mlag = 4,nw = 100)
 
