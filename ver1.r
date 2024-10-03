# Team Members:
# - Member1 ()
# - Member2 ()
# - Member3 (s2663848)


# Description of Contributions:
# - Member1: Implemented file reading and tokenization (30%)
# - Member2: Developed Markov model and simulation logic (40%)
# - Member3: Handled pre-processing, punctuation splitting, and code documentation (30%)


#resource preparation
#setwd("C:/Users/Steven/Desktop/R/test1/R-test1") ## comment out of submitted
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8")
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("


#split punction
## 加入空格
## 对空格进行分割
split_punct = function(words,punct){
  punct_escaped <- gsub("([\\.^$|?*+(){}\\[\\]])", "\\\\\\1", punct) 
  words <- gsub(paste0("([", punct_escaped, "])"), " \\1", words)
  split_words <- unlist(strsplit(words, " "))
  split_words <- split_words[split_words != ""]
  return(split_words)
}
output <- split_punct(a, ",.?!;:-_")


# Convert text to lowercase and get the vector of unique words
output <- tolower(output)
unique_words <- unique(output)

# Create a vector of indices for each word in the text
index_vector <- match(output, unique_words)

# Count the frequency of each unique word using the tabulate function
word_freq <- tabulate(index_vector)

# Determine the threshold to keep ~1000 most common words
m_threshold <- sort(word_freq, decreasing = TRUE)[1000]
common_words <- unique_words[word_freq >= m_threshold]


# Create a matrix of word sequences for the Markov model
mlag <- 4  # max lag (you can adjust this)
n <- length(output)
M <- matrix(NA, nrow = n - mlag, ncol = mlag + 1)

for (i in 1:(mlag + 1)) {
  M[, i] <- c(rep(NA, mlag - i + 1 ), index_vector[1:(n + i - 1 - 2 * mlag)])
}


# Function to simulate text using the Markov model
simulate_text <- function(M, common_words, nw = 50, mlag = 4) {
  text <- c()
  first_word <- sample(common_words, 1)
  text <- append(text, first_word)
  
  # 确保 token_seq 和 common_words 都是字符向量
  common_words <- as.character(common_words)
  # 去除前后空格
  common_words <- trimws(common_words)
  
  for (i in 2:nw) {
    for (j in mlag:1) {
      if (length(text) >= j) {
        token_seq <- tail(text, j)  # 取最后 j 个单词
      } else {
        token_seq <- text  # 如果生成的单词还不足 j 个，使用所有已生成的单词
      }
      
      #测试部分
      print(j)
      print("text:")
      print(text)
      token_seq <- trimws(token_seq)
      token_seq <- as.character(token_seq)
      token_seq_numeric <- match(token_seq, common_words)
      print("token_seq_num")
      print(token_seq_numeric)
      
      # Find matching rows in M（弃用版本）
      #matching_rows <- which(rowSums(M[, 1:j, drop = FALSE] == token_seq & !is.na(M[,1:j,drop = FALSE])) == j)
      #matching_rows <- which(rowSums(M[, 1:j, drop = FALSE] == token_seq) == j)
      #matching_rows <- which(apply(M[, 1:j, drop = FALSE], 1, function(row) all(row == token_seq, na.rm = FALSE)))
      
      matching_rows <- c()

      for (row_index in 5:nrow(M)) {
        current_row <- M[row_index, 1:j]

        # 如果当前行与 token_seq 匹配
        if (all(current_row == token_seq_numeric, na.rm = TRUE)) {
          matching_rows <- c(matching_rows, row_index)  # 添加匹配行的索引
        }
      }
      
      #测试
      print("Current token_seq:")
      print(token_seq) 
      found_match <- TRUE
      
      
      #input new word
      if (length(matching_rows) > 0) {
        next_word <- na.omit(sample(common_words[M[matching_rows, j+1]], 1))
        text <- append(text, next_word)
        found_match <- FALSE
        break
      }
      }
  }
  
  
  # 如果没有找到匹配的行，则随机选择一个单词
  if (found_match) {
    next_word <- sample(common_words, 1)  # 随机选择一个单词
    text <- append(text, next_word)
    print("No matching rows found, selecting a random word")
  }
  
  # Return generated text
  return(paste(text, collapse = " "))
}

# Simulate 50-word section from the Markov model
simulated_text <- simulate_text(M, common_words, nw = 50, mlag = 4)
cat(simulated_text)

# Generate text based on word frequencies alone for comparison
common_word_freq <- word_freq[match(common_words, unique_words)]

# 计算每个单词的出现概率 (基于其频率)
word_prob <- common_word_freq / sum(common_word_freq)

# 使用 word_prob 进行加权抽样生成 nw 个单词
generated_text <- sample(common_words, size = 50, replace = TRUE, prob = word_prob)

# 返回生成的文本
cat(paste(generated_text, collapse = " "))
