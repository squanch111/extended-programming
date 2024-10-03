# Team Members:
# - Member1 (S)
# - Member2 (S)
# - Member3 (S)

# Member 1 was primarily responsible for the split_function and simulating the text.
# Member 2 was primarily responsible for generating common_words and matrix M, as well as developing the Markov model for the simulation.
# Member 3 was primarily responsible for identifying capital letters in the final step and testing and debugging the overall code.


#######################################################################
# resource preparation
setwd("D:/rstudio/clone/extended-programming")
a <- scan("4300-0.txt",what="character",skip=73,nlines=32858-73,
          fileEncoding="UTF-8") 
a <- gsub("_(","",a,fixed=TRUE) ## remove "_("


# split punctuation
split_punct = function(words,punct){
  # add spaces around punctuation
  punct_escaped <- gsub("([\\.^$|?*+(){}\\[\\]])", "\\\\\\1", punct) 
  words <- gsub(paste0("([", punct_escaped, "])"), " \\1", words) 
  # then split the string by spaces and remove empty elements
  split_words <- unlist(strsplit(words, " "))
  split_words <- split_words[split_words != ""]
  return(split_words)
}
a_split <- split_punct(a, ",.?!;:-_")


# convert all words to lowercase 
a_lower = tolower(a_split)
# get unique words to form a dictionary, and find their position
b = unique(a_lower)
index = match(a_lower,b)
# count the frequency of each word
frequency = tabulate(index)
# sort words by frequency in descending order
sort_index = order(frequency,decreasing = TRUE)
sort_freq = frequency[sort_index]
sort_words = dictionary[sort_index]
# take the first 1000 most common words
threshold = sort_freq[1000]
b = sort_words[sort_freq>=threshold]
index_vector = match(a_lower,b)
freq = tabulate(index_vector)

# build matrix M to store sequences of words
mlag = 4
n = length(index_vector)
M = matrix(NA,nrow = n-mlag,ncol = mlag+1)
for (i in 1:(mlag + 1)) {
  M[, i] <- index_vector[i:(n - mlag + i - 1)]
}


# # 判断是否是字母
# is_word <- function(word) {
#   grepl("^[a-z]+$", word)
# }

# check whether word is punctuation
is_punctuation <- function(word) {
  grepl("^[^a-z]+$", word)
}


#simulate function
text_simulator <- function(b, M, mlag, nw) {
  # choose a non-empty word to be the first word randomly
  generator <- c()
  current_word <- sample(na.omit(M[, 1]), 1)  # 随机选择起始词
  generator <- c(generator, b[current_word])
  
  
  # set a loop to generate remaining words
  for (i in 2:nw) {
    word_found <- FALSE
    # gradually decrease the lag to find matching words in M
    for (j in mlag:1) {
      if (i > j) {
        # find matching rows of current sequences
        current_seq <- match(generator[(i - j):(i - 1)], b)
        match_rows <- which(rowSums(M[, 1:j, drop = FALSE] == current_seq) == j)
        
        # generate the next word
        if (length(match_rows) > 0) {
          valid_next_words <- na.omit(M[match_rows, j + 1])
          
          # avoid choosing words that are the same as the previous word
          valid_next_words <- valid_next_words[b[valid_next_words] != generator[i - 1]]
          
          # avoid continuous generation of punctuation                                                                                            
          if (is_punctuation(generator[i - 1])) {
            valid_next_words <- valid_next_words[grepl("^[a-z]+$", b[valid_next_words])]
          } 
          
          # choose next words
          if (length(valid_next_words) > 0) {
            next_word_index <- sample(valid_next_words, 1)
            generator[i] <- b[next_word_index]
            word_found <- TRUE
            break
          }
        }
      }
    }
    
    # if can't find next word, choose it in b randomly 
    if (!word_found) 
      generator[i] <- sample(na.omit(b), 1)
  }
  
  
  # remove space before punctuation
  formatted_output <- generator[1]
  for (i in 2:nw) {
    formatted_word <- generator[i]
    # if word is punctuation, paste directly; if not, add space before it
    if (is_punctuation(formatted_word))
      formatted_output <- paste0(formatted_output, formatted_word)
    else
      formatted_output <- paste(formatted_output, formatted_word)
  }

  # print text
  cat("simulated text:","\n",formatted_output,"\n")
  
}


#Q10#########
# search capital words
indicies = match(b,a_lower)
b2 = a_split[indicies]

# check whether word starts with capital letter 
is_capitalized <- function(word) {
  grepl("^[A-Z]", word)
}


 # change b to new_b which remains capital format of often capital words
 mark_capitalized = function(b2){
   # count appearance number of capital words
   lowercase_b <- tolower(b)
   total_count <- table(lowercase_b)
   capitalized_words <- b[is_capitalized(b)]
   capitalized_count <- table(tolower(capitalized_words))
   new_b <- b2
   
   # handle words that have a capitalization rate of more than 50%
   for (word in names(total_count)) {
     total_freq <- total_count[word] 
     capitalized_freq <- capitalized_count[word]  
     
     # handle NA
     if (is.na(capitalized_freq)) {
       capitalized_freq <- 0
     }
     
     
     if (!is.na(total_freq) && total_freq > 0 && capitalized_freq / total_freq > 0.5) {
       # convert eligible words into capital format
       new_b[lowercase_b == word] <- str_to_title(word)
     }
   }

   return(new_b)  
 }
 
new_b = mark_capitalized(b2)



#####Q9#######
# generate text based on word frequencies alone for comparison
unique_words = unique(a_lower)
common_word_freq <- frequency[match(b, unique_words)]
word_prob <- common_word_freq / sum(common_word_freq)
generated_text <- sample(b, size = 50, replace = TRUE, prob = word_prob)
cat("random text:","\n",paste(generated_text, collapse = " "))



################
# finally simulate
text_simulator(new_b,M,mlag = 4,nw = 50)
 
