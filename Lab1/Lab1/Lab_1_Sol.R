# 1
df <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6), c=c(7, 8, 9))
new_vec <- df$c

# 2
new_vec <- df[1,2]

# 3
new_vec <- df[1,] * df[,2]

# 4
df$d <- as.numeric(new_vec)

# 5
df[1,] <- df[1,] * c(2,3)


# 6
write.table(df, "Lab1", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

