# This implemenation of the max-min-hill-climbing algorithm has its Copyright by Michael Bauer
# For any questions send me an email via michael1.bauer@stud.uni-r.de
# File description comment, including purpose of program, inputs and outputs

Example <- function(rows, cols = 5, char = TRUE) {
  binary <- c()
  trinary <- c()
  if(char) {
    binary <- c("LOW", "HIGH")
    trinary <- c("A", "B", "C")
  } else {
    binary <- c(1, 2)
    trinary <- c(1, 2, 3)
  }
  dimnames <- list(c(), c("difficulty", "intelligence", "SAT", "grade", "letter"))#, "prop"))
  studentMatrix <- matrix(, rows, cols, dimnames = dimnames)
  student <- data.frame(studentMatrix, check.names = FALSE)
  # cardinality <- c(2, 2, 2, 3, 2, sample(0, (rows - 5), replace = TRUE))
  # student$prop <- cardinality
  student$difficulty <- sample(binary, rows, replace = TRUE, prob = c(0.6, 0.4))
  student$intelligence <- sample(binary, rows, replace = TRUE, prob = c(0.7, 0.3))
  n <- 1
  for (i in student$intelligence) {
    if (i == binary[1])
      student$SAT[n] <- sample(binary, 1, replace = TRUE, prob = c(0.95, 0.05))

    if (i == binary[2])
      student$SAT[n] <- sample(binary, 1, replace = TRUE,prob = c(0.2, 0.8))

    if (i == binary[2] && student$difficulty[n] == binary[2])
      student$grade[n] <- sample(trinary, 1, replace = TRUE,prob = c(0.5, 0.3, 0.2))

    if (i == binary[2] && student$difficulty[n] == binary[1])
      student$grade[n] <- sample(trinary, 1, replace = TRUE,prob = c(0.9, 0.08, 0.02))

    if (i == binary[1] && student$difficulty[n] == binary[2])
      student$grade[n] <- sample(trinary, 1, replace = TRUE,prob = c(0.05, 0.25, 0.7))

    if (i == binary[1] && student$difficulty[n] == binary[1])
      student$grade[n] <- sample(trinary, 1, replace = TRUE,prob = c(0.3, 0.4, 0.3))

    if (student$grade[n] == trinary[1])
      student$letter[n] <- sample(binary, 1,replace = TRUE,prob = c(0.1, 0.9))

    if (student$grade[n] == trinary[2])
      student$letter[n] <- sample(binary, 1,replace = TRUE,prob = c(0.4, 0.6))

    if (student$grade[n] == trinary[3])
      student$letter[n] <- sample(binary, 1,replace = TRUE,prob = c(0.99, 0.01))

    n <- n+1
  }
  return (student)
}

Rainy <- function(dim, cols = 3) {
  binary <- c(1, 2)
  dimnames <- list(c(), c("sprinkler", "rain", "grassWet"))
  testMatrix <- matrix(, dim, cols, dimnames = dimnames)
  df <- data.frame(testMatrix, check.names = FALSE)
  df$rain <- sample(binary, dim, replace = TRUE, prob = c(0.2, 0.8))
  n <- 1
  for (i in df$rain) {

    if (i == binary[1])
      df$sprinkler[n] <- sample(binary, 1, replace = TRUE, prob = c(0.01, 0.99))

    if (i == binary[2])
      df$sprinkler[n] <- sample(binary, 1, replace = TRUE, prob = c(0.40, 0.60))

    if (i == binary[1] && df$sprinkler[n] == 1)
      df$grassWet[n] <- sample(binary, 1, replace = TRUE, prob = c(0.99, 0.01))

    if (i == binary[1] && df$sprinkler[n] == 2)
      df$grassWet[n] <- sample(binary, 1, replace = TRUE, prob = c(0.80, 0.20))

    if (i == binary[2] && df$sprinkler[n] == 1)
      df$grassWet[n] <- sample(binary, 1, replace = TRUE, prob = c(0.90, 0.10))

    if (i == binary[2] && df$sprinkler[n] == 2)
      df$grassWet[n] <- sample(binary, 1, replace = TRUE, prob = c(0.00, 1.00))


    n <- n+1
  }

  return (df)
}

# Allocating the example of the book. Make it global and allocate it as a matrix for testing
# purpose. Later on it should be a data frame.
MyExample <- Example(5000, char = FALSE)