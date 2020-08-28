# Create sample vector
X <- c(1:100); print(X)

# Create sample matrix
M <- matrix(c(1:100),nrow=10); print(M)

# Set limits
minV <- 15; maxV <- 85;

# Limit vector
sapply(X, function(y) min(max(y,minV),maxV))

# Limit matrix
apply(M, c(1, 2), function(x) min(max(x,minV),maxV))
