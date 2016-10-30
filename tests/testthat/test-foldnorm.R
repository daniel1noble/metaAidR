
context("Check the folded normal function")

test_that("Check folded normal", {
	set.seed(60)
	x <- rnorm(1000, 1, 1)
	# Calculate the mean of the folded norm from raw x
	y <- round(foldnorm(mu = mean(x), sd = sd(x), type = "raw"), digits = 2)
	#Should be close to: 
	absX <- mean(abs(x))

	expect_is(y, "data.frame", info = "x is not a data frame")
	expect_equal(absX, y[,1], tolerance = 0.01)
})