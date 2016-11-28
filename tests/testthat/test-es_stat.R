context("Check that effect size statistics are being calculated correctly")

test_that("Checking effect size statistics...", {
	
	# Example data from Borenstein et al. 2009.	
	d <- es_stat(m1 = 103, m2 = 100, sd1 = 5.5, sd2 = 4.5, n1 = 50, n2 = 50, type = "d")
	g <- es_stat(m1 = 103, m2 = 100, sd1 = 5.5, sd2 = 4.5, n1 = 50, n2 = 50, type = "g")

	expect_equal(as.numeric(round(d[,1], digits =4)), 0.5970,  info = "d failed")
	expect_equal(as.numeric(round(d[,2], digits =4)), 0.0418,  info = "Vd failed")
	expect_equal(as.numeric(round(g[,1], digits =4)), 0.5924, info = "g failed")
	expect_equal(as.numeric(round(g[,2], digits =4)), 0.0412,  info = "Vg failed")

	lnRR <- es_ratio(m1 = 103, m2 = 100, sd1 = 5.5, sd2 = 4.5, n1 = 50, n2 = 50,  type = "lnRR")
	lnVR <- es_ratio(m1 = 169.76, m2 = 110.91, sd1 = 49.29, sd2 = 59.162, n1 = 75, n2 = 57,  type = "lnVR")
	lnCVR <- es_ratio(m1 = 169.76, m2 = 110.91, sd1 = 49.29, sd2 = 59.162, n1 = 75, n2 = 57,  cor = 0.5, type = "lnCVR")

	expect_equal(as.numeric(round(lnRR[,1], digits =4)), 0.0296, info = "lnRR failed")
	expect_equal(as.numeric(round(lnRR[,2], digits =4)), 0.0025, info = "VlnRR failed")
	expect_equal(as.numeric(round(lnCVR[,1], digits =4)), -0.6104, info = "lnCVR failed")
	expect_equal(as.numeric(round(lnCVR[,2], digits =4)), 0.0124, info = "VlnCVR failed")
	expect_equal(as.numeric(round(lnVR[,1], digits =4)), -0.1847, info = "lnVR failed")
	expect_equal(as.numeric(round(lnVR[,2], digits =4)), 0.0157, info = "VlnVR failed")
})