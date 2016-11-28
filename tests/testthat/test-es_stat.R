context("Check that effect size statistics are being calculated correctly")

test_that("Checking effect size statistics...", {
	
	# Example data from Borenstein et al. 2009.	
	d <- es_stat(m1 = 103, m2 = 100, sd1 = 5.5, sd2 = 4.5, n1 = 50, n2 = 50, type = "d")
	g <- es_stat(m1 = 103, m2 = 100, sd1 = 5.5, sd2 = 4.5, n1 = 50, n2 = 50, type = "g")
	lnRR <- es_stat(m1 = 103, m2 = 100, sd1 = 5.5, sd2 = 4.5, n1 = 50, n2 = 50,  type = "lnRR")

	expect_equal(as.numeric(round(d[,1], digits =4)), 0.5970,  info = "d failed")
	expect_equal(as.numeric(round(d[,2], digits =4)), 0.0418,  info = "Vd failed")
	expect_equal(as.numeric(round(g[,1], digits =4)), 0.5924, info = "g failed")
	expect_equal(as.numeric(round(g[,2], digits =4)), 0.0412,  info = "Vg failed")
	expect_equal(as.numeric(round(lnRR[,1], digits =4)), 0.0296, info = "lnRR failed")
	expect_equal(as.numeric(round(lnRR[,1], digits =4)), 0.0296, info = "VlnRR failed")
})