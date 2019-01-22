context('Test summaryFD function')

testthat::test_that('Output of summaryStatistics() is the same as output of
                    AutoFD/pft_JC_FDStatistics.m', {
            fd_pft_JC <- read.table("./summaryFD_matlab.csv", sep=",",
                                    header=TRUE, row.names=1)
            summary_pft_JC <- as.matrix(fd_pft_JC[,1:6])
            summary_summaryStatistics <- t(apply(fd_pft_JC[,7:26], 1,
                                                 summaryStatistics))
            colnames(summary_summaryStatistics) <- colnames(summary_pft_JC)
            testthat::expect_equal(summary_pft_JC, summary_summaryStatistics,
                                   tolerance = 0.001)
})
