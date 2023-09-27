test_that("basic estimation is consistent", {
    set.seed(302019)
    df <- read.csv("testcase.csv")
    expect_snapshot(cran = TRUE, pgc::estimate(T_seq = df$T_seq, N_seq = df$N_seq,
        n_seq = df$n_seq))
})
