# Create a probability vector
x <- rdirichlet(1L, runif(20))

test_that("random dirichlet distribution sums to 1", {

    expect_true(sum(x) == 1)
})

# randomly shuffle x to be 10 aitchison away
y <- perturb_by_relab(x, by = 10)

d <- rbind(x, y) |>
    .clr(cols_as_features = TRUE) |>
    t() |>
    dist(method = "euclidean")

test_that("perturbation works", {

    expect_true(0.1 > ((d - 10)^2)  )
})
