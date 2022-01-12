context("transfactor correctly estimates transcription factor activities.")

set.seed(3)
n <- 400
G <- 100
nTF <- 20
means <- matrix(rep(rlnorm(n = G, meanlog = 4, sdlog = 1), n),
                nrow = G, ncol = n, byrow = FALSE
)
dispersions <- matrix(rep(runif(n = G, min = 0.8, max = 3), n),
                      nrow = G, ncol = n, byrow = FALSE
)
# simulate NB counts
counts <- matrix(rnbinom(n = G * n, mu = means, size = 1 / dispersions),
                 nrow = G, ncol = n)
X <- matrix(0, nrow=G, ncol=nTF)
X[cbind(sample(G, size=450, replace=TRUE),
        sample(nTF, size=450, replace=TRUE))] <- 1
rownames(X) <- rownames(counts) <- paste0("gene",1:G)

test_that("Default activity estimation doesn't error", {
  resPoisson <- estimateActivity(counts, X, model="poisson")
  resDirMult <- estimateActivity(counts, X, model="dirMult")
  resDirMultAlpha <- estimateActivity(counts, X, model="dirMultAlpha")
})

test_that("Sequencing depth of results is equal to input", {
  resPoisson <- estimateActivity(counts, X, model="poisson")
  YPoisson <- tfCounts(mu_gtc = resPoisson$mu_gtc,
                       counts = counts)
  expect_equal(colSums(counts), colSums(YPoisson))

  resDirMult <- estimateActivity(counts, X, model="dirMult")
  YDirMult <- tfCounts(mu_gtc = resDirMult$mu_gtc,
                       counts = counts)
  expect_equal(colSums(counts), colSums(YDirMult))

  resDirMultAlpha <- estimateActivity(counts, X, model="dirMultAlpha")
  YDirMultAlpha <- tfCounts(mu_gtc = resDirMultAlpha$mu_gtc,
                       counts = counts)
  expect_equal(colSums(counts), colSums(YDirMultAlpha))
})


test_that("DirMult without prior equals Poisson result in simple setting", {
  resPoisson <- estimateActivity(counts, X, model="poisson")
  resDirMultNoAlpha <- estimateActivity(counts, X, model="dirMult", alphaScale = "none")
  resDirMultNoAlpha2 <- estimateActivity(counts, X, model="dirMultAlpha", alphaScale = "none")
  expect_equal(resPoisson$mu_tc, resDirMultNoAlpha$mu_tc)
  expect_equal(c(resPoisson$mu_tc), c(resDirMultNoAlpha2$mu_tc), ignore_attr = TRUE)
})
