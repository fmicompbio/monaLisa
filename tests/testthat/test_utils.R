test_that(".assertScalar works", {
    expect_error(.assertScalar(numeric(0)), "scalar value")
    expect_error(.assertScalar(1:3),        "scalar value")
    expect_true(.assertScalar(1))
    expect_true(.assertScalar("test"))
    
    expect_error(.assertScalar(1, type = "character"), "must be of type")
    expect_error(.assertScalar("1", type = "numeric"), "must be of type")
    expect_error(.assertScalar(1.5, type = "integer"), "must be of type")
    expect_error(.assertScalar("0.5", rngIncl = 0:1),  "must be of type")
    expect_true(.assertScalar(1, type = "numeric"))
    expect_true(.assertScalar(TRUE, type = "logical"))
    
    expect_error(.assertScalar( 2, rngIncl = 0:1), "inclusive")
    expect_error(.assertScalar(-1, rngIncl = 0:1), "inclusive")
    expect_true(.assertScalar( 4, rngIncl = c(1,10)))
    expect_true(.assertScalar(-4, rngIncl = c(-4,0)))
    
    expect_error(.assertScalar( 2, rngExcl = 1:2),    "exclusive")
    expect_error(.assertScalar(-1, rngExcl = c(-1,1), "exclusive"))
    expect_true(.assertScalar(1.5, rngExcl = 1:2))
    expect_true(.assertScalar(0, rngExcl = c(-1,1)))
    
    expect_error(.assertScalar("a", validValues = c("b", "c")))
    expect_error(.assertScalar(2, validValues = c(1, 3)))
    expect_true(.assertScalar("b", validValues = c("b", "c")))
    expect_true(.assertScalar(1, validValues = c(1, 3)))
})

test_that(".assertVector works", {
    expect_true(.assertVector(1))
    expect_true(.assertVector(1:3))
    expect_true(.assertVector("test"))
    expect_true(.assertVector(list()))
    
    expect_error(.assertVector(1, type = "character"),    "must be of class")
    expect_error(.assertVector("test", type = "numeric"), "must be of class")
    expect_error(.assertVector("0.5", rngIncl = 0:1),     "must be of class")
    expect_true(.assertVector(1:3, type = "numeric"))
    expect_true(.assertVector(letters, type = "character"))
    
    expect_error(.assertVector(1:4, rngIncl = c(2,10)), "inclusive")
    expect_error(.assertVector(1:4, rngIncl = c(-2,3)), "inclusive")
    expect_true(.assertVector(1:4, rngIncl = c(1,4)))
    expect_true(.assertVector(-3, rngIncl = c(-3,-3)))
    
    expect_error(.assertVector(1:4, rngExcl = c(1,5)),    "exclusive")
    expect_error(.assertVector(-1:4, rngExcl = c(10,12)), "exclusive")
    expect_true(.assertVector(1:4, rngExcl = c(.5,4.5)))
    expect_true(.assertVector(-1:4, rngExcl = c(-2,5)))
    
    expect_error(.assertVector(1:4, len = 2),    "length")
    expect_error(.assertVector("test", len = 2), "length")
    expect_true(.assertVector("test", len = 1))
    expect_true(.assertVector(1:3, len = 3))
})