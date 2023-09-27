library("BioCor")
context("Testing pathSim")

test_that("pathSim", {
  test <- pathSim("194840", "156580", info)
  expect_equal(test, 0.4)
  expect_error(pathSim(NULL, "156580", info), "several")
  expect_error(pathSim(1, "156580", info), "characters")
  expect_error(pathSim("1", "156580", "a"), "list")
  expect_true(is.na(pathSim("1", "156580", info)))
  expect_equal(pathSim("109582", "114608", info), 1)
})

test_that("pathSim with GSC", {
  Info <- new("GeneSetCollection",
              .Data = list(
                new("GeneSet",
                    geneIdType = new("NullIdentifier", type = "Null", annotation = ""),
                    geneIds = c("10", "9"),
                    setName = "1430728",
                    setIdentifier = "linux:7312:Fri Mar 30 13:41:18 2018:36",
                    shortDescription = "", longDescription = "", organism = "",
                    pubMedIds = character(0), urls = character(0),
                    contributor = character(0),
                    creationDate = character(0),
                    collectionType = new("NullCollection", type = "Null")),
                new("GeneSet",
                    geneIdType = new("NullIdentifier", type = "Null", annotation = ""),
                    geneIds = c("10", "3", "9"), setName = "156580",
                    setIdentifier = "linux:7312:Fri Mar 30 13:41:18 2018:37",
                    shortDescription = "", longDescription = "", organism = "",
                    pubMedIds = character(0), urls = character(0),
                    contributor = character(0), creationDate = character(0),
                    collectionType = new("NullCollection", type = "Null")),
                new("GeneSet",
                    geneIdType = new("NullIdentifier", type = "Null", annotation = ""),
                    geneIds = c("10", "9"), setName = "156582",
                    setIdentifier = "linux:7312:Fri Mar 30 13:41:18 2018:38",
                    shortDescription = "", longDescription = "", organism = "",
                    pubMedIds = character(0), urls = character(0),
                    contributor = character(0),

                    collectionType = new("NullCollection", type = "Null")),
                new("GeneSet",
                    geneIdType = new("NullIdentifier", type = "Null", annotation = ""),
                    geneIds = c("2", "3"), setName = "194840",
                    setIdentifier = "linux:7312:Fri Mar 30 13:41:18 2018:39",
                    shortDescription = "", longDescription = "", organism = "",
                    pubMedIds = character(0), urls = character(0),
                    contributor = character(0),  creationDate = character(0),
                    collectionType = new("NullCollection", type = "Null")),
                new("GeneSet",
                    geneIdType = new("NullIdentifier", type = "Null", annotation = ""),
                    geneIds = c("10", "9"), setName = "211859",
                    setIdentifier = "linux:7312:Fri Mar 30 13:41:18 2018:40",
                    shortDescription = "", longDescription = "", organism = "",
                    pubMedIds = character(0), urls = character(0),
                    contributor = character(0),
                    creationDate = character(0),
                    collectionType = new("NullCollection", type = "Null"))))
  test <- pathSim("194840", "156580", Info)
  expect_equal(test, 0.4)
  expect_error(pathSim(NULL, "156580", Info), "several")
  expect_error(pathSim(1, "156580", Info), "characters")
  expect_error(pathSim("1", "156580", "a"), "list")
  expect_true(is.na(pathSim("1", "156580", Info)))
  expect_true(is.na(pathSim("109582", "114608", Info)))
})
