# Sample data for test

fl <- system.file("extdata", "Broad.xml", package = "GSEABase")
gss <- getBroadSets(fl) # GeneSetCollection of 2 sets


info <- structure(list(
  `1` = NA_character_,
  `10` = c("1430728", "156580", "156582", "211859"),
  `2` = c(
    "109582", "114608", "140837", "140877",
    "1474228", "1474244", "162582", "194315",
    "194840", "76002", "76005"
  ),
  `3` = c("194840", "156580"), `4` = NA_character_,
  `5` = NA_character_, `6` = NA_character_,
  `7` = NA_character_, `8` = NA_character_,
  `9` = c("1430728", "156580", "156582", "211859")
),
.Names = c("1", "10", "2", "3", "4", "5", "6", "7", "8", "9")
)

clusters <- structure(list(
  cluster1 = c("10", "3"),
  cluster2 = c("10", "2", "9"),
  cluster3 = c("2", "9", "3", "4")
),
.Names = c("cluster1", "cluster2", "cluster3")
)

genes.id <- c(
  "52", "11342", "80895", "57654", "58493", "1164", "1163", "4150",
  "2130", "159"
)
genes.symbol <- c(
  "ACP1", "RNF13", "ILKAP", "UVSSA", "INIP", "CKS2", "CKS1B",
  "MAZ", "EWSR1", "ADSS"
)

d <- structure(c(
  0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
  0.285714285714286
),
.Dim = c(3L, 3L),
.Dimnames = list(c("a", "b", "c"), c("d", "e", "f"))
)

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
