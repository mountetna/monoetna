library('RUnit')

source('server.R')

test.suite <- defineTestSuite(
  "rtemis",
  dirs = file.path("test"),
  testFileRegexp = '^.*\\.R'
)

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)
