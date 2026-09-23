# Run at the end of the final analysis environment and commit the generated files.

dir.create("environment", showWarnings=FALSE)
writeLines(capture.output(sessionInfo()), "environment/R_sessionInfo.txt")

ip <- as.data.frame(installed.packages()[,c("Package","Version")])
write.csv(ip, "environment/package_versions.csv", row.names=FALSE)
