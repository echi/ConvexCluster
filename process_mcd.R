#
# Process McDonald's data
#
mcd = read.table('McDonaldsNutritionFacts.csv',sep=",",header=TRUE)
#
# Remove duplicates
#
mcd = mcd[which(!duplicated(mcd$Name)),]
#
# Remove Daily Rx & Serving Size
#
mcd = mcd[,-c(2,3,7,9,12,14,16,18)]

for (i in 1:nrow(mcd)) {
  if (mcd[i,2] > 0)
    mcd[i,-1] = mcd[i,-1]/mcd[i,2]
}
#
# Remove Large, Medium, and Child sized drinks
#
mcd = mcd[-grep(pattern="(Large)", x=mcd$Name),]
mcd = mcd[-grep(pattern="(Medium)", x=mcd$Name),]
mcd = mcd[-grep(pattern="(Child)", x=mcd$Name),]

write.table(scale(mcd[,-c(1:3)]),"mcd.csv",sep=",",row.names=FALSE,col.names=FALSE)
write.table(mcd[,-1],"mcd.csv",sep=",",row.names=FALSE,col.names=FALSE)

write.table(mcd[,1],"mcd_names.csv",sep=",",row.names=FALSE,col.names=FALSE)