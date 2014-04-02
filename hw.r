newM <- newmat - matrix(rep(rsum,3),3,3)
FinalM / sqrt(rowSums(newM*newM) / 3)


##Color ramp def.
colors<-c('white','black')
cus_col<-colorRampPalette(colors=colors)
 
## Plot the average image of each digit
par(mfrow=c(4,3),pty='s',mar=c(1,1,1,1),xaxt='n',yaxt='n')
all_img<-array(dim=c(10,28*28))
for(di in 0:9)
{
print(di)
all_img[di+1,]<-apply(train[train[,1]==di,-1],2,sum)
all_img[di+1,]<-all_img[di+1,]/max(all_img[di+1,])*255
 
z<-array(all_img[di+1,],dim=c(28,28))
z<-z[,28:1] ##right side up
image(1:28,1:28,z,main=di,col=cus_col(256))
}