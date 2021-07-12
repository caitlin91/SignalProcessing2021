# # # # #
# Day 1 #
# # # # #

### install packages
install.packages(c("RCurl","sound","phonTools","raster","ggplot2"))


### function to download data
dl_from_dropbox <- function(x, key) {
  bin <- RCurl::getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),
                             ssl.verifypeer = FALSE)
  con <- file(x, open = "wb")
  writeBin(bin, con)
  close(con)
}


### clean up memory
rm(list=setdiff(ls(), "dl_from_dropbox"))

gc()

.rs.restartR()


### set up audio player
# Linux
sound::setWavPlayer('/bin/aplay')
sound::setWavPlayer('/bin/mplayer')

# Windows
sound::setWavPlayer('C:/Program Files/Window Media Player/wmplayer.exe')

# Mac
sound::setWavPlayer('/Applications/"QuickTime Player.app"/Contents/MacOS/"QuickTime Player"') 


### download first data set
dl_from_dropbox('EMA-EMG_data.Rda', 'cdzc6vcjnn8ppqw')
load('EMA-EMG_data.Rda')

summary(EMA.EMG)


### sampling rate
EMA.EMG$SR

sr <- EMA.EMG$SR$audio


### play audio 
audio <- sound::as.Sample(as.vector(EMA.EMG$audio), fs=sr)
sound::play(audio)


### view spectrogram
phonTools::spectrogram(EMA.EMG$audio, fs=sr)

phonTools::spectrogram(EMA.EMG$audio, fs=sr, colors=F)


### choose time points
EMA.EMG$segments

t1 <- EMA.EMG$segments$start[1]
t2 <- EMA.EMG$segments$end[5]

t1; t2


### time to samples
t1 * sr

t1 <- round(EMA.EMG$segments$start[1] * sr)
t2 <- round(EMA.EMG$segments$end[5] * sr)

t1; t2


### view spectrogram of audio segment
phonTools::spectrogram(EMA.EMG$audio[t1:t2], fs=sr, colors=F, dynamicrange=65)

phonTools::spectrogram(EMA.EMG$audio[t1:t2], fs=sr, colors=F, dynamicrange=65, maxfreq=8000)



### segment vowel by hand
coords <- locator(2)

coords

t1 <- round(coords$x[1]*sr/1000)
t2 <- round(coords$x[2]*sr/1000)


### spectrogram of segmented vowel
phonTools::spectrogram(EMA.EMG$audio[t1:t2], fs=sr, colors=F, dynamicrange=65, maxfreq=8000)


### add offset
offset <- EMA.EMG$segments$start[1]
offset

coords

t1 <- round(
  (offset + coords$x[1]/1000)*sr
)

t2 <- round(
  (offset + coords$x[2]/1000)*sr
)


### spectrogram of corrected segmented vowel
phonTools::spectrogram(EMA.EMG$audio[t1:t2], fs=sr, colors=F, dynamicrange=65, maxfreq=8000)


### spectral slice formant analysis
phonTools::findformants(EMA.EMG$audio[t1:t2], fs=sr)

phonTools::findformants(EMA.EMG$audio[t1:t2], fs=sr, showbws=T, showrejected=F)


### save formant results
formants <- phonTools::findformants(EMA.EMG$audio[t1:t2], fs=sr, verify=F)
formants


### acoustic vowel space
EMA.EMG$segments

vowels <- c("a","e","i","o","{","@")

vowelnum <- sum(EMA.EMG$segments$phone %in% vowels)
vowelnum

length(unique(EMA.EMG$segments$word))

form.dat <- data.frame(
  vowel = character(vowelnum),
  f1 = numeric(vowelnum),
  f2 = numeric(vowelnum)
)

x <- 1
for (phone in 1:nrow(EMA.EMG$segments)) {
  if (EMA.EMG$segments$phone[phone] %in% vowels) {
    t1 <- round(EMA.EMG$segments$start[phone]*sr)
    t2 <- round(EMA.EMG$segments$end[phone]*sr)
    
    dur <- t2 - t1
    
    t20 <- t1 + round(dur*0.2)
    t80 <- t1 + round(dur*0.8)
    
    formants <- phonTools::findformants(EMA.EMG$audio[t20:t80], fs=sr, verify=F)
    form.dat$vowel[x] <- EMA.EMG$segments$phone[phone]
    form.dat[x, c("f1","f2")] <- formants$formant[1:2]
    x <- x+1
  }
}

par(mfrow=c(1,1))

ggplot(form.dat, aes(x=f2,y=f1,col=vowel,shape=vowel)) + geom_point(cex=3) + 
  scale_y_reverse() + scale_x_reverse() + theme_bw()


### multiple slices: formant tracking
winlen  <- 0.01

winlen <- round(sr*winlen)
bins   <- floor(length(EMA.EMG$audio)/winlen)

chunk <- EMA.EMG$audio[1:winlen]
plot(chunk, type='l')


### windowing
hannwin <- ( 0.5 - (0.5 * cos(2*pi*(0:(winlen-1))/(winlen-1))) )
plot(hannwin)

hannwin <- phonTools::windowfunc(winlen, type="hanning")
plot(hannwin*chunk, type='l')


### formant tracking
formant.tracks <- data.frame(time=numeric(bins),
                             f1=numeric(bins),
                             f2=numeric(bins),
                             f3=numeric(bins))
for (bin in 1:bins) {
  t1 <- (bin-1)*winlen + 1
  t2 <- bin*winlen
  
  chunk <- EMA.EMG$audio[t1:t2]
  
  formants <- phonTools::findformants(hannwin*chunk, fs=sr, verify=F)
  
  midpoint <- mean(c(t1,t2))/sr
  
  formant.tracks[bin,] <- c(midpoint, formants$formant[1:3])
}


### manual clean up
formant.tracks$f1[formant.tracks$f1 < 250 | formant.tracks$f1 > 1200] <- NA
formant.tracks$f2[formant.tracks$f2 < 800 | formant.tracks$f2 > 2800] <- NA
formant.tracks$f3[formant.tracks$f3 < 2500 | formant.tracks$f3 > 4000] <- NA


### view results of formant tracking
t1 <- 5.3
t2 <- 6.3

formant.chunk <- formant.tracks[formant.tracks$time>=t1 & formant.tracks$time<=t2,]


### compare to phonTools::spectrogram()
phonTools::spectrogram(EMA.EMG$audio[round(t1*sr):round(t2*sr)], fs=sr, colors=F)

par(new=T,mfrow=c(1,1))

plot(formant.chunk$time - t1, formant.chunk$f1,
     col=1, pch=16, cex=0.75, 
     xaxt='n', yaxt='n', xlab=NA, ylab=NA, ylim=c(180,4820))
points(formant.chunk$time - t1, formant.chunk$f2,
       col=2, pch=16, cex=0.75, 
       xaxt='n', yaxt='n', xlab=NA, ylab=NA)
points(formant.chunk$time - t1, formant.chunk$f3,
       col=3, pch=16, cex=0.75, 
       xaxt='n', yaxt='n', xlab=NA, ylab=NA)

phonTools::formanttrack(EMA.EMG$audio[round(t1*sr):round(t2*sr)], formants=3,
                        fs=sr, windowlength=10, timestep=10)


### compare to phonTools::spectrogram() using default parameters
phonTools::formanttrack(EMA.EMG$audio[round(t1*sr):round(t2*sr)], formants=3, fs=sr)


### save formant tracks
formants <- phonTools::formanttrack(EMA.EMG$audio[round(t1*sr):round(t2*sr)], formants=3, fs=sr, show=F)
head(formants)


### download second data set
dl_from_dropbox('ultrasound_data.Rda', '68nfc5yie653kf7')
load('ultrasound_data.Rda')
summary(US.data)


### listen to audio
sr    <- US.data$SR$audio
audio <- sound::as.Sample(as.vector(US.data$data1$audio), sr)
sound::play(audio)


### view first ultrasound frame
grayscale <- gray(seq(0, 1, length = 256))
image(US.data$data1$images[,,1], xaxt='n', yaxt='n', col=grayscale)


### plot tongue contour
plot(US.data$data1$contours$coords$x50, US.data$data1$contours$coords$y50, type='l')


### breakout group exercise!


### Cartesian to polar coordinates
r  <- sqrt(x^2 + y^2)
th <- atan2(y, x)


### radians to degrees
th <- th*180/pi


### Cartesian to polar function
cart2polar <- function(x, y, degrees=F) {
  r  <- sqrt(x^2 + y^2)
  th <- atan2(y, x)
  if (degrees) {
    th <- th*180/pi
  }
  return(data.frame(r, th))
}


### polar to Cartesian coordinates
x <- r*cos(th)
y <- r*sin(th)


### radians to degrees
x <- r*cos(th*pi/180)
y <- r*sin(th*pi/180)


### polar to Cartesian function
polar2cart <- function(r, th, degrees=F) {
  if (degrees) {
    x <- r*cos(th*pi/180)
    y <- r*sin(th*pi/180)
  } else {
    x <- r*cos(th)
    y <- r*sin(th)
  }
  return(data.frame(x, y))
}


### convert contour data to polar coordinates
US.data$data1$contours$polar <- c()

pol25 <- cart2polar(US.data$data1$contours$coords$x25,
                    US.data$data1$contours$coords$y25,
                    degrees=T)
US.data$data1$contours$polar$r25  <- pol25$r
US.data$data1$contours$polar$th25 <- pol25$th

pol50 <- cart2polar(US.data$data1$contours$coords$x50,
                    US.data$data1$contours$coords$y50,
                    degrees=T)
US.data$data1$contours$polar$r50  <- pol50$r
US.data$data1$contours$polar$th50 <- pol50$th

pol75 <- cart2polar(US.data$data1$contours$coords$x75,
                    US.data$data1$contours$coords$y75,
                    degrees=T)
US.data$data1$contours$polar$r75  <- pol75$r
US.data$data1$contours$polar$th75 <- pol75$th


polpal <- cart2polar(US.data$data1$contours$palate$x,
                     US.data$data1$contours$palate$y,
                     degrees=T)
US.data$data1$contours$polpalate <- polpal


### visualize the results
xrange <- range(c(
  US.data$data1$contours$polpalate$th,
  US.data$data1$contours$polar$th25,
  US.data$data1$contours$polar$th50,
  US.data$data1$contours$polar$th75
))
yrange <- range(c(
  US.data$data1$contours$polpalate$r,
  US.data$data1$contours$polar$r25,
  US.data$data1$contours$polar$r50,
  US.data$data1$contours$polar$r75
))

plot(NA, xlim=rev(xrange), ylim=yrange,
     xlab="Angular coordinate (theta)",
     ylab="Radial coordinate (r)")
lines(US.data$data1$contours$polpalate$th, US.data$data1$contours$polpalate$r)
lines(US.data$data1$contours$polar$th50, US.data$data1$contours$polar$r25, col='blue')
lines(US.data$data1$contours$polar$th25, US.data$data1$contours$polar$r50, col='purple')
lines(US.data$data1$contours$polar$th75, US.data$data1$contours$polar$r75, col='red')


### add theta lines
pol.lines <- list(th=seq(0,180,by=10),
                  r1=rep(10,19),
                  r2=rep(100,19))

for (line in 1:19){
  lines(
    c(pol.lines$th[line],pol.lines$th[line]),
    c(pol.lines$r1[line],pol.lines$r2[line]),
    lty=2)
}


### visualize the results
xrange <- range(c(
  US.data$data1$contours$palate$x,
  US.data$data1$contours$coords$x25,
  US.data$data1$contours$coords$x50,
  US.data$data1$contours$coords$x75
))
yrange <- range(c(
  US.data$data1$contours$palate$y,
  US.data$data1$contours$coords$y25,
  US.data$data1$contours$coords$y50,
  US.data$data1$contours$coords$y75
))

plot(NA, xlim=xrange, ylim=yrange,
     xlab="<-- Posterior (mm) | Anterior (mm) -->",
     ylab="<-- Inferior (mm) | Superior (mm) -->")
lines(US.data$data1$contours$palate$x, US.data$data1$contours$palate$y)
lines(US.data$data1$contours$coords$x25, US.data$data1$contours$coords$y25, col='blue')
lines(US.data$data1$contours$coords$x50, US.data$data1$contours$coords$y50, col='purple')
lines(US.data$data1$contours$coords$x75, US.data$data1$contours$coords$y75, col='red')

for (x in 1:19){
  lines(
    c(cart.lines1$x[x],cart.lines2$x[x]),
    c(cart.lines1$y[x],cart.lines2$y[x]),
    lty=2)
}


### convert theta lines to Cartesian coordinates
cart.lines1 <- polar2cart(pol.lines$r1,
                          pol.lines$th,
                          degrees=T)
cart.lines2 <- polar2cart(pol.lines$r2,
                          pol.lines$th,
                          degrees=T)


### visualize the results
xrange <- range(c(
  US.data$data1$contours$palate$x,
  US.data$data1$contours$coords$x25,
  US.data$data1$contours$coords$x50,
  US.data$data1$contours$coords$x75
))
yrange <- range(c(
  US.data$data1$contours$palate$y,
  US.data$data1$contours$coords$y25,
  US.data$data1$contours$coords$y50,
  US.data$data1$contours$coords$y75
))

plot(NA, xlim=xrange, ylim=yrange,
     xlab="<-- Posterior (mm) | Anterior (mm) -->",
     ylab="<-- Inferior (mm) | Superior (mm) -->")
lines(US.data$data1$contours$palate$x, US.data$data1$contours$palate$y)
lines(US.data$data1$contours$coords$x25, US.data$data1$contours$coords$y25, col='blue')
lines(US.data$data1$contours$coords$x50, US.data$data1$contours$coords$y50, col='purple')
lines(US.data$data1$contours$coords$x75, US.data$data1$contours$coords$y75, col='red')

for (x in 1:19){
  lines(
    c(cart.lines1$x[x],cart.lines2$x[x]),
    c(cart.lines1$y[x],cart.lines2$y[x]),
    lty=2)
}



# # # # #
# Day 2 #
# # # # #

### listen to audio
sr    <- US.data$SR$audio
audio <- sound::as.Sample(as.vector(US.data$data2$audio), sr)
sound::play(audio)


### create variance image
varimg <- apply(US.data$data2$images, 1:2, var)
image(varimg, col=grayscale)


### select a line of interest
coords <- locator(2)

lines(coords$x, coords$y, col='red', lwd=2)
points(coords$x[1], coords$y[1], col='cyan', pch=1, cex=2)
points(coords$x[2], coords$y[2], col='orange', pch=1, cex=2)


### beware of image origin
? image


### convert coordinates to pixels
coords2 <- c()

coords2$x <- round(coords$x * dim(varimg)[1])
coords2$y <- round(coords$y * dim(varimg)[2])

coords2


### Euclidean distance (in pixels)
edist <- round(
  sqrt(
    (coords2$x[1] - coords2$x[2])^2 + (coords2$y[1] - coords2$y[2])^2
    )
)
edist


### generate pixel values along the line
xx <- seq(coords2$x[1], coords2$x[2], length.out=edist)
yy <- seq(coords2$y[1], coords2$y[2], length.out=edist)


head(cbind(round(xx),round(yy))); tail(cbind(round(xx),round(yy)))


### view a test image
image(US.data$data2$images[,,1], col=grayscale)

lines(coords$x, coords$y, col='red', lwd=2)
points(coords$x[1], coords$y[1], col='cyan', pch=1, cex=2)
points(coords$x[2], coords$y[2], col='orange', pch=1, cex=2)


### view pixel values along the line of interest
vals <- numeric(edist)

for (val in 1:edist) {
  vals[val] <- US.data$data2$images[round(xx)[val],round(yy)[val],1]
}

plot(vals, type='l')


### center of gravity
cog <- sum(seq(1,edist) * vals) / sum(vals)
cog


### calculate COG along line of interest for all images
imgnum <- dim(US.data$data2$images)[3]

cogs  <- numeric(imgnum)

for (img in 1:imgnum){
  this.img <- US.data$data2$images[,,img]
  
  vals <- numeric(edist)
  
  for (val in 1:edist) {
    vals[val] <- this.img[round(xx)[val],round(yy)[val]]
  }
  
  cogs[img] <- sum(seq(1,edist) * vals) / sum(vals)
}

plot(cogs, type='l')


### convert COG values to percentages
cogs <- 100 * (cogs - min(cogs)) / (max(cogs) - min(cogs))


### check signals for front and back vowels
US.data$data2$segments$word


t1 <- round(US.data$data2$segments$wstart[3] * US.data$SR$US)
t2 <- round(US.data$data2$segments$wend[3] * US.data$SR$US)

chunk <- cogs[t1:t2]
wtime <- seq(0, 100, length.out=length(chunk))

plot(wtime, chunk, type='l', col='red', ylim=range(cogs),
     xlab="Normalized time (%)", ylab="Tongue retraction (%)")

t1 <- round(US.data$data2$segments$wstart[11] * US.data$SR$US)
t2 <- round(US.data$data2$segments$wend[11] * US.data$SR$US)

chunk <- cogs[t1:t2]
wtime <- seq(0, 100, length.out=length(chunk))

points(wtime, chunk, type='l', col='red')


t1 <- round(US.data$data2$segments$wstart[7] * US.data$SR$US)
t2 <- round(US.data$data2$segments$wend[7] * US.data$SR$US)

chunk <- cogs[t1:t2]
wtime <- seq(0, 100, length.out=length(chunk))

points(wtime, chunk, type='l', col='blue')

t1 <- round(US.data$data2$segments$wstart[8] * US.data$SR$US)
t2 <- round(US.data$data2$segments$wend[8] * US.data$SR$US)

chunk <- cogs[t1:t2]
wtime <- seq(0, 100, length.out=length(chunk))

points(wtime, chunk, type='l', col='blue')

legend("topright", legend=c("front Vs", "back Vs"), col=c("red","blue"), lty=1)


### play audio of data to be used in PCA
audio <- sound::as.Sample(as.vector(US.data$data2$audio), sr)
sound::play(audio)


### down-sampling images for PCA model
nframes     <- dim(US.data$data2$images)[3]
imfact      <- 37
tongue.dat  <- matrix(nrow=nframes,
                      ncol=round(nrow(US.data$data2$images[,,1])/imfact) *
                        round(ncol(US.data$data2$images[,,1])/imfact))
ncol(tongue.dat)


### apply down-sampling to all images
for (frame in 1:nframes){
  this.frame <- US.data$data2$images[,,frame]
  s <- raster::raster(nrow=round(nrow(this.frame)/imfact), 
                      ncol=round(ncol(this.frame)/imfact))
  r <- raster::raster(this.frame)
  raster::extent(r) <- raster::extent(c(-180, 180, -90, 90))
  new.frame <- raster::resample(r, s)
  tongue.dat[frame, ] <- raster::as.matrix(new.frame, nrow=1, byrow=T)
}


### visualize the results
image(this.frame, xaxt='n', yaxt='n', col=grayscale)
image(raster::as.matrix(new.frame), xaxt='n', yaxt='n', col=grayscale)


### perform the PCA
pca <- prcomp(tongue.dat)


### determine number of PCs to retain (not necessary, but useful info)
vars      <- apply(pca$x, 2, var)
props     <- vars / sum(vars)
var.exp   <- cumsum(props)
to.keep   <- var.exp[var.exp < 0.8]
to.keep   <- length(to.keep) + 1

PCs <- rep(paste0('PC', seq(1, to.keep)))


### PC loading heatmaps
coeff.maps <- c()
for (pc in PCs){
  coeff.maps[[pc]] <- t(matrix(pca$rotation[,pc], 
                               nrow = round(ncol(US.data$data2$images[,,1])/imfact),
                               byrow = TRUE))
}

image(coeff.maps$PC1, xaxt='n', yaxt='n', col=grayscale)
image(coeff.maps$PC2, xaxt='n', yaxt='n', col=grayscale)


### articulatory vowel space using PC1 and PC2 scores
midpoints <- rowMeans(US.data$data2$segments[,c('vstart','vend')])
midframes <- round(midpoints * US.data$SR$US)

PC1.mid   <- pca$x[midframes, 1]
PC2.mid   <- pca$x[midframes, 2]

plot(NA, xlim=range(PC1.mid), ylim=range(PC2.mid), xlab='PC1 score', ylab='PC2 score')
text(PC1.mid, PC2.mid, labels=US.data$data2$segments[,1])


### time varying PC1 score
word    <- "who'd"

word.num<- which(US.data$data2$segments==word)

vstart  <- US.data$data2$segments$vstart[word.num]
vend    <- US.data$data2$segments$vend[word.num]
wstart  <- US.data$data2$segments$wstart[word.num]
wend    <- US.data$data2$segments$wend[word.num]

int     <- round(seq(wstart*US.data$SR$US, wend*US.data$SR$US))
time    <- seq(wstart, wend, length.out=length(int))

scores  <- pca$x[int,]
PC      <- 'PC1'

plot(time, scores[,c(PC)], type='l', xlab='Time (s)', ylab=paste(PC,'score'), 
     main=paste0('Australian English /hVd/ word: "',word,'"'))
abline(v=vstart, col='blue'); abline(v=vend, col='red')


### breakout group exercise!