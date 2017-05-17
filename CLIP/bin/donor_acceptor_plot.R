library(ggplot2)
library(RColorBrewer)

max_dist <- 10000
sample_count <- 1000

df <- read.table('results-bam/EWS.filtered.reads.dist.txt', header=F, sep='\t', comment.char='', stringsAsFactors=F)
colnames(df) <- c('read', 'class', 'dist')
df$sample<-'enriched'
df$repl<-0

df_d_density <- density(df[df$class=='DONOR' & abs(df$dist) < max_dist,]$dist, from=-max_dist, to=max_dist, n=2*max_dist)$y
df_a_density <- density(df[df$class=='ACCEPTOR' & abs(df$dist) < max_dist,]$dist, from=-max_dist, to=max_dist, n=2*max_dist)$y

df0<-df
i<-1

acceptor_density <- c()
donor_density <- c()

while (i <= sample_count) {
    cat(paste0('Reading file: results-bam/readlists/EWS.failed.reads.',i,'.txt.bed.dist\n'))
    df1 <- read.table(paste0('results-bam/readlists/EWS.failed.reads.',i,'.txt.bed.dist'), header=F, sep='\t', comment.char='', stringsAsFactors=F)
    colnames(df1) <- c('read', 'class', 'dist')
    df1$sample<-'bg'
    df1$repl<-i
    
    dd <- density(df1[df1$class=='DONOR' & abs(df1$dist) < max_dist,]$dist, from=-max_dist, to=max_dist, n=2*max_dist)
    ad <- density(df1[df1$class=='ACCEPTOR' & abs(df1$dist) < max_dist,]$dist, from=-max_dist, to=max_dist, n=2*max_dist)

    donor_density<-c(donor_density,dd$y)
    acceptor_density<-c(acceptor_density,ad$y)
    
    df0 <- data.frame(read=c(df0$read, df1$read), 
                      class=c(df0$class, df1$class), 
                      dist=c(df0$dist, df1$dist), 
                      sample=c(df0$sample, df1$sample),
                      repl=c(df0$repl, df1$repl),
                      stringsAsFactors=F
    )

    i <- i + 1
}

df0$sample <- factor(df0$sample, levels=c('enriched', 'bg'))
df0$class <- factor(df0$class, levels=c('DONOR', 'ACCEPTOR'))

ddensity_m <- matrix(donor_density, nrow=2*max_dist)
adensity_m <- matrix(acceptor_density, nrow=2*max_dist)

d_ci_lo <- c()
d_ci_hi <- c()
d_dens_mean <- c()

a_ci_lo <- c()
a_ci_hi <- c()
a_dens_mean <- c()


i<-1
while (i <= 2*max_dist) {
    d_dens_mean<- c(d_dens_mean, mean(ddensity_m[i,]))
    ci <- t.test(ddensity_m[i,])$conf.int
    d_ci_lo <- c(d_ci_lo,ci[1])
    d_ci_hi <- c(d_ci_hi,ci[2])

    a_dens_mean<- c(a_dens_mean, mean(adensity_m[i,]))
    ci <- t.test(adensity_m[i,])$conf.int
    a_ci_lo <- c(a_ci_lo,ci[1])
    a_ci_hi <- c(a_ci_hi,ci[2])
    
    i<-i+1
}

rm(df1, dd, ad,ci, acceptor_density, donor_density)

donor_df <- data.frame(x=c(seq(-max_dist,max_dist-1), seq(-max_dist,max_dist-1)), val=c(df_d_density, d_dens_mean), class=c(rep('enriched', 2*max_dist), rep('bg', 2*max_dist)))
acceptor_df <- data.frame(x=c(seq(-max_dist,max_dist-1), seq(-max_dist,max_dist-1)), val=c(df_a_density, a_dens_mean), class=c(rep('enriched', 2*max_dist), rep('bg', 2*max_dist)))



#plot(x=seq(-max_dist,max_dist-1), df_d_density, type='l', xlim=c(-500,2000), ylim=c(0,1e-3))
#lines(density(df0[df0$class=='DONOR' & abs(df0$dist) < max_dist & df0$sample=='bg',]$dist), type='l', xlim=c(-500,2000), col='red')
#lines(x=seq(-max_dist,max_dist-1), d_dens_mean, type='l', xlim=c(-500,2000), lty=2, col='blue')
#abline(v=0)


#lines(x=seq(-10000,10000-1), d_ci_hi, type='l', xlim=c(-500,2000), lty=2, col='grey')
#lines(x=seq(-10000,10000-1), d_ci_lo, type='l', xlim=c(-500,2000), lty=2, col='grey')

#plot(x=seq(-max_dist,max_dist-1), df_a_density, type='l', xlim=c(-2000,500), ylim=c(0,1e-3))

#lines(density(df0[df0$class=='ACCEPTOR' & abs(df0$dist) < max_dist & df0$sample=='bg',]$dist), type='l', xlim=c(-2000,500), col='red')
#lines(x=seq(-max_dist,max_dist-1), a_dens_mean, type='l', xlim=c(-2000,500), lty=2, col='blue')
#abline(v=0)
#lines(x=seq(-10000,10000-1), a_ci_hi, type='l', xlim=c(-500,2000), lty=2, col='grey')
#lines(x=seq(-10000,10000-1), a_ci_lo, type='l', xlim=c(-500,2000), lty=2, col='grey')



# restrict to w/in 10kb or a splice site

p <- ggplot(donor_df, aes(x=x, y=val, linetype=class)) + 
    geom_line() +
    scale_linetype_manual(values=c(2,1), labels=c('enriched', 'bg'), breaks=c('enriched', 'bg')) +
    coord_cartesian(xlim=c(-500, 2000), ylim=c(0,1.2e-3)) + 
    geom_vline(xintercept=0)+
    labs(title="Distance to splice donor", x="Distance (bp)", y="Density") + 
    theme_bw() + 
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank())
p

ggsave("figures/donor-dist.pdf", plot=p, height=4, width=4)


p <- ggplot(acceptor_df, aes(x=x, y=val, linetype=class)) + 
    geom_line() + 
    scale_linetype_manual(values=c(2,1), labels=c('enriched', 'bg'), breaks=c('enriched', 'bg')) +
    coord_cartesian(xlim=c(-2000, 500), ylim=c(0,1.2e-3)) + 
    geom_vline(xintercept=0)+
    labs(title="Distance to splice acceptor", x="Distance (bp)", y="Density") + 
    theme_bw() + 
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank())
p
ggsave("figures/acceptor-dist.pdf", plot=p, height=4, width=4)

# p <- ggplot(tmp, aes(x=dist, linetype=sample)) + 
#     geom_line(stat="density") + 
#     coord_cartesian(xlim=c(-1000, 1000)) + 
#     geom_vline(xintercept=0)+
#     labs(title="Distance to splice site", x="Distance (bp)", y="Density") + 
#     theme_bw() + 
#     theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank()) +
#     facet_grid(. ~ class)
# p

# ggsave("figures/donor-acceptor-dist.pdf", plot=p, height=4, width=4)

#ks.test(donors[donors$sample=='test',]$dist, donors[donors$sample=='bg',]$dist, )
#ks.test(acceptors[acceptors$sample=='test',]$dist, acceptors[acceptors$sample=='bg',]$dist, )
