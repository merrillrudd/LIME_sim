print.letter <- function(label="(a)",xy=c(0.1,0.925),...) { 
            tmp <- par("usr") 
            text.x <- tmp[1]+xy[1]*diff(tmp[1:2]) #x position, diff=difference 
            text.y <- tmp[3]+xy[2]*diff(tmp[3:4]) #y position 
            text(x=text.x, y=text.y, labels=label,...) 
}