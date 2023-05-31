errbar = function(x, y, bar, double = T, width=0.01,col="black")
{
  segments(x0 = x, y0 = y,
           x1 = x, y1 = y + bar, col=col)
  segments(x0 = x - width, y0 = y + bar,
           x1 = x + width, y1 = y + bar,col=col)
  if(double)
  {
    segments(x0 = x, y0 = y,
             x1 = x, y1 = y - bar,col=col)
    segments(x0 = x - width, y0 = y - bar,
             x1 = x + width, y1 = y - bar,col=col)
  }
}