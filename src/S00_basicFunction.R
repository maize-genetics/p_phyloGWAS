#basic function
cont_table=function(query,background,classifyer){
  p1=length(intersect(query,classifyer))
  q1=length(query)-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}