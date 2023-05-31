asterisk.default = function(pvalue)
{
	if(pvalue <= 0.05 & pvalue > 0.001)
	return("*")
	if(pvalue <= 0.001 & pvalue > 0.0001)
	return("**")
	if(pvalue < 0.0001)
	return("***")
	if(pvalue >= 0.05)
	return("NS")
}

asterisk=function(pvalue)
{
	sapply(pvalue, asterisk.default)
}