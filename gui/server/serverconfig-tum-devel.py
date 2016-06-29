###########################################
# Change the following directories
website = "http://www.attract.ph.tum.de"
localserverdir = "/home/server/services"
servicename = "ATTRACT-devel"
###########################################

webdir = website + "/services/%s/" % servicename
cgidir = website + "/cgi/services/%s/" % servicename
webresultdir = website + "/results/services/%s/" % servicename
localdir = localserverdir + "/%s/html/" % servicename
localresultdir = localserverdir + "/%s/results/" % servicename
