###########################################
# Change the following directory
localserverdir = "/home/sjoerd/server"
servicename = "ATTRACT"
###########################################

website = "http://localhost"
webdir = website + "/server/%s/" % servicename
cgidir = website + "/cgi/server/%s/" % servicename
webresultdir = website + "/results/server/%s/" % servicename
localdir = localserverdir + "/%s/html/" % servicename
localresultdir = localserverdir + "/%s/results/" % servicename
