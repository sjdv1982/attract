How to set up an ATTRACT web server

Step 1. Develop your web interfaces

While developing your web forms, don't test them in your browser but use Silk instead.
Only set up a real web server once it is ready.

Step 2. Setting up Apache

- Install Apache (e.g. "sudo apt-get install apache2")
- Create a directory for your web server(s), e.g. /home/user/server
- Add the following lines to your Apache configuration file 
(e.g. /etc/apache2/sites-enabled/000-default )

<Directory /home/user/server>
        Options Indexes FollowSymLinks
        AllowOverride All
        Order allow,deny
        Allow from all
</Directory>
AliasMatch ^/server/([^/]*)/(.*) /home/user/server/$1/html/$2
AliasMatch ^/results/server/([^/]*)/(.*) /home/user/server/$1/results/$2
ScriptAliasMatch ^/cgi/server/([^/]*)/(.*) /home/user/server/$1/cgi/$2


- Restart Apache, e.g. "sudo /etc/init.d/apache2 restart"

Step 3. Installing the ATTRACT server

- Go to $ATTRACTGUI,
- Choose one of the generate-html*.sh scripts.  
generate-html-tum.sh sets up a TUM web server (www.attract.ph.tum.de)
generate-html-local.sh sets up a web server for local use (localhost). 
- Edit serverconfig-tum.py or serverconfig-local.py 
- Run your chosen generate-html scripts
- Create a directory <servicename> inside your web server directory, e.g /home/user/server/ATTRACT
- Go to this directory
- Create (or softlink) a results directory, e.g /home/user/server/ATTRACT/results
  
  !!! Make sure that this directory is writeable to everyone !!!  

- Finally, create the following links: 
  ln -s $ATTRACTGUI/server/cgi-local cgi
  ln -s $ATTRACTGUI/server/html-local html  (generate-html-local.sh)
 OR
  ln -s $ATTRACTGUI/server cgi
  ln -s $ATTRACTGUI/server/html html (generate-html-tum.sh)

Done! 
If you set it up for local use, the web interfaces will be now available under http://localhost/server/ATTRACT/XXX
  XXX can be easy.html, full.html or cryo.html
  
Step 4. Maintaining the ATTRACT server

- Result directories will accumulate in the results directory (by default in /home/user/server/ATTRACT/results)
These originate from the generation of parameter files and shell scripts.
You may want to delete these directories after a few weeks
- Temporary file directories will accumulate under /tmp/attractrunXXXX . 
These originate from the editing of existing parameter files (embedded resources)
If you don't restart your server, you will want to delete these periodically

One way to periodically run cleaning jobs is with CRON. An example cron job "clean.cron" has been provided.
To run this script periodically, add it to your crontab using the syntax shown in "clean.crontab". If you don't have a crontab,
(type "crontab -l" to verify this), simply set "clean.crontab" as your crontab with "crontab clean.crontab"