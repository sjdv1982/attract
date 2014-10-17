header = """<!DOCTYPE html>
<!--[if lt IE 7]>      <html class="no-js lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>         <html class="no-js lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>         <html class="no-js lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class="no-js"> <!--<![endif]-->
    
    <head>
        <title>ATTRACT NARefine Online</title>
        
        <meta charset="utf-8" />
        <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1" />
        <meta name="description" content="" />
        <meta name="viewport" content="width=device-width"/>
        
        <!-- Always load normalize CSS first -->
        <link rel="stylesheet" href="css/normalize.min.css">
        <link href='http://fonts.googleapis.com/css?family=Open+Sans:400,300,600' rel='stylesheet' type='text/css'>
        <link href='http://fonts.googleapis.com/css?family=Lobster' rel='stylesheet' type='text/css'>
        <link rel="stylesheet" href="css/attract.min.css">
    </head>
    
    <body>
        <!--[if lt IE 7]>
            <p class="chromeframe">You are using an <strong>outdated</strong> browser. Please <a href="http://browsehappy.com/">upgrade your browser</a> or <a href="http://www.google.com/chromeframe/?redirect=true">activate Google Chrome Frame</a> to improve your experience.</p>
        <![endif]-->

        <main id="container" class="col row">
          
          <section id="sidebar" class="col">
            <div id="title"><span class="t1">Attract</span> <span class="t2">online</span></div>
            <div id="logo"></div>
            
            <nav id="form-category-navigation" class="row">
              <ul>
                <li id="view-page1"><a id="nav-block-setup" class="puzzle-icon">Setup</a></li>
                <li id="view-page2"><a id="nav-block-simulation" class="sampling-icon">Simulation</a></li>
              </ul>
            </nav>
            
            <div id="download-button" class="row download-icon" onClick="submitForm();">
              <p>Get configuration</p>
            </div>  
            
          </section>
            
          <section id="content" class="col">
            
            <header id="header">
              
              <div id="show-hide-sidepanel" class="button menu-icon float-left"><div class="header-tooltip"><p>Hide sidepanel</p></div></div>
              <div id="reload-form" class="button reload-icon-active float-left"><div class="header-tooltip"><p>Restore default values</p></div></div>
              <div id="unfold-all" class="button double-arrow-icon float-left"><div class="header-tooltip"><p>Unfold all form blocks</p></div></div>
              <div id="header-nav" class="button-text down-icon-active float-left">
                <p>Menu</p>
                <nav class="header-tooltip">
                  <ul>
                    <li><a href="#" class="message-icon">Messages<span id="message-counter">2</span></a></li>
                    <li><a href="#" class="download-icon" onClick="submitForm();">Save</a></li>
                    <li><a href="#" class="contact-icon">Contact</a></li>
                    <li><a href="documentation.html" target="_blank" class="help-icon">Help</a></li>
                  </ul>
                </nav> 
              </div>
              <div id="close-app" class="button close-icon-active float-right"></div>
              
            </header>
"""

footer = """           
            
        </section>
      </main>  
        
      <script src="js/jquery-1.10.2.min.js"></script>
      <script src="js/main.min.js"></script>

    </body>
</html>"""

def webform(f, model=None):
  import copy
  f = copy.deepcopy(f)  
  f.resourcefilevar = "_tempresource"
  f.arraymarker = "_clonearraymarker"

  ### START setup category
  c = f.new_group("c_setup", "category")
  c.page = 1
  c.icon = "puzzle-icon"
  c.title = "Setup"
  c.categoryname = "setup"
  c.description = """Setup"""
  c.members.append("runname")
  c.members.append("pdbfile")
  c.members.append("offset")
  c.members.append("use_simulation")
  c.members.append("min2_cycles")
  c.members.append("minc_cycles")
  c.members.append("nparallel")
  c.members.append("use_cuda")
  c.members.append("use_5PO")

  ### START simulation category
  c = f.new_group("c_simulation", "category")
  c.page = 2
  c.icon = "sampling-icon"
  c.title = "Setup"
  c.categoryname = "simulation"
  c.description = """Simulation"""
  
  b = f.new_group("b_md1", "block")
  b.title = "Simulation 1"
  b.has_switch = False 
  b.members.append("md1")
  c.members.append("b_md1")

  b = f.new_group("b_md2", "block")
  b.title = "Simulation 2"
  b.has_switch = False 
  b.members.append("md2")
  c.members.append("b_md2")
  
  return f

def webserverform(webdict, form=None, spydertype=None):
  if spydertype is not None: form = spydertype._form()
  f = webform(form)
  return f
  
def html(form, cgi, spyderobj, newtab=False):
  from form_model import html
  return html(form, cgi, spyderobj, newtab, header=header)