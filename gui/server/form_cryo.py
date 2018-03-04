header = """<!DOCTYPE html>
<!--[if lt IE 7]>      <html class="no-js lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>         <html class="no-js lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>         <html class="no-js lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class="no-js"> <!--<![endif]-->
    
    <head>
        <title>ATTRACT-EM Online</title>
        
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
                <li id="view-page1"><a id="nav-block-partners" class="puzzle-icon">Partners</a></li>
                <li id="view-page2"><a id="nav-block-symmetry" class="symmetry-icon">Symmetry</a></li>
                <li id="view-page3"><a id="nav-block-cryo" class="cryo-icon">Cryo-EM data</a></li>  
                <li id="view-page4"><a id="nav-block-sampling" class="sampling-icon">Sampling</a></li>
                <li id="view-page5"><a id="nav-block-computation" class="computation-icon">Computation</a></li>
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
                    <li><a href="#" class="download-icon" onClick="submitForm();">Save</a></li>
                    <li><a href="documentation.html#contact" target="_blank" class="contact-icon">Contact</a></li>
                    <li><a href="documentation.html" target="_blank" class="help-icon">Help and documentation</a></li>
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

def _assign_category(f, category, groupname, span = False):
  for member in f._membernames:
    ff = getattr(f, member)
    if not hasattr(ff, "group"): continue
    if ff.group != groupname: continue
    category.members.append(member)
    ff.group = None
    if span: ff.span = True

def webform(f, model=None):
  import copy
  f = copy.deepcopy(f)  
  f.resourcefilevar = "_tempresource"
  f.arraymarker = "_clonearraymarker"

  if model is not None:
      partnerslength = max(1,len(model.partners))
      symmetrieslength = max(1,len(model.axsymmetry))
  else:
      partnerslength = 1
      symmetrieslength = 1
  f.partners.length = partnerslength
  f.axsymmetry.length = symmetrieslength
  f.partners.clonebutton = "Add partner"
  f.partners.clonelength = 10

  ### START partners category
  c = f.new_group("c_partners", "category")
  c.page = 1
  c.icon = "puzzle-icon"
  c.title = "Docking partners"
  c.categoryname = "partners"
  c.description = """Define docking partners by uploading a PDB structure file.""" 
  c.members.append("partners")   
  f.partners.clonelength = 6
  f.partners.controltitle = "Docking partners"  
  for fpnr in range(f.partners.length):
    fp = f.partners[fpnr]
    fp.group = None
    fp.multi_active = True

    ### START b_struc block
    b = fp.new_group("b_struc", "block")
    b.title = "Structure Sources"
    b.has_switch = False
    b.members.append("pdbfile") 
    ff = fp.pdbfile
    ff.name = "Structure file"
    ff.tooltip = "Upload PDB structure file"
    ff.tooltip_doc = "documentation.html#partners-structure_file"
    ff.span = True
    ### END b_struc block


  ### START symmetry category
  c = f.new_group("c_symmetry", "category")
  c.page = 2
  c.icon = "symmetry-icon"
  c.title = "Generative symmetry"
  c.description = c.title
  c.always_active = True
  c.categoryname = "symmetry"  
  c.members.append("axsymmetry")
  f.axsymmetry.clonebutton = "Add symmetry"
  f.axsymmetry.clonelength = 6
  f.axsymmetry.controltitle = "Symmetry"    
  ### END symmetry category  

  ### START cryo category
  c = f.new_group("c_cryo", "category")
  c.page = 3
  c.icon = "cryo-icon"
  c.title = "Cryo-EM data"
  c.description = c.title
  c.always_active = True
  c.categoryname = "cryo"  
  _assign_category(f, c, c.title, span = True)
  f.harmonic_restraints_file.type = None 
  ### END cryo category  

  ### START sampling category
  c = f.new_group("c_sampling", "category")
  c.page = 4
  c.icon = "sampling-icon"
  c.title = "Sampling"
  c.description = c.title
  c.always_active = True
  c.categoryname = "sampling"  
  _assign_category(f, c, c.title, span = True)
  ### END sampling category  

  ### START computation category
  c = f.new_group("c_computation", "category")
  c.page = 5
  c.icon = "computation-icon"
  c.title = "Computation"
  c.description = c.title
  c.always_active = True
  c.categoryname = "computation"  
  _assign_category(f, c, c.title, span = True)
  ### END computation category  
  
  return f

def webserverform(webdict, form=None, spydertype=None):
  if spydertype is not None: form = spydertype._form()
  f = webform(form)
  return f
  
def html(form, cgi, spyderobj, newtab=False):
  from form_model import html
  return html(form, cgi, spyderobj, newtab, header=header)
