header = """<!DOCTYPE html>
<!--[if lt IE 7]>      <html class="no-js lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>         <html class="no-js lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>         <html class="no-js lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class="no-js"> <!--<![endif]-->
    
    <head>
        <title>ATTRACT Online</title>
        
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
                <li id="view-page1"><a id="nav-block-sequences" class="grid-icon">Sequences</a></li>              
                <li id="view-page2"><a id="nav-block-partners" class="puzzle-icon">Partners</a></li>
                <li id="view-page3"><a id="nav-block-symmetry" class="symmetry-icon">Symmetry</a></li>
                <li id="view-page4"><a id="nav-block-cryoem" class="cryo-icon">Cryo-EM density maps</a></li>
                <li id="view-page5"><a id="nav-block-reference" class="analysis-icon">Reference structure</a></li>                
                <li id="view-page6"><a id="nav-block-protocol" class="iteration-icon">Assembly protocol</a></li>
                <li id="view-page7"><a id="nav-block-computation" class="computation-icon">Computation</a></li>
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

def _assign_category(f, category, groupname, span = False):
  for member in f._membernames:
    ff = getattr(f, member)
    if not hasattr(ff, "group"): continue
    if ff.group != groupname: continue
    category.members.append(member)
    ff.group = None
    if span: ff.span = True

def webform(f, model=None,
 partnerslength=None, sequenceslength=None,
):
  if model is not None:
    if partnerslength is None:
      partnerslength = max(1,len(model.partners))
    if sequenceslength is None:
      sequenceslength = max(1,len(model.sequences))
  else:  
    if partnerslength is None: partnerslength = 1
    if sequenceslength is None: sequenceslength = 1  
  import copy
  f = copy.deepcopy(f)
        
  ### START sequences category
  f.sequences.length = sequenceslength      
  c = f.new_group("c_sequences", "category")
  c.page = 1
  c.icon = "grid-icon"
  c.title = "Protein sequences"
  c.categoryname = "sequences"
  c.description = "" #TODO
  c.members.append("sequences")   
  f.sequences.clonebutton = "Add sequence"
  f.sequences.clonelength = 50
  f.sequences.controltitle = "Assembly sequences"  
  for fpnr in range(f.sequences.length):
    fp = f.sequences[fpnr]
    fp.group = None
        
    ### START b_sequence block
    b = fp.new_group("b_sequence", "block")
    b.title = "Sequence origin"
    b.has_switch = False
    b.members.append("mode")
    b.members.append("proteincode")
    b.members.append("sequence")
    b.members.append("chain")
    ### END b_sequence block
    
    fp.subsequences.length = 5
    
  ### START partners category
  f.partners.length = partnerslength      
  c = f.new_group("c_partners", "category")
  c.page = 2
  c.icon = "puzzle-icon"
  c.title = "Docking partners"
  c.categoryname = "partners"
  c.description = """High-level definition of a component that is to be fitted
A CryoPartner has a single PDB, and/or a sequence mapping, and/or a DSSP assignment  
The PDB consists of one or more chains  
A CryoPartner has different options for how it is to be converted into CryoBodies  
- mode "single": 
  The PDB consists of a single chain. If a sequence range has been defined, the PDB must align to that range, else it will be auto-aligned
  The CryoPartner will be mapped onto a single CryoBody. 
- mode "multi"
  The PDB consists of multiple chains. If a sequence range has been defined, all chains together must align to that range, else it will be auto-aligned
  -relplace "fixed": The relative positions of the chains is fixed. The partner maps to a single CryoBody
  -relplace "guess": (guess) The relative positions of the chains are an initial estimate, they can re-adjust during refinement. The partner maps to one CryoBody per chain
  -relplace "estimate": (good estimate) Same as above, but a restraining force will enforce resemblance to the initial positioning
  -relplace "free": The CryoPartner will be split into independent CryoBodies, one per chain
- mode "SSE":  
  The PDB consists of a single chain. If a sequence range has been defined, the PDB must align to that range, else it will be auto-aligned
  The PDB will be auto-fragmented into secondary structure element bodies, based on the DSSP assignment if provided, or else automatically
- mode "ab_initio"  
  No PDB has been provided. Sequence and DSSP alignment must be provided. The sequence will be fragmented into secondary structure element bodies, based on the DSSP assignment.
  Idealized SSE fragment PDBs will be built ab initio.    
"""
  c.members.append("partners")   
  f.partners.clonebutton = "Add partner"
  f.partners.clonelength = 50
  f.partners.controltitle = "Assembly partners"  
  for fpnr in range(f.partners.length):
    fp = f.partners[fpnr]
    fp.group = None
    
    ### START b_struc block
    b = fp.new_group("b_struc", "block")
    b.title =  "Structure sources"
    b.has_switch = False
    b.members.append("mode")
    b.members.append("pdb")
    b.members.append("dssp")
    ### END b_struc block
        
    ### START b_mappings block
    for n in range(fp.sequencemappings.length):
      b = fp.new_group("b_mappings-%d" % n, "block")
      b.title =  "Sequence mappings %d" % (n+1)      
      b.has_switch = False          
      fp.sequencemappings[n].subsequences[None].span = True
      for nn in range(fp.sequencemappings[0].subsequences.length): #TODO: swap
	fp.sequencemappings[n].subsequences[0].name = "1" #TODO: fix
        b.members.append("sequencemappings[%d].subsequences[%d]" % (n, nn))
      b.members.append("sequencemappings[%d]" % n) #TODO: swap
    ### END b_mappings block
    
    ### START b_advanced block    
    #insert placeholder boolean to determine the state of the Advanced switch
    b = fp.new_group("b_advanced", "block")
    b.title =  "Advanced settings"
    fp._membernames.append("use_advanced")
    sw = f.cryozoom._members["pre_mini"].get_copy()        
    sw.name = "Use advanced settings"
    sw.type = "switch"
    fp._members["use_advanced"] = sw
    b.has_switch = True
    b.members.append("use_advanced")
    b.members.append("absplace")
    b.members.append("relplace")
    b.members.append("missing_loops")
    b.members.append("floppiness")
    ### END b_advanced block
  
  ### START symmetry category
  c = f.new_group("c_symmetry", "category")
  c.page = 3
  c.title = "Symmetry"
  c.icon = "symmetry-icon"
  c.categoryname = "symmetry"
  c.description = "" #TODO
  b = f.new_group("b_symmetry", "block")
  b.title = c.title
  b.members.append("symmetry")
  c.members.append("b_symmetry")
  ### END symmetry category
  
  ### START cryo-EM category
  c = f.new_group("c_cryoem", "category")
  c.page = 4
  c.title = "Cryo-EM data"
  c.icon = "cryoem-icon"
  c.categoryname = "cryoem"
  c.description = "" #TODO
  b = f.new_group("b_cryodata", "block")
  b.title = f.cryodata.name
  f.cryodata.name = ""
  b.members.append("cryodata")
  c.members.append("b_cryodata")
  b = f.new_group("b_simdata", "block")
  b.title = f.simdata.name
  f.simdata.name = ""
  b.members.append("simdata")
  c.members.append("b_simdata")  
  ### END cryo-EM category

  ### START reference category
  c = f.new_group("c_reference", "category")
  c.page = 5
  c.title = f.reference.name
  c.icon = "analysis-icon"
  c.categoryname = "reference"
  c.description = "If you want to compare the assembly result against a known reference structure, you can define it here"
  
  ### START sequencemapping blocks #TODO: swap 1
  for n in range(f.reference.sequencemappings.length):
    f.reference.sequencemappings[n].subsequences[None].span = True
    bname = "b_reference-%d" % n
    b = f.new_group(bname, "block")
    b.title = f.reference.sequencemappings[n].name
    b.members.append("reference.sequencemappings[%d].subsequences" % n) #TODO: swap 2
    b.members.append("reference.sequencemappings[%d]" % n) #TODO: swap 2    
    f.reference.sequencemappings[n].subsequences[0].name = "1" #TODO: fix
    c.members.append(bname)
  ### END sequencemapping blocks
  ### START reference block #TODO: swap 1
  b = f.new_group("b_reference", "block") #TODO: swap 1
  b.title = "Reference"
  b.members.append("reference")
  c.members.append("b_reference")      
  ### END reference block
  
  ### END reference category

  ### START protocol category
  c = f.new_group("c_protocol", "category")
  c.page = 6
  c.title = f.cryozoom.name
  c.icon = "iteration-icon"
  c.categoryname = "protocol"
  c.description = "Protocol settings"
  b = f.new_group("b_protocol", "block")
  b.title = "General settings"
  b.members.append("cryozoom.struc_initial")
  b.members.append("cryozoom.pre_mini")
  b.members.append("cryozoom.pre_frac")
  b.members.append("cryozoom.select")
  b.members.append("cryozoom.clone")
  b.members.append("cryozoom.ori")  
  b.members.append("cryozoom.trans")
  b.members.append("cryozoom.iterations")  
  c.members.append("b_protocol")  
  b = f.new_group("b_initial", "block")
  b.title = "Initial iteration"
  b.members.append("cryozoom.initial")
  c.members.append("b_initial")
  b = f.new_group("b_subsequent", "block")
  b.title = "Subsequent iterations"
  b.members.append("cryozoom.subsequent")  
  c.members.append("b_subsequent")  
  for fi in (f.cryozoom.initial, f.cryozoom.subsequent):
    fi.group = None
    b = fi.new_group("b_mc", "block")
    b.insert_at_member = True
    b.has_switch = True
    _assign_category(fi, b, "Monte Carlo", span = True)
    b.members.remove("mc")
    b.members.insert(0, "mc")
    ff = fi.mc
    ff.type = "switch"
  
  ### END protocol category
  
  ### START computation category
  c = f.new_group("c_computation", "category")
  c.page = 7
  c.icon = "computation-icon"
  c.title = "Computation"
  c.categoryname = "computation"
  c.description = ""
  _assign_category(f, c, "Computing and parallelization parameters", span = True)
  ### END computation category
    
  return f

import spyder.htmlform
def webserverform(webdict, form=None, spydertype=None):
  if spydertype is not None: form = spydertype._form()
  f = webform(
   form,
   partnerslength = 10,
   sequenceslength = 10,
  )  
  return f
  
def html(form, cgi,newtab=False):
  import attracthtmlform 
  html = attracthtmlform.htmlform(
   form=form, cgi=cgi, 
   header=header, footer=footer, header_indentation = 12, 
   newtab=newtab
  )
  return html
  